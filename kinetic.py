# import packges
import pandas as pd
import numpy as np
import math as m
import thermodynamic as thermo

R = 8.314472             # ideal gas universal constant (J/mol.K)
T_ref = 300+273.15       # reference temperature (K)

def reaction_constant(T, kinetic_info):

    """calculates the reaction's overall kinetic constant using the reparametrized
    Arrhenius model

    inputs: temperature in Kelvin, dictionary of data containing kinetic
             information
    outputs: overall reaction kinetic constants
    """

    # create empty list to store reaction constants
    overall_k = []

    # loop through each reaction 
    for reaction in kinetic_info['kref']:

        # extract reference constant and activation energy
        k_ref = kinetic_info['kref'][reaction]
        Ea = kinetic_info['E_a'][reaction]*1000

        # calculate the constant corrected by temperature
        overall_k.append(k_ref*np.exp((Ea/R)*((1/T_ref)-(1/T))))

    # return calculated values
    return overall_k

def adsorption_K(T, comp_dict):

    """calculates the corrected adsorption equilibrium constant using the
    Vant'Hoff equation

    inputs: temperature in Kelvin, dictionary containing chemical compound
             information
    outputs: dictionary containing the corrected adsorption equilibrium
              constant for each compound
    """

    # creat empty dictionary to store equilibrium constants
    ads_K = dict()

    # loop through each chemical compound's information
    for comp in comp_dict['Kref_ads']:
        
        # extract the reference equilibrium constants and the variation of
        # enthaply due to adsorption
        Kref = comp_dict['Kref_ads'][comp]
        deltaHads = comp_dict['deltaHads'][comp]

        # calculate the corrected equilibrium constant
        ads_K[comp] = Kref*np.exp((deltaHads/R)*((1/T_ref)-(1/T)))

    # return calculated values
    return ads_K

def reaction_rates(T, P, F_dict, initial_cond, comp_dict, kinetic_dict, single = True):

    """calculates the reaction rates for all chemical compounds in the system

    inputs: temperature in Kelvin, pressure in bar, dictionary of number of
             mols of each component of the mixture, initial conditions list,
             dictionary with chemical compounds information, 
             dictionary with kinetic information, flag to indicate which model
             (single site or dual site) will be applied
    outputs: dictionary with reaction rates for all chemical components of 
              the system
    """

    # calculate the total molar flow rate
    Ft = 0
    for comp in F_dict:
        Ft += F_dict[comp]

    # calculate the correction to volumetric flow rate
    v0 = initial_cond[0]             # initial volumetric flow rate (m3/s)
    T0 = initial_cond[1]             # initial temperature (K)
    P0 = initial_cond[2]             # initial pressure (Pa)
    Ft0 = initial_cond[3]            # initial total molar flow rate (mol/s)

    v = v0*(P0/P)*(T/T0)*(Ft/Ft0)    # corrected volumetric flow rate (m3/s)

    # calculate the partial pressures of each component
    Pi = dict()

    for comp in F_dict:
        Pi[comp] = (F_dict[comp]*R*T)/v

    # calculate reaction constants
    k_list = reaction_constant(T, kinetic_info=kinetic_dict)

    # calculate adsorption equilibrium constants
    K_ads = adsorption_K(T, comp_dict=comp_dict)

    # calculate thermodynamic chemical equilibrium constants
    K_eq = dict()
    for i in range(len(k_list)):
        K_eq[i+1] = thermo.calculate_K_and_Hr(comp_list=comp_dict['comp_list'],
        stoic_list=kinetic_dict['stoichiometry'][str(i+1)], T = T, prop_dict=comp_dict)[0]

    # calculate the reaction rates based on the Langmuir-Hinshelwood mechanism
    r_rates = dict()

    # based on site model choice, calculate inhibition factor
    if single:
        inh_factor = (1+(K_ads['CO2']*Pi['CO2'])+np.sqrt(K_ads['H2']*Pi['H2']))**2
    else:
        inh_factor = (1+(K_ads['CO2']*Pi['CO2']))*(1+np.sqrt(K_ads['H2']*Pi['H2']))
    
    # methanol formation reaction rate
    r_rates[1] = ((k_list[0]*((Pi['CO2']*(Pi['H2']**3))-((Pi['MeOH']*Pi['H2O'])/K_eq[1])))/(Pi['H2']**2))/inh_factor

    # reverse water gas shift reaction rate
    r_rates[2] = ((k_list[1]*((Pi['CO2']*Pi['H2'])-((Pi['CO']*Pi['H2O'])/K_eq[2])))/(np.sqrt(Pi['H2'])))/inh_factor

    # CO2 methanation reaction reaction rate
    r_rates[3] = ((k_list[2]*np.sqrt(Pi['CO2'])*np.sqrt(Pi['H2']))*(1-((Pi['CH4']*(Pi['H2O']**2)/(Pi['CO2']*(Pi['H2']**4)*K_eq[3])))))/inh_factor

    # create dataframe with stoichiometry for all components
    reaction_df = pd.DataFrame(kinetic_dict['stoichiometry']['1'], columns=['methanol'],
    index = [comp_dict['comp_list']])
    reaction_df['rwgs'] = kinetic_dict['stoichiometry']['2']
    reaction_df['methanation'] = kinetic_dict['stoichiometry']['3']

    # create columns at reaction datafram with the reaction rates
    reaction_df['methanol_rate'] = r_rates[1]
    reaction_df['rwgs_rate'] = r_rates[2]
    reaction_df['methanation_rate'] = r_rates[3]

    # calculate global rate for each compound
    reaction_df['global_rate'] = 0
    for i in range(reaction_df.shape[0]):
        for j in range(3):
            reaction_df.iloc[i,6] += reaction_df.iloc[i,j]*reaction_df.iloc[i,j+3]

    # return calculated reaction rates
    return reaction_df