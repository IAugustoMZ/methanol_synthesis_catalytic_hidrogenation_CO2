from chemical_properties import IdealGasCp
import pandas as pd
import numpy as np
import math as m
import thermodynamic as thermo
import chemical_properties as cp
import kinetic as kn

# constant definition
phi = 0.92                   # catalyst bed void fraction
Dp = 16.37*(10**(-6))         # catalyst average particle diameter (m)
rhoc = 7180                   # catalyst density (kg/m³)
R = 8.314472                  # ideal gas universal constant (J/mol.K)

def ergun(rho, mf, mu, R):

    """
    calculates the pressure drop from a axial flow rate using Ergun's equation
    
    inputs: fluid density, total mass flow rate and the fluid dynamic viscosity, reactor radius
    outputs: the pressure drop in the interval
    """

    # calculate the axial flow velocity
    u = mf/(rho*(m.pi*(R**2)))

    # calculate the three terms of Ergun's correlation
    term1 = (150*(1-phi)*mu)/Dp
    term2 = 1.75*rho*u
    term3 = -(u/Dp)*((1-phi)/(phi**3))

    print(rho, mu)

    # calculte the pressure drop
    dP = term3*(term1 + term2)

    print(dP)

    return dP

def energy_balance(data_dict, react_info, T, P, F_dict, initial_condition,single_model=True):

    """
    calculates the energy balance of the reaction system

    inputs: dictionary containing chemical data, dicitionary containing reaction data, temperature of the system
            total pressure of the system, dictionary containing the molar flow rates of each component
            initial conditions list, flag to indicate which site model is to be used
    outputs: the temperature variation in the evaluated interval
    """

    # extract the list of components of the system
    comp_list = data_dict['comp_list']

    # calculate the heat of reaction based on the temperature and equilibria constants
    for reaction in react_info['stoichiometry']:

        # extracting the stoichiometric information from each reaction
        nu_list = react_info['stoichiometry'][reaction]
        K, Hr = thermo.calculate_K_and_Hr(comp_list=comp_list,stoic_list=nu_list, T=T, prop_dict=data_dict)

        # inserting calculated values at dictionary data
        react_info['equilibrium_constant'][reaction] = K
        react_info['heat_reaction'][reaction] = Hr

    # calculate reaction kinetic constants
    rate = kn.reaction_rates(T=T, P=P, F_dict=F_dict, 
                            initial_cond=initial_condition, 
                            comp_dict=data_dict,
                            kinetic_dict=react_info,
                            single=single_model)

    # calculate the energy transport term sum for all compounds
    energy_flux = 0
    for i in range(len(F_dict)):

        # calculate the Cp for each compound and convert to J/mol.K
        Cp = cp.IdealGasCp(comp_list[i], T, data_dict)/1000

        # calculate the energy flux (J/s.K)
        energy_flux += F_dict[comp_list[i]]*Cp

    # calculate the generation term
    methanol_rate = rate['methanol_rate'].head(1).values
    rwgs_rate = rate['rwgs_rate'].head(1).values
    methanation_rate = rate['methanation_rate'].head(1).values
    energy_gen = -(react_info['heat_reaction']['1']*methanol_rate + \
                    react_info['heat_reaction']['2']*rwgs_rate + \
                    react_info['heat_reaction']['3']*methanation_rate)

    # calculate the temperature variation
    return energy_gen[0]/energy_flux    


def prepare_data(T0, P0, WHSV, reactor_setup, h2toco2, comp_dict):

    """
    prepares the data to initialize reactor simulation
    inputs: feed temperature, feed stream pressure,
            weight hourly space velocity, dictionary containing reactor setup,
            h2toco2 molar ratio at feed stream, dictionary containing 
            chemical data
    outputs: initial conditions list and initial flow rates for all chemical
            components
    """

    # extract the volume of the fixed bed (m3)
    V = reactor_setup['volume']

    # calculate the mass of catalyst (kg)
    W = rhoc*V*(1-phi)

    # calculate the total inlet flow rate (m³/s) in each tube
    v0 = WHSV*W*(1000/(3600*1000000))/reactor_setup['N_tubes']

    # calculate the total molar flow rate (mol/s) using ideal gas law
    P0 = P0*(10**5)                  # convert P in bar to Pa
    T0 = T0 +273.15                  # convert T in degC to K   
    F0 = (P0*v0)/(R*T0)              

    # calculate the initial flow rate (mol/s) for all chemical compounds
    F_dict = {}

    for comp in comp_dict['comp_list']:
        if comp == 'CO2':
            F_dict[comp] = (1/(1+h2toco2))*F0
        elif comp == 'H2':
            F_dict[comp] = (h2toco2/(1+h2toco2))*F0
        else:
            F_dict[comp] = 0

    # return the initial conditions and initial flow rates
    return [v0, T0, P0, F0], F_dict

def material_balance(T, P, F_dict, initial_condition, data_dict, react_info, single_model):

    """calculates the material balance in the reactor
    inputs: temperature (K), pressure (bar), initial conditions list, 
            dictionary containing chemical properties, 
            dictionary containing kinetic information
            option to which site model is to use
    outputs: the molar flux variation for each component
    """

    # calculate reaction kinetic constants
    rate = kn.reaction_rates(T=T, P=P, F_dict=F_dict, 
                            initial_cond=initial_condition, 
                            comp_dict=data_dict,
                            kinetic_dict=react_info,
                            single=single_model)

    # extract the global rates of reaction
    rates_j = rate['global_rate']

    return rates_j

def calculate_k(F_dict, comp_dict, T, P, reactor_setup, m_in, react_dict, initial_condition, single_model=True):

    """calculates the Runge-Kutta k's
    inputs: dictionary containing molar flow rates of each compounds,
            dictionary containing chemical and physical properties,
            temperature (K), pressure (bar), dictionary containing
            reactor information, feed mass volumetric rate (kg/s),
            dictionary containin reaction information,
            list containing initial conditions information
    outputs: list of Runge-Kutta k's
    """

    # step 1 - total molar flow rate
    Ft = 0 
    for comp in F_dict:
        Ft += F_dict[comp]
    
    # step 2 - molar fraction
    y = []
    for comp in F_dict:
        y.append(F_dict[comp]/Ft)

    # step 3 - physical and chemical properties
    # 3.1. mixture density (kg/m³)
    rho_mix = cp.mixture_density(comp_list=comp_dict['comp_list'],
                                mol_fraction=y,
                                T = T, P = P, data_dict=comp_dict)

    # 3.2. mixture vapor viscosity (Pa.s)
    mu_list = []
    for comp in comp_dict['comp_list']:
        mu_list.append(cp.vapor_viscosity(comp, T, comp_dict))
    MM_list = [comp_dict['molar_mass'][comp] for comp in comp_dict['molar_mass']]
    mu_mix = cp.mixture_viscosity(mu_list, MM_list=MM_list, molar_fraction=y)

    # step 3 - pressure drop - Ergun's equation (bar)
    # dP = ergun(rho_mix, m_in, mu_mix, reactor_setup['D']/2)/100000
    dP = 0

    # step 4 - temperature variation by energy balance (K)
    dT = energy_balance(data_dict=comp_dict, react_info=react_dict, T = T, 
                        P = P, F_dict=F_dict, initial_condition=initial_condition)

    # step 5 - material balance (mol/s)
    dF = material_balance(T=T, P=P, F_dict=F_dict, initial_condition=initial_condition,
                          data_dict=comp_dict, react_info=react_dict, single_model=single_model)

    return [dP, dT, dF['CO2'], dF['H2'], dF['MeOH'], dF['H2O'], dF['H2O'], dF['CO'], dF['CH4']]

def runge_kutta(initial_conditions, F_dict, m_rate, reactor_setup, comp_dict, react_info, N = 1000, single_model = True):

    """
    executes the integration method of 4th order Runge Kutta
    inputs: list containining initial conditions, dictionary containing the molar flow rate of each component, 
            mass flow rate, dictionary containing the reactor setup information, dictionary containing chemical data,
            dictionary containing kinetic information, number of integration points in the grid,
            flag selecting the inhibition model to be used
    # outputs: dataframe containing pressure, temperature and composition
    """

    # extract information about reactor lenght and diameter
    L = reactor_setup['L']
    D = reactor_setup['D']

    # create integration mesh
    l = np.linspace(0, L, N)
    
    # calculate integration step
    h = l[1] - l[0]

    # calculate conversion factor from dw to dx
    dw_to_dl = (m.pi*(D**2)*rhoc*(1-phi))/4

    # calculate the mass flow rate per tube
    m_in = m_rate

    # output dataframe definition
    df_cols = ['x', 'T_K', 'P_bar']
    for comp in comp_dict['comp_list']:
        df_cols.append('F_'+comp)
    profile_df = pd.DataFrame(columns = df_cols)

    # set initial conditions
    T = initial_conditions[1]
    P = initial_conditions[2]/100000 # convert from Pa to bar

    # loop through the integration grid
    for i in range(N):

        # capture coordinate
        x = l[i]

        # create copies of the dictionaries of molar flow rates
        F1_dict = F2_dict = F3_dict = F_dict

        # store values in the dataframe
        profile_df.loc[i, df_cols] = x, T, P, F_dict['CO2'], F_dict['H2'], F_dict['MeOH'], F_dict['H2O'], F_dict['CO'], F_dict['CH4']

        # calculate k1 ---------------------------------------------------------------------------------------------------------------------
        k_list = calculate_k(F_dict, comp_dict, T, P, reactor_setup, m_in, react_info, initial_conditions, single_model)
        k1P = k_list[0]
        k1T = k_list[1]
        k1CO2 = k_list[2]
        k1H2 = k_list[3]
        k1MeOH = k_list[4]
        k1H2O = k_list[5]
        k1CO = k_list[6]
        k1CH4 = k_list[7]

        # update variables
        P1 = P + k1P*(h/2)
        T1 = T + k1T*(h/2)
        F1_dict['CO2'] += k1CO2*(h/2)
        F1_dict['H2'] += k1H2*(h/2)
        F1_dict['MeOH'] += k1MeOH*(h/2)
        F1_dict['H2O'] += k1H2O*(h/2)
        F1_dict['CO'] += k1CO*(h/2)
        F1_dict['CH4'] += k1CH4*(h/2)

        # calculate k2 ---------------------------------------------------------------------------------------------------------------------
        k_list = calculate_k(F1_dict, comp_dict, T1, P1, reactor_setup, m_in, react_info, initial_conditions, single_model)
        k2P = k_list[0]
        k2T = k_list[1]
        k2CO2 = k_list[2]
        k2H2 = k_list[3]
        k2MeOH = k_list[4]
        k2H2O = k_list[5]
        k2CO = k_list[6]
        k2CH4 = k_list[7]

        # update variables
        P2 = P + k2P*(h/2)
        T2 = T + k2T*(h/2)
        F2_dict['CO2'] += k2CO2*(h/2)
        F2_dict['H2'] += k2H2*(h/2)
        F2_dict['MeOH'] += k2MeOH*(h/2)
        F2_dict['H2O'] += k2H2O*(h/2)
        F2_dict['CO'] += k2CO*(h/2)
        F2_dict['CH4'] += k2CH4*(h/2)

        # calculate k3 ---------------------------------------------------------------------------------------------------------------------
        k_list = calculate_k(F2_dict, comp_dict, T2, P2, reactor_setup, m_in, react_info, initial_conditions, single_model)
        k3P = k_list[0]
        k3T = k_list[1]
        k3CO2 = k_list[2]
        k3H2 = k_list[3]
        k3MeOH = k_list[4]
        k3H2O = k_list[5]
        k3CO = k_list[6]
        k3CH4 = k_list[7]

        # update variables
        P3 = P + k3P*h
        T3 = T + k3T*h
        F3_dict['CO2'] += k3CO2*h
        F3_dict['H2'] += k3H2*h
        F3_dict['MeOH'] += k3MeOH*h
        F3_dict['H2O'] += k3H2O*h
        F3_dict['CO'] += k3CO*h
        F3_dict['CH4'] += k3CH4*h

        # calculate k4 ---------------------------------------------------------------------------------------------------------------------
        k_list = calculate_k(F3_dict, comp_dict, T3, P3, reactor_setup, m_in, react_info, initial_conditions, single_model)
        k4P = k_list[0]
        k4T = k_list[1]
        k4CO2 = k_list[2]
        k4H2 = k_list[3]
        k4MeOH = k_list[4]
        k4H2O = k_list[5]
        k4CO = k_list[6]
        k4CH4 = k_list[7]

        # update functions
        T += ((h/6)*(k1T + (2*k2T) + (2*k3T) + k4T))
        P += ((h/6)*(k1P + (2*k2P) + (2*k3P) + k4P))
        F_dict['CO2'] += dw_to_dl*((h/6)*(k1CO2 + (2*k2CO2) + (2*k3CO2) + k4CO2))
        F_dict['H2'] += dw_to_dl*((h/6)*(k1H2 + (2*k2H2) + (2*k3H2) + k4H2))
        F_dict['MeOH'] += dw_to_dl*((h/6)*(k1MeOH + (2*k2MeOH) + (2*k3MeOH) + k4MeOH))
        F_dict['H2O'] += dw_to_dl*((h/6)*(k1H2O + (2*k2H2O) + (2*k3H2O) + k4H2O))
        F_dict['CO'] += dw_to_dl*((h/6)*(k1CO + (2*k2CO) + (2*k3CO) + k4CO))
        F_dict['CH4'] += dw_to_dl*((h/6)*(k1CH4 + (2*k2CH4) + (2*k3CH4) + k4CH4))

    return profile_df


def reactor_dimensions(reactor_setup):

    """calculate the reactor dimensions based on the reactor setup
    inputs: dictionary containing reactor setup
    output: dictionary containing reactor setup
    """

    # extract necessary parameters
    V = reactor_setup['volume']/reactor_setup['N_tubes']
    k = reactor_setup['LtoD']

    # calculate reactor dimensions
    reactor_setup['D'] = ((4*V)/(m.pi*k))**(1/3) # reactor diameter
    reactor_setup['L'] = k*reactor_setup['D']