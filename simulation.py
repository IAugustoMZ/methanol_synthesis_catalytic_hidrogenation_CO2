from chemical_properties import IdealGasCp
import pandas as pd
import numpy as np
import math as m
import thermodynamic as thermo
import chemical_properties as cp
import kinetic as kn

# constant definition
phi = 0.92                    # catalyst bed void fraction
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

    print(term1)
    print(term2)
    print(term3)

    # calculte the pressure drop
    dP = term3*(term1 + term2)

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

    # calculate the total inlet flow rate (m³/s)
    v0 = WHSV*W*(1000/(3600*1000000))

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

def calculate_k(F_dict, comp_dict, T, P, reactor_setup, m_in):

    """calculates the Runge-Kutta k's
    inputs: dictionary containing molar flow rates of each compounds,
            dictionary containing chemical and physical properties,
            temperature (K), pressure (bar), dictionary containing
            reactor information, feed mass volumetric rate (kg/s)
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

    # step 3 - pressure drop - Ergun's equation 
    dP = ergun(rho_mix, m_in, mu_mix, reactor_setup['D']/2)

    return dP

def reactor_dimensions(reactor_setup):

    """calculate the reactor dimensions based on the reactor setup
    inputs: dictionary containing reactor setup
    output: dictionary containing reactor setup
    """

    # extract necessary parameters
    V = reactor_setup['volume']
    k = reactor_setup['LtoD']

    # calculate reactor dimensions
    reactor_setup['D'] = ((4*V)/(m.pi*k))**(1/3) # reactor diameter
    reactor_setup['L'] = k*reactor_setup['D']