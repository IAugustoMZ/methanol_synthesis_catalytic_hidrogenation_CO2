# import libraries
import pandas as pd
import numpy as np

R = 8.314472 # ideal gas universal constant (J/mol.K)

def IdealGasCp(comp_name, T, data_dict):

    # calculates the ideal gas heat capacity

    # input: name of component, temperature in Kelvin, dictionary of physical and chemical properties
    # output: ideal gas heat capacity (J/kmol.K)

    # extracting correlation coefficients from dictionary
    C1 = data_dict['IdealGasCp'][comp_name]['C1']
    C2 = data_dict['IdealGasCp'][comp_name]['C2']
    C3 = data_dict['IdealGasCp'][comp_name]['C3']
    C4 = data_dict['IdealGasCp'][comp_name]['C4']
    C5 = data_dict['IdealGasCp'][comp_name]['C5']

    # calculate terms for correlation
    Cp = C1 + (C2*(((C3/T)/np.sinh(C3/T))**2))+(C4*(((C5/T)/np.cosh(C5/T))**2))

    return Cp

def vapor_viscosity(comp_name, T, data_dict):

    # calculates the ideal gas vapor viscosity

    # input: name of component, temperature in Kelvin, dictionary of physical and chemical properties
    # output: ideal gas vapor viscosity (Pa.s)

    # extracting correlation coefficients from dictionary
    C1 = data_dict['vapor_viscosity'][comp_name]['C1']
    C2 = data_dict['vapor_viscosity'][comp_name]['C2']
    C3 = data_dict['vapor_viscosity'][comp_name]['C3']
    C4 = data_dict['vapor_viscosity'][comp_name]['C4']

    # calculate vapor viscosity correlation
    mu = (C1*(T**C2))/(1+(C3/T)+(C4/(T**2)))

    return mu

def mixture_density(comp_list, mol_fraction, T, P, data_dict):

    # calculates the ideal gas mixture density

    # input: list of mixture components, list of each component mol fraction (same order)
    #        temperature (K), pressure (Pa), dictionary of physical and chemical properties
    # output: the ideal gas mixture (kg/m³)

    # extract molar mass of each component and convert to kg/mol
    MM_list = []
    for comp in comp_list:
        MM_list.append(data_dict['molar_mass'][comp]/1000)

    # calculate the mixture's molar mass
    MM_mix = 0
    for i in range(len(MM_list)):
        MM_mix += mol_fraction[i]*MM_list[i]

    # calculate mixture density using ideal gas law
    rho = (P*MM_mix)/(R*T)

    return rho