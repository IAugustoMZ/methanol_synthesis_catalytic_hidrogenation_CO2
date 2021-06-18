# import packages
import pandas as pd
import numpy as np
import chemical_properties as cp

R = 8.314472             # ideal gas universal constant (J/mol.K)
T_std = 298.15           # reference state temperature (K) 
prog_din = {}            # hash table dictionary to store already calculated integrals of heat capacity

def integrate_cp(component, T_list, therm_func, prop_dict):

    """calculates the integration of heat capacity expression for a specific component
    using the Simpson's composite rule

    inputs: symbol of component, list with initital and end limits of temperature range
             desired thermodynamic function ('h' for enthalpy or 's' for entropy)
             dictionary containing the correlations coefficients for heat capacity expression
    """

    # dynamic programming approach to faster computations
    global prog_din

    # key to hash table containing already calculated values
    dic_key = component + str(np.round(T_list[0], 2)) + str(np.round(T_list[1], 2)) + \
        therm_func

    # try to look if the inputed value was already calculated
    try:
        return prog_din[dic_key]
    except:

        # if the value is not found at hash table, calculate it and save it

        # number of points for integration
        N = 1000

        # create mesh for integration
        T = np.linspace(T_list[0], T_list[1], N)

        # calculate integration step
        h = T[1] - T[0]

        # loop through temperature range to integrate (composite Simpson's rule)
        integral = 0
        for i in range(N):

            # compute function values for the desired thermodynamic function
            if therm_func == 'h':
                f = cp.IdealGasCp(component, T[i], prop_dict)
            else:
                f = cp.IdealGasCp(component, T[i], prop_dict)/T[i]
            
            # add the function evaluations to the integral, according to Simpson's rule
            if(i == 0 | i == N-1):
                integral += f
            elif(i%2 == 0):
                integral += 2*f
            else:
                integral += 4*f

        # after the integration, correct the value using integration step
        integral *= (h/3)

        # convert the integral to J/mol (enthalpy) / J/mol.K (entropy)
        integral /= 1000

        # save the computed integral to hash table
        prog_din[dic_key] = integral

        return integral

def calculate_K_and_Hr(comp_list, stoic_list, T, prop_dict):

    """calculates the chemical equilibrium reaction and the heat of reaction at temperature of reaction

    inputs: list of components, list of stoichiometric coefficients, temperature of reaction
             dictionary containing the heat capacity correlation coefficients
     output: chemical equilibrium constant (adm.) and enthalpy of reaction at desired temperature (J/mol)
     """

    # global variables
    global T_std

    # empty lists to store individual components values
    hf0 = []
    sf0 = []
    deltaH = []
    deltaS = []
    hf = []
    sf = []

    # loop through components to obtain thermal properties
    for i in range(len(comp_list)):

        # select standard enthalpy of formation and standard entropy of formation and convert to J/mol or J/mol.K
        hf0.append(prop_dict['H0f'][comp_list[i]] / 1000)
        sf0.append(prop_dict['S0f'][comp_list[i]] / 1000)

        # integrating the heat capacity from standard temperature to desired temperature
        deltaH.append(integrate_cp(comp_list[i], [T_std, T], 'h', prop_dict))
        deltaS.append(integrate_cp(comp_list[i], [T_std, T], 's', prop_dict))

        # calculating the enthalpy and entropy of each component at reaction temperature, multiplied by stoichiometric coefficient
        hf.append((hf0[i] + deltaH[i])*stoic_list[i])
        sf.append((sf0[i] + deltaS[i])*stoic_list[i])

    # calculate reaction properties
    hr = np.sum(hf)                            # enthalpy of reaction
    sr = np.sum(sf)                            # entropy of reaction
    gr = hr - (T*sr)                           # gibbs energy of reaction

    # logarithm of chemical equilibrium constant
    logK = (-gr/(R*T))

    return [np.exp(logK), hr]