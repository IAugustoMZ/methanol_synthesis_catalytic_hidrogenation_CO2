# this script is applied to validate the implementation of 
# chemical properties correlation

# validation data collected at Perry's Chemical Engineering Handbook
validation_data = {
    "IdealGasCp":{
        "T_lim":{
            "CO2": [50, 5000],
            "H2": [250, 1500],
            "MeOH": [200, 1500],
            "H2O": [100, 2273.15],
            "CO": [60, 1500],
            "CH4": [50, 1500]
        },
        "values":{
            "CO2": [0.2937*(10**5), 0.6335*(10**5)],
            "H2": [0.2843*(10**5), 0.3225*(10**5)],
            "MeOH": [0.3980*(10**5), 1.0533*(10**5)],
            "H2O": [0.3336*(10**5), 0.5276*(10**5)],
            "CO": [0.2911*(10**5), 0.3521*(10**5)],
            "CH4": [0.3330*(10**5), 0.8890*(10**5)]
        }
    },
    "vapor_viscosity":{
        "T_lim":{
            "CO2": [194.67, 1500],
            "H2": [13.95, 3000],
            "MeOH": [240, 1000],
            "H2O": [273.16, 1073.15],
            "CO": [68.15, 1250],
            "CH4": [90.69, 1000]
        },
        "values":{
            "CO2": [9.749*(10**-6), 5.203*(10**-5)],
            "H2": [6.517*(10**-7), 4.330*(10**-5)],
            "MeOH": [7.523*(10**-6), 3.128*(10**-5)],
            "H2O": [8.882*(10**-6), 4.082*(10**-5)],
            "CO": [4.434*(10**-6), 4.654*(10**-5)],
            "CH4": [3.468*(10**-6), 2.800*(10**-5)]
        }
    }
}

# importing libraries
import numpy as np
import pandas as pd
import json
import chemical_properties as cp

# import physical and chemical data from json
f = open('data_dict.json')
data_dict = json.load(f)

# validation of ideal gas heat capacity
for comp in validation_data['IdealGasCp']['T_lim']:
    for i in range(len(validation_data['IdealGasCp']['T_lim'][comp])):
        ref_value = validation_data['IdealGasCp']['values'][comp][i]
        T = validation_data['IdealGasCp']['T_lim'][comp][i]
        print('Erro Absoluto de Cp para %s em %.2f K: %.2f%%'%(comp, T, (cp.IdealGasCp(comp, T, data_dict) - ref_value)*(100/ref_value)))
    print('\n')

# validation of ideal gas vapor viscosity
for comp in validation_data['vapor_viscosity']['T_lim']:
    for i in range(len(validation_data['vapor_viscosity']['T_lim'][comp])):
        ref_value = validation_data['vapor_viscosity']['values'][comp][i]
        T = validation_data['vapor_viscosity']['T_lim'][comp][i]
        print('Erro Absoluto de mu para %s em %.2f K: %.2f%%'%(comp, T, (cp.vapor_viscosity(comp, T, data_dict) - ref_value)*(100/ref_value)))
    print('\n')