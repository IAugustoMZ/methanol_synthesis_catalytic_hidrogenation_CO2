# **Methanol Synthesis by Catalytic Hydrogenation of CO2**
<br>

## **1. Introduction**

This script aims to simulate a catalytic reactor design to perform the synthesis of methanol by hydrogenation of CO2. This simulations is based on the work of Gosh *et al.*, 2021.

## **2. Hypothesis about the system**

- gaseous phase is ideal
- stationary operation
- reactions follow a Langmuir-Hinshelwood mechanism, where the reactant must adsorve over the catalyst surface before reaction
- reference state is 25 °C and 1 bar
- methanol adsorption obver catalyst is negligible
- both single site and dual site approaches are tested
- enthalpy of adsorption is constant over the temperature range
- pressure in the bed length is given by Ergun's equation (only for PBR approach)
- mass transfer resistances are negligible
- no radial dispersion (PBR approach)
- the catalyst suffers no desactivation over time
- both isothermal and adiabatic approaches are tested
- a $In_2O_3$ catalyst is applied
