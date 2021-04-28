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

## **3. Reactions considered in this system**

<br>

$$CO_{2(g)} + 3 H_{2(g)}\rightleftharpoons CH_3OH_{(g)} + H_2O_{(g)}$$

$$CO_{2(g)} + H_{2(g)} \rightleftharpoons CO_{(g)} + H_2O_{(g)}$$

$$CO_{2(g)} + 4 H_{2(g)}\rightleftharpoons CH_4{(g)} + 2 H_2O_{(g)}$$

## **4. Rate laws for Langmuir-Hinshelwood mechanism**

By following the Langmuir-Hinshelwood mechanism, the rate laws are given by equations below:

$$r_{CH_3OH}=\frac{k_1 \frac{\Bigg(P_{CO_2}P_{H_2}^3 - \frac{P_{CH_3OH}P_{H_2O}}{K_{eq, CH_3OH}}\Bigg)}{P_{H_2}^2}}{\text{Inhibition Term}}$$

<br>

$$r_{RWGS} = \frac{k_2 \frac{\Bigg(P_{CO_2}P_{H_2} - \frac{P_{CO} P_{H_2O}}{K_{eq, RWGS}}\Bigg)}{\sqrt{P_{H_2}}}}{\text{Inhibition Term}}$$

<br>

$$r_{CH4} = k_3\sqrt{P_{CO_2}} \sqrt{P_{H2}}\frac{\Bigg(1-\frac{P_{CH_4}P_{H_2O}^2}{P_{CO_2}P_{H_2}^4 K_{eq, CH_4}}\Bigg)}{\text{Inhibition Term}}$$

<br>

The expression for the inhibition term will depend on the approach, as follows:

<br>

$$\text{Inhibition Term} = (1+K_{CO_2}P_{CO_2}+\sqrt{K_{H2}P_{H2}})^2\text{ for the single-site approach}$$

<br>

$$\text{Inhibition Term} = (1+K_{CO_2}P_{CO_2})(1+\sqrt{K_{H2}P_{H2}})\text{ for the single-site approach}$$

## 7. Physical and chemical properties correlations

<br>

### 7.1 Ideal Gas Heat Capacity (J/kmol.K)

<br>

$$C_{P} = C_1+C_2\Bigg[\frac{\frac{C_3}{T}}{sinh\frac{C_3}{T}}\Bigg]^2+C_4\Bigg[\frac{\frac{C_5}{T}}{cosh \frac{C_5}{T}}\Bigg]^2$$

<br>

### 7.2. Vapor Viscosity (Pa.s)

<br>

$$\mu = \frac{C_1 T^{C_2}}{1+\frac{C_3}{T}+\frac{C_4}{T^2}}$$

<br>

### 7.3. Ideal Gas Mixture Vapor Viscosity (Pa.s)

<br>

$$\mu_{mix} = \sum_{i=1}^N \frac{y_i \mu_i}{\sum_{j=1}^N y_j \phi_{ij}}$$

<br>

$$\phi_{ij} = \frac{\Bigg[1+\sqrt{\frac{\mu_i}{\mu_j}}\Bigg(\frac{M_j}{M_i}\Bigg)^4\Bigg]^2}{\frac{4}{\sqrt{2}}\sqrt{1+\frac{M_i}{M_j}}}$$

<br>

### 7.4. Mixture Density (kg/m³)

<br>

$$\rho_{mix} = \frac{PM_{mix}}{RT}$$

<br>

$$M_{mix} = \sum_{i=1}^N y_i M_i$$

<br>

Where: $M_{mix}$ is the mixture molar mass and $y_i$ is the molar fraction of component $i$
