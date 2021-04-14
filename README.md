# glfMB - Gas-Liquid Flow model of Mukherjee and Brill

R package for calculation of gas-liquid two-phase flow in a circular pipe with the model of Mukherjee &amp; Brill (1985) 


## About this 

This R package provides functions to calculate flow regime, liquid holdup, and pressure drop with the empirical model of Mukherjee &amp; Brill (1985). In the model, there are four types of flow regime: stratified, annular, slug, and bubble. 

* Mukherjee, H., and J. P. Brill. 1985. Pressure Drop Correlations for Inclined Two-Phase Flow. Journal of Energy Resources Technology, Transactions of the ASME 107 (4)


## Installing 

```
# install.packages("remotes")
remotes::install_github("sshunsuke/glfMB")
```

## Examples

Prediction of flow regime, holdup, and pressure drop. 

```
library(glfMB)

# Flow conditions (SI unit is used for all properties)
vsG <- 5                 # m/s   - superficial velocity of gas
vsL <- 1                 # m/s   - superficial velocity of liquid
ID <- 0.1                # m     - pipe inner diameter
densityG <- 1            # kg/m3 - density of gas
densityL <- 1000         # kg/m3 - density of liquid 
viscosityG <- 10^(-5)    # Pa-s  - viscosity of gas
viscosityL <- 10^(-3)    # Pa-s  - viscosity of liquid
surfaceTension <- 0.072  # N/m
angle <- pi/2            # rad (positive means upward) - Pipe inclination 
roughness <- 0.0001      # m  - Pipe roughness

# Predict flow regime, holdup, and pressure drop. 
call_MB(vsG, vsL, ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle, roughness)
#   fr        hl     dPdL   dPdL_H   dPdL_F dPdL_A
# 1  3 0.3272534 4403.032 3213.676 1189.355      0
```

Create a flow regime map. 

```
vector_vsG <- 1:10
vector_vsL <- 1:10
frm <- generate_frm_MB(vector_vsG, vector_vsL, ID, densityG, densityL,
                       viscosityG, viscosityL, surfaceTension, angle)
# Plot a flow regime map
#   s (black): Stratified, A (red): Annular, S (green): Slug, B (blue): Bubble
plot(frm)
```

You can fine more examples here. 


# License 



## Related documents


