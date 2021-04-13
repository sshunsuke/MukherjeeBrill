# glfMB - Gas-Liquid Flow model of Mukherjee and Brill

R package for calculation of gas-liquid two-phase flow in a circular pipe with the model of Mukherjee &amp; Brill (1985) 

* Mukherjee, H., and J. P. Brill. 1985. Pressure Drop Correlations for Inclined Two-Phase Flow. Journal of Energy Resources Technology, Transactions of the ASME 107 (4)


## About this 

This R package provides functions to calculate flow regime, liquid holdup, and pressure drop with the empirical model of Mukherjee &amp; Brill (1985). In the model, there are four types of flow regime: stratified, annular, slug, and bubble. 



## Installing 

```
# install.packages("remotes")
remotes::install_github("sshunsuke/glfMB")
```

## Examples

Prediction of flow regime, holdup, and pressure drop. 

```
library(glfMB)

# Flow conditions
vsG <- 5               # m/s   - superficial velocity of gas
vsL <- 1               # m/s   - superficial velocity of liquid
ID <- 0.1              # m     - pipe inner diameter
densityG <- 5          # kg/m3 - density of gas
densityL <- 1000       # kg/m3 - density of liquid 
viscosityG <- 10^(-5)  # Pa-s  - viscosity of gas
viscosityL <- 10^(-3)  # Pa-s  - viscosity of liquid
surfaceTension <- 8.41 / 1000   # 8.41 dynes/cm
angle <- pi/2          # 90 deg (upward)
roughness <- 0.0001      # m

# Predict flow regime, holdup, and pressure drop. 
call_MB(vsG, vsL, ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle, roughness)
#   fr        hl     dPdL   dPdL_H   dPdL_F dPdL_A
# 1  3 0.2911177 3955.877 2887.688 1068.189      0
```

Create a flow regime map. 

```
vs_range = glfMB::vs_vector_MB(0.1, 10, 20, TRUE)
frm <- generate_frm_MB(vs_range, vs_range, 0.1,
                       40, 1002, 1.1E-05, 1.6E-03, 0.0695, pi/2)
plot(frm)
```
