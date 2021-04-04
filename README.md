# glfMB - Gas-Liquid Flow model of Mukherjee and Brill

R package for calculation of gas-liquid two-phase flow in a circular pipe with the model of Mukherjee &amp; Brill (1985) 

* Mukherjee, H., and J. P. Brill. 1985. Pressure Drop Correlations for Inclined Two-Phase Flow. Journal of Energy Resources Technology, Transactions of the ASME 107 (4)


## About this 

This R package provides functions to calculate flow regime, liquid holdup, and pressure drop with the empirical model of Mukherjee &amp; Brill (1985). In the model, there are four types of flow regime: stratified, annular, slug, and bubble. 



## Installing 

will be written later.

## Examples

```
library(glfMB)

vsG <- 3.86 * 0.3048    # 3.86 ft/s
vsL <- 3.97 * 0.3048    # 3.97 ft/s
D   <- 6 * 0.0254       # 6 inch
densityG <- 5.88 * 16.01845     # 5.88 lbm/ft3  - (1 lbm/ft3 = 16.01845 kg/m3)
densityL <- 47.61 * 16.01845    # 47.61 lbm/ft3 - (1 lbm/ft3 = 16.01845 kg/m3)
viscosityG <- 0.016 / 1000      # 0.016 cp
viscosityL <- 0.97  / 1000      # 0.970 cp
surfaceTension <- 8.41 / 1000   # 8.41 dynes/cm
angle <- pi/2                   # 90 deg (upward)

# NLv=11.87, NGv=11.54, NGvSM=350.8, NLvBS_up=18.40, NL=0.0118
l_dlns_MB(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
```

