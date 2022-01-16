# MukherjeeBrill - Gas-Liquid Flow model of Mukherjee and Brill

R package for calculation of gas-liquid two-phase flow in a circular pipe with the model of Mukherjee &amp; Brill (1985) 


## About this 

This R package provides functions to calculate flow regime, liquid holdup, and pressure drop with the empirical model of Mukherjee &amp; Brill (1985). In the model, there are four types of flow regime: stratified, annular, slug, and bubble. 

* Mukherjee, H., and J. P. Brill. 1985. Pressure Drop Correlations for Inclined Two-Phase Flow. Journal of Energy Resources Technology, Transactions of the ASME 107 (4)


## Installing 

```r
# install.packages("remotes")
remotes::install_github("sshunsuke/MukherjeeBrill")
```

## Examples

Prediction of flow regime, holdup, and pressure drop. 

```r
library(MukherjeeBrill)

# Flow conditions (SI unit is used for all properties)
vsG <- 5                 # Superficial velocity of gas [m/s]
vsL <- 1                 # Superficial velocity of liquid [m/s]
ID <- 0.1                # Pipe inner diameter [m]
densityG <- 1            # Density of gas [kg/m3]
densityL <- 1000         # Density of liquid [kg/m3]
viscosityG <- 10^(-5)    # Viscosity of gas [Pa-s]
viscosityL <- 10^(-3)    # Viscosity of liquid [Pa-s]
surfaceTension <- 0.072  # Surface tension [N/m]
angle <- pi/2            # Pipe inclination [rad] (positive means upward)
roughness <- 0.0001      # Pipe roughness [m] 

# Predict flow regime, holdup, and pressure drop. 
# call_MB() returns a data frame including six columns:
#   - fr: Flow regime (1: Stratified, 2: Annular, 3: Slug, and 4: Bubbly)
#   - hl: Holdup
#   - dPdL: Pressure drop per unit length [Pa/m]
#   - dPdL_H: Pressure drop per unit length due to hydrostatic [Pa/m]
#   - dPdL_F: Pressure drop per unit length due to friction [Pa/m]
#   - dPdL_A: Pressure drop per unit length due to acceleration [Pa/m]
call_MB(vsG, vsL, ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle, roughness)

#   fr        hl     dPdL   dPdL_H   dPdL_F dPdL_A
# 1  3 0.3272534 4403.032 3213.676 1189.355      0
```

Create a flow regime map. 

```r
vector_vsG <- 1:10
vector_vsL <- 1:10
angle2 <- - pi/ 6

# generate_frm_MB() returns a S3 class object, the class name of which is 'frm_MB'. 
frm <- generate_frm_MB(vector_vsG, vector_vsL, ID, densityG, densityL,
                       viscosityG, viscosityL, surfaceTension, angle2)

# Plot a flow regime map
#   s (black): Stratified, A (red): Annular, S (green): Slug, B (blue): Bubble
plot(frm)
```

If you use log-scale for axis, `util_MB_vs_vector()` may be helpful.

```r
# util_MB_vs_vector(from, to, num_points, log_scale)
vector_vsG2 <- util_MB_vs_vector(0.1, 100, 10, TRUE)
vector_vsL2 <- util_MB_vs_vector(0.1, 100, 10, TRUE)
frm_wider <- generate_frm_MB(vector_vsG2, vector_vsL2, ID,
                             densityG, densityL,
                             viscosityG, viscosityL,
                             surfaceTension, angle2)

# s (black): Stratified, A (red): Annular, S (green): Slug, B (blue): Bubble
plot(frm_wider, log="xy", main="downward flow (-30 deg)", las=1)
vector_vsG2
```

Show help. 

```r
help("MukherjeeBrill")
```

You can find more examples [here](EXAMPLES.md). 


## License 

[MIT License](LICENSE.md)


## Related documents

* Mukherjee, H. and Brill J. P. 1985. Empirical Equations to Predict Flow Patterns in Two-Phase Inclined Flow. International Journal of Multiphase Flow 11 (3)
* Mukherjee, H. and Brill J. P. 1983. Liquid Holdup Correlations for Inclined Two-Phase Flow. JPT, Journal of Petroleum Technology 35(5):1003–8.
* Brill J. P. and Mukherjee, H. 1999. Multiphase Flow in Wells. Society of Petroleum Engineers.

