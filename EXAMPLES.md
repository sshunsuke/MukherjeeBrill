# Examples

## Basic

### Estimate flow regime, holdup, and pressure drop

```
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



```
vs_range = glfMB::vs_vector_MB(0.1, 10, 20, TRUE)
frm <- generate_frm_MB(vs_range, vs_range, 0.1,
                       40, 1002, 1.1E-05, 1.6E-03, 0.0695, pi/2)
plot(frm)
```

