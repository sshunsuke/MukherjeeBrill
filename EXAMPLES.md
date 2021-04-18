# Examples

## Setup

```r
library(MukherjeeBrill)
```

----

## Basic

### Estimate flow regime, holdup, and pressure drop


```r
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

### Flow regime map

```r
vector_vsG <- 1:10
vector_vsL <- 1:10
angle2 <- - pi/ 6
frm <- generate_frm_MB(vector_vsG, vector_vsL, ID,
                       densityG, densityL,
                       viscosityG, viscosityL,
                       surfaceTension, angle2)

# s (black): Stratified, A (red): Annular, S (green): Slug, B (blue): Bubble
plot(frm, main="downward flow (-30 deg)", las=1)
```

If you use log-scale for axis, `util_MB_vs_vector()` may be helpful.

```r
# vs_vector_MB(from, to, num_points, log_scale)
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

----

## Advanced

### Very simple flow simulation

Air and water two-phase flow in an inclined pipe. 

```r
# Properties of air and water
density_air <- function(temperature, pressure) {
  testdata_MB$ideal_gas_density(temperature, pressure, testdata_MB$Mair)
}
viscosity_air <- function(temperature) {
  testdata_MB$Sutherland_viscosityG(temperature,
                                    testdata_MB$Sutherland_vis0_air,
                                    testdata_MB$Sutherland_T0_air,
                                    testdata_MB$Sutherland_C_air)
}
density_water <- 1000          # kg/m3
viscosity_water <- 0.001       # Pa-s
surfaceTension_water <- 0.072  # N/m
```

```r
# Flow conditions
temperature <- 20 + 273.15     # 20 degC
inlet_pressure <- 20 * 100000  # 20 bar
ID <- 0.1                      # 10 cm
angle_up <- pi/6               # 30 deg (upward)
roughness <- 0.0001            # m  - Pipe roughness
pipe_length <- 1000            # m
num_section <- 20

# Mass flow rates of air and water
massflowrate_air <- 1          # kg/s
massflowrate_water <- 2        # kg/s

# Length of each section and flow area
section_length <- pipe_length / num_section
area <- (ID/2)^2 * pi

# Prepare an empty data frame for result
NAs <- rep(NA, num_section)
df <- data.frame(L=seq(section_length, pipe_length, by=section_length),
                 pressure=NAs, fr=NAs, hl=NAs, dPdL=NAs,
                 dPdL_H=NAs, dPdL_F=NAs, dPdL_A=NAs,
                 vsG=NAs, vsL=NAs)
```

```r
# Flow simulaion
pressure_ <- inlet_pressure

for (i in 1:num_section) {
  densityG_ <- density_air(temperature, pressure_)
  densityL_ <- density_water
  viscosityG_ <- viscosity_air(temperature)
  viscosityL_ <- viscosity_water
  surfaceTension_ <- surfaceTension_water
  vsG_ <- massflowrate_air / densityG_ / area
  vsL_ <- massflowrate_water / densityL_ / area
  
  ret <- call_MB(vsG_, vsL_, ID, densityG_, densityL_,
                 viscosityG_, viscosityL_, surfaceTension_,
                 angle_up, roughness, pressure_)
  pressure_ <- pressure_ - (ret$dPdL * section_length)
  
  df$pressure[i] <- pressure_
  df$fr[i] <- ret$fr
  df$hl[i] <- ret$hl
  df$dPdL[i] <- ret$dPdL
  df$dPdL_H[i] <- ret$dPdL_H
  df$dPdL_F[i] <- ret$dPdL_F
  df$dPdL_A[i] <- ret$dPdL_A
  df$vsG[i] <- vsG_
  df$vsL[i] <- vsL_
}

plot(df$L, df$pressure/100/1000, type="l", xlab='L [m]', ylab='P [bar]', las=1)
# 3: Slug, 2: Annular
plot(df$L, df$fr, type="l", col="red", xlab='L [m]', ylab='Flow regime', las=1)

df
```



