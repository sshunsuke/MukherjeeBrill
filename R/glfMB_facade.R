# - - - - - - - - - - - - - - - - - - - - - - - - -
# Facade ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - -
# Main functions ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

#' Calculate flow regime, holdup, and pressure drop with the model of Mukherjee and Brill (1985)
#'
#' Call the core functions of Mukherjee & Brill model: `dlns_MB()`, `flow_regime_MB()`,
#' `holdup_MB()`, and `dPdL_MB()` to calculate flow regime, holdup, and pressure drop.
#' The last argument (`pressure`) is for consideration of accelaration, which you can ignore.
#'
#' @usage call_MB(vsG, vsL, D, densityG, densityL,
#'         viscosityG, viscosityL, surfaceTension,
#'         angle, roughness, pressure)
#'
#' @param vsG Superficial velocity of gas - m/s
#' @param vsL Superficial velocity of liquid - m/s
#' @param D Pipe diameter - m
#' @param densityG Density of gas - kg/m3
#' @param densityL Density of liquid - kg/m3
#' @param viscosityG Visosity of gas - Pa-s
#' @param viscosityL Visosity of liquid - Pa-s
#' @param surfaceTension Surface tension - N/m
#' @param angle Pipe angle (0 is horizontal flow) - radian
#' @param roughness Pipe roughness - m
#' @param pressure Pressure (optional) - Pa
#'
#' @return Data frame including following data:
#' * `fr`: Flow regime (1: Stratified, 2: Annular, 3: Slug, and 4: Bubbly)
#' * `hl`: Holdup
#' * `dPdL`: Pressure drop per unit length - Pa/m
#'
#' @references
#' * Mukherjee, H., and J. P. Brill. 1985. Pressure Drop Correlations for Inclined Two-Phase Flow. Journal of Energy Resources Technology, Transactions of the ASME 107 (4)
#' * Mukherjee, Hemanta, and James P. Brill. 1985. Empirical Equations to Predict Flow Patterns in Two-Phase Inclined Flow. International Journal of Multiphase Flow 11 (3)
#'
#' @examples
#' # This example is from Brill and Mukherjee (1999) "Multiphase Flow in Wells"
#' vsG <- 3.86 * 0.3048    # 3.86 ft/s
#' vsL <- 3.97 * 0.3048    # 3.97 ft/s
#' D   <- 6 * 0.0254       # 6 inch
#' densityG <- 5.88 * 16.01845     # 5.88 lbm/ft3  - (1 lbm/ft3 = 16.01845 kg/m3)
#' densityL <- 47.61 * 16.01845    # 47.61 lbm/ft3 - (1 lbm/ft3 = 16.01845 kg/m3)
#' viscosityG <- 0.016 / 1000      # 0.016 cp
#' viscosityL <- 0.97  / 1000      # 0.970 cp
#' surfaceTension <- 8.41 / 1000   # 8.41 dynes/cm
#' angle <- pi/2                   # 90 deg
#' roughness <- 0.00006 * 0.3048   # 0.00006 ft
#' pressure <- 1700 * 6894.76      # 1700 psia (Example 3.2)
#'
#' # Results should be 3 (slug), 0.560, and about 4727 Pa/m (= 0.209 psi/ft)
#' call_MB(vsG, vsL, D, densityG, densityL,
#'         viscosityG, viscosityL, surfaceTension,
#'         angle, roughness, pressure)
#'
#' @note You can execute the calculation step by step without this function if want. Below is an example.
#' ```
#' dlns <- dlns_MB(vsG, vsL, D, densityG, densityL,
#'                 viscosityG, viscosityL, surfaceTension, angle)
#' fr   <- flow_regime_MB(dlns)
#' hl <- holdup_MB(dlns, fr)
#' dPdL <- dPdL_MB(dlns, fr, hl, roughness, pressure, debug=FALSE)
#' ```
#'
#' @export
#' @md
call_MB <- function(vsG, vsL, D, densityG, densityL,
                    viscosityG, viscosityL, surfaceTension, angle, roughness, pressure) {
  dlns <- dlns_MB(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr   <- flow_regime_MB(dlns)
  hl <- holdup_MB(dlns, fr)
  dPdL <- dPdL_MB(dlns, fr, hl, roughness, pressure, debug=FALSE)
  
  cbind("fr" = fr, "hl" = hl, "dPdL" = as.vector(dPdL))
}



#' Generate flow regime map (frm) data
#'
#' Generate data (matrix) for flow regime map. Range of the map is specified by `vector_vsG` and `vector_vsL`.
#'
#' @param vector_vsG Vector of Superficial velocity of gas - m/s
#' @param vector_vsL Vector of superficial velocity of liquid - m/s
#' @param D Pipe diameter - m
#' @param densityG Density of gas - kg/m3
#' @param densityL Density of liquid - kg/m3
#' @param viscosityG Visosity of gas - Pa-s
#' @param viscosityL Visosity of liquid - Pa-s
#' @param surfaceTension Surface tension - N/m
#' @param angle Pipe angle (0 is horizontal flow) - radian
#'
#' @return A matrix (class="frm_MB") including the information of flow regime at
#' specified superfical velocities of gas and liquid.
#' * `vsG`: Superficial velocity of gas
#' * `vsL`: Superficial velocity of liquid
#' * `fr`: Flow regime (1: Stratified, 2: Annular, 3: Slug, and 4: Bubbly)
#' * `NGv`, `NLv`, `NL`, `NGvSM`, `NGvBS`, `NLvBS_up`, `NLvST`: Properties used in model
#'
#' @examples
#' vs_range = glfMB:::frm_vs_range(0.1, 10, 20, TRUE)
#' frm <- generate_frm_MB(vs_range, vs_range, 0.1, 40, 1002, 1.1E-05, 1.6E-03, 0.0695, pi/2)
#' glfMB:::plot.frm_MB(frm)
#'
#' @export
#' @md
generate_frm_MB <- function(vector_vsG, vector_vsL, D, densityG, densityL,
                            viscosityG, viscosityL, surfaceTension, angle) {
  pairs_vs <- expand.grid(vector_vsG, vector_vsL)
  frm <- glfMB:::frm0_MB(pairs_vs[,1], pairs_vs[,2], D, densityG, densityL,
                         viscosityG, viscosityL, surfaceTension, angle)
  
  #if (plotting == TRUE) {
  #  glfMB:::plot_flow_regime_map(frm)
  #}
  
  frm
}


#' Plot a flow regime map
#'
#' @param x Flow regime map data created by `generate_frm_MB()`
#' @param xlab label x (optional)
#' @param ylab label y (optional)
#' @param xval 'vsG' or 'NGv' (optional)
#' @param yval 'vsL' or 'NLv' (optional)
#' @param ... graphical parameters to plot
#'
#' @examples # plot(a)
#'
#' @importFrom graphics plot points
#'
#' @rdname plot.frm_MB
#' @export
plot.frm_MB <- function(x, xlab, ylab, xval='vsG', yval='vsL', ...) {
  #axis_log <- ifelse(missing(axis_log), FALSE, axis_log)
  
  xlab <- ifelse(missing(xlab), xval, xlab)
  ylab <- ifelse(missing(ylab), yval, ylab)
  
  plot(c(min(x[,xval]), max(x[,xval])), c(min(x[,yval]), max(x[,yval])),
       type="n", xlab=xlab, ylab=ylab, ...)
  
  st <- which(x[,'fr'] == 1)
  an <- which(x[,'fr'] == 2)
  sl <- which(x[,'fr'] == 3)
  bl <- which(x[,'fr'] == 4)
  
  points(x[st, xval], x[st, yval], pch='s', col='black')
  points(x[an, xval], x[an, yval], pch='A', col='red')
  points(x[sl, xval], x[sl, yval], pch='S', col='green')
  points(x[bl, xval], x[bl, yval], pch='B', col='blue')
}


# - - - - - - - - - - - - - - - - - - - - - - - - -
# Darcy friction factor ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

#' Blasius correlation for the Darcy friction factor
#'
#' Calculate the Darcy friction factor with assuming a smooth pipe
#'
#' @param Re Reynold number
#' @return Darcy friction factor
#' @export
Blasius <- function(Re) {
  0.3164 / (Re ^ 0.25)
}


#' Colebrook correlation for the Darcy friction factor
#'
#' Calculate the Darcy friction factor with the Colebrook correlation.
#' As the correlation cannot be resolved explicitly, Newton-Raphson is used.
#'
#' @param Re Reynold number
#' @param roughness Pipe roughness
#' @param D Pipe diameter
#' @param tol Tolerance in Newton-Raphson method (optional)
#' @param itMax Maximum number of iteration  (optional)
#' @param warn show warnings when Re <= 4000  (optional)
#'
#' @return Darcy friction factor
#' @export
Colebrook <- function(Re, roughness, D, tol=1e-8, itMax=10, warn=TRUE) {
  mapply(glfMB:::Colebrook_core, Re, roughness, D, tol, itMax, warn)
}


# - - - - - - - - - - - - - - - - - - - - - - - - -
# Utilities ----
# - - - - - - - - - - - - - - - - - - - - - - - - -


frm_vs_range <- function(min_, max_, np, log_scale) {
  log_scale <- ifelse(missing(log_scale), FALSE, log_scale)
  
  if (log_scale == TRUE) {
    ret <- seq(log10(min_), log10(max_), length.out = np)
    ret <- 10^(ret)
  } else {
    ret <- seq(min_, max_, length.out = np)
  }
  ret
}




#' Get a list of functions for test (and example)
#' 
#' This function returns a list of functions for fluid properties (see also Value).
#' 
#' @return a list of functions for fluid properties. SI unit is used in the functions. 
#' * Kerosene_surfaceTension(temperature) - K
#' * Kerosene_density(temperature) - K
#' * Kerosene_viscosity(temperature) - K
#' * Lubricating_surfaceTension(temperature) - K
#' * Lubricating_density(temperature) - K
#' * Lubricating_viscosity(temperature) - K
#' * Air_density(temperature, pressure) - K, Pa
#' * Air_viscosity(temperature, pressure) - K, Pa
#' 
#' @examples td <- testdata_MB()
#' td$Kerosene_viscosity(283.15)     # Kerosene Viscosity at 10 degC
#' td$Air_density(293.15, 5*100000)  # Air density at 20 degC and 5 bar
#'
#' @export
#' @md
testdata_MB <- function() {
  
  list(
    # Kerosene (Mukherjee & Brill, 1985)
    Kerosene_surfaceTension = function(temperature) {
      (27.6 - 0.09 * (temperature-273.15)) / 1000           # N/m
    },
    Kerosene_density = function(temperature) {
      832.34 - 0.8333 * (temperature-273.15)                # kg/m3
    },
    Kerosene_viscosity = function(temperature) {
      1e-03 * exp(1.0664 - 0.0207 * (temperature-273.15))   # Pa
    },
    
    # Lubricating Oil (Mukherjee & Brill, 1985)
    Lubricating_surfaceTension = function(temperature) {
      (36.6094 - 0.117 * (temperature-273.15)) / 1000       # N/m
    },
    Lubricating_density = function(temperature) {
      861.22 - 0.7585 * (temperature-273.15)                # kg/m3
    },
    Lubricating_viscosity = function(temperature) {
      1e-03 * exp(3.99 - 0.0412 * (temperature-273.15))     # Pa
    },
    
    # Air (regarded as ideal gas)
    Air_density = function(temperature, pressure) {
      M <- 28.96 / 1000  # Molar mass of air (kg/mol)
      R <- 8.314         # Gas constant
      pressure * M / (R * temperature)
    },
    Air_viscosity = function(temperature, pressure) { 1.81 * 10^(-5)}
  )
}



