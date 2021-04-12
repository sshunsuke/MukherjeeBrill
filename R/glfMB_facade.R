# - - - - - - - - - - - - - - - - - - - - - - - - -
# Facade ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - -
# Main functions ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

#' Calculate flow regime, holdup, and pressure drop with the model of Mukherjee and Brill (1985)
#'
#' Call core functions of Mukherjee & Brill model: `l_dlns_MB()`, `l_flow_regime_MB()`,
#' `l_holdup_MB()`, and `l_dPdL_MB()` to calculate flow regime, holdup, and pressure drop.
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
#' @param angle Pipe angle (a positive value means upward) - radian
#' @param roughness Pipe roughness - m
#' @param pressure Pressure (optional) - Pa
#'
#' @return Data frame including following data:
#' * `fr`: Flow regime (1: Stratified, 2: Annular, 3: Slug, and 4: Bubbly)
#' * `hl`: Holdup
#' * `dPdL`: Pressure drop per unit length - Pa/m
#' * `dPdL_H`: Pressure drop per unit length due to hydrostatic - Pa/m
#' * `dPdL_F`: Pressure drop per unit length due to friction - Pa/m
#' * `dPdL_A`: Pressure drop per unit length due to acceleration - Pa/m
#'
#' @references
#' * Mukherjee, H., and J. P. Brill. 1985. Pressure Drop Correlations for Inclined Two-Phase Flow. Journal of Energy Resources Technology, Transactions of the ASME 107 (4)
#' * Mukherjee, H., and J. P. Brill. 1985. Empirical Equations to Predict Flow Patterns in Two-Phase Inclined Flow. International Journal of Multiphase Flow 11 (3)
#'
#' @examples
#' # This example is from Brill and Mukherjee (1999)
#' # "Multiphase Flow in Wells"
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
#' @note You can execute the calculation step by step by calling low-level functions if you need. Below is an example.
#' ```
#' dlns <- l_dlns_MB(vsG, vsL, D, densityG, densityL,
#'                   viscosityG, viscosityL, surfaceTension, angle)
#' fr <- l_flow_regime_MB(dlns)
#' hl <- l_holdup_MB(dlns, fr)
#' dPdL <- l_dPdL_MB(dlns, fr, hl, roughness, pressure, debug=FALSE)
#' ```
#'
#' @export
#' @md
call_MB <- function(vsG, vsL, D, densityG, densityL,
                    viscosityG, viscosityL, surfaceTension, angle, roughness, pressure) {
  dlns <- l_dlns_MB(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr   <- l_flow_regime_MB(dlns)
  hl <- l_holdup_MB(dlns, fr)
  dPdL <- l_dPdL_MB(dlns, fr, hl, roughness, pressure, debug=FALSE)
  
  if (any(hl < 0) || any(hl > 1)) {
    stop(sprintf("holdup was not calculated correctly (hl=%.3f). call_MB() may not support viscous liquid.", hl))
  }
  
  data.frame("fr" = fr, "hl" = hl, "dPdL" = as.vector(dPdL[,'dPdL']),
             "dPdL_H" = as.vector(dPdL[,'dPdL_H']),
             "dPdL_F" = as.vector(dPdL[,'dPdL_F']),
             "dPdL_A" = as.vector(dPdL[,'dPdL_A']))
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
#' @return A matrix (class="frm_MB") including the information of flow regime (columns) at
#' specified superfical velocities of gas and liquid.
#' * `vsG`: Superficial velocity of gas
#' * `vsL`: Superficial velocity of liquid
#' * `fr`: Flow regime (1: Stratified, 2: Annular, 3: Slug, and 4: Bubbly)
#' * `NGv`, `NLv`, `NL`, `NGvSM`, `NGvBS`, `NLvBS_up`, `NLvST`: Properties used in the prediction of flow regime
#'
#' @examples
#' vs_range = glfMB::vs_vector_MB(0.1, 10, 20, TRUE)
#' frm <- generate_frm_MB(vs_range, vs_range, 0.1,
#'                        40, 1002, 1.1E-05, 1.6E-03, 0.0695, pi/2)
#' plot(frm)    # glfMB:::plot.frm_MB(frm)
#'
#' @references
#' Mukherjee, H., and J. P. Brill. 1985. Empirical Equations to Predict Flow Patterns in Two-Phase Inclined Flow. International Journal of Multiphase Flow 11 (3)
#'
#' @export
#' @md
generate_frm_MB <- function(vector_vsG, vector_vsL, D, densityG, densityL,
                            viscosityG, viscosityL, surfaceTension, angle) {
  pairs_vs <- expand.grid(vector_vsG, vector_vsL)
  frm <- glfMB:::frm0_MB(pairs_vs[,1], pairs_vs[,2], D, densityG, densityL,
                         viscosityG, viscosityL, surfaceTension, angle)
  frm
}

# internal function
frm0_MB <- function(vsG, vsL, D, densityG, densityL,
                    viscosityG, viscosityL, surfaceTension, angle) {
  dlns <- l_dlns_MB(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr   <- l_flow_regime_MB(dlns)
  frm <- cbind('vsG'=vsG, 'vsL'=vsL, 'fr'=fr,
               'NGv'=dlns$NGv, 'NLv'=dlns$NLv, "NL"=dlns$NL,
               'NGvSM'=dlns$NGvSM, 'NGvBS'=dlns$NGvBS, 'NLvBS_up'=dlns$NLvBS_up, 'NLvST'=dlns$NLvST)
  class(frm) <- "frm_MB"
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
#' @examples
#' \dontrun{
#' frm <- generate_frm_MB(vector_vsG, vector_vsL, D, densityG, densityL,
#'                       viscosityG, viscosityL, surfaceTension, angle)
#' plot(frm)
#' }
#' @importFrom graphics plot points
#'
#' @rdname plot.frm_MB
#' @export
plot.frm_MB <- function(x, xlab, ylab, xval='vsG', yval='vsL', ...) {
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
#' Calculate the Darcy friction factor with assuming a smooth pipe.
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
#' As the correlation cannot be resolved explicitly, Newton-Raphson method is used.
#'
#' @param Re Reynold number
#' @param roughness Pipe roughness
#' @param D Pipe diameter
#' @param tol Tolerance in Newton-Raphson method (optional)
#' @param itMax Maximum number of iteration  (optional)
#' @param warn If FALSE, not show warnings when Re <= 4000  (optional)
#'
#' @return Darcy friction factor
#' @export
Colebrook <- function(Re, roughness, D, tol=1e-8, itMax=10, warn=TRUE) {
  mapply(glfMB:::Colebrook_core, Re, roughness, D, tol, itMax, warn)
}


# - - - - - - - - - - - - - - - - - - - - - - - - -
# Utilities ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

#' Generate a vector of superficial velocities (vs)
#' 
#' @param from minimum value of superficial velocities
#' @param to maximum value of superficial velocities
#' @param num_points number of data points of the returned vector
#' @param log_scale a
#'
#' @export
vs_vector_MB <- function(from, to, num_points, log_scale) {
  log_scale <- ifelse(missing(log_scale), FALSE, log_scale)
  
  stopifnot(from < to)
  
  if (log_scale == TRUE) {
    ret <- seq(log10(from), log10(to), length.out = num_points)
    ret <- 10^(ret)
  } else {
    ret <- seq(from, to, length.out = num_points)
  }
  ret
}


# testdata_MB$ ----

#' List of functions and constants for test (and example)
#' 
#' A list containing functions and constants for fluid properties used in tests (see also Value).
#' SI unit is used.
#' * Mair: Molar mass of air - kg/mol
#' * Sutherland_vis0_air: Viscosity of air at 273 K - Pa-s
#' * Sutherland_T0_air - K
#' * Sutherland_C_air - K
#' * ideal_gas_density(temperature, pressure, M) - kg/m3 (K, Pa)
#' * Sutherland_viscosityG(temperature, viscosity0, t0, const) - Pa-s (K, Pa-s, K, K)
#' * Kerosene_surfaceTension(temperature): K
#' * Kerosene_density(temperature) - kg/m3 (K, Pa)
#' * Kerosene_viscosity(temperature): K
#' * Lubricating_surfaceTension(temperature) - K
#' * Lubricating_density(temperature) - kg/m3 (K, Pa)
#' * Lubricating_viscosity(temperature) - K
#' 
#' @name testdata_MB
#' @docType data
#' 
#' @examples
#' testdata_MB$Kerosene_viscosity(283.15)     # Kerosene Viscosity at 10 degC
#'
#' density_air <- function(temperature, pressure) {
#'   testdata_MB$ideal_gas_density(temperature, pressure, testdata_MB$Mair)
#' } 
#' density_air(293.15, 5*100000)  # Air density at 20 degC and 5 bar
#'
#' @export
#' @md
testdata_MB <- list(
  
  # Molar mass of air (kg/mol)
  Mair = 28.96 / 1000,  
  
  # Constants for Sutherland's formula
  Sutherland_vis0_air = 1.71E-5,    # Pa-s
  Sutherland_T0_air = 273,          # K
  Sutherland_C_air = 110.4,         # K
  
  # Density of ideal gas
  ideal_gas_density = function(temperature, pressure, M) {
    pressure * M / (glfMB:::R * temperature)
  },
    
  # Sutherland's formula
  Sutherland_viscosityG = function(temperature, viscosity0, t0, const) {
    viscosity0 * (temperature / t0)^(3/2) * (t0 + const) / (temperature + const)
  },
    
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
  }
    
)




