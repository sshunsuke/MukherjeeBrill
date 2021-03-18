#' glfMB - Gas-Liquid Flow model of Mukherjee and Brill
#'
#' A package for calculation of gas-liquid two-phase flow in a circular pipe
#' with the model of Mukherjee and Brill (1985), which predicts flow regime, liquid holdup, and pressure drop.
#'
#' @references
#' * Mukherjee, H., and J. P. Brill. 1985. Pressure Drop Correlations for Inclined Two-Phase Flow. Journal of Energy Resources Technology, Transactions of the ASME 107 (4)
#' * Mukherjee, Hemanta, and James P. Brill. 1985. Empirical Equations to Predict Flow Patterns in Two-Phase Inclined Flow. International Journal of Multiphase Flow 11 (3)
#'
#' @name glfMB
#' @docType package
#' @import stats
#' @md
NULL

# Note:
# The numbers in comments, such as "(4.130)", indicate the equation numbers in Brill and Mukherjee ().
#
#



# - - - - - - - - - - - - - - - - - - - - - - - - -
# Constants ----
# * pi ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

g <- 9.8    # Gravitational acceleration (m/s2)



# - - - - - - - - - - - - - - - - - - - - - - - - -
# Facade ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

#' Call core functions of Mukherjee & Brill model
#'
#' Calculate flow regime, holdup, and pressure drop with the model of Mukherjee and Brill (1985).
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
#' @param roughness Pipe roughness
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
#' # Results should be 3 (slug), 0.560, and 4727 Pa/m (= 0.209 psi/ft)
#' call_MB(vsG, vsL, D, densityG, densityL,
#'         viscosityG, viscosityL, surfaceTension,
#'         angle, roughness, pressure)
#'
#' @note You can execute the calculation step by step without this function if want. Below is an example.
#' ```
#'   dlns <- dlns_MB(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
#'   fr   <- flow_regime_MB(dlns)
#'   hl <- holdup_MB(dlns, fr)
#'   dPdL <- dPdL_MB(dlns, fr, hl, roughness, pressure, debug=FALSE)
#' ```
#'
#' @export
#' @md
call_MB <- function(vsG, vsL, D, densityG, densityL,
                    viscosityG, viscosityL, surfaceTension, angle, roughness, pressure) {
  dlns <- dlns_MB(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr   <- flow_regime_MB(dlns)
  hl <- holdup_MB(dlns, fr)
  dPdL <- dPdL_MB(dlns, fr, hl, roughness, pressure)

  cbind("fr" = fr, "hl" = hl, "dPdL" = dPdL)
}


# developing function

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
#' @example
#' vs_range = glfMB:::frm_vs_range(0.1, 10, 20, TRUE)
#' frm <- generate_frm_MB(vs_range, vs_range, 0.1, 40, 1002, 1.1E-05, 1.6E-03, 0.0695, pi/2)
#' plot.frm_MB(frm)
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


# developing function
frm0_MB <- function(vsG, vsL, D, densityG, densityL,
                    viscosityG, viscosityL, surfaceTension, angle) {
  dlns <- dlns_MB(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr   <- flow_regime_MB(dlns)
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
# Utilities ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

circle_area = function(r) { r * r * pi }

Reynolds <- function(density, v, D, viscosity) {
  density * v * D / viscosity
}

# Converter from radian to degree
rad2deg = function(rad) { rad * 180 / pi }

# Convertor from degree to radian
deg2rad = function(deg) { deg * pi / 180 }




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
#' Calculate the Darcy friction factor with assuming a smooth pipe
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
  core_ <- function(Re, roughness, D) {
    if (Re <= 4000 && (warn == TRUE)) { warning("Re <= 4000 !") }

    fun <- function(fD) {
      (1 / sqrt(fD)) + 2 * log10( roughness / D / 3.71 + 2.51 / Re / sqrt(fD))
    }

    # Derivative of fun().
    dFun <- function(fD) {
      - fD^(-3/2) * (1/2 + 2.51 / log(10) / (2.51 / Re / sqrt(fD) + roughness / 3.71 / D) / Re)
    }

    fD_n <- glfMB:::Blasius(Re)
    d <- fun(fD_n) / dFun(fD_n)

    # Newton-Raphson method
    it <- 0
    while (abs(d) >= tol) {
      it <- it + 1
      if (it > itMax) {
        stop("Calculation did not converge.")
      }

      d <- fun(fD_n) / dFun(fD_n)
      fD_n <- fD_n - d
    }

    fD_n
  }

  mapply(core_, Re, roughness, D)
}


laminar <- function(Re) {
  64 / Re
}

transition <- function(Re, roughness, D, tol=1e-8, itMax=10, warn=TRUE) {
  fl <- glfMB:::laminar(Re)
  ft <- glfMB:::Colebrook(Re, roughness, D)

  (fl * (4000 - Re) + ft * (Re - 2000)) / 2000
}


#' Calculate the Darcy friction factor
#'
#' @param Re Reynold number
#' @param roughness Pipe roughness
#' @param D Pipe diameter
#' @param tol Tolerance in `Colebrook()` (optional)
#' @param itMax Maximum number of iteration in `Colebrook()` (optional)
#'
#' @return Darcy friction factor
#' @export
Darcy_friction_factor <- function(Re, roughness, D, tol=1e-8, itMax=10) {
  if (Re >= 4000) {
    ret <- glfMB:::Colebrook(Re, roughness, D, warn=FALSE)
  } else if (Re <= 2000) {
    ret <- glfMB:::laminar(Re, roughness, D)
  } else {
    ret <- glfMB:::transition(Re, roughness, D)
  }
  ret
}






# - - - - - - - - - - - - - - - - - - - - - - - - -
# Dimensionless groups proposed by Duns & Ros ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

#' Calculate dimensionless numbers
#'
#' Return the dimensionless groups proposed by Duns & Ros:
#' * `NLv`: Liquid velocity number
#' * `NGv`: Gas velocity number
#' * `Nd`: Pipe diameter number
#' * `NL`: Liquid viscosity number
#' * `NGvSM`: Slug/(Annular Mist) transition
#' * `NGvBS`: Bubble/Slug transition (Downflow)
#' * `NLvBS_up`: Bubble/Slug transition (Upflow)
#' * `NLvST`: Stratified (Downflow)
#'
#' @param vsG Superficial velocity of gas
#' @param vsL Superficial velocity of liquid
#' @param D Pipe diameter
#' @param densityG Density of gas
#' @param densityL Density of liquid
#' @param viscosityG Visosity of gas
#' @param viscosityL Visosity of liquid
#' @param surfaceTension Surface tension
#' @param angle Pipe angle in radian (0 is horizontal flow)
#'
#' @return Data frame of dimensionless numbers (`NLv`, `NGv`, `Nd`, `NL`, `NGvSM`, `NLvBS`, `NGvBS_up`, and `NLvST`) and input values
#'
#' @usage dlns_MB(vsG, vsL, D, densityG, densityL,
#'         viscosityG, viscosityL, surfaceTension, angle)
#'
#' @examples
#' \dontrun{
#' vsG <- 3.86 * 0.3048    # 3.86 ft/s
#' vsL <- 3.97 * 0.3048    # 3.97 ft/s
#' D   <- 6 * 0.0254       # 6 inch
#' densityG <- 5.88 * 16.01845     # 5.88 lbm/ft3  - (1 lbm/ft3 = 16.01845 kg/m3)
#' densityL <- 47.61 * 16.01845    # 47.61 lbm/ft3 - (1 lbm/ft3 = 16.01845 kg/m3)
#' viscosityG <- 0.016 / 1000      # 0.016 cp
#' viscosityL <- 0.97  / 1000      # 0.970 cp
#' surfaceTension <- 8.41 / 1000   # 8.41 dynes/cm
#' angle <- pi/2                   # 90 deg
#'
#' #  NLv=11.84, NGv=11.54, NGvSM=350.8, NLvBS_up=18.40, NL=0.0118
#' dlns_MB(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
#' }
#' @export
#' @md

dlns_MB <- function(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle) {
  g <- glfMB:::g

  NLv <- vsL * (densityL / g / surfaceTension)^(0.25)
  NGv <- vsG * (densityL / g / surfaceTension)^(0.25)
  Nd <- D * sqrt(densityL * g / surfaceTension)
  NL <- viscosityL * (g / densityL / surfaceTension^3)^(0.25)

  # Gas velocity number for Slug/(Annular Mist) transition (4.130)
  NGvSM <- 10 ^ (1.401 - 2.694 * NL + 0.521 * NLv ^ 0.329)

  # Upflow: Bubble/Slug transition (4.128)
  x <- log10(NGv) + 0.940 + (0.074 * sin(angle)) - (0.855 * sin(angle)^2) + (3.695 * NL)
  NLvBS_up <- 10^x

  # Downflow: Bubble/Slug transition (4.131)
  y <- 0.431 - (3.003 * NL) - (1.138 * log10(NLv) * sin(angle)) - (0.429 * log10(NLv)^2 * sin(angle)) + (1.132 * sin(angle))
  NGvBS <- 10^y

  # Downflow: Stratified (4.133)
  z <- 0.321 - (0.017 * NGv) - (4.267 * sin(angle)) - (2.972 * NL) - (0.033 * log10(NGv)^2) - (3.925 * sin(angle)^2)
  NLvST <- 10^z

  data.frame(
    NLv = NLv,      # Liquid velocity number  (4.3)
    NGv = NGv,      # Gas velocity number     (4.4)
    Nd = Nd,        # Pipe diameter number    (4.5)
    NL = NL,        # Liquid viscosity number (4.6)

    NGvSM = NGvSM,            # Slug/(Annular Mist) transition (4.130)
    NGvBS = NGvBS,            # Downflow: Bubble/Slug transition (4.131)
    NLvBS_up = NLvBS_up,      # Upflow: Bubble/Slug transition (4.128)
    NLvST = NLvST,            # Downflow: Stratified (4.133)

    # Input values (these are used in the calculations of holdup and dPdL).
    vsG = vsG,
    vsL = vsL,
    D = D,
    densityG = densityG,
    densityL = densityL,
    viscosityG = viscosityG,
    viscosityL = viscosityL,
    surfaceTension = surfaceTension,
    angle = angle   # [rad]
  )
}


# - - - - - - - - - - - - - - - - - - - - - - - - -
# Flow regime ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

#' Determine flow regime
#'
#' There are four types of defined flow regime:
#' * 1: Stratified
#' * 2: Annular
#' * 3: Slug
#' * 4: Bubbly
#'
#' @param DLNs Dimensionless numbers calculated by `dlns_MB()`
#' @return a vector of numbers indicating flow regime (1: Stratified, 2: Annular, 3: Slug, and 4: Bubbly)
#' @export
#' @md
flow_regime_MB <- function(DLNs) {
  mapply(glfMB:::flow_regime_MB_core,
         DLNs$NGv, DLNs$NLv, DLNs$angle, DLNs$NGvSM, DLNs$NGvBS, DLNs$NLvBS_up, DLNs$NLvST)
}

# Core logic to determine flow regime
#
# Internal function, called from `flow_regime_MB()`.
#
# @param NGv Gas velocity number
# @param NLv Liquid velocity number
# @param angle Pipe angle
# @param NGvSM Slug/(Annular Mist) transition
# @param NGvBS Bubble/Slug transition (Downflow)
# @param NLvBS_up Bubble/Slug transition (Upflow)
# @param NLvST Stratified (Downflow)
#
# @return a number indicating flow regime
flow_regime_MB_core <- function(NGv, NLv, angle, NGvSM, NGvBS, NLvBS_up, NLvST) {
  flowRegime <- 2  # annular

  if (NGv > NGvSM) {
    return (flowRegime)  # Annular
  }

  if (angle > 0) {
    # Upfill
    if (NLv > NLvBS_up) {
      flowRegime <- 4  # bubbly
    } else {
      flowRegime <- 3  # slug
    }
  } else if (abs(angle) > glfMB:::deg2rad(30)) {
    # Downhill
    if (NGv > NGvBS) {
      if (NLv > NLvST) {
        flowRegime <- 3  # Slug
      } else {
        flowRegime <- 1  # Stratified
      }
    } else {
      flowRegime <- 4    # bubbly
    }
  } else {
    # DownStratified
    if (NLv > NLvST) {
      if (NGv > NGvBS) {
        flowRegime <- 3    # Slug
      } else {
        flowRegime <- 4    # bubbly
      }
    } else {
      flowRegime <- 1  # Stratified
    }
  }

  return (flowRegime)
}


# - - - - - - - - - - - - - - - - - - - - - - - - -
# Liquid holdup    ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

#' Calculate holdup
#'
#' @param DLNs dimensionless numbers calculated by `dlns_MB()`
#' @param flowRegime flow regime estimated by `flow_regime_MB()`
#' @return holdup
#' @export
#' @md
holdup_MB <- function(DLNs, flowRegime) {
  mapply(glfMB:::holdup_MB_core,
         DLNs$angle, DLNs$NL, DLNs$NGv, DLNs$NLv, flowRegime)
}

# Core logic to calculate holdup
#
# Internal function, called from `holdup_MB()`
#
# @param angle Pipe angle (radian)
# @param NL a
# @param NGv b
# @param NLv c
# @param flowRegime d
# @return holdup
holdup_MB_core <- function(angle, NL, NGv, NLv, flowRegime) {
  # coefficients
  co <- cbind(
    c(-0.380113, 0.129875, -0.119788,  2.343227, 0.475686, 0.288657),   # Up
    c(-1.330282, 4.808139,  4.171584, 56.262268, 0.079951, 0.504887),   # DownStratified
    c(-0.516644, 0.789805,  0.551627, 15.519214, 0.371771, 0.393952)    # DownOther
  )

  # 1: Up, 2: DownStratified, 3: DownOther
  j <- ifelse((angle > 0), 1, ifelse((flowRegime == 1), 2, 3))

  t1 <- co[1,j] + (co[2,j] * sin(angle)) + (co[3,j] * sin(angle)^2) + (co[4,j] * NL^2)
  t2 <- NGv^co[5,j] / NLv^co[6,j]
  exp(t1 * t2)
}



# - - - - - - - - - - - - - - - - - - - - - - - - -
# dPdL    ----
# - - - - - - - - - - - - - - - - - - - - - - - - -


#' Calculate pressure drop (dP/L)
#'
#' @param DLNs dimensionless numbers calculated by `dlns_MB()`
#' @param flowRegime flow regime estimated by `flow_regime_MB()`
#' @param HL Holdup
#' @param roughness Pipe roughness
#' @param pressure Pressure (optional)
#' @param debug If TRUE, print parameters (usually set to FALSE)
#'
#' @return pressure drop (Pa/m)
#'
#' @export
#' @md
dPdL_MB <- function(DLNs, flowRegime, HL, roughness, pressure, debug=FALSE) {
  if (missing(pressure) == TRUE) {
    pressure <- NA
  }

  mapply(glfMB:::dPdL_core_MB,
         DLNs$D, DLNs$vsG, DLNs$vsL, DLNs$densityG, DLNs$densityL,
         DLNs$viscosityG, DLNs$viscosityL, DLNs$angle,
         flowRegime, HL, roughness, pressure, debug)
}


# Usage: MukherjeeBrill$ffRatio(HR)
ffRatio <- (function(){
  # Table 4.5
  correspondence <- cbind(
    c(1.00, 0.98, 1.20, 1.25, 1.30, 1.25, 1.00,  1.00),    # fR
    c(0.01, 0.20, 0.30, 0.40, 0.50, 0.70, 1.00, 10.00)     # HR
  )
  colnames(correspondence) <- c("fR", "HR")
  approxfun(correspondence[,"HR"], correspondence[,"fR"])
})()


dPdL_core_MB <- function(D, vsG, vsL, densityG, densityL, viscosityG, viscosityL, angle,
                         flowRegime, HL, roughness, pressure, debug) {
  g <- glfMB:::g

  if (flowRegime == 1) {
    # Stratified

    # delta: angle related to liquid level (shown in Fig 4.20)
    fun <- function(delta) { 1/(2*pi) * (delta - sin(delta)) - HL }  # (4.147)
    f <- stats::uniroot(fun, c(0,2*pi))
    delta <- f$root

    # Flow area of each phase
    AL <- glfMB:::circle_area(D / 2) * HL         # (4.147)
    AG <- glfMB:::circle_area(D / 2) * (1 - HL)

    # 4.146 for hL/d
    #   4.150 and 4.151
    dhG <- D * (2*pi - (delta - sin(delta))) / (2*pi - delta + 2 * sin(delta/2))    # (4.150)
    dhL <- D * (delta - sin(delta)) / (delta + 2 * sin(delta / 2))                  # (4.151)

    # Perimeter
    P <- D * pi
    PG <- (1 - delta / (2*pi)) * P    # 4.149
    PL <- P - PG                      # 4.148

    # Actual velocity
    vG <- vsG / (1-HL)           # 4.157
    vL <- vsL / HL               # 4.156

    ReG <- glfMB:::Reynolds(densityG, vG, dhG, viscosityG)   # (4.155)
    ReL <- glfMB:::Reynolds(densityL, vL, dhL, viscosityL)   # (4.154)

    fDG <- glfMB:::Darcy_friction_factor(ReG, roughness, D)
    fDL <- glfMB:::Darcy_friction_factor(ReL, roughness, D)

    shearStressG <- fDG * densityG * vG^2 / (2*g)    # 4.153
    shearStressL <- fDL * densityL * vL^2 / (2*g)    # 4.152

    if (debug == TRUE) {
      cat(sprintf("delta: %.2f, dhG: %.2f, dhL: %.2f", delta, dhG, dhL))
      cat(sprintf("PG: %.2f, PL: %.2f, ReG: %.1f, ReL: %.1f"), PG, PL, ReG, ReL)
    }

    # dPdL (4.144)
    dPdL <- - (shearStressL * PL + shearStressG * PG) - (viscosityL * AL + viscosityG * AG) * g * sin(angle)
  } else {
    vmix <- vsG + vsL
    HLnoslip <- vsL / vmix

    densityMixS <- densityG * (1 - HL) + densityL * HL                    # (3.22)
    densityMixN <- densityG * (1 - HLnoslip) + densityL * HLnoslip        # (3.23)
    viscosityMixS <- viscosityG * (1 - HL) + viscosityL * HL              # (3.19)
    viscosityMixN <- viscosityG * (1 - HLnoslip) + viscosityL * HLnoslip  # (3.21)

    ReN <- glfMB:::Reynolds(densityMixN, vmix, D, viscosityMixN)

    if (is.na(pressure) == TRUE) {
      Ek <- 0
    } else {
      Ek <- densityMixS * vmix * vsG / pressure    # (4.53)  (4.137)
    }

    if (flowRegime == 2) {
      # Annular
      fn <- glfMB:::Darcy_friction_factor(ReN, roughness, D)    # no-slip friction factor
      HR <- HLnoslip / HL                    # (4.140)
      fR <- glfMB:::ffRatio(HR)              # friction factor ratio
      fD <- fn * fR                          # (4.141)

      dPdL <- (fD * densityMixS * vmix^2 / (2 * D) + densityMixS * g * sin(angle)) / (1 - Ek)  # (4.139)
    } else if (flowRegime == 3 | flowRegime == 4) {
      # Slug or Bubble
      fD <- glfMB:::Darcy_friction_factor(ReN, roughness, D)
      dPdL <- (fD * densityMixS * vmix^2 / (2 * D) + densityMixS * g * sin(angle)) / (1 - Ek)  # (4.136)
    } else {
      # error!
      stop("Undefined flow regime.")
    }
  }

  ret <- dPdL

  if (debug == TRUE) {
    if (flowRegime == 1) {
      # Stratified
      ret <- c('dPdL'=ret, 'delta'=delta, 'dhG'=dhG, 'dhL'=dhL, 'PG'=PG, 'PL'=PL,
               'ReG'=ReG, 'ReL'=ReL, 'fDG'=fDG, 'fDL'=fDL,
               'shearStressG'=shearStressG, 'shearStressL'=shearStressL)
    } else if (flowRegime == 2) {
      # Annular
      ret <- c('dPdL'=ret, 'densityMixS'=densityMixS, 'densityMixN'=densityMixN,
               'viscosityMixS'=viscosityMixS, 'viscosityMixN'=viscosityMixN,
               'ReN'=ReN, 'Ek'=Ek, 'fn'=fn, 'HR'=HR, 'fR'=fR, 'fD'=fD)
    } else {
      # Slug or Bubble
      ret <- c('dPdL'=ret, 'densityMixS'=densityMixS, 'densityMixN'=densityMixN,
               'viscosityMixS'=viscosityMixS, 'viscosityMixN'=viscosityMixN,
               'ReN'=ReN, 'Ek'=Ek, 'fD'=fD)
    }
  }

  ret
}

