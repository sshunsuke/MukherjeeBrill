# - - - - - - - - - - - - - - - - - - - - - - - - -
# Constants ----
# * pi ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

g <- 9.8              # Gravitational acceleration (m/s2)
R <- 8.314            # Gas constant





# - - - - - - - - - - - - - - - - - - - - - - - - -
# Utilities (not exported) ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

circle_area = function(D) { (D/2)^2 * pi }

Reynolds <- function(density, v, D, viscosity) {
  density * v * D / viscosity
}

# Converter from radian to degree
rad2deg = function(rad) { rad * 180 / pi }

# Convertor from degree to radian
deg2rad = function(deg) { deg * pi / 180 }


ft2m <- function(ft) { ft * 0.3048 }
m2ft <- function(m) { m / 0.3048 }

psi2Pa <- function(psi) { psi * 6894.76}






