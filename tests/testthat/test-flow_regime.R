temperature <- 10 + 273.15   # 10 degC
pressure <- 8.6 * 100000     # 8.6 bar
ID <- 0.0381                 # 3.81 cm

td <- testdata_MB()

densityG <- td$Air_density(temperature, pressure)
densityL <- td$Kerosene_density(temperature)
viscosityG <- td$Air_viscosity(temperature, pressure)
viscosityL <- td$Kerosene_viscosity(temperature)
surfaceTension <- td$Kerosene_surfaceTension(temperature)

coefficient <- (densityL / glfMB:::g / surfaceTension)^(0.25)



# - - - - - - - - - - - - - - - - - - - - - - - - -
# Kerosene, 70 deg downflow ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

test_that("Kerosene, -70 deg: Stratifired", {
  angle <- glfMB:::rad2deg(-70)
  NGv <- c(0.5,   5,  10, 1, 10, 2)
  NLv <- c(0.5, 0.5, 0.5, 1,  1, 2)
  dlns <- l_dlns_MB(NGv/coefficient, NLv/coefficient,
                  ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr <- l_flow_regime_MB(dlns)
  expect_equal(fr, rep(1, 6))
})

test_that("Kerosene, -70 deg: Annular", {
  angle <- glfMB:::rad2deg(-70)
  NGv <- c( 50, 100, 200, 300)
  NLv <- c(0.1,  1,    5,  8)
  dlns <- l_dlns_MB(NGv/coefficient, NLv/coefficient,
                  ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr <- l_flow_regime_MB(dlns)
  expect_equal(fr, rep(2, 4))
})

test_that("Kerosene, -70 deg: Slug", {
  angle <- glfMB:::rad2deg(-70)
  NGv <- c(10, 20, 20, 100)
  NLv <- c(10,  5, 10,   2)
  dlns <- l_dlns_MB(NGv/coefficient, NLv/coefficient,
                    ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr <- l_flow_regime_MB(dlns)
  expect_equal(fr, rep(3, 4))
})

test_that("Kerosene, -70 deg: Bubble", {
  angle <- glfMB:::rad2deg(-70)
  NGv <- c(0.1, 0.1, 1, 10, 100)
  NLv <- c(0.5, 100, 5, 20,  50)
  dlns <- l_dlns_MB(NGv/coefficient, NLv/coefficient,
                  ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr <- l_flow_regime_MB(dlns)
  expect_equal(fr, rep(4, 5))
})

test_that("Kerosene, -70 deg: generate_frm_MB()", {
  angle <- glfMB:::rad2deg(-70)
  NGv <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
  NLv <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)

  frm <- generate_frm_MB(NGv/coefficient, NLv/coefficient, ID, densityG, densityL,
                        viscosityG, viscosityL, surfaceTension, angle)
  
  f_get_fr <- function(frm_, NGv_, NLv_) {
    rn_NGv <- which(round(frm_[,'NGv'], 1) == NGv_)
    rn_NLv <- which(round(frm_[,'NLv'], 1) == NLv_)
    rn <- intersect(rn_NGv, rn_NLv)
    frm_[rn, 'fr']
  }
  
  plot(frm, log="xy", xval='NGv', yval='NLv', xlim=c(0.1, 200), ylim=c(0.1, 200), main="Kerosene, -70 deg")

  expect_equal(f_get_fr(frm, 0.1, 0.2), c('fr'=1))
  expect_equal(f_get_fr(frm, 0.1, 0.5), c('fr'=4))
  expect_equal(f_get_fr(frm, 1, 2), c('fr'=1))
  expect_equal(f_get_fr(frm, 1, 5), c('fr'=4))
  expect_equal(f_get_fr(frm, 10, 5), c('fr'=1))
  expect_equal(f_get_fr(frm, 10, 10), c('fr'=3))
  expect_equal(f_get_fr(frm, 10, 20), c('fr'=4))
  expect_equal(f_get_fr(frm, 100, 1), c('fr'=2))
  expect_equal(f_get_fr(frm, 100, 2), c('fr'=3))
  expect_equal(f_get_fr(frm, 100, 20), c('fr'=3))
  expect_equal(f_get_fr(frm, 100, 50), c('fr'=4))
})


# - - - - - - - - - - - - - - - - - - - - - - - - -
# Kerosene, 0 deg ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

test_that("Kerosene, 0 deg: Stratifired", {
  angle <- 0
  NGv <- c(0.1, 1, 10,  20)
  NLv <- c(  1, 1,  1, 0.5)
  dlns <- l_dlns_MB(NGv/coefficient, NLv/coefficient,
                  ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr <- l_flow_regime_MB(dlns)
  expect_equal(fr, rep(1, 4))
})

test_that("Kerosene, 0 deg: Annular", {
  angle <- 0
  NGv <- c( 50, 100, 200)
  NLv <- c(0.2,   1,   5)
  dlns <- l_dlns_MB(NGv/coefficient, NLv/coefficient,
                  ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr <- l_flow_regime_MB(dlns)
  expect_equal(fr, rep(2, 3))
})

test_that("Kerosene, 0 deg: Slug", {
  angle <- 0
  NGv <- c(5,   5, 20, 50, 100)
  NLv <- c(2, 100,  1,  1,   2)
  dlns <- l_dlns_MB(NGv/coefficient, NLv/coefficient,
                  ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr <- l_flow_regime_MB(dlns)
  expect_equal(fr, rep(3, 5))
})

test_that("Kerosene, 0 deg: generate_frm_MB()", {
  angle <- glfMB:::rad2deg(0)
  NGv <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
  NLv <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
  
  frm <- generate_frm_MB(NGv/coefficient, NLv/coefficient, ID, densityG, densityL,
                         viscosityG, viscosityL, surfaceTension, angle)
  
  f_get_fr <- function(frm_, NGv_, NLv_) {
    rn_NGv <- which(round(frm_[,'NGv'], 1) == NGv_)
    rn_NLv <- which(round(frm_[,'NLv'], 1) == NLv_)
    rn <- intersect(rn_NGv, rn_NLv)
    frm_[rn, 'fr']
  }
  
  plot(frm, log="xy", xval='NGv', yval='NLv', xlim=c(0.1, 200), ylim=c(0.1, 200), main="Kerosene, 0 deg")
  
  expect_equal(f_get_fr(frm, 0.1, 1), c('fr'=1))
  expect_equal(f_get_fr(frm, 0.1, 2), c('fr'=4))
  expect_equal(f_get_fr(frm, 1, 1), c('fr'=1))
  expect_equal(f_get_fr(frm, 1, 2), c('fr'=4))
  expect_equal(f_get_fr(frm, 10, 1), c('fr'=1))
  expect_equal(f_get_fr(frm, 10, 2), c('fr'=3))
  expect_equal(f_get_fr(frm, 100, 1), c('fr'=2))
  expect_equal(f_get_fr(frm, 100, 2), c('fr'=3))
})


# - - - - - - - - - - - - - - - - - - - - - - - - -
# Kerosene, 70 deg upflow ----
# - - - - - - - - - - - - - - - - - - - - - - - - -

test_that("flow_regime_MB(): Kerosene, 70 deg, Annular", {
  angle <- glfMB:::rad2deg(70)
  NGv <- c( 50, 100, 200)
  NLv <- c(0.2,   1,   5)
  dlns <- l_dlns_MB(NGv/coefficient, NLv/coefficient,
                  ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr <- l_flow_regime_MB(dlns)
  expect_equal(fr, rep(2, 3))
})

test_that("flow_regime_MB(): Kerosene, 70 deg, Slug", {
  angle <- glfMB:::rad2deg(70)
  NGv <- c(0.1, 1, 10,  50,  50)
  NLv <- c(0.2, 2, 20, 0.5, 100)
  dlns <- l_dlns_MB(NGv/coefficient, NLv/coefficient,
                    ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr <- l_flow_regime_MB(dlns)
  expect_equal(fr, rep(3, 5))
})

test_that("flow_regime_MB(): Kerosene, 70 deg, Bubble", {
  angle <- glfMB:::rad2deg(-70)
  NGv <- c(0.1, 1, 10, 50)
  NLv <- c(0.5, 5, 50, 200)
  dlns <- l_dlns_MB(NGv/coefficient, NLv/coefficient,
                  ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr <- l_flow_regime_MB(dlns)
  expect_equal(fr, rep(4, 4))
})

test_that("Kerosene, 70 deg upflow: generate_frm_MB()", {
  angle <- glfMB:::rad2deg(70)
  NGv <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
  NLv <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
  
  frm <- generate_frm_MB(NGv/coefficient, NLv/coefficient, ID, densityG, densityL,
                         viscosityG, viscosityL, surfaceTension, angle)
  
  f_get_fr <- function(frm_, NGv_, NLv_) {
    rn_NGv <- which(round(frm_[,'NGv'], 1) == NGv_)
    rn_NLv <- which(round(frm_[,'NLv'], 1) == NLv_)
    rn <- intersect(rn_NGv, rn_NLv)
    frm_[rn, 'fr']
  }
  
  plot(frm, log="xy", xval='NGv', yval='NLv', xlim=c(0.1, 200), ylim=c(0.1, 200), main="Kerosene: upflow 70 deg")
  
  expect_equal(f_get_fr(frm, 0.1, 0.2), c('fr'=3))
  expect_equal(f_get_fr(frm, 0.1, 0.5), c('fr'=4))
  expect_equal(f_get_fr(frm, 1, 2), c('fr'=3))
  expect_equal(f_get_fr(frm, 1, 5), c('fr'=4))
  expect_equal(f_get_fr(frm, 10, 20), c('fr'=3))
  expect_equal(f_get_fr(frm, 10, 50), c('fr'=4))
  expect_equal(f_get_fr(frm, 100, 1), c('fr'=2))
  expect_equal(f_get_fr(frm, 100, 5), c('fr'=3))
})

