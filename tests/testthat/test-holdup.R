# Mukherjee, Hemanta, and James P. Brill. 1983. “Liquid Holdup Correlations for Inclined Two-Phase Flow.” JPT, Journal of Petroleum Technology 35(5):1003–8.

temperature <- 10 + 273.15   # 10 degC
pressure <- 8.6 * 100000     # 8.6 bar
ID <- 0.0381                 # 3.81 cm

td <- testdata_MB()

densityG <- td$Air_density(temperature, pressure)
densityL <- td$Kerosene_density(temperature)
viscosityG <- td$Air_viscosity(temperature, pressure)
viscosityL <- td$Kerosene_viscosity(temperature)
surfaceTension <- td$Kerosene_surfaceTension(temperature)


check_range <- function(val, calc, tol) {
  min_ <- calc - tol
  max_ <- calc + tol
  
  mapply(expect_gt, val, min_)
  mapply(expect_lt, val, max_)
}



# Note: Calculation at 0 deg is not explicitly defined
test_that("Fig.7: 0 deg ", {
  angle <- deg2rad(0)
  vsL <- ft2m( c(0.094, 0.363, 3.9, 7.3, 12) ) 
  vsG <- ft2m( c(5, 20, 40, 60) )
  
  dlns_1 <- l_dlns_MB(vsG, vsL[1], 
                      ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr_1 <- l_flow_regime_MB(dlns_1)
  hl_1 <- l_holdup_MB(dlns_1, fr_1)

  dlns_2 <- l_dlns_MB(vsG, vsL[2], 
                      ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr_2 <- l_flow_regime_MB(dlns_2)
  hl_2 <- l_holdup_MB(dlns_2, fr_2)
  
  dlns_3 <- l_dlns_MB(vsG, vsL[3], 
                      ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr_3 <- l_flow_regime_MB(dlns_3)
  hl_3 <- l_holdup_MB(dlns_3, fr_3)
 
  
  dlns_4 <- l_dlns_MB(vsG, vsL[4], 
                      ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr_4 <- l_flow_regime_MB(dlns_4)
  hl_4 <- l_holdup_MB(dlns_4, fr_4)
  
  dlns_5 <- l_dlns_MB(vsG, vsL[5], 
                      ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr_5 <- l_flow_regime_MB(dlns_5)
  hl_5 <- l_holdup_MB(dlns_5, fr_5)
  
  
  print(densityL)
  print(viscosityL)
  print(surfaceTension)
  
  print("Fig.7: 0 deg")
  print(sprintf("%f ", 1-hl_1))
  print(sprintf("%f ", 1-hl_2))
  print(sprintf("%f ", 1-hl_3))
  print(sprintf("%f ", 1-hl_4))
  print(sprintf("%f ", 1-hl_5))
  
  check_range(1-hl_1, c(0.8, 0.95, 0.95, 0.95), 0.05)
  check_range(1-hl_2, c(0.7, 0.9, 0.95, 0.95), 0.05)
  check_range(1-hl_3, c(0.45, 0.7, 0.8, 0.85), 0.05)
  check_range(1-hl_4, c(0.4, 0.6, 0.75, 0.8), 0.05)
  check_range(1-hl_5, c(0.35, 0.55, 0.75, 0.8), 0.05)
})

test_that("Fig.8: 30 deg ", {
  angle <- deg2rad(30)
  vsL <- ft2m( c(0.094, 3.9, 12) ) 
  vsG <- ft2m( c(5, 20, 40, 60) )
  
  dlns_1 <- l_dlns_MB(vsG, vsL[1], 
                      ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr_1 <- l_flow_regime_MB(dlns_1)
  hl_1 <- l_holdup_MB(dlns_1, fr_1)
  check_range(1-hl_1, c(0.85, 0.95, 0.95, 0.95), 0.05)
  
  dlns_2 <- l_dlns_MB(vsG, vsL[2], 
                      ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr_2 <- l_flow_regime_MB(dlns_2)
  hl_2 <- l_holdup_MB(dlns_2, fr_2)
  check_range(1-hl_2, c(0.45, 0.65, 0.75, 0.85), 0.05)
  
  dlns_3 <- l_dlns_MB(vsG, vsL[3], 
                      ID, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
  fr_3 <- l_flow_regime_MB(dlns_3)
  hl_3 <- l_holdup_MB(dlns_3, fr_3)
  check_range(1-hl_3, c(0.35, 0.55, 0.65, 0.7), 0.05)
  
  print("Fig.8: 30 deg")
  print(sprintf("%f ", 1-hl_1))  
  print(sprintf("%f ", 1-hl_2))
  print(sprintf("%f ", 1-hl_3)) 
})





