
P = 100    # 100 psi


#denstiyG = td$ideal_gas_density(temperature, pressure, glfMB:::Mair)
densityL = 1000
#viscosityG = td$Sutherland_viscosityG(temperature, s$vis0_air, s$T0_air, s$vis0_air)
viscosityL = 0.0001



test_that("multiplication works", {
  
  #dlns <- l_dlns_MB(vsG, vsL, D, densityG, densityL,
  #                 viscosityG, viscosityL, surfaceTension, angle)
  #fr   <- l_flow_regime_MB(dlns)
  #hl <- l_holdup_MB(dlns, fr)
  #dPdL <- l_dPdL_MB(dlns, fr, hl, roughness, pressure, debug=FALSE)
  
  
  
  expect_equal(2 * 2, 4)
})
