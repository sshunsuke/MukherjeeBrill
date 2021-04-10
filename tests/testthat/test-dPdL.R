
P = 100    # 100 psi

ID = 0.052  # m


#denstiyG = td$ideal_gas_density(temperature, pressure, glfMB:::Mair)
densityL = 1000
#viscosityG = td$Sutherland_viscosityG(temperature, s$vis0_air, s$T0_air, s$vis0_air)
viscosityL = 0.0001



test_that("Stratified flow", {
  # Dummy data
  ID = 0.1
  vsG <- 1
  vsL <- 1
  densityG <- 1
  densityL <- 1000
  viscosityG <- 10^(-5)
  viscosityL <- 10^(-3)
  angle <- 0
  fr <- 1
  holdup <- 0.5
  roughness <- 0.0001
  pressure <- 100 * 1000    # 100 kPa
  
  # stratified
  ret0 <- glfMB:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, angle,
                                 fr, 0.99, roughness, pressure, TRUE)
  
  
  
  
  
  ret_st <- glfMB:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, angle,
                                   fr, holdup, roughness, pressure, TRUE)
  
  expect_equal(as.numeric(ret_st['delta']), pi)     # angle related to liquid level
  expect_equal(as.numeric(ret_st['PG']), ID*pi/2)
  expect_equal(as.numeric(ret_st['PL']), ID*pi/2)
  expect_equal(as.numeric(ret_st['dhG']), ID*(2 * pi - pi + sin(pi))/(2 * pi - pi + 2*sin(pi/2)))
  expect_equal(as.numeric(ret_st['dhL']), ID*(pi - sin(pi))/(pi + 2*sin(pi/2)))
  expect_equal(as.numeric(ret_st['ReG']), ID*(pi + sin(pi))/(pi + 2*sin(pi/2)) * (vsG/(1-holdup)) * densityG /  viscosityG)
  expect_equal(as.numeric(ret_st['ReL']), ID*(pi - sin(pi))/(pi + 2*sin(pi/2)) * (vsL/holdup) * densityL /  viscosityL)
  
  
  
  
  ret_st2 <- glfMB:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, angle,
                                   fr, 3/4 + 1/(2*pi), roughness, pressure, TRUE)
  expect_equal(as.numeric(ret_st2['delta']), pi*3/2, tolerance = 1e-5)     # angle related to liquid level
  
  
  expect_equal(2 * 2, 4)
})
