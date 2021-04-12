
# Dummy data for tests
ID = 0.1
densityG <- 1
densityL <- 1000
viscosityG <- 10^(-5)
viscosityL <- 10^(-3)
roughness <- 0.0001
pressure <- 100 * 1000    # 100 kPa


test_that("Stratified flow", {
  # Dummy data
  vsG <- 1
  vsL <- 1
  fr <- 1
  
  # stratified (holdup = 0.5, angle = 0)
  ret_st <- glfMB:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, 0,
                                   fr, 0.5, roughness, pressure, TRUE)
  
  exp_dh  <- (pi*ID/2) / (1+pi/2)                                # Dh = 4A / P
  exp_ReG <- exp_dh * (vsG/(1-holdup)) * densityG / viscosityG   # Reynolds number
  exp_ReL <- exp_dh * (vsL/holdup) * densityL /  viscosityL
  
  expect_equal(ret_st['delta'], c('delta' = pi))      # angle related to liquid level
  expect_equal(ret_st['PG'], c('PG' = ID*pi/2))
  expect_equal(ret_st['PL'], c('PL' = ID*pi/2))
  expect_equal(ret_st['dhG'], c('dhG' = exp_dh))    # = 4AG / PG
  expect_equal(ret_st['dhL'], c('dhL' = exp_dh))    # = 4AL / PL
  expect_equal(ret_st['ReG'], c('ReG' = exp_ReG))
  expect_equal(ret_st['ReL'], c('ReL' = exp_ReL))
  expect_equal(ret_st['fDG'], c('fDG' = Colebrook(exp_ReG, roughness, ID)))
  expect_equal(ret_st['fDL'], c('fDL' = Colebrook(exp_ReL, roughness, ID)))
  
  dPdL_FG <- ret_st['fDG'] * densityG * (vsG/0.5)^2 / 2 / ID
  dPdL_FL <- ret_st['fDL'] * densityL * (vsL/0.5)^2 / 2 / ID 
  expect_equal(as.numeric(ret_st['dPdL']), as.numeric((dPdL_FG + dPdL_FL) / 2))
  
  
  # stratified (holdup = 3/4 + 1/(2*pi), angle = 0)
  holdup2 <- 3/4 + 1/(2*pi)
  ret_st2 <- glfMB:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, angle,
                                    fr, holdup2, roughness, pressure, TRUE)
  
  exp_dhG2 <- 4*(1-holdup2)*(ID/2)^2*pi / (ID*pi*1/4 + sqrt(2)*ID/2)    # = 4AG / PG
  exp_dhL2 <- 4*holdup2*(ID/2)^2*pi     / (ID*pi*3/4 + sqrt(2)*ID/2)    # = 4AL / PL
  exp_ReG2 <- exp_dhG2 * (vsG/(1-holdup2)) * densityG / viscosityG
  exp_ReL2 <- exp_dhL2 * (vsL/holdup2) * densityL /  viscosityL
  
  expect_equal(ret_st2['delta'], c('delta' = pi*3/2))     # angle related to liquid level
  expect_equal(ret_st2['PG'], c('PG' = ID*pi*1/4))
  expect_equal(ret_st2['PL'], c('PL' = ID*pi*3/4), tolerance = 1e-7)
  expect_equal(ret_st2['dhG'], c('dhG' = exp_dhG2))
  expect_equal(ret_st2['dhL'], c('dhL' = exp_dhL2))
  expect_equal(ret_st2['ReG'], c('ReG' = exp_ReG2))
  expect_equal(ret_st2['ReL'], c('ReL' = exp_ReL2))
  expect_equal(ret_st2['fDG'], c('fDG' = Colebrook(exp_ReG2, roughness, ID)))
  expect_equal(ret_st2['fDL'], c('fDL' = Colebrook(exp_ReL2, roughness, ID)))
  
  dPdL_FG2 <- ret_st2['fDG'] * densityG * (vsG/(1-holdup2))^2 / 2 / ID
  dPdL_FL2 <- ret_st2['fDL'] * densityL * (vsL/holdup2)^2     / 2 / ID 
  expect_equal(as.numeric(ret_st2['dPdL']),
               as.numeric(dPdL_FG2 * (1/4) + dPdL_FL2 * (3/4)))
  
  
  # stratified (holdup = 0.5, angle = -90 deg)
  ret_st3 <- glfMB:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, -pi/2,
                                    fr, 0.5, roughness, pressure, TRUE)
  
  expect_equal(ret_st3['delta'], c('delta' = pi))      # angle related to liquid level
  expect_equal(ret_st3['PG'], c('PG' = ID*pi/2))
  expect_equal(ret_st3['PL'], c('PL' = ID*pi/2))
  expect_equal(ret_st3['dhG'], c('dhG' = exp_dh))    # = 4AG / PG
  expect_equal(ret_st3['dhL'], c('dhL' = exp_dh))    # = 4AL / PL
  expect_equal(ret_st3['ReG'], c('ReG' = exp_ReG))
  expect_equal(ret_st3['ReL'], c('ReL' = exp_ReL))
  expect_equal(ret_st3['fDG'], c('fDG' = Colebrook(exp_ReG, roughness, ID)))
  expect_equal(ret_st3['fDL'], c('fDL' = Colebrook(exp_ReL, roughness, ID)))
  expect_equal(as.numeric(ret_st3['dPdL']), as.numeric(dPdL_FG/2 + dPdL_FL/2 - (densityG + densityL)/2 * glfMB:::g))
  
  # stratified (holdup = 0.5, angle = -30 deg)
  ret_st4 <- glfMB:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, -pi/6,
                                    fr, 0.5, roughness, pressure, TRUE)
  expect_equal(as.numeric(ret_st4['dPdL']), as.numeric(dPdL_FG/2 + dPdL_FL/2 + sin(-pi/6) * (densityG + densityL)/2 * glfMB:::g))
  
  
  
})
