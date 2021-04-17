
# Dummy data for tests
ID = 0.1
densityG <- 1
densityL <- 1000
viscosityG <- 10^(-5)
viscosityL <- 10^(-3)
roughness <- 0.0001
pressure <- 10 * 1000 * 1000    # 10 MPa


test_that("l_dPdL_core_MB(): Stratified flow", {
  # Dummy data
  vsG <- 1
  vsL <- 1
  fr <- 1
  
  # stratified (holdup = 0.5, angle = 0)
  holdup <- 0.5
  ret_st <- MukherjeeBrill:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, 0,
                                   fr, holdup, roughness, pressure, TRUE)
  
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
  expect_equal(ret_st['fDG'], c('fDG' = util_MB_Colebrook(exp_ReG, roughness, ID)))
  expect_equal(ret_st['fDL'], c('fDL' = util_MB_Colebrook(exp_ReL, roughness, ID)))
  
  dPdL_FG <- ret_st['fDG'] * densityG * (vsG/0.5)^2 / 2 / ID
  dPdL_FL <- ret_st['fDL'] * densityL * (vsL/0.5)^2 / 2 / ID 
  expect_equal(as.numeric(ret_st['dPdL']), as.numeric((dPdL_FG + dPdL_FL) / 2))
  expect_equal(as.numeric(ret_st['dPdL_H']), 0)
  expect_equal(as.numeric(ret_st['dPdL_A']), 0)
  
  
  # stratified (holdup = 3/4 + 1/(2*pi), angle = 0)
  holdup2 <- 3/4 + 1/(2*pi)
  ret_st2 <- MukherjeeBrill:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, 0,
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
  expect_equal(ret_st2['fDG'], c('fDG' = util_MB_Colebrook(exp_ReG2, roughness, ID)))
  expect_equal(ret_st2['fDL'], c('fDL' = util_MB_Colebrook(exp_ReL2, roughness, ID)))
  
  dPdL_FG2 <- ret_st2['fDG'] * densityG * (vsG/(1-holdup2))^2 / 2 / ID
  dPdL_FL2 <- ret_st2['fDL'] * densityL * (vsL/holdup2)^2     / 2 / ID 
  expect_equal(as.numeric(ret_st2['dPdL']),
               as.numeric(dPdL_FG2 * (1/4) + dPdL_FL2 * (3/4)))
  
  
  # stratified (holdup = 0.5, angle = -90 deg)
  ret_st3 <- MukherjeeBrill:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, -pi/2,
                                    fr, 0.5, roughness, pressure, TRUE)
  
  expect_equal(ret_st3['delta'], c('delta' = pi))      # angle related to liquid level
  expect_equal(ret_st3['PG'], c('PG' = ID*pi/2))
  expect_equal(ret_st3['PL'], c('PL' = ID*pi/2))
  expect_equal(ret_st3['dhG'], c('dhG' = exp_dh))    # = 4AG / PG
  expect_equal(ret_st3['dhL'], c('dhL' = exp_dh))    # = 4AL / PL
  expect_equal(ret_st3['ReG'], c('ReG' = exp_ReG))
  expect_equal(ret_st3['ReL'], c('ReL' = exp_ReL))
  expect_equal(ret_st3['fDG'], c('fDG' = util_MB_Colebrook(exp_ReG, roughness, ID)))
  expect_equal(ret_st3['fDL'], c('fDL' = util_MB_Colebrook(exp_ReL, roughness, ID)))
  expect_equal(as.numeric(ret_st3['dPdL']), as.numeric(dPdL_FG/2 + dPdL_FL/2 - (densityG + densityL)/2 * MukherjeeBrill:::g))
  
  # stratified (holdup = 0.5, angle = -30 deg)
  ret_st4 <- MukherjeeBrill:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, -pi/6,
                                    fr, 0.5, roughness, pressure, TRUE)
  expect_equal(as.numeric(ret_st4['dPdL']), as.numeric(dPdL_FG/2 + dPdL_FL/2 + sin(-pi/6) * (densityG + densityL)/2 * MukherjeeBrill:::g))
})


test_that("l_dPdL_core_MB(): Annular flow", {
  # Dummy data
  vsG <- 3
  vsL <- 1
  fr <- 2
  
  # annular (holdup = 0.5, angle = 0)
  holdup <- 0.5
  ret_an <- MukherjeeBrill:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, 0,
                                   fr, holdup, roughness, pressure, TRUE)
  
  holdup_noslip <- vsL / (vsL+vsG)    # 1 / (1+3) = 0.25
  exp_ReN <- as.numeric(ID * ret_an['densityMixN'] * (vsG+vsL) / ret_an['viscosityMixN'])
  exp_densityMixS <- densityG*0.5+densityL*0.5
  exp_Ek <- exp_densityMixS * (vsG+vsL) * vsG / pressure
  
  expect_equal(ret_an['HR'], c('HR'=0.5))    # 0.25 / 0.5 = 0.5
  expect_equal(ret_an['fR'], c('fR'=1.3))
  expect_equal(ret_an['densityMixN'], c('densityMixN'=densityG*0.75+densityL*0.25))
  expect_equal(ret_an['viscosityMixN'], c('viscosityMixN'=viscosityG*0.75+viscosityL*0.25))
  expect_equal(ret_an['ReN'], c('ReN'=exp_ReN))
  expect_equal(ret_an['fn'], c('fn' = util_MB_Colebrook(exp_ReN, roughness, ID)))
  expect_equal(ret_an['fD'], c('fD' = util_MB_Colebrook(exp_ReN, roughness, ID) * 1.3))
  expect_equal(ret_an['Ek'], c('Ek' = exp_Ek))
  
  exp_dPdL_F <- ret_an['fD'] * ret_an['densityMixN'] * (vsL+vsG)^2 / 2 / ID
  names(exp_dPdL_F) <- NULL
  
  expect_equal(ret_an['dPdL_F'], c('dPdL_F' = exp_dPdL_F))
  expect_equal(ret_an['dPdL'], c('dPdL' = exp_dPdL_F / (1 - exp_Ek)))
  expect_equal(ret_an['dPdL_A'], c('dPdL_A' = exp_dPdL_F * exp_Ek / (1-exp_Ek)))
  expect_equal(as.numeric(ret_an['dPdL_H']), 0)
  
  
  # annular (holdup = 0.5, angle = -30 deg)
  ret_an2 <- MukherjeeBrill:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL,  -pi/6,
                                    fr, holdup, roughness, pressure, TRUE)
  expect_equal(ret_an2['HR'], c('HR'=0.5))    # 0.25 / 0.5 = 0.5
  expect_equal(ret_an2['fR'], c('fR'=1.3))
  expect_equal(ret_an2['densityMixN'], c('densityMixN'=densityG*0.75+densityL*0.25))
  expect_equal(ret_an2['viscosityMixN'], c('viscosityMixN'=viscosityG*0.75+viscosityL*0.25))
  expect_equal(ret_an2['ReN'], c('ReN'=exp_ReN))
  expect_equal(ret_an2['fn'], c('fn' = util_MB_Colebrook(exp_ReN, roughness, ID)))
  expect_equal(ret_an2['fD'], c('fD' = util_MB_Colebrook(exp_ReN, roughness, ID) * 1.3))
  expect_equal(ret_an2['Ek'], c('Ek' = exp_Ek))
  
  exp_dPdL_H <- sin(-pi/6) * (densityG + densityL)/2 * MukherjeeBrill:::g
  
  expect_equal(ret_an2['dPdL_F'], c('dPdL_F' = exp_dPdL_F))
  expect_equal(ret_an2['dPdL_H'], c('dPdL_H' = exp_dPdL_H))
  expect_equal(ret_an2['dPdL_A'], c('dPdL_A' = (exp_dPdL_F + exp_dPdL_H) * exp_Ek / (1-exp_Ek)))
  expect_equal(ret_an2['dPdL'], c('dPdL' = (exp_dPdL_F + exp_dPdL_H) / (1-exp_Ek)))
})


test_that("l_dPdL_core_MB(): Bubble flow", {
  # Dummy data
  vsG <- 1
  vsL <- 3
  fr <- 4
  
  # Bubble (holdup = 0.5, angle = 0)
  holdup <- 0.5
  ret_b <- MukherjeeBrill:::l_dPdL_core_MB(ID, vsG, vsL, densityG, densityL, viscosityG, viscosityL, 0,
                                  fr, holdup, roughness, pressure, TRUE)
  exp_ReN <- as.numeric(ID * ret_b['densityMixN'] * (vsG+vsL) / ret_b['viscosityMixN'])
  exp_densityMixS <- densityG*0.5+densityL*0.5
  exp_Ek <- exp_densityMixS * (vsG+vsL) * vsG / pressure
  
  expect_equal(ret_b['densityMixN'], c('densityMixN'=densityG*0.25+densityL*0.75))
  expect_equal(ret_b['viscosityMixN'], c('viscosityMixN'=viscosityG*0.25+viscosityL*0.75))
  expect_equal(ret_b['ReN'], c('ReN'=exp_ReN))
  expect_equal(ret_b['fD'], c('fD' = util_MB_Colebrook(exp_ReN, roughness, ID)))
  expect_equal(ret_b['Ek'], c('Ek' = exp_Ek))
  
  
  
})


