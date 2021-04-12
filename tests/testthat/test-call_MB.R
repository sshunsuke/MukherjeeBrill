# Test code for call_MB()
#
# * call_MB()
# * l_dlns_MB(), l_holdup_MB(), l_dPdL_MB



test_that("call_MB(): Slug flow", {
  vsG <- 3.86 * 0.3048            # 3.86 ft/s
  vsL <- 3.97 * 0.3048            # 3.97 ft/s
  ID  <- 6 * 0.0254               # 6 inch
  densityG <- 5.88 * 16.01845     # 5.88 lbm/ft3  - (1 lbm/ft3 = 16.01845 kg/m3)
  densityL <- 47.61 * 16.01845    # 47.61 lbm/ft3 - (1 lbm/ft3 = 16.01845 kg/m3)
  viscosityG <- 0.016 / 1000      # 0.016 cp
  viscosityL <- 0.97  / 1000      # 0.970 cp
  surfaceTension <- 8.41 / 1000   # 8.41 dynes/cm
  angle <- pi/2                   # 90 deg
  roughness <- 0.00006 * 0.3048   # 0.00006 ft
  pressure <- 1700 * 6894.76      # 1700 psia (Example 3.2)

  ret <- call_MB(vsG, vsL, ID, densityG, densityL, viscosityG, viscosityL,
                 surfaceTension, angle, roughness, pressure)
  
  expect_equal(ret$fr, 3)
  expect_equal(ret$hl, 0.560, tolerance = 0.0005)  # diff <= 0.005%
  expect_equal(ret$dPdL, 4727, tolerance = 0.005)  # diff <= 0.05%
  expect_equal(ret$dPdL_H, 29.249*0.006944*6894.76 / 0.3048, tolerance = 0.005)
  expect_equal(ret$dPdL_F, 0.864*0.006944*6894.76 / 0.3048, tolerance = 0.005)
})


test_that("l_dlns_MB(), l_holdup_MB(), l_dPdL_MB: Slug flow", {
  vsG <- 3.86 * 0.3048            # 3.86 ft/s
  vsL <- 3.97 * 0.3048            # 3.97 ft/s
  ID  <- 6 * 0.0254               # 6 inch
  densityG <- 5.88 * 16.01845     # 5.88 lbm/ft3  - (1 lbm/ft3 = 16.01845 kg/m3)
  densityL <- 47.61 * 16.01845    # 47.61 lbm/ft3 - (1 lbm/ft3 = 16.01845 kg/m3)
  viscosityG <- 0.016 / 1000      # 0.016 cp
  viscosityL <- 0.97  / 1000      # 0.970 cp
  surfaceTension <- 8.41 / 1000   # 8.41 dynes/cm
  angle <- pi/2                   # 90 deg
  roughness <- 0.00006 * 0.3048   # 0.00006 ft
  pressure <- 1700 * 6894.76      # 1700 psia (Example 3.2)
  
  dlns <- l_dlns_MB(vsG, vsL, ID, densityG, densityL,
                  viscosityG, viscosityL, surfaceTension, angle)
  fr <- l_flow_regime_MB(dlns)
  hl <- l_holdup_MB(dlns, fr)
  dPdL <- l_dPdL_MB(dlns, fr, hl, roughness, pressure, debug=TRUE)
  
  # Check dlns
  #   NLv=11.87, NGv=11.54, NGvSM=350.8, NLvBS_up=18.40, NL=0.0118
  expect_equal(dlns$NGv, 11.54, tolerance=0.001)       # diff <= 0.1% of 11.54
  expect_equal(dlns$NLv, 11.87, tolerance=0.001)       # diff <= 0.1% of 11.84
  expect_equal(dlns$NGvSM, 350.8, tolerance=0.0001)    # diff <= 0.01% of 350.8
  expect_equal(dlns$NLvBS_up, 18.40, tolerance=0.001)  # diff <= 0.1%
  expect_equal(dlns$NL, 0.0118, tolerance=0.004)       # diff <- 0.4%
  
  expect_equal(fr, 3)
  expect_equal(hl, 0.560, tolerance=0.001)            # diff <= 0.1%
  
  expect_equal(dPdL[1,'dPdL'], c(dPdL=4727), tolerance = 0.005)  # diff <= 0.05%
  expect_equal(dPdL[1,'dPdL_H'], c('dPdL_H'=29.249*0.006944*6894.76 / 0.3048), tolerance = 0.005)
  expect_equal(dPdL[1,'dPdL_F'], c('dPdL_F'=0.864*0.006944*6894.76 / 0.3048), tolerance = 0.005)
  
  expect_equal(dPdL[1,'ReN'], c(ReN=315000), tolerance = 0.001)  # diff <= 0.01%
  expect_equal(dPdL[1,'fD'], c(fD=0.0155), tolerance = 0.0005)   # diff <= 0.005%
  expect_equal(dPdL[1,'densityMixS'], c(densityMixS=densityG*(1-hl)+densityL*hl))
})
