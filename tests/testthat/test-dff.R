






test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})




test_that("dlns_MB(), holdup_MB(): Slug flow", {
  vsG <- 3.86 * 0.3048    # 3.86 ft/s
  vsL <- 3.97 * 0.3048    # 3.97 ft/s
  ID  <- 6 * 0.0254       # 6 inch
  densityG <- 5.88 * 16.01845     # 5.88 lbm/ft3  - (1 lbm/ft3 = 16.01845 kg/m3)
  densityL <- 47.61 * 16.01845    # 47.61 lbm/ft3 - (1 lbm/ft3 = 16.01845 kg/m3)
  viscosityG <- 0.016 / 1000      # 0.016 cp
  viscosityL <- 0.97  / 1000      # 0.970 cp
  surfaceTension <- 8.41 / 1000   # 8.41 dynes/cm
  angle <- pi/2                   # 90 deg
  
  dlns <- dlns_MB(vsG, vsL, ID, densityG, densityL,
                  viscosityG, viscosityL, surfaceTension, angle)
  fr <- flow_regime_MB(dlns)
  hl <- holdup_MB(dlns, fr)
  #dPdL <- dPdL_MB(dlns, fr, hl, roughness, pressure, debug=TRUE)
  
  # NLv=11.84, NGv=11.54, NGvSM=350.8, NLvBS_up=18.40, NL=0.0118
  #expect_equal(dlns$NLv, 11.84)
  expect_equal(round(dlns$NGv, 2), 11.54)
  expect_equal(round(dlns$NGvSM, 1), 350.8)
  expect_equal(round(dlns$NLvBS_up, 1), 18.4)
  expect_equal(round(dlns$NL, 4), 0.0118)
  
  expect_equal(fr, 3)
  expect_equal(round(hl,2), 0.56)
  
  #expect_equal(dPdL[1,'dPdL'], 4730)
  
})
