


test_that("laminar()", {
  # laminar flow
  Re_lam <- c(500, 1000, 1500, 2000)
  expect_equal(MukherjeeBrill:::laminar(Re_lam), 64/Re_lam)
})

# Compare return values with https://www.engineeringtoolbox.com/colebrook-equation-d_1031.html
test_that("Colebrook()", {
  ID <- 0.1
  Re <- c(5000, 10^4, 10^5, 10^6)
  
  # relative roughness: 0.001
  expect_equal(round(util_MB_Colebrook(Re, 0.0001, ID), 4), c(0.0385, 0.0324, 0.0222, 0.0199))
  # relative roughness: 0.01
  expect_equal(round(util_MB_Colebrook(Re, 0.001, ID), 4), c(0.0472, 0.0431, 0.0384, 0.0379), tolerance=0.005)
  # relative roughness: 0.05
  expect_equal(round(util_MB_Colebrook(Re, 0.005, ID), 4), c(0.0758, 0.0736, 0.0716, 0.0716), tolerance=0.005)
})






# - - - - - - - - - - 

test_that("l_Darcy_friction_factor()", {
  ID <- 0.1
  
  # laminar flow
  Re_lam <- c(500, 1000, 1500, 2000)
  expect_equal(MukherjeeBrill:::l_Darcy_friction_factor(Re_lam, 0.0001, ID), 64/Re_lam)
  expect_equal(MukherjeeBrill:::l_Darcy_friction_factor_core(Re_lam[1], 0.0001, ID), 64/Re_lam[1])
  expect_equal(MukherjeeBrill:::l_Darcy_friction_factor_core(Re_lam[2], 0.0001, ID), 64/Re_lam[2])
  
  
  # Turbulent flow
  Re_tur <- c(5000, 10^4, 10^5, 10^6)
  expect_equal(round(MukherjeeBrill:::l_Darcy_friction_factor(Re_tur, 0.0001, ID), 4), c(0.0385, 0.0324, 0.0222, 0.0199))
})


