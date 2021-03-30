




# Compare 
# https://www.engineeringtoolbox.com/colebrook-equation-d_1031.html

test_that("Colebrook()", {
  ID <- 0.1
  Re <- c(5000, 10^4, 10^5, 10^6)
  
  # relative roughness: 0.001
  expect_equal(round(Colebrook(Re, 0.0001, ID), 4), c(0.0385, 0.0324, 0.0222, 0.0199))
  # relative roughness: 0.01
  expect_equal(round(Colebrook(Re, 0.001, ID), 4), c(0.0472, 0.0431, 0.0384, 0.0379), tolerance=0.005)
  # relative roughness: 0.05
  expect_equal(round(Colebrook(Re, 0.005, ID), 4), c(0.0758, 0.0736, 0.0716, 0.0716), tolerance=0.005)
})


