library(bacon)
context("bacon estimates")

test_that("bacon estimates reproducible", {
              set.seed(12345)
              n <- 1000
              tol <- 0.05
              y <- rnormmix(n, c(0.9, 0, 1, 0, 4, 1))
              bc <- bacon(y)
              est <- estimates(bc)
              print(est[c(1,4,7)])
              expect_equal(est[1], 0.9, tolerance=tol)
              expect_equal(est[4], 0.0,  tolerance = tol)
              expect_equal(est[7], 1.0,tolerance = tol)
          })
