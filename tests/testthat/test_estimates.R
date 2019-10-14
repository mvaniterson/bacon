library(bacon)
context("bacon estimates")

test_that("bacon estimates reproducible", {
  set.seed(12345)
  n <- 1000
  tol <- 0.05
  y <- rnormmix(n, c(0.9, 0, 1, 0, 4, 1))
  bc <- bacon(y)
  est <- estimates(bc)
  print(est[c(1, 4, 7)])
  expect_equal(est[1], 0.9, tolerance = tol)
  expect_equal(est[4], 0.0,  tolerance = tol)
  expect_equal(est[7], 1.0, tolerance = tol)
})

test_that("Handling of NA's in test-statistics", {
  set.seed(12345)
  n <- 1000
  tol <- 0.05
  y <- rnormmix(n, c(0.9, 0, 1, 0, 4, 1))
  y[1] <- NA
  bc <- bacon(y, na.exclude = TRUE)
  est <- estimates(bc)
  expect_equal(est[1], 0.91, tolerance = tol)
  expect_equal(est[4], 0.001,  tolerance = tol)
  expect_equal(est[7], 0.96, tolerance = tol)
  expect_equal(length(pval(bc)), length(y))
})



test_that("Handling of NA's in meta-analysis", {
  set.seed(12345)
  biases <- runif(6,-0.2, 0.2)
  inflations <- runif(6, 1, 1.3)
  es <- matrix(nrow = 5000, ncol = 6)
  for (i in 1:6)
    es[, i] <-
    rnormmix(5000, c(0.9, biases[i], inflations[i], 0, 4, 1), shuffle = FALSE)
  se <- replicate(6, 0.8 * sqrt(4 / rchisq(5000, df = 4)))
  colnames(es) <- colnames(se) <- LETTERS[1:ncol(se)]
  rownames(es) <- rownames(se) <- 1:5000
  head(rownames(es))
  
  idx <- sample(nrow(se), size = 10)
  
  es[idx, ] <- se[idx, ] <- NA
  
  bc <- bacon(NULL, es, se, na.exclude = TRUE)
  
  bcm <- meta(bc)
  expect_equal(dim(pval(bcm)), c(5000, 6+1))
})
