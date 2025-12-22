# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

test_that("infer_orientation_signs returns correct length", {
  set.seed(123)
  Y <- rnorm(20)
  obs <- as.data.frame(matrix(sample(c(0,1), 20*3, replace=TRUE), ncol=3))
  signs <- infer_orientation_signs(Y, obs)
  expect_length(signs, 3)
  expect_true(all(signs %in% c(1, -1)))
})

test_that("infer_orientation_signs correlates with Y", {
  set.seed(456)
  Y <- rnorm(30)
  obs <- as.data.frame(matrix(sample(c(0,1), 30*2, replace=TRUE), ncol=2))
  signs <- infer_orientation_signs(Y, obs)
  oriented <- sweep(obs, 2, signs, FUN = "*")
  cors <- sapply(oriented, function(x) cor(x, Y))
  expect_true(all(cors >= 0 | is.na(cors)))
})

test_that("infer_orientation_signs works with PC1 method", {
  set.seed(789)
  obs <- as.data.frame(matrix(sample(c(0,1), 50*4, replace=TRUE), ncol=4))
  signs <- infer_orientation_signs(observables = obs, method = "PC1")
  expect_length(signs, 4)
  expect_true(all(signs %in% c(1, -1)))
})

test_that("infer_orientation_signs PC1 method produces positive correlations with PC1", {
  set.seed(321)
  obs <- as.data.frame(matrix(sample(c(0,1), 100*5, replace=TRUE), ncol=5))
  signs <- infer_orientation_signs(observables = obs, method = "PC1")
  pc1 <- prcomp(obs, scale. = TRUE)$x[, 1]
  oriented <- sweep(obs, 2, signs, FUN = "*")
  cors <- sapply(oriented, function(x) cor(x, pc1))
  expect_true(all(cors >= 0 | is.na(cors)))
})

test_that("infer_orientation_signs errors when Y missing for Y method", {
  obs <- as.data.frame(matrix(sample(c(0,1), 20*3, replace=TRUE), ncol=3))
  expect_error(
    infer_orientation_signs(observables = obs, method = "Y"),
    "Y must be provided"
  )
})

test_that("infer_orientation_signs errors for non-binary observables", {
  Y <- rnorm(20)
  obs <- as.data.frame(matrix(rnorm(20*3), ncol=3))
  expect_error(
    infer_orientation_signs(Y, obs),
    "binary observables"
  )
})
