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
