# Skip all tests on CRAN to avoid timeouts
skip_on_cran()

test_that("lpmec_onerun basic functionality", {
  set.seed(123)
  Y <- rnorm(20)
  obs <- as.data.frame(matrix(sample(c(0,1), 20*4, replace=TRUE), ncol=4))
  res <- lpmec_onerun(Y, obs, estimation_method = "pca")
  expect_s3_class(res, "lpmec_onerun")
  expect_true(all(c("ols_coef", "x_est1", "x_est2") %in% names(res)))
})

test_that("lpmec_onerun S3 methods work", {
  set.seed(123)
  Y <- rnorm(20)
  obs <- as.data.frame(matrix(sample(c(0,1), 20*4, replace=TRUE), ncol=4))
  res <- lpmec_onerun(Y, obs, estimation_method = "pca")
  sum_df <- summary(res)
  expect_s3_class(sum_df, "data.frame")
  expect_equal(row.names(sum_df), c("OLS", "IV", "Corrected IV", "Corrected OLS"))
  expect_type(capture.output(print(res)), "character")
  expect_silent(plot(res))
})
