library(ccar3)

test_that("cca_group_rrr_cv returns correct output structure", {
  set.seed(123)
  r = 2
  gen = generate_example_group(n=10000, p1=100, p2=5, 
                               r_pca = 3,
                               nnzeros = 5,
                               theta = diag(c(0.9,  0.8)),
                               lambda_pca = 1,
                               r = r, overlapping_amount = 1,
                               normalize_diagonal = TRUE,
                               n_new = 5000)
  
  X = gen$X
  Y = gen$Y
  groups = gen$groups

  result <- cca_group_rrr_cv(X, Y, groups=groups, r = 2, kfolds=3, parallelize = FALSE)
  
  expect_type(result, "list")
  expect_true(all(c("U", "V", "rmse", "cor") %in% names(result)))
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  expect_equal(dim(result$U)[1], 100)
  expect_equal(dim(result$V)[1], 5)
})



test_that("cca_group_rrr_cv can run in parallel", {
  r = 2
  gen = generate_example_group(n=500, p1=100, p2=5, 
                               r_pca = 3,
                               nnzeros = 5,
                               theta = diag(c(0.9,  0.8)),
                               lambda_pca = 1,
                               r = r, overlapping_amount = 1,
                               normalize_diagonal = TRUE,
                               n_new = 5000)
  
  X = gen$X
  Y = gen$Y
  groups = gen$groups
  result <- cca_group_rrr_cv(X, Y, groups=groups, r = 2, kfolds=3,  parallelize = TRUE)
  
  expect_type(result, "list")
  expect_true(subdistance(result$U, gen$u) < 0.35)  
  expect_true(subdistance(result$V, gen$v) < 0.35)
})


test_that("cca_group_rrr_cv returns the correct answer", {

  ##### Generate toy example data
  set.seed(123)
  r = 2
  gen = generate_example_group(n=500, p1=100, p2=5, 
                               r_pca = 3,
                               nnzeros = 5,
                               theta = diag(c(0.9,  0.8)),
                               lambda_pca = 1,
                               r = r, overlapping_amount = 1,
                               normalize_diagonal = TRUE,
                               n_new = 5000)
  
  X = gen$X
  Y = gen$Y
  groups = gen$groups
  result <- cca_rrr_cv(X, Y, r = r,
                      solver = "ADMM", param_lambda=10^seq(-5, 1.5, length.out = 100),
                      kfolds=3, Kx=NULL, r, parallelize = TRUE,
                      LW_Sy = TRUE)
  
  expect_type(result, "list")
  expect_true(subdistance(result$U, gen$u) < 0.35)  
  expect_true(subdistance(result$V, gen$v) < 0.35)
})


