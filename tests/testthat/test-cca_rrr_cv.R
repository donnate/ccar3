library(ccar3)

test_that("cca_rrr_cv returns correct output structure", {
  set.seed(123)
  r = 2
  gen = generate_example_sparse_U(n=100, p1=50, p2=10,
                                  r_pca = 5,
                                  nnzeros = 10,
                                  theta = diag(seq(0.9, 0.85, length.out = r)),
                                  lambda_pca = 1,
                                  r = r,  overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000) 
  X = gen$X
  Y = gen$Y  

  result <- cca_rrr_cv(X, Y, r = r, kfolds=3, parallelize = FALSE)
  
  expect_type(result, "list")
  expect_true(all(c("U", "V", "rmse", "lambda", "cor") %in% names(result)))
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  expect_equal(dim(result$U)[1], 50)
  expect_equal(dim(result$V)[1], 10)
})

test_that("cca_rrr_cv can run in parallel", {
    skip_if_not_installed("doParallel")
  set.seed(123)
  r = 3
  gen = generate_example_sparse_U(n=100, p1=50, p2=10,
                                  r_pca = 5,
                                  nnzeros = 10,
                                  theta = diag(seq(0.9, 0.85, length.out = r)),
                                  lambda_pca = 1,
                                  r = r,  overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000) 
  X = gen$X
  Y = gen$Y  

  result <- cca_rrr_cv(X, Y, r = r, kfolds=3, parallelize = TRUE)
  
  expect_type(result, "list")
  expect_true(all(c("U", "V", "rmse", "lambda", "cor") %in% names(result)))
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  expect_equal(dim(result$U)[1], 50)
  expect_equal(dim(result$V)[1], 10)
})


test_that("cca_rrr_cv returns the correct answer", {

  ##### Generate toy example data
  set.seed(123)
  r = 2
  gen = generate_example_sparse_U(n=5000, p1=30, p2=5,
                                  r_pca = 5,
                                  nnzeros = 10,
                                  theta = diag(seq(0.9, 0.85, length.out = r)),
                                  lambda_pca = 1,
                                  r = r,  overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000) 
  X = gen$X
  Y = gen$Y
  result <- cca_rrr_cv(X, Y, r = r,
                      solver = "ADMM", lambdas=10^seq(-5, 1.5, length.out = 10),
                      kfolds=3, parallelize = TRUE,
                      LW_Sy = FALSE)
  
  expect_type(result, "list")
  expect_true(subdistance(result$U, gen$u) < 0.2)  
  expect_true(subdistance(result$V, gen$v) < 0.2)
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
})




