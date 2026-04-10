library(ccar3)

test_that("cca_rrr returns correct output structure", {
  X <- matrix(rnorm(100*10), 100, 10)
  Y <- matrix(rnorm(100*5), 100, 5)
  result <- cca_rrr(X, Y, r = 2)
  
  expect_type(result, "list")
  expect_true(all(c("U", "V", "loss", "cor") %in% names(result)))
  expect_equal(dim(result$U)[2], 2)
  expect_equal(dim(result$V)[2], 2)
  expect_equal(dim(result$U)[1], 10)
  expect_equal(dim(result$V)[1], 5)
})


test_that("cca_rrr returns the correct answer", {

  ##### Generate toy example data
  set.seed(123)
  r = 3
  gen = generate_example_sparse_U(n=10000, p1=100, p2=10,
                                  r_pca = 5,
                                  nnzeros = 10,
                                  theta = diag(seq(0.9, 0.85, length.out = r)),
                                  lambda_pca = 1,
                                  r = r,  overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000) 
  X = gen$X
  Y = gen$Y

  result <- cca_rrr(X, Y, r = r, Sx=NULL,  Sy=NULL, lambda=0, highdim=TRUE,
                    LW_Sy = TRUE, solver = "ADMM")
  
  expect_type(result, "list")
  expect_true(subdistance(result$U, gen$u) < 0.15)  
  expect_true(subdistance(result$V, gen$v) < 0.15)
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
})





test_that("cca_rrr computes the same solutions across solvers", {
  set.seed(123)
  skip_if_not_installed("CVXR")


  r = 3
  gen = generate_example_sparse_U(n=10000, p1=100, p2=10,
                                  r_pca = 5,
                                  nnzeros = 10,
                                  theta = diag(seq(0.9, 0.85, length.out = r)),
                                  lambda_pca = 1,
                                  r = r,  overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000) 
  X = gen$X
  Y = gen$Y  

  result2 <- cca_rrr(X, Y, r = r, Sx=NULL,  Sy=NULL, lambda=0.1, rho=10,
                     highdim=TRUE,
                     LW_Sy = TRUE, solver="ADMM")
  print("Done ADMM")
  result3 <- cca_rrr(X, Y, r = r, Sx=NULL,  Sy=NULL, lambda=0.1,
                    highdim=TRUE,
                    LW_Sy = TRUE, solver="CVX")
  print("Done CVXR")
  
  expect_true(subdistance(result2$U, gen$u) <0.2 )
  expect_true(subdistance(result3$U, gen$u) <0.2 )
})


test_that("cca_rrr exposes matrix-free ADMM tuning arguments", {
  set.seed(42)
  X <- matrix(rnorm(80 * 12), 80, 12)
  Y <- matrix(rnorm(80 * 4), 80, 4)

  result <- cca_rrr(
    X, Y, r = 2,
    highdim = TRUE,
    solver = "ADMM",
    matrix_free_threshold = 1L,
    cg_tol = 1e-5,
    cg_maxiter = 50
  )

  expect_type(result, "list")
  expect_equal(dim(result$U), c(12, 2))
  expect_equal(dim(result$V), c(4, 2))
})

test_that("cca_rrr can skip internal preprocessing", {
  set.seed(21)
  X <- scale(matrix(rnorm(90 * 10), 90, 10), center = TRUE, scale = FALSE)
  Y <- scale(matrix(rnorm(90 * 4), 90, 4), center = TRUE, scale = FALSE)

  result <- cca_rrr(
    X, Y, r = 2,
    highdim = TRUE,
    solver = "ADMM",
    standardize = FALSE,
    preprocess = FALSE
  )

  expect_type(result, "list")
  expect_equal(dim(result$U), c(10, 2))
  expect_equal(dim(result$V), c(4, 2))
})

test_that("cca_rrr supports preprocess modes and logical backward compatibility", {
  set.seed(33)
  X <- matrix(rnorm(120 * 8), 120, 8)
  Y <- matrix(rnorm(120 * 4), 120, 4)

  fit_center_mode <- cca_rrr(
    X, Y, r = 2,
    highdim = TRUE,
    solver = "ADMM",
    preprocess = "center",
    LW_Sy = FALSE,
    niter = 200
  )
  fit_center_logical <- cca_rrr(
    X, Y, r = 2,
    highdim = TRUE,
    solver = "ADMM",
    standardize = FALSE,
    preprocess = TRUE,
    LW_Sy = FALSE,
    niter = 200
  )

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  Y_centered <- scale(Y, center = TRUE, scale = FALSE)
  fit_none_char <- cca_rrr(
    X_centered, Y_centered, r = 2,
    highdim = TRUE,
    solver = "ADMM",
    preprocess = "none",
    LW_Sy = FALSE,
    niter = 200
  )
  fit_none_logical <- cca_rrr(
    X_centered, Y_centered, r = 2,
    highdim = TRUE,
    solver = "ADMM",
    standardize = FALSE,
    preprocess = FALSE,
    LW_Sy = FALSE,
    niter = 200
  )

  expect_equal(fit_center_mode$B_opt, fit_center_logical$B_opt, tolerance = 1e-8)
  expect_equal(fit_none_char$B_opt, fit_none_logical$B_opt, tolerance = 1e-8)
})
