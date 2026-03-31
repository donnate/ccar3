library(ccar3)

test_that("cca_group returns correct output structure", {
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

  result <- cca_group_rrr(X, Y, lambda=0, groups=groups, r = r)
  expect_type(result, "list")
  expect_true(all(c("U", "V", "loss", "cor") %in% names(result)))
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  expect_equal(dim(result$U)[1], 100)
  expect_equal(dim(result$V)[1], 5)
})


test_that("cca_group returns the correct answer", {

  ##### Generate toy example data
  set.seed(123)
  r = 3
  gen = generate_example_group(n=10000, p1=100, p2=5, 
                                  r_pca = 3,
                                  nnzeros = 5,
                                  theta = diag(c(0.9,  0.85,  0.8)),
                                 lambda_pca = 1,
                                  r = r, overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                 n_new = 5000)
  X = gen$X
  Y = gen$Y
  result <- cca_group_rrr(X, Y, r = r, group=gen$groups, Sx=NULL,  Sy=NULL, lambda=0.01,
                    LW_Sy = TRUE)
  
  expect_type(result, "list")
  expect_true(subdistance(result$U, gen$u) < 0.15)  
  expect_true(subdistance(result$V, gen$v) < 0.15)
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  

})





test_that("cca_group computes the same solutions across solvers", {
  skip_if_not_installed("CVXR")
  set.seed(123)
  r = 2
  gen = generate_example_group(n=500, p1=100, p2=5, 
                               r_pca = 3,
                               nnzeros = 5,
                               theta = diag(c(0.9,  0.85,  0.8)[1:r]),
                               lambda_pca = 1,
                               r = r, overlapping_amount = 1,
                               normalize_diagonal = TRUE,
                               n_new = 5000)
  X = gen$X
  Y = gen$Y                          
  result1 <- cca_group_rrr(X, Y, r = r, groups=gen$groups, Sx=NULL,  Sy=NULL, 
                          lambda=0.1, LW_Sy = TRUE, solver="ADMM")
  result2 <- cca_group_rrr(X, Y, r = r, 
                          groups=gen$groups, 
                          Sx=NULL,  Sy=NULL, lambda=0.1,
                          LW_Sy = TRUE, solver="CVXR")
  
  expect_true(subdistance(result1$U, gen$u) <0.35)
  expect_true(subdistance(result2$U, gen$u) <0.35)
  
  expect_true(subdistance(result1$V, gen$v) <0.35)
  expect_true(subdistance(result2$V, gen$v) <0.35)
})

test_that("cca_group_rrr supports preprocess modes and p > n", {
  set.seed(42)
  n <- 40
  p <- 80
  q <- 5
  r <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  groups <- split(seq_len(p), ceiling(seq_len(p) / 5))

  fit_none <- cca_group_rrr(X, Y, groups = groups, r = r, lambda = 0.01,
                            preprocess = "none", LW_Sy = FALSE, niter = 300)
  fit_center <- cca_group_rrr(X, Y, groups = groups, r = r, lambda = 0.01,
                              preprocess = "center", LW_Sy = FALSE, niter = 300)
  fit_scale <- cca_group_rrr(X, Y, groups = groups, r = r, lambda = 0.01,
                             preprocess = "scale", LW_Sy = FALSE, niter = 300)

  expect_equal(dim(fit_none$U), c(p, r))
  expect_equal(dim(fit_center$U), c(p, r))
  expect_equal(dim(fit_scale$U), c(p, r))
})

test_that("cca_group_rrr supports matrix-free ADMM solves", {
  set.seed(7)
  n <- 50
  p <- 70
  q <- 4
  r <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  groups <- split(seq_len(p), ceiling(seq_len(p) / 5))

  fit <- cca_group_rrr(
    X, Y, groups = groups, r = r, lambda = 0.01,
    preprocess = "center", LW_Sy = FALSE, niter = 200,
    matrix_free_threshold = 1L,
    cg_tol = 1e-5,
    cg_maxiter = 50
  )

  expect_equal(dim(fit$U), c(p, r))
  expect_equal(dim(fit$V), c(q, r))
  expect_true(all(is.finite(fit$cor)))
})
