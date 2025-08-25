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

  result <- cca_rrr(X, Y, r = r, Sx=NULL,  Sy=NULL, lambda=0, r, highdim=TRUE,
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