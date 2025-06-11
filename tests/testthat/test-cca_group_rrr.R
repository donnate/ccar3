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