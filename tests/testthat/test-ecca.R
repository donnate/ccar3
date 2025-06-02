library(ccar3)

test_that("ecca returns correct output structure", {
  
  
  set.seed(123)
  gen = generate_example_none_trivial_pca(n=100, p1=50, p2=50, r_pca = 3, nnzeros = 5, 
                                           theta = diag(c(.9, .9)), overlapping_amount = 1,
                                           lambda_pca = 1, r = 2)
  result <- ecca(gen$X, gen$Y, r = 2)
  
  expect_type(result, "list")
  expect_true(all(c("U", "V", "loss", "cor") %in% names(result)))
  expect_equal(dim(result$U)[2], 2)
  expect_equal(dim(result$V)[2], 2)
  expect_equal(dim(result$U)[1], dim(gen$X)[2])
  expect_equal(dim(result$V)[1], dim(gen$Y)[2])
})


test_that("ecca returns the correct answer", {
  
  ##### Generate toy example data
  set.seed(123)
  r = 2
  gen = generate_example_none_trivial_pca(n=700, p1=50, p2=50, r_pca = 3, nnzeros = 5, 
                                          theta = diag(c(.9, .9)), overlapping_amount = 1,
                                          lambda_pca = 1, r = r)
  
  
  result <- ecca(gen$X, gen$Y, r = r, lambda = 0.01)
  
  expect_type(result, "list")
  expect_true(subdistance(result$U, gen$u) < 0.3)  
  expect_true(subdistance(result$V, gen$v) < 0.3)
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
})





test_that("ecca.cv returns the correct answer", {
  
  ##### Generate toy example data
  set.seed(123)
  r = 2
  gen = generate_example_none_trivial_pca(n=700, p1=50, p2=50, r_pca = 3, nnzeros = 5, 
                                          theta = diag(c(.9, .9)), overlapping_amount = 1,
                                          lambda_pca = 1, r = r)
  
  
  result <- ecca.cv(gen$X, gen$Y, r = r, lambdas = c(10^(-4), 10^(-3), 0.01, 0.05, 0.075, 0.1, 0.5),, parallel=F)
  
  expect_type(result, "list")
  expect_true(subdistance(result$U, gen$u) < 0.3)  
  expect_true(subdistance(result$V, gen$v) < 0.3)
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
})
  


test_that("ecca.cv parallelization works", {
  
  ##### Generate toy example data
  set.seed(123)
  r = 2
  gen = generate_example_none_trivial_pca(n=700, p1=50, p2=50, r_pca = 3, nnzeros = 5, 
                                          theta = diag(c(.9, .9)), overlapping_amount = 1,
                                          lambda_pca = 1, r = r)
  
  
  result <- ecca.cv(gen$X, gen$Y, r = r, lambdas = c(10^(-4), 10^(-3), 0.01, 0.05, 0.075, 0.1, 0.5), parallel=T)
  
  expect_type(result, "list")
  expect_true(subdistance(result$U, gen$u) < 0.3)  
  expect_true(subdistance(result$V, gen$v) < 0.3)
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
})


