library(ccar3)

test_that("cca_graph returns correct output structure", {
  set.seed(123)
  r = 2
  gen = generate_example_graph(n=500, p1=10, p2=5, 
                               r_pca = 3,
                               nnzeros = 5,
                               theta = diag(c(0.9,  0.8)),
                               lambda_pca = 1,
                               r = r, overlapping_amount = 1,
                               normalize_diagonal = TRUE,
                               n_new = 5000)
  
  X = gen$X
  Y = gen$Y
  Gamma = gen$Gamma
  result <- cca_graph_rrr(X, Y, Gamma = Gamma, r = r)
  
  expect_type(result, "list")
  expect_true(all(c("U", "V", "loss", "cor") %in% names(result)))
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  expect_equal(dim(result$U)[1], dim(X)[2])
  expect_equal(dim(result$V)[1], dim(Y)[2])
})


test_that("cca_graph returns the correct answer", {

  ##### Generate toy example data
  set.seed(123)
  r = 2
  gen = generate_example_graph(n=500, p1=10, p2=5, 
                               r_pca = 3,
                               nnzeros = 5,
                               theta = diag(c(0.9,  0.8)),
                               lambda_pca = 1,
                               r = r, overlapping_amount = 1,
                               normalize_diagonal = TRUE,
                               n_new = 5000)
  
  X = gen$X
  Y = gen$Y
  Gamma = gen$Gamma
  result <- cca_graph_rrr(X, Y, Gamma = Gamma, lambda = 0.025, r = r)
  expect_type(result, "list")
  expect_true(subdistance(result$U, gen$u) < 0.5)  
  expect_true(subdistance(result$V, gen$v) < 0.5)
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  
})

