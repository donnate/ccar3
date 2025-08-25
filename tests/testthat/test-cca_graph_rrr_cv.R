library(ccar3)

test_that("cca_graph_rrr_cv returns correct output structure", {
  set.seed(123)
  skip_if_not_installed("igraph")

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
  
  result <- cca_graph_rrr_cv(X, Y, Gamma=Gamma, r = r, kfolds=3, parallelize = FALSE)
  
  expect_type(result, "list")
  expect_true(all(c("U", "V", "rmse", "cor") %in% names(result)))
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  expect_equal(dim(result$U)[1], dim(X)[2])
  expect_equal(dim(result$V)[1], dim(Y)[2])
})



test_that("cca_graph_rrr_cv can run in parallel", {
  skip_if_not_installed("igraph")

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
  result <- cca_graph_rrr_cv(X, Y, Gamma=Gamma, r = 2, kfolds=3,  parallelize = TRUE)
  
  expect_type(result, "list")
  expect_true(subdistance(result$U, gen$u) < 0.7)  
  expect_true(subdistance(result$V, gen$v) < 0.7)
})


test_that("cca_graph_rrr_cv returns the correct answer", {
  skip_if_not_installed("igraph")

  
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
  
  result <- cca_graph_rrr_cv(X, Y, Gamma=Gamma, r = r, kfolds=3, parallelize =TRUE)
  
  expect_type(result, "list")
  expect_true(subdistance(result$U, gen$u) < 0.5)  
  expect_true(subdistance(result$V, gen$v) < 0.5)
})
