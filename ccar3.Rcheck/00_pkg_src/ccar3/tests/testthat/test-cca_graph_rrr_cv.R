library(ccar3)

expected_cv_names <- c(
  "U", "V", "lambda", "rmse", "cor",
  "lambda_x", "lambda_x_se", "lambda_y", "lambda_y_se",
  "resultsx", "cv_summary", "cv_folds", "Lambda", "B", "fit"
)

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
  expect_identical(names(result), expected_cv_names)
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

test_that("cca_graph_rrr_cv accepts preprocessing mode and sparse Gamma", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("Matrix")

  set.seed(42)
  r <- 2
  gen <- generate_example_graph(
    n = 250, p1 = 10, p2 = 5,
    r_pca = 3,
    nnzeros = 5,
    theta = diag(c(0.9, 0.8)),
    lambda_pca = 1,
    r = r, overlapping_amount = 1,
    normalize_diagonal = TRUE,
    n_new = 1000
  )

  Gamma_sparse <- Matrix::Matrix(gen$Gamma, sparse = TRUE)
  result <- cca_graph_rrr_cv(
    gen$X, gen$Y, Gamma = Gamma_sparse, r = r,
    kfolds = 3, parallelize = FALSE, preprocess = "center"
  )

  expect_type(result, "list")
  expect_identical(names(result), expected_cv_names)
})

test_that("cca_graph_rrr_cv_folds honors preprocessing modes on raw inputs", {
  set.seed(88)
  X <- matrix(rnorm(60 * 8), 60, 8)
  Y <- matrix(rnorm(60 * 4), 60, 4)
  Gamma <- diag(ncol(X) - 1)
  Gamma <- cbind(Gamma, 0) - cbind(0, diag(ncol(X) - 1))
  folds <- ccar3:::.create_cv_folds(nrow(X), 3)

  fold_values <- ccar3:::cca_graph_rrr_cv_folds(
    X, Y, Gamma = Gamma,
    Sx = NULL, Sy = NULL,
    kfolds = 3, lambda = 0.01, r = 2,
    preprocess = "center",
    LW_Sy = FALSE,
    niter = 200,
    folds = folds,
    return_fold_values = TRUE
  )

  expect_length(fold_values, length(folds))
  expect_false(all(is.na(fold_values)))
})
