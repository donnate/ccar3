library(ccar3)
test_that("SAR works", {
  set.seed(123)
  r = 3
  gen = generate_example_sparse_U(n=1000, p1=100, p2=10,
                                  r_pca = 5,
                                  nnzeros = 10,
                                  theta = diag(seq(0.9, 0.85, length.out = r)),
                                  lambda_pca = 1,
                                  r = r,  overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000) 
  X = gen$X
  Y = gen$Y
  
  result <- sparse_CCA_benchmarks(X,
                              Y, S=NULL, 
                              rank=r, kfolds=5, method.type = "FIT_SAR_CV",
                              lambdax = 10^seq(from=-3,to=,length=10),
                              lambday = 10^seq(from=-8,to=-1,length=8))
  
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  expect_equal(dim(result$U)[1], dim(X)[2])
  expect_equal(dim(result$V)[1], dim(Y)[2])
  
  expect_true(subdistance(result$U, gen$u) <0.2 )
  
})


test_that("SCCA_Parkhomenko works", {
  set.seed(123)
  r = 3
  gen = generate_example_sparse_U(n=1000, p1=100, p2=10,
                                  r_pca = 5,
                                  nnzeros = 10,
                                  theta = diag(seq(0.9, 0.85, length.out = r)),
                                  lambda_pca = 1,
                                  r = r,  overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000) 
  X = gen$X
  Y = gen$Y
  
  result <- sparse_CCA_benchmarks(X,
                                  Y, S=NULL, 
                                  rank=r, kfolds=5, method.type = "SCCA_Parkhomenko",
                                  lambdax = 10^seq(from=-3,to=,length=10),
                                  lambday = 10^seq(from=-8,to=-1,length=8))
  
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  expect_equal(dim(result$U)[1], dim(X)[2])
  expect_equal(dim(result$V)[1], dim(Y)[2])
  
  
})


test_that("Waaijenborg-CV works", {
  set.seed(123)
  r = 3
  gen = generate_example_sparse_U(n=1000, p1=100, p2=10,
                                  r_pca = 5,
                                  nnzeros = 10,
                                  theta = diag(seq(0.9, 0.85, length.out = r)),
                                  lambda_pca = 1,
                                  r = r,  overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000) 
  X = gen$X
  Y = gen$Y
  
  result <- sparse_CCA_benchmarks(X,
                                  Y, S=NULL, 
                                  rank=r, kfolds=5, method.type = "Waaijenborg-CV",
                                  lambdax = 10^seq(from=-3,to=,length=10),
                                  lambday = 10^seq(from=-8,to=-1,length=8))
  
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  expect_equal(dim(result$U)[1], dim(X)[2])
  expect_equal(dim(result$V)[1], dim(Y)[2])
  
  
})



test_that("Witten.CV works", {
  set.seed(123)
  r = 3
  gen = generate_example_sparse_U(n=1000, p1=100, p2=10,
                                  r_pca = 5,
                                  nnzeros = 10,
                                  theta = diag(seq(0.9, 0.85, length.out = r)),
                                  lambda_pca = 1,
                                  r = r,  overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000) 
  X = gen$X
  Y = gen$Y
  
  result <- sparse_CCA_benchmarks(X,
                                  Y, S=NULL, 
                                  rank=r, kfolds=5, method.type = "Witten.CV",
                                  lambdax = 10^seq(from=-3,to=,length=10),
                                  lambday = 10^seq(from=-8,to=-1,length=8))
  
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  expect_equal(dim(result$U)[1], dim(X)[2])
  expect_equal(dim(result$V)[1], dim(Y)[2])
  
  
})


test_that("Witten_Perm works", {
  set.seed(123)
  r = 3
  gen = generate_example_sparse_U(n=1000, p1=100, p2=10,
                                  r_pca = 5,
                                  nnzeros = 10,
                                  theta = diag(seq(0.9, 0.85, length.out = r)),
                                  lambda_pca = 1,
                                  r = r,  overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000) 
  X = gen$X
  Y = gen$Y
  
  result <- sparse_CCA_benchmarks(X,
                                  Y, S=NULL, 
                                  rank=r, kfolds=5, method.type = "Witten_Perm",
                                  lambdax = 10^seq(from=-3,to=0,length=10),
                                  lambday = 10^seq(from=-8,to=0,length=10))
  
  expect_equal(dim(result$U)[2], r)
  expect_equal(dim(result$V)[2], r)
  expect_equal(dim(result$U)[1], dim(X)[2])
  expect_equal(dim(result$V)[1], dim(Y)[2])
  
  
})

