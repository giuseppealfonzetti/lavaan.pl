#### test compute_frequencies() ####
# compute_frequencies() output
local({
  set.seed(123)
  p <- 10; q <- 2; n <- 10
  y <- matrix(as.numeric(rbinom(n*p, 1, .5)), n, p)
  f <- compute_frequencies(Y = y, C_VEC = rep(2, p))
  expect_equal(ncol(f), 4*p*(p-1)/2)
  expect_equal(nrow(f), 5)
  expect_identical(sum(f[5,]<0), 0L)
  expect_identical(sum(f[5,]>n), 0L)
  expect_identical(sum(!is.finite(f)), 0L)
})

# compute_frequencies() errors
local({
  set.seed(123)
  p <- 10; q <- 2; n <- 10
  y <- matrix(as.numeric(rbinom(n*p, 1, .5)), n, p)
  y1 <- y; y1[3,5] <- NaN
  expect_error(compute_frequencies(Y = y,  C_VEC = c(rep(2, p),NA)), 'C_VEC is not a numeric vector.')
  expect_error(compute_frequencies(Y = y,  C_VEC = matrix(2, 3, 4)), 'C_VEC is not a numeric vector.')
  expect_error(compute_frequencies(Y = y,  C_VEC = rep(2, p-1)), 'Y and C_VEC dimensions do not match.')
})

#### test build_constrMat() ####
# build_constrMat() output
local({
  p <- 10; q <- 2
  A <- is.na(build_constrMat(P = p, Q = q, STRUCT = 'triangular'))
  expect_equal(sum(A), (p*q)-(q*(q-1)/2))
  expect_equal(sum(A[lower.tri(A, diag = T)]), (p*q)-(q*(q-1)/2))
  expect_equal(sum(A[upper.tri(A, diag = F)]), 0)
  expect_equal(sum(A%in%c(0,1)), p*q)
  A <- is.na(build_constrMat(P = p, Q = q, STRUCT = 'simple'))
  expect_equal(sum(A), p)
  expect_equal(sum(rowSums(A)==1), p)
})

#### test gen_loadings() ####

# gen_loadings() errors
local({
  p <- 10; q <- 2
  A <- build_constrMat(P = p, Q = q, STRUCT = 'triangular')
  expect_error(gen_loadings(CONSTRMAT = 2))
  expect_error(gen_loadings(CONSTRMAT = NA))
})
