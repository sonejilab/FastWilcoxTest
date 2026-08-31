test_that("historical FC direction is preserved", {
  expect_equal(FC(2, 8), 4)
  expect_equal(logFC(2, 8), log(4))
  expect_equal(FC(c(8, 8), c(2, 2)), 0.25)
})

test_that("vector Wilcoxon agrees with stats normal approximation", {
  x <- c(0, 0, 1, 2, 4, 4)
  y <- c(0, 1, 1, 1, 2, 3, 6)

  got <- cppWilcoxTest(x, y, type = 2)
  expected <- stats::wilcox.test(x, y, exact = FALSE, correct = TRUE)

  expect_equal(unname(got[["p.value"]]), unname(expected$p.value), tolerance = 1e-10)
})

test_that("StatTest works directly on sparse matrices", {
  dense <- cbind(
    gene_a = c(8, 9, 10, 0, 0, 0),
    gene_b = c(0, 0, 0, 8, 9, 10),
    gene_c = c(1, 0, 1, 1, 0, 1)
  )
  X <- Matrix::Matrix(dense, sparse = TRUE)

  got <- StatTest(
    X,
    interest = 1:3,
    background = 4:6,
    logFCcut = 0,
    minPct = 0
  )

  expect_true(is.matrix(got))
  expect_equal(colnames(got), c(
    "colID", "logFC", "fracExprIN", "fracExprOUT", "rank.sum", "p.value", "FC"
  ))
  expect_true(all(got[, "colID"] %in% seq_len(ncol(X))))
})

test_that("StatTest p-values agree feature-wise with stats", {
  dense <- cbind(
    a = c(0, 0, 1, 2, 4, 4, 0, 1, 1, 1, 2, 3, 6),
    b = c(2, 3, 3, 4, 5, 8, 1, 1, 2, 2, 2, 3, 3)
  )
  X <- Matrix::Matrix(dense, sparse = TRUE)
  interest <- 1:6
  background <- 7:13

  got <- StatTest(X, interest, background, logFCcut = 0, minPct = 0)

  for (row in seq_len(nrow(got))) {
    j <- got[row, "colID"]
    expected <- stats::wilcox.test(
      dense[interest, j],
      dense[background, j],
      exact = FALSE,
      correct = TRUE
    )$p.value
    expect_equal(unname(got[row, "p.value"]), expected, tolerance = 1e-10)
  }
})

test_that("invalid identities are rejected", {
  X <- Matrix::Matrix(matrix(1:12, nrow = 4), sparse = TRUE)
  expect_error(StatTest(X, integer(), 2:4), "interest")
  expect_error(StatTest(X, 1:2, integer()), "background")
  expect_error(StatTest(X, c(1, 2), c(2, 3)), "overlap")
  expect_error(StatTest(X, 0L, 2:4), "outside")
})
