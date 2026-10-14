test_that("FindMarkersSeurat agrees with StatTest on a Seurat object", {
  skip_if_not_installed("SeuratObject")

  counts <- Matrix::Matrix(
    matrix(
      c(
        10, 9, 8, 0, 0, 0,
         0, 0, 1, 8, 9, 10,
         4, 3, 5, 4, 3, 5
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(
        c("geneA", "geneB", "geneC"),
        paste0("cell", 1:6)
      )
    ),
    sparse = TRUE
  )

  object <- SeuratObject::CreateSeuratObject(counts = counts)
  SeuratObject::Idents(object) <- rep(c("A", "B"), each = 3)

  got <- FindMarkersSeurat(
    object,
    ident.1 = "A",
    assay = "RNA",
    layer = "counts",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE
  )

  direct <- StatTest(
    methods::as(Matrix::t(counts), "dgCMatrix"),
    interest = 1:3,
    background = 4:6,
    logFCcut = 0,
    minPct = 0,
    onlyPos = FALSE
  )
  direct <- as.data.frame(direct)

  expect_equal(
    got$p_val[match(rownames(counts)[direct$colID], got$gene)],
    unname(direct$p.value),
    tolerance = 1e-12
  )

  expect_true(got$avg_log2FC[got$gene == "geneA"] > 0)
  expect_true(got$avg_log2FC[got$gene == "geneB"] < 0)
})

test_that("FindAllMarkersSeurat agrees with repeated FindMarkersSeurat calls", {
  skip_if_not_installed("SeuratObject")

  counts <- Matrix::Matrix(
    matrix(
      c(
        8, 7, 0, 0, 0, 0,
        0, 0, 7, 8, 0, 0,
        0, 0, 0, 0, 8, 7,
        2, 2, 2, 2, 2, 2
      ),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(
        paste0("gene", 1:4),
        paste0("cell", 1:6)
      )
    ),
    sparse = TRUE
  )

  object <- SeuratObject::CreateSeuratObject(counts = counts)
  SeuratObject::Idents(object) <- rep(c("A", "B", "C"), each = 2)

  all <- FindAllMarkersSeurat(
    object,
    layer = "counts",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE
  )

  for (identity in c("A", "B", "C")) {
    one <- FindMarkersSeurat(
      object,
      ident.1 = identity,
      layer = "counts",
      logfc.threshold = 0,
      min.pct = 0,
      only.pos = FALSE
    )
    got <- all[all$cluster == identity, setdiff(names(all), "cluster"), drop = FALSE]
    rownames(got) <- NULL
    expect_equal(got, one)
  }
})


test_that("FindAllMarkersSeurat returns markers for all identities", {
  skip_if_not_installed("SeuratObject")

  counts <- Matrix::Matrix(
    matrix(
      c(
        10, 9, 8, 0, 0, 0,
         0, 0, 1, 8, 9, 10,
         4, 3, 5, 4, 3, 5
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(
        c("geneA", "geneB", "geneC"),
        paste0("cell", 1:6)
      )
    ),
    sparse = TRUE
  )

  object <- SeuratObject::CreateSeuratObject(counts = counts)
  SeuratObject::Idents(object) <- rep(c("A", "B"), each = 3)

  got <- FindAllMarkersSeurat(
    object,
    assay = "RNA",
    layer = "counts",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

  rownames(got) <- NULL

  expect_setequal(unique(got$cluster), c("A", "B"))

  expect_true(
    any(got$cluster == "A" & got$gene == "geneA")
  )

  expect_true(
    any(got$cluster == "B" & got$gene == "geneB")
  )

  expect_true(
    got$avg_log2FC[
      got$cluster == "A" & got$gene == "geneA"
    ] > 0
  )

  expect_true(
    got$avg_log2FC[
      got$cluster == "B" & got$gene == "geneB"
    ] > 0
  )
})