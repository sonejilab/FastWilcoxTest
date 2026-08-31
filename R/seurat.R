#' Fast marker testing for a Seurat object
#'
#' Thin Seurat interface around [StatTest()]. Expression is extracted from a
#' Seurat assay layer, kept sparse, transposed to cells-by-features, and passed
#' directly to the FastWilcoxTest C++ implementation.
#'
#' @param object A Seurat object.
#' @param ident.1 Identity or vector of identities defining the population of
#'   interest.
#' @param ident.2 Optional identity or vector of identities defining the
#'   background. If `NULL`, every cell not in `ident.1` is used.
#' @param group.by Optional metadata column defining the identities. If `NULL`,
#'   `SeuratObject::Idents(object)` is used.
#' @param assay Assay containing the expression matrix. Default is `"RNA"`.
#' @param layer Assay layer to test. Default is `"data"`.
#' @param logfc.threshold Minimum absolute log2 fold change required before a
#'   feature is tested. Default is `0.25`.
#' @param min.pct Minimum fraction of cells expressing a feature in either
#'   population. Default is `0.1`.
#' @param only.pos If `TRUE`, return only features with higher mean expression
#'   in `ident.1` than in the background.
#' @param p.adjust.method Method passed to [stats::p.adjust()].
#'
#' @return A data frame with feature names, raw and adjusted p-values, log2 fold
#'   change, detection fractions, rank sum, and fold change. Positive
#'   `avg_log2FC` means higher expression in `ident.1`.
#'
#' @details
#' The historical FastWilcoxTest C++ API reports fold change as
#' `mean(background) / mean(interest)`. This wrapper reverses that convention
#' on output so that its direction matches the usual Seurat interpretation.
#'
#' The wrapper deliberately depends only on `SeuratObject`, not the full Seurat
#' package. The requested assay layer must exist as a single matrix layer.
#'
#' @examples
#' \dontrun{
#' markers <- FindMarkersSeurat(
#'   object,
#'   ident.1 = "0",
#'   only.pos = TRUE
#' )
#'
#' markers_0_vs_1 <- FindMarkersSeurat(
#'   object,
#'   ident.1 = "0",
#'   ident.2 = "1"
#' )
#' }
#'
#' @export
FindMarkersSeurat <- function(
    object,
    ident.1,
    ident.2 = NULL,
    group.by = NULL,
    assay = "RNA",
    layer = "data",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = FALSE,
    p.adjust.method = "bonferroni") {

  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop(
      "Package 'SeuratObject' is required for FindMarkersSeurat()",
      call. = FALSE
    )
  }

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object", call. = FALSE)
  }

  if (length(assay) != 1L || !is.character(assay) || !nzchar(assay)) {
    stop("assay must be one non-empty character string", call. = FALSE)
  }
  if (length(layer) != 1L || !is.character(layer) || !nzchar(layer)) {
    stop("layer must be one non-empty character string", call. = FALSE)
  }
  if (length(logfc.threshold) != 1L || !is.finite(logfc.threshold) ||
      logfc.threshold < 0) {
    stop("logfc.threshold must be one finite non-negative number", call. = FALSE)
  }
  if (length(min.pct) != 1L || !is.finite(min.pct) ||
      min.pct < 0 || min.pct > 1) {
    stop("min.pct must be between 0 and 1", call. = FALSE)
  }

  metadata <- object[[]]

  if (is.null(group.by)) {
    groups <- as.character(SeuratObject::Idents(object))
  } else {
    if (length(group.by) != 1L || !is.character(group.by) || !nzchar(group.by)) {
      stop("group.by must be NULL or one metadata column name", call. = FALSE)
    }
    if (!group.by %in% colnames(metadata)) {
      stop("group.by metadata column not found: ", group.by, call. = FALSE)
    }
    groups <- as.character(metadata[[group.by]])
  }

  ident.1 <- as.character(ident.1)
  if (!length(ident.1) || anyNA(ident.1)) {
    stop("ident.1 must contain at least one non-NA identity", call. = FALSE)
  }

  cells.1 <- which(groups %in% ident.1)
  if (!length(cells.1)) {
    stop(
      "No cells found for ident.1: ",
      paste(ident.1, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(ident.2)) {
    cells.2 <- which(!(groups %in% ident.1))
  } else {
    ident.2 <- as.character(ident.2)
    if (!length(ident.2) || anyNA(ident.2)) {
      stop("ident.2 must contain at least one non-NA identity", call. = FALSE)
    }
    if (length(intersect(ident.1, ident.2))) {
      stop("ident.1 and ident.2 must not overlap", call. = FALSE)
    }
    cells.2 <- which(groups %in% ident.2)
  }

  if (!length(cells.2)) {
    stop("No cells found for the background population", call. = FALSE)
  }

  expression <- tryCatch(
    SeuratObject::LayerData(object, assay = assay, layer = layer),
    error = function(e) {
      stop(
        "Unable to read assay '", assay, "' layer '", layer, "': ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!inherits(expression, "sparseMatrix")) {
    expression <- Matrix::Matrix(expression, sparse = TRUE)
  }

  if (ncol(expression) != length(groups)) {
    stop(
      "Expression matrix has ", ncol(expression),
      " cells but the identity vector has ", length(groups),
      call. = FALSE
    )
  }
  if (is.null(rownames(expression))) {
    stop("Expression layer has no feature names", call. = FALSE)
  }

  test.matrix <- methods::as(Matrix::t(expression), "dgCMatrix")

  # StatTest uses natural-log thresholds internally.
  cpp.logfc.cut <- logfc.threshold * log(2)

  # Always disable the historical directional filter here because its positive
  # direction is background / interest. Direction is normalized below.
  raw <- StatTest(
    test.matrix,
    interest = cells.1,
    background = cells.2,
    logFCcut = cpp.logfc.cut,
    minPct = min.pct,
    onlyPos = FALSE
  )

  raw <- as.data.frame(raw)

  if (!nrow(raw)) {
    return(data.frame(
      gene = character(),
      p_val = numeric(),
      p_val_adj = numeric(),
      avg_log2FC = numeric(),
      pct.1 = numeric(),
      pct.2 = numeric(),
      rank_sum = numeric(),
      FC = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  feature.id <- as.integer(raw$colID)

  result <- data.frame(
    gene = rownames(expression)[feature.id],
    p_val = raw$p.value,
    avg_log2FC = -raw$logFC / log(2),
    pct.1 = raw$fracExprIN,
    pct.2 = raw$fracExprOUT,
    rank_sum = raw$rank.sum,
    FC = 1 / raw$FC,
    stringsAsFactors = FALSE
  )

  result$p_val_adj <- stats::p.adjust(
    result$p_val,
    method = p.adjust.method,
    n = nrow(expression)
  )

  if (isTRUE(only.pos)) {
    result <- result[result$avg_log2FC > 0, , drop = FALSE]
  }

  result <- result[
    order(result$p_val, -abs(result$avg_log2FC)),
    ,
    drop = FALSE
  ]
  rownames(result) <- NULL

  result[
    ,
    c(
      "gene",
      "p_val",
      "p_val_adj",
      "avg_log2FC",
      "pct.1",
      "pct.2",
      "rank_sum",
      "FC"
    )
  ]
}
