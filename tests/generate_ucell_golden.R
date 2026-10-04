#!/usr/bin/env Rscript
# Golden-value generator for the h5ad_ucell_score plugin node.
#
# Scores every fixture produced by tests/make_ucell_fixtures.py with the
# OFFICIAL UCell R package (ScoreSignatures_UCell) and writes the expected
# per-cell scores to a JSON file the pytest suite compares against. The
# plugin implementation must never be used to produce or adjust these
# values -- that would be circular validation.
#
# Usage (from the plugin root):
#   R_LIBS=~/R/library Rscript tests/generate_ucell_golden.R \
#       [fixtures_dir] [output_json]
# Defaults: tests/golden_inputs  tests/ucell_golden.json

args <- commandArgs(trailingOnly = TRUE)
fixtures_dir <- if (length(args) >= 1) args[1] else "tests/golden_inputs"
output_json <- if (length(args) >= 2) args[2] else "tests/ucell_golden.json"

suppressPackageStartupMessages({
  library(UCell)
  library(jsonlite)
})

ucell_version <- as.character(packageVersion("UCell"))
cat("UCell version:", ucell_version, "\n")

spec <- jsonlite::fromJSON(file.path(fixtures_dir, "fixtures.json"),
                           simplifyVector = FALSE)

golden <- list()
for (fixture in spec) {
  params <- fixture$params
  max_rank <- if (is.null(params$max_rank)) 1500 else as.integer(params$max_rank)
  ties_method <- if (is.null(params$ties_method)) "average" else params$ties_method
  missing_genes <- if (is.null(params$missing_genes)) "impute" else params$missing_genes

  frame <- read.csv(file.path(fixtures_dir, fixture$csv),
                    row.names = 1, check.names = FALSE)
  # UCell expects genes x cells; read.csv gives cells x genes.
  matrix <- t(as.matrix(frame))
  features <- lapply(fixture$gene_sets, function(genes) unlist(genes))

  scores <- ScoreSignatures_UCell(
    matrix = matrix,
    features = features,
    maxRank = max_rank,
    ties.method = ties_method,
    missing_genes = missing_genes
  )
  stopifnot(identical(rownames(scores), unlist(fixture$cells)))

  per_signature <- list()
  signature_names <- names(fixture$gene_sets)
  for (name in signature_names) {
    column <- paste0(name, "_UCell")
    # I() keeps length-1 score vectors encoded as JSON arrays.
    per_signature[[name]] <- I(as.numeric(scores[, column]))
  }
  golden[[fixture$name]] <- list(
    params = list(max_rank = max_rank, ties_method = ties_method,
                  missing_genes = missing_genes),
    cells = I(unlist(fixture$cells)),
    scores = per_signature
  )
  cat(sprintf("scored %-24s (%d cells x %d genes, %d signatures)\n",
              fixture$name, ncol(matrix), nrow(matrix), length(features)))
}

payload <- list(
  `_provenance` = list(
    package = "UCell",
    package_version = ucell_version,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    entrypoint = "ScoreSignatures_UCell",
    generated = format(Sys.time(), tz = "UTC", usetz = TRUE),
    note = paste("Expected values produced solely by the official UCell R",
                 "package; regenerate with tests/generate_ucell_golden.R",
                 "after changing fixtures, never from plugin output.")
  ),
  fixtures = golden
)

jsonlite::write_json(payload, output_json, auto_unbox = TRUE, digits = 17)
cat("wrote", output_json, "\n")
