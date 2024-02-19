#' @useDynLib RcppPileup, .registration = TRUE
#' @importFrom Rcpp sourceCpp

#' @export
sc_mutations <- function(bam_path, seqnames, positions, indel, barcodes, threads = 1) {
  stopifnot(
    "bam_path must be a single character string" =
      is.character(bam_path) && length(bam_path) == 1
  )
  stopifnot(
    "seqnames not the same length as positions" =
      length(seqnames) == length(positions)
  )
  stopifnot(
    "barcodes must be a character vector" =
      is.character(barcodes)
  )

  variants <- parallel::mcmapply(
    FUN = function(seqname, pos) {
      variant_count_matrix_cpp(
        bam_path = bam_path,
        seqname = seqname, pos = pos, indel = indel, barcodes = barcodes
      ) |>
        tibble::as_tibble(rownames = "allele") |>
        tidyr::pivot_longer(
          cols = -tidyselect::matches("allele"),
          values_to = "allele_count", names_to = "barcode"
        ) |>
        dplyr::group_by(barcode) |>
        dplyr::mutate(
          cell_total_reads = sum(allele_count),
          pct = allele_count / cell_total_reads,
          pos = pos, seqname = seqname
        ) |>
        dplyr::ungroup()
    },
    seqname = seqnames, pos = positions, SIMPLIFY = FALSE, mc.cores = threads
  ) |>
    dplyr::bind_rows()

  return(variants)
}
