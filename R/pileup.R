#' @useDynLib RcppPileup, .registration = TRUE
#' @importFrom Rcpp sourceCpp

variant_count_tb <- function(bam_path, seqname, pos, indel, barcodes) {
  # allele by barcode matrix (value: read count)
  variant_count_matrix_cpp(
    bam_path = bam_path,
    seqname = seqname, pos = pos, indel = indel, barcodes = barcodes
  ) |>
    tibble::as_tibble(rownames = "allele") |>
    # pivot to long format: allele, barcode, allele_count
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
}

#' @export
sc_mutations <- function(bam_path, seqnames, positions, indel, barcodes, threads = 1) {
  stopifnot(
    "seqnames not the same length as positions" =
      length(seqnames) == length(positions)
  )

  if (length(bam_path) == 1) {
    # single bam file, parallelize over positions
    stopifnot(
      "barcodes must be a character vector" =
        is.character(barcodes)
    )
    variants <- parallel::mcmapply(
      FUN = function(seqname, pos) {
        variant_count_tb(bam_path, seqname, pos, indel, barcodes)
      },
      seqname = seqnames, pos = positions, SIMPLIFY = FALSE, mc.cores = threads
    ) |>
      dplyr::bind_rows()
  } else {
    # multiple bam files, parallelize over bam files
    stopifnot(
      "barcodes must be a list of character vectors, same length as bam_path" =
        is.list(barcodes) & length(barcodes) == length(bam_path)
    )
    # data frame of all combinations between (seqname, pos) and (bam_path, barcodes)
    args_grid <- expand.grid(
      mutation_index = seq_along(positions),
      bam_index = seq_along(bam_path),
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(
        seqname = seqnames[mutation_index],
        pos = positions[mutation_index],
        sample_bam = bam_path[bam_index],
        sample_barcodes = barcodes[bam_index]
      ) |>
      dplyr::select(-mutation_index, -bam_index)

    variants <- parallel::mcmapply(
      FUN = function(sample_bam, seqname, pos, sample_barcodes) {
        variant_count_tb(sample_bam, seqname, pos, indel, sample_barcodes) |>
          dplyr::mutate(bam_file = sample_bam)
      },
      sample_bam = args_grid$sample_bam, seqname = args_grid$seqname,
      pos = args_grid$pos, sample_barcodes = args_grid$sample_barcodes,
      SIMPLIFY = FALSE, mc.cores = threads
    ) |>
      dplyr::bind_rows()
  }
  return(variants)
}
