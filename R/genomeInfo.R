#' Genome helpers for the R scripts: the R side of bin/genomeInfo.sh.
#'
#' A supported genome is a row of lib/genomes/GENOMES.tsv. The genome of
#' a run is detected once by pipe.sh from the BAM headers and recorded as
#' GENOME= in 00.RUNSTATUS.txt; post-run scripts read it from there and
#' take no genome argument.

#' Get one field of a genome from lib/genomes/GENOMES.tsv
#'
#' @param tag Character genome tag (b37, hg38, mm10)
#' @param field Character column name of GENOMES.tsv (TXDB, ORGDB, ...)
#' @param sdir Character path of the ATAC-seq checkout
#' @return Character value of the field
genomeInfo <- function(tag, field, sdir) {
  table_file <- file.path(sdir, "lib", "genomes", "GENOMES.tsv")
  genomes <- readr::read_tsv(table_file,
                             col_types = readr::cols(.default = "c"))

  if (!field %in% names(genomes)) {
    stop(stringr::str_glue(
      "unknown genome field [{field}]; fields: ",
      "{stringr::str_c(names(genomes), collapse = ' ')}"
    ), call. = FALSE)
  }

  row <- dplyr::filter(genomes, .data$TAG == tag)
  if (nrow(row) != 1) {
    stop(stringr::str_glue(
      "unsupported genome [{tag}]; supported: ",
      "{stringr::str_c(genomes$TAG, collapse = ' ')}"
    ), call. = FALSE)
  }

  row[[field]]
}

#' Get the genome of the run in the current directory
#'
#' @param status_file Character path of the run status file
#' @return Character genome tag from the GENOME= line
runGenome <- function(status_file = "00.RUNSTATUS.txt") {
  if (!fs::file_exists(status_file)) {
    stop(stringr::str_glue(
      "no {status_file} in {getwd()}; run this from the analysis ",
      "directory of a pipe.sh run"
    ), call. = FALSE)
  }

  genome <- readr::read_lines(status_file) |>
    stringr::str_subset("^GENOME=") |>
    stringr::str_remove("^GENOME=")

  if (length(genome) != 1 || genome == "") {
    stop(stringr::str_glue("no GENOME= value in {status_file}"),
         call. = FALSE)
  }

  genome
}
