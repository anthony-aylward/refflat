#===============================================================================
# refflat.R
#===============================================================================

# plot refFlat data




# Functions ====================================================================

#' Fetch refflat data
#'
#' Fetch refSeq refFlat data from the UCSC genome browser
#'
#' @param tmpdir character. A non-empty character vector giving the temporary
#'   directory name (provided to tempfile()).
#' @param method character. Method to be used for downloading the data file.
#' @return A data frame containing the refFlat data (see also ?refflat_data)
#' @export
fetch_refflat_data <- function(tmpdir = tempdir(), method = "auto") {
  temp_file_name <- tempfile(tmpdir = tmpdir)
  download.file(
    "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz",
    temp_file_name,
    method
  )
  refflat_data <- setNames(
    read.table(temp_file_name, stringsAsFactors = FALSE),
    c(
      "geneName",
      "name",
      "chrom",
      "strand",
      "txStart",
      "txEnd",
      "cdsStart",
      "cdsEnd",
      "exonCount",
      "exonStarts",
      "exonEnds"
    )
  )
  file.remove(temp_file_name)
  refflat_data
}

directional_gene_coordinates <- function(start, end, strand = c("+", "-")) {
  if (strand == "+") {
    c(start, end)
  } else if (strand == "-") {
    c(end, start)
  }
}

#' Plot a gene
#'
#' Plot a gene
#'
#' @param name character. The name of the gene to plot.
#' @param refflat data frame. The refflat dataset.
#' @export
plot_gene <- function(
  name,
  refflat = refflat_data,
  y_coord = 0,
  flatten = FALSE
) {
  gene_data <- refflat[
    refflat[["geneName"]] == name | refflat[["name"]] == name,
  ]
  gene_length <- gene_data[["cdsEnd"]] - gene_data[["cdsStart"]]
  gene_center <- (gene_data[["cdsStart"]] + gene_data[["cdsEnd"]]) / 2
  plot(
    c(gene_center[[1]] - gene_length[[1]], gene_center[[1]] + gene_length[[1]]),
    c(y_coord, y_coord + 1),
    col = "white"
  )
  if (strand == "+") {
    arrows(gene_data[1, "cdsStart"], 0.5, x1 = gene_data[1, "cdsEnd"], lwd = 4)
  } else if (strand == "-") {
    arrows(gene_data[1, "cdsEnd"], 0.5, gene_data[1, "cdsStart"], lwd = 4)
  }
}