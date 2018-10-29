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
#' @param y_coord numeric. Y-coordinate for the arrow.
#' @param flatten logical. Flatten multiple transcripts if present.
#' @param angle numeric. angle from the shaft of the arrow to the edge of the
#'   arrow head
#' @param lwd numeric. Weight of lines.
#' @export
plot_gene <- function(
  name,
  refflat = refflat_data,
  y_coord = 0,
  flatten = TRUE,
  angle = 15,
  lwd = 2
) {
  gene_data <- refflat[
    refflat[["geneName"]] == name | refflat[["name"]] == name,
  ]
  if (flatten) {
    gene_data <- data.frame(
      geneName = gene_data[1, "geneName"],
      chrom = gene_data[1, "chrom"],
      strand = gene_data[1, "strand"],
      cdsStart = min(gene_data[["cdsStart"]]),
      cdsEnd = max(gene_data[["cdsEnd"]]),
      stringsAsFactors = FALSE
    )
  }
  gene_length <- gene_data[["cdsEnd"]] - gene_data[["cdsStart"]]
  gene_center <- (gene_data[["cdsStart"]] + gene_data[["cdsEnd"]]) / 2
  plot(
    c(
      gene_center[[1]] - gene_length[[1]],
      gene_center[[1]] + gene_length[[1]]
    ),
    c(y_coord, y_coord + 1),
    col = "white",
    ann = FALSE,
    yaxt = "n"
  )
  title(
    xlab = paste("Chromosome", sub("chr", "", gene_data[1, "chrom"]), "(Mb)")
  )
  if (gene_data[1, "strand"] == "+") {
    arrowhead_code = 2
  } else if (gene_data[1, "strand"] == "-") {
    arrowhead_code = 1
  }
  arrows(
    gene_data[1, "cdsStart"],
    y_coord + 0.5,
    x1 = gene_data[1, "cdsEnd"],
    code = arrowhead_code,
    angle = angle,
    lwd = lwd
  )
  text(gene_center[[1]], y_coord + 0.5, labels = name, font = 3, pos = 3)
}