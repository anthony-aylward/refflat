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

#' Slice the refFlat dataset
#'
#' Slice the refFlat dataset
#'
#' @param chrom character. The chromosome of the slice.
#' @param start integer. The start position of the slice.
#' @param end integer. The end position of the slice.
#' @param refflat data frame. The refflat dataset.
#' @return A data frame representing the requested slice of the refflat data.
#' @export
slice_refflat <- function(chrom, start, end, refflat = refflat_data) {
  refflat[
    (
      refflat[["chrom"]] == chrom 
      & refflat[["cdsStart"]] < end 
      & refflat[["cdsEnd"]] > start
    ),
  ]
}

#' Flatten refFlat data
#'
#' Flatten refFlat data
#'
#' @param chrom character. The chromosome of the slice.
#' @param start integer. The start position of the slice.
#' @param end integer. The end position of the slice.
#' @param refflat data frame. The refflat dataset.
#' @return A data frame representing the requested slice of the refflat data.
#' @export
flatten_refflat <- function(refflat) {
  as.data.frame(
    do.call(
      "rbind",
      lapply(
        unique(refflat[["geneName"]]),
        function(gene_name) {
          gene_data <- refflat[refflat[["geneName"]] == gene_name,]
          list(
            geneName = gene_data[1, "geneName"],
            chrom = gene_data[1, "chrom"],
            strand = gene_data[1, "strand"],
            cdsStart = min(gene_data[["cdsStart"]]),
            cdsEnd = max(gene_data[["cdsEnd"]])
          )
        }
      )
    )
  )
}


#' Draw a gene on a plot
#'
#' Draw a gene on a plot
#'
#' @param gene_data list. List representing the gene.
#' @param y numeric. The y coordinate at which to draw the gene.
#' @param arrowhead_length numeric. Length of the edges of the arrow head (in
#'   inches).
#' @param angle numeric. Angle from the shaft of the arrow to the edge of the
#'   arrow head.
#' @param lwd numeric. Weight of lines.
#' @export
draw_gene <- function(
  gene_data,
  y,
  arrowhead_length = 0.125,
  angle = 15,
  lwd = 2
) {
  gene_center <- (gene_data[["cdsStart"]] + gene_data[["cdsEnd"]]) / 2
  if (gene_data[["strand"]] == "+") {
    arrowhead_code = 2
  } else if (gene_data[["strand"]] == "-") {
    arrowhead_code = 1
  }
  arrows(
    gene_data[["cdsStart"]],
    y_coord + 0.5,
    x1 = gene_data[["cdsEnd"]],
    length = arrowhead_length,
    angle = angle,
    code = arrowhead_code,
    lwd = lwd
  )
  text(
    gene_center,
    y_coord + 0.5,
    labels = gene_data[["geneName"]],
    font = 3,
    pos = 3
  )
}

#' Plot a gene
#'
#' Plot a gene
#'
#' @param name character. The name of the gene to plot.
#' @param refflat data frame. The refflat dataset.
#' @param y_coord numeric. Y-coordinate for the arrow.
#' @param flatten logical. Flatten multiple transcripts if present.
#' @param arrowhead_length numeric. Length of the edges of the arrow head (in
#'   inches).
#' @param angle numeric. Angle from the shaft of the arrow to the edge of the
#'   arrow head.
#' @param lwd numeric. Weight of lines.
#' @export
plot_gene <- function(
  name,
  refflat = refflat_data,
  y_coord = 0,
  flatten = TRUE,
  arrowhead_length = 0.125,
  angle = 15,
  lwd = 2
) {
  gene_data <- refflat[
    refflat[["geneName"]] == name | refflat[["name"]] == name,
  ]
  if (flatten) gene_data <- flatten_refflat(gene_data)
  print(gene_data)
  print(str(gene_data))
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
    xlab = paste("Chromosome", sub("chr", "", gene_data[1, "chrom"]))
  )
  draw_gene(
    gene_data,
    y_coord + 0.5,
    arrowhead_length = arrowhead_length,
    angle = angle,
    lwd = lwd
  )
}
