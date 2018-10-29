#===============================================================================
# refflat.R
#===============================================================================

# plot refFlat data




# Functions ====================================================================

fetch_refflat_data <- function(tmpdir = tempdir(), method = "auto") {
  temp_file_name <- tempfile(tmpdir = tmpdir)
  download.file(
    "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz",
    temp_file_name,
    method
  )
  rf <- setNames(
    read.table(temp_file_name, stringsAsFactors = FALSE)
    c(
      "geneName"
      "name"
      "chrom"
      "strand"
      "txStart"
      "txEnd"
      "cdsStart"
      "cdsEnd"
      "exonCount"
      "exonStarts"
      "exonEnds"
    )
  )
}