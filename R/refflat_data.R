#' RefFlat gene annotations from RefSeq
#'
#' RefFlat gene annotations from RefSeq
#'
#' @format A data frame with 72153 rows and 11 variables:
#' \describe{
#'   \item{geneName}{name of gene as it appears in genome browser}
#'   \item{name}{refFlat name of gene}
#'   \item{chrom}{reference sequence chromosome or scaffold}
#'   \item{strand}{+ or - for strand}
#'   \item{txStart}{transcription start position}
#'   \item{txEnd}{transcription end position}
#'   \item{cdsStart}{coding region start}
#'   \item{cdsEnd}{coding region end}
#'   \item{exonCount}{Number of exons}
#'   \item{exonStarts}{exon start positions}
#'   \item{exonEnds}{exon end positions}
#' }
#' @source \url{http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz}
"refflat_data"