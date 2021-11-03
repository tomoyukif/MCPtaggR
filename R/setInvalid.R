#' Set invalid interval for genotyping
#'
#' @param object A SInG object.
#' @param interval A GRanges object of intervals to set invalid.
#' @param genome A string of either "ref" or "alt to indicate which genome
#'  the specified intervals should be set as invalid in.
#'
#' @return A SInG object.
#'
#' @importFrom GenomicRanges reduce
#'
#' @examples
#' # Load a sample SInG object.
#' sing_object <- system.file("extdata", "sample.Rdata", package = "SInG")
#' load(sing_object)
#'
#' #' # import the GFF file of repeat elements in O. sativa cv. Nipponbare
#' # identified by RepeatMasker.
#' filepath <- system.file("extdata", "sample_ref.gff.gz", package = "SInG")
#' ref_rep <- rtracklayer::readGFFAsGRanges(filepath = filepath)
#'
#' # Extract repeat elements in the chromosomes not in the unanchored fragments.
#' chrs <- grepl("chr", tolower(as.character(GenomicRanges::seqnames(ref_rep))))
#' ref_rep <- ref_rep[chrs]
#' GenomeInfoDb::seqlevels(ref_rep) <- GenomeInfoDb::seqlevelsInUse(ref_rep)
#'
#' # Change and sort chromosome names
#' GenomeInfoDb::seqlevels(ref_rep) <- sub("*[^0-9]", "", sub("[^0-9]*", "",
#' as.character(GenomeInfoDb::seqlevels(ref_rep))))
#' ref_rep <- GenomeInfoDb::sortSeqlevels(ref_rep)
#' GenomeInfoDb::seqlevels(ref_rep) <- GenomeInfoDb::seqlevels(object$sin_ref)
#'
#' # import the GFF file of repeat elements in O. glaberrima acc. WK21
#' # identified by RepeatMasker.
#' filepath <- system.file("extdata", "sample_alt.gff.gz", package = "SInG")
#' alt_rep <- rtracklayer::readGFFAsGRanges(filepath = filepath)
#'
#' # Extract repeat elements in the chromosomes not in the unanchored fragments.
#' chrs <- grepl("chr", tolower(as.character(GenomicRanges::seqnames(alt_rep))))
#' alt_rep <- alt_rep[chrs]
#' GenomeInfoDb::seqlevels(alt_rep) <- GenomeInfoDb::seqlevelsInUse(alt_rep)
#'
#' # Change and sort chromosome names
#' GenomeInfoDb::seqlevels(alt_rep) <- sub("*[^0-9]", "", sub("[^0-9]*0*", "",
#' as.character(GenomeInfoDb::seqlevels(alt_rep))))
#' alt_rep <- GenomeInfoDb::sortSeqlevels(alt_rep)
#' GenomeInfoDb::seqlevels(alt_rep) <- GenomeInfoDb::seqlevels(object$sin_alt)
#'
#' @export
#'
#'
setInvalid <- function(object, interval, genome) {
stopifnot(inherits(object, "SInG"))
if (genome == "ref") {
    if(is.null(object$ints$ref_inv)){
        object$ints$ref_inv <- interval

    } else {
        object$ints$ref_inv <- GenomicRanges::reduce(c(object$ints$ref_inv, interval))
    }
} else if (genome == "alt") {
    if(is.null(object$ints$alt_inv)){
        object$ints$alt_inv <- interval

    } else {
        object$ints$alt_inv <- GenomicRanges::reduce(c(object$ints$alt_inv, interval))
    }
}
return(object)
}
