#' Write a CSV file.
#'
#' Output SInG results in a CSV file with interval position information,
#' read counts, and read ratio.
#'
#' @param object A SInG object.
#' @param file_path A string to indicate the file path of
#' the output CSV file.
#' @param dummy_parents A logical value to indicate whether insert dummy
#' data for parental samples or not.
#'
#' @return Invisibly returns TRUE.
#'
#' @importFrom GenomicRanges start end
#' @importFrom GenomeInfoDb seqnames
#'
#' @examples
#' # Load a sample SInG object.
#' sing_object <- system.file("extdata", "sample.Rdata", package = "SInG")
#' load(sing_object)
#'
#' file_path <- tempfile("_SInG", fileext = "vcf")
#'
#' sing2VCF(object, file_path)
#'
#' @export
#'
#'
#'
sing2CSV <- function(object, file_path, dummy_parents = TRUE) {
    options(scipen = 10)

    chr <- as.numeric(seqnames(object$varInt$ref_vs))
    pos_start <- start(object$varInt$ref_vs)
    pos_end <- end(object$varInt$ref_vs)

    df <- lapply(object$counts, function(x){
        return(data.frame(x$sin, ratio = x$ratio, smooth = x$smoothed))
    })

    df <- data.frame(chr = chr,
                     pos_start = pos_start,
                     pos_end = pos_end,
                     df)

    if (file_path == "") {
        file_path <- tempfile(pattern = "singOutput", fileext = ".vcf")
    }
    write.table(x = df,
                file = file_path,
                sep = ",",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)

    invisible(TRUE)
}
