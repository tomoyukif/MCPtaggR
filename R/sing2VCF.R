#' Write a VCF file.
#'
#' Output SInG results in a VCF file with interval position information,
#' read counts, and read ratio.
#'
#' @param object A SInG object.
#' @param file_path A string to indicate the file path of
#' the output VCF file.
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
sing2VCF <- function(object, file_path, dummy_parents = TRUE) {
    options(scipen = 10)

    format <- lapply(object$counts, function(x){
        ref <- x$sin$ref_raw > 0
        alt <- x$sin$alt_raw > 0
        gt <- rep("./.", nrow(x$sin))
        gt[ref & !alt] <- "0/0"
        gt[ref & alt] <- "0/1"
        gt[!ref & alt] <- "1/1"
        ad <- paste(ref, alt, sep = ",")
        nad <- paste(round(x$sin$ref_norm, 5), round(x$sin$alt_norm, 5), sep = ",")
        rt <- round(x$ratio, 5)
        srt <- round(x$smoothed, 5)
        out <- paste(gt, ad, nad, rt, srt, sep = ":")
        return(out)
    })
    format <- do.call("cbind", format)


    chr <- as.numeric(seqnames(object$varInt$ref_vs))
    pos <- start(object$varInt$ref_vs) + round((end(object$varInt$ref_vs) - start(object$varInt$ref_vs)) / 2)
    n_mar <- length(chr)
    n_chr <- length(unique(chr))
    chr <- paste0("chr", sprintf(fmt = paste0("%0", nchar(n_chr), "d"), chr))
    id <- paste0(chr, "_", pos)

    ref_allele <- rep("A", n_mar)
    alt_allele <- rep("G", n_mar)

    sample_id <- names(object$counts)

    if(dummy_parents){
        p1 <- paste(rep("0/0", n_mar),
                    rep("50,0", n_mar),
                    rep("50,0", n_mar),
                    rep("1,0", n_mar),
                    rep("1,0", n_mar), sep = ":")
        p2 <- paste(rep("1/1", n_mar),
                    rep("0,50", n_mar),
                    rep("0,50", n_mar),
                    rep("0,1", n_mar),
                    rep("0,1", n_mar), sep = ":")
        format <- cbind(p1, p2, format)
        sample_id <- c("Ref", "Alt", sample_id)
    }

    out_matrix <- cbind(chr,
                        pos,
                        id,
                        ref_allele,
                        alt_allele,
                        ".",
                        "PASS",
                        ".",
                        "GT:AD:NAD:RT:SRT",
                        format)

    header <- c("#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                sample_id)

    out_matrix <- rbind(header, out_matrix)

    if (file_path == "") {
        file_path <- tempfile(pattern = "singOutput", fileext = ".vcf")
    }
    meta <- c('##fileformat=VCFv4.0',
              '##FILTER=<ID=PASS,Description="All filters passed">',
              '##FILTER=<ID=PASS,Description="All filters passed">',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
              '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the reference and alternate alleles in the order listed">',
              '##FORMAT=<ID=NAD,Number=.,Type=Float,Description="Normalized allelic depths for the reference and alternate alleles in the order listed">',
              '##FORMAT=<ID=RT,Number=.,Type=Float,Description="Reference read ratio">',
              '##FORMAT=<ID=SRT,Number=.,Type=Float,Description="Smoothed reference read ratio">')
    write.table(
        x = meta,
        file = file_path,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        append = FALSE
    )

    write.table(
        x = out_matrix,
        file = file_path,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        append = TRUE
    )

    invisible(TRUE)
}
