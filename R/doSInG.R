#' Execute syntenic interval genotyping
#'
#' Import bam file(s) and count reads in the syntenic intervals.
#'
#' @param object A SInG obejct.
#' @param bam_files A string vector of paths to BAM files.
#'
#' @return A SInG object with the `counts` element which stores a list of
#' read counts in each interval
#'
#' @importFrom Rsamtools scanBamFlag ScanBamParam BamFile indexBam
#' @importClassesFrom Rsamtools ScanBamParam BamFile
#' @importFrom GenomicAlignments readGAlignmentPairs first second cigar
#' @importClassesFrom GenomicAlignments GAlignmentPairs
#' @importFrom GenomicRanges gaps start setdiff findOverlaps
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom IRanges width
#' @importFrom S4Vectors subjectHits
#'
#' @examples
#' # Load a sample SInG object.
#' sing_object <- system.file("extdata", "sample.Rdata", package = "SInG")
#' load(sing_object)
#'
#' # Get path to a sample BAM file.
#' bam_files <- system.file("extdata", "sample.bam", package = "SInG")
#'
#' object <- doSInG(object, bam_files)
#'
#' @export
#'
#'
doSInG <- function(object, bam_files) {
    stopifnot(inherits(object, "SInG"))
    object <- .setVariableRegion(object)

    bamflag <- Rsamtools::scanBamFlag(isPaired = TRUE,
                                      isProperPair = TRUE,
                                      isUnmappedQuery = FALSE,
                                      hasUnmappedMate = FALSE,
                                      isMinusStrand = NA,
                                      isMateMinusStrand = NA,
                                      isFirstMateRead = NA,
                                      isSecondMateRead = NA,
                                      isSecondaryAlignment = FALSE,
                                      isNotPassingQualityControls = FALSE,
                                      isDuplicate = FALSE,
                                      isSupplementaryAlignment = FALSE)

    param <- Rsamtools::ScanBamParam(flag = bamflag,
                                     what = c("rname", "pos", "cigar", "isize", "groupid"))

    for(bam_file in bam_files){
        read <- .getComplReads(bam_file = bam_file, param = param)
        counts <- .getCounts(object = object, read = read)
        if(is.null(object$counts)){
            object$counts <- list(counts)
        } else {
            object$counts <- c(object$counts, list(counts))
        }
        names(object$counts)[length(object$counts)] <- sub("\\..*", "", sub(".*/", "", bam_file))
    }
    object <- .countIntRead(object)

    return(object)
}


.setVariableRegion <- function(object) {
    n <- length(object$ints$ref_sin)
    ref_vr <- setdiff(object$ints$ref_sin, object$ints$ref_inv, ignore.strand = TRUE)
    parent_ref <- findOverlaps(ref_vr, object$ints$ref_sin)
    ref_vr$parent <- subjectHits(parent_ref)

    alt_vr <- setdiff(object$ints$alt_sin, object$ints$alt_inv, ignore.strand = TRUE)
    parent_alt <- findOverlaps(alt_vr, object$ints$alt_sin)
    alt_vr$parent <- subjectHits(parent_alt)

    ref_parents <- unique(ref_vr$parent)
    alt_parents <- unique(alt_vr$parent)
    useless_ints <- c(setdiff(1:n, ref_parents),
                      setdiff(1:n, alt_parents))
    if(length(useless_ints) != 0){
        ref_vr <- ref_vr[!ref_vr$parent %in% useless_ints]
        alt_vr <- alt_vr[!alt_vr$parent %in% useless_ints]
        ref_vs <- object$ints$ref_sin[-useless_ints]
        alt_vs <- object$ints$alt_sin[-useless_ints]
    }

    ref_vr$parent <- as.numeric(as.factor(ref_vr$parent))
    alt_vr$parent <- as.numeric(as.factor(alt_vr$parent))

    object$varInt <- list(ref_vr = ref_vr,
                          alt_vr = alt_vr,
                          ref_vs = ref_vs,
                          alt_vs = alt_vs)
    return(object)
}


.getComplReads <- function(bam_file = bam_file, param = param){
    ga <- GenomicAlignments::readGAlignmentPairs(file = bam_file, param = param)
    f <- GenomicAlignments::first(ga)
    s <- GenomicAlignments::second(ga)
    comple_reads <- grepl("^[0-9]*M$", GenomicAlignments::cigar(f)) & grepl("^[0-9]*M$", GenomicAlignments::cigar(s))
    return(c(f[comple_reads], s[comple_reads]))
}


.getCounts <- function(object = object, read = read){
    varInt <- object$varInt

    ref_count <- findOverlaps(read, varInt$ref_vr)
    ref_count <- as.data.frame(ref_count)
    ref_count <- factor(ref_count$subjectHits, levels = 1:length(varInt$ref_vr))
    ref_count <- as.numeric(table(ref_count))
    ref_len <- IRanges::width(varInt$ref_vr)
    ref_norm_count <- ref_count / ref_len
    ref <- data.frame(raw = ref_count,
                      norm = ref_norm_count,
                      len = ref_len,
                      parents = varInt$ref_vr$parent)

    alt_count <- findOverlaps(read, varInt$alt_vr)
    alt_count <- as.data.frame(alt_count)
    alt_count <- factor(alt_count$subjectHits, levels = 1:length(varInt$alt_vr))
    alt_count <- as.numeric(table(alt_count))
    alt_len <- IRanges::width(varInt$alt_vr)
    alt_norm_count <- alt_count / alt_len
    alt <- data.frame(raw = alt_count,
                      norm = alt_norm_count,
                      len = alt_len,
                      parents = varInt$alt_vr$parent)

    return(list(raw = list(ref = ref, alt = alt)))
}

.countIntRead <- function(object){
    for(i in seq_along(object$counts)){
        i_count <- object$counts[[i]]$raw
        ref_int <- tapply(seq_along(object$varInt$ref_vr$parent), object$varInt$ref_vr$parent, function(x){
            sum_count <- sum(i_count$ref$raw[x])
            sum_len <- sum(i_count$ref$len[x])
            sum_norm_count <- sum_count / sum_len
            return(c(sum_count, sum_norm_count))
        })
        ref_int <- do.call("rbind", ref_int)

        alt_int <- tapply(seq_along(object$varInt$alt_vr$parent), object$varInt$alt_vr$parent, function(x){
            sum_count <- sum(i_count$alt$raw[x])
            sum_len <- sum(i_count$alt$len[x])
            sum_norm_count <- sum_count / sum_len
            return(c(sum_count, sum_norm_count))
        })
        alt_int <- do.call("rbind", alt_int)

        object$counts[[i]]$sin <- data.frame(ref_raw = as.integer(ref_int[, 1]),
                                             alt_raw = as.integer(alt_int[, 1]),
                                             ref_norm = as.numeric(ref_int[, 2]),
                                             alt_norm = as.numeric(alt_int[, 2]))
    }
    return(object)
}
