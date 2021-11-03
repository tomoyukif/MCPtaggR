#' Load Mummer's output file#' Load Mummer's output file#' Load Mummer's o#' Load Mummer's output file
#'
#' @param mummer_fn A string to specify the path to output file of
#'  Mummer to load.
#'
#'
#' @return A list with two elements which are [GRanges-class]
#'  objects storing stretches described in Mummer's output for
#'  reference genome and alternative genome.
#'
#' @importFrom GenomicRanges start end strand reduce findOverlaps GRanges
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomeInfoDb seqlevels seqnames seqlevels<-
#'
#' @examples
#' mummer_fn <- system.file("extdata", "sample.mummer", package = "SInG")
#' object <- mum2SInG(mummer_fn = mummer_fn)
#'
#' @export
#'
mum2SInG <- function(mummer_fn) {
    mummer <- unlist(read.table(mummer_fn, sep = "\t", fill = TRUE))
    chr_label_pos <- grep(">", mummer)
    n <- length(mummer)
    group  <- cut(1:n,
                  breaks = c(chr_label_pos, n),
                  include.lowest = TRUE,
                  right = FALSE)

    object <- do.call("rbind", tapply(mummer, group, function(x) {
        lab <- sub(">\ +", "", x[1])
        df <-
            do.call("rbind", strsplit(sub("\ ", "", gsub(
                "\ +", "\ ", x[-1]
            )), "\ "))
        df <- data.frame(
            query_chr = sub("\ ", "_", lab),
            query_pos = as.numeric(df[, 3]),
            ref_chr = sub("\ ", "_", df[, 1]),
            ref_pos = as.numeric(df[, 2]),
            width = as.numeric(df[, 4]),
            row.names = NULL
        )
        df <- df[order(df$query_pos, df$ref_chr, df$ref_pos),]
        return(df)
    }))
    rownames(object) <- NULL
    object <- .df2Grange(object)
    object <- .findInt(object)
    class(object) <- c("SInG", class(object))
    return(object)
}


.df2Grange <- function(object) {
    ref <- subset(object, select = grepl("ref", names(object)))
    alt <- subset(object, select = grepl("query", names(object)))
    width <- object$width
    ref_strand <- rep("+", nrow(ref))
    ref_strand[grepl("_Reverse$", ref$ref_chr)] <- "-"
    ref_chr <- sub("_Reverse$", "", ref$ref_chr)
    ref_start <- ref$ref_pos
    ref_start[ref_strand == "-"] <-
        (ref_start - width + 1)[ref_strand == "-"]
    ref_end <- ref$ref_pos + width - 1
    ref_end[ref_strand == "-"] <- ref$ref_pos[ref_strand == "-"]
    ref <- GRanges(
        seqnames = ref_chr,
        ranges = IRanges(start = ref_start,
                         end = ref_end),
        strand = ref_strand
    )
    alt_strand <- rep("+", nrow(alt))
    alt_strand[grepl("_Reverse$", alt$query_chr)] <- "-"
    alt_chr <- sub("_Reverse$", "", alt$query_chr)
    alt_start <- alt$query_pos
    alt_start[alt_strand == "-"] <-
        (alt_start - width + 1)[alt_strand == "-"]
    alt_end <- alt$query_pos + width - 1
    alt_end[alt_strand == "-"] <- alt$query_pos[alt_strand == "-"]
    alt <- GRanges(
        seqnames = alt_chr,
        ranges = IRanges(start = alt_start,
                         end = alt_end),
        strand = alt_strand
    )
    return(list(ref = ref, alt = alt))
}

.findInt <- function(object) {
    ref <- object$ref
    ref_chr <- grepl("chr", tolower(as.character(seqnames(ref))))
    ref <- ref[ref_chr]
    num_levels <- as.numeric(sub("*[^0-9]", "", sub("[^0-9]*", "", seqlevels(ref))))
    seqlevels(ref) <- seqlevels(ref)[order(num_levels)]
    unique_ref_i <- table(queryHits(findOverlaps(ref, ref))) == 1

    alt <- object$alt
    alt_chr <- grepl("chr", tolower(as.character(seqnames(alt))))
    alt <- alt[alt_chr]
    num_levels <- as.numeric(sub("*[^0-9]", "", sub("[^0-9]*", "", seqlevels(alt))))
    seqlevels(alt) <- seqlevels(alt)[order(num_levels)]
    unique_alt_i <- table(queryHits(findOverlaps(alt, alt))) == 1

    forward_alt <- as.vector(strand(alt) == "+")

    same_chr <-
        as.numeric(seqnames(ref)) == as.numeric(seqnames(alt))

    valid_record <-
        unique_ref_i & unique_alt_i & forward_alt & same_chr
    valid_ref <- ref[as.vector(valid_record)]
    valid_alt <- alt[as.vector(valid_record)]

    n <- length(valid_ref)
    ref_order <- order(as.integer(seqnames(valid_ref)), start(valid_ref))
    valid_ref <- valid_ref[ref_order]
    valid_alt <- valid_alt[ref_order]
    alt_order <- order(as.integer(seqnames(valid_alt)), start(valid_alt))
    alt_order_i <- (1:n)[alt_order]

    index <- seq_along(valid_ref)
    index_gaps <- cbind(head(index, -1), tail(index, -1))

    not_chr_border <- diff(as.numeric(seqnames(valid_ref))) == 0
    index_gaps <- index_gaps[not_chr_border, ]

    alt_ordered <- abs(alt_order_i[index_gaps[, 1]] - alt_order_i[index_gaps[, 2]]) == 1
    index_gaps <- index_gaps[alt_ordered, ]

    alt_index_gaps <- cbind(alt_order_i[index_gaps[, 1]], alt_order_i[index_gaps[, 2]])

    ref_start <- start(valid_ref)
    ref_end <- end(valid_ref)
    ref_seqnames <- as.character(seqnames(valid_ref))
    alt_start <- start(valid_alt)
    alt_end <- end(valid_alt)
    alt_seqnames <- as.character(seqnames(valid_alt))

    ref_sin <- GRanges(seqnames = ref_seqnames[index_gaps[, 1]],
                   ranges = IRanges(start = ref_end[index_gaps[, 1]] + 1,
                                    end = ref_start[index_gaps[, 2]] - 1))
    alt_sin <- GRanges(seqnames = alt_seqnames[alt_index_gaps[, 1]],
                   ranges = IRanges(start = alt_end[alt_index_gaps[, 1]] + 1,
                                    end = alt_start[alt_index_gaps[, 2]] - 1))

    ref_inv <- ref[queryHits(findOverlaps(ref, ref_sin))]
    ref_inv <- ref_inv[order(as.numeric(seqnames(ref_inv)), start(ref_inv))]
    alt_inv <- alt[queryHits(findOverlaps(alt, alt_sin))]
    alt_inv <- alt_inv[order(as.numeric(seqnames(alt_inv)), start(alt_inv))]

    object <- list(ints = list(ref_sin = ref_sin,
                               alt_sin = alt_sin,
                               ref_inv = ref_inv,
                               alt_inv = alt_inv))
    return(object)
}
