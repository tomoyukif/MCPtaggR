library(ggplot2)

# Define the function to convert a mummer output to the data.frame.


# Convert the data.frame to a list of Ranges objects for reference and alternative genomes.
df2gff <- function(df, only = NULL) {
    ref <- subset(df, select = grepl("ref", names(df)))
    alt <- subset(df, select = grepl("query", names(df)))
    width <- df$width
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
        ranges = IRanges::IRanges(start = ref_start,
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
        ranges = IRanges::IRanges(start = alt_start,
                                  end = alt_end),
        strand = alt_strand
    )

    if (is.null(only)) {
        return(list(ref = ref, alt = alt))
    } else if (only == "ref") {
        return(list(ref = c(ref, alt)))
    } else if (only == "alt") {
        return(list(alt = c(ref, alt)))
    } else {
        return(NULL)
    }
}

# Define the function to find valid syntenic intervals.
.rmShortInterval <- function(ref, alt, len) {
    while (TRUE) {
        ref_end <- end(ref)
        ref_start <- start(ref)
        ref_gap_len <- ref_start[-1] - 1 - head(ref_end,-1)
        chr_border <-
            which(diff(as.numeric(as.factor(
                seqnames(ref)
            ))) != 0)
        ref_gap_len[chr_border] <- Inf
        alt_end <- end(alt)
        alt_start <- start(alt)
        alt_gap_len <- alt_start[-1] - 1 - head(alt_end,-1)
        alt_gap_len[chr_border] <- Inf

        invalid_intervals <-
            which(ref_gap_len < len | alt_gap_len < len)
        if (length(invalid_intervals) == 0) {
            break
        } else {
            invalid_intervals <- min(invalid_intervals)
            invalid_intervals <- invalid_intervals + 1
            rm_ref <- ref[invalid_intervals]
            rm_alt <- alt[invalid_intervals]
            ref <- ref[-invalid_intervals]
            alt <- alt[-invalid_intervals]
        }
    }

    return(list(
        ref = ref,
        alt = alt,
        rm_ref = rm_ref,
        rm_alt = rm_alt
    ))
}

getSin <- function(mummer) {
    ref <- mummer$ref
    seqlevels(ref) <- sort(seqlevels(ref))
    unique_ref_i <- table(queryHits(findOverlaps(ref, ref))) == 1

    alt <- mummer$alt
    seqlevels(alt) <- sort(seqlevels(alt))
    unique_alt_i <- table(queryHits(findOverlaps(alt, alt))) == 1

    forward_alt <- as.vector(strand(alt) == "+")

    same_chr <-
        as.character(seqnames(ref)) == as.character(seqnames(alt))

    valid_record <-
        unique_ref_i & unique_alt_i & forward_alt & same_chr
    valid_ref <- ref[as.vector(valid_record)]
    valid_alt <- alt[as.vector(valid_record)]


    n <- length(valid_ref)
    alt_order <- order(seqnames(valid_alt), start(valid_alt))
    valid_ref <- valid_ref[alt_order]
    valid_alt <- valid_alt[alt_order]
    valid_ref$order_alt <- 1:n
    ref_order <- order(seqnames(valid_ref), start(valid_ref))
    valid_ref <- valid_ref[ref_order]
    valid_alt <- valid_alt[ref_order]
    valid_alt$order_ref <- 1:n
    ref_order <- order(seqnames(valid_ref), start(valid_ref))
    valid_ref <- valid_ref[ref_order]
    alt_order <- order(seqnames(valid_alt), start(valid_alt))
    valid_alt <- valid_alt[alt_order]

    ref_i <- 1
    alt_i <- 1
    while (TRUE) {
        check1 <- valid_ref$order_alt[ref_i] - alt_i
        check2 <- valid_alt$order_ref[alt_i] - ref_i
        if (check1 == 0 & check2 == 0) {
            ref_i <- ref_i + 1
            alt_i <- alt_i + 1
        } else if (check1 == check2) {
            valid_ref$order_alt[ref_i] <- NA
            valid_alt$order_ref[alt_i] <- NA
            ref_i <- ref_i + 1
            alt_i <- alt_i + 1
        } else if ({
            check1 > 0 & check1 < check2
        } | check2 < 0) {
            valid_alt$order_ref[alt_i] <- NA
            alt_i <- alt_i + 1
        } else if ({
            check2 > 0 & check1 > check2
        } | check1 < 0) {
            valid_ref$order_alt[ref_i] <- NA
            ref_i <- ref_i + 1
        }
        if (ref_i > n | alt_i > n) {
            break
        }
    }
    valid_ref <- valid_ref[!is.na(valid_ref$order_alt)]
    valid_alt <- valid_alt[!is.na(valid_alt$order_ref)]

    out <-
        .rmShortInterval(ref = valid_ref, alt = valid_alt, len = 1)
    valid_ref <- out$ref
    valid_alt <- out$alt

    invalid_ref <-
        reduce(ref[-subjectHits(findOverlaps(valid_ref, ref))])
    invalid_alt <-
        reduce(alt[-subjectHits(findOverlaps(valid_alt, alt))])
    map <-
        list(
            sin_ref = valid_ref,
            sin_alt = valid_alt,
            invs_ref = invalid_ref,
            invs_alt = invalid_alt
        )
    class(map) <- c("sinmap", class(map))
    return(map)
}

# Combine GRanges for repeat elements and invalid stretches.
mergeInv <- function(map, inv, allele) {
    stopifnot(inherits(map, "sinmap"))
    if (allele == "ref") {
        sinmap$invs_ref <- reduce(c(sinmap$invs_ref, inv))
    } else if (allele == "alt") {
        sinmap$invs_alt <- reduce(c(sinmap$invs_alt, inv))
    }
    return(map)
}

# Build the map of valid variable stretches for genotyping.
buildVvsMap <- function(map) {
    stopifnot(inherits(map, "sinmap"))

    ref_vs <- gaps(map$sin_ref)
    ref_vs <- ref_vs[start(ref_vs) != 1]
    ref <- setdiff(ref_vs, map$invs_ref)
    parent_ref <- findOverlaps(ref, ref_vs)
    ref$parent <- subjectHits(parent_ref)

    alt_vs <- gaps(map$sin_alt)
    alt_vs <- alt_vs[start(alt_vs) != 1]
    alt <- setdiff(alt_vs, map$invs_alt)
    parent_alt <- findOverlaps(alt, alt_vs)
    alt$parent <- subjectHits(parent_alt)

    ref_parents <- unique(ref$parent)
    alt_parents <- unique(alt$parent)
    useless_ints <-
        c(setdiff(ref_parents, alt_parents),
          setdiff(alt_parents, ref_parents))
    ref$parent[ref$parent %in% useless_ints] <-
        ref$parent[ref$parent %in% useless_ints] - 1
    alt$parent[alt$parent %in% useless_ints] <-
        alt$parent[alt$parent %in% useless_ints] - 1

    out_msp <- list(ref = ref, alt = alt)
    class(out_msp) <- c("vvsmap", class(out_msp))
    return(out_msp)
}

# Widen intervals with too short lengths.
widenVarInt <- function(map, len) {
    stopifnot(inherits(map, "vvsmap"))
    int_len_ref <- tapply(width(map$ref), map$ref$parent, sum)
    int_len_alt <- tapply(width(map$alt), map$alt$parent, sum)
    n <- length(int_len_ref)
    cum_len_ref <- 0
    cum_len_alt <- 0
    int_group <- NULL
    int_count <- 1
    for (i in 1:n) {
        cum_len_ref <- cum_len_ref + int_len_ref[i]
        cum_len_alt <- cum_len_alt + int_len_alt[i]
        if (cum_len_ref > len & cum_len_alt > len) {
            cum_len_ref <- 0
            cum_len_alt <- 0
            int_group <- c(int_group, int_count)
            int_count <- int_count + 1
        } else {
            int_group <- c(int_group, int_count)
        }
    }
    map$ref$parent <- int_group[map$ref$parent]
    map$alt$parent <- int_group[map$alt$parent]

    return(map)
}

# Execute SInG.
sing <- function(bam_file, map) {
    bamflag <-
        scanBamFlag(
            isPaired = TRUE,
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
            isSupplementaryAlignment = FALSE
        )
    param <-
        ScanBamParam(flag = bamflag,
                     what = c("rname", "pos", "cigar", "isize", "groupid"))

    ga <- readGAlignmentPairs(file = bam_file, param = param)
    f <- first(ga)
    s <- second(ga)
    comple_reads <- cigar(f) == "150M" & cigar(s) == "150M"
    f <- f[comple_reads]
    s <- s[comple_reads]
    read <- c(f, s)

    parent_levels <- unique(map$ref$parent)

    ref_count <- findOverlaps(read, map$ref)
    ref_count <- as.data.frame(ref_count)
    ref_count$parents <- map$ref$parent[ref_count$subjectHits]
    ref_count <- subset(ref_count, select = -subjectHits)
    ref_count <- unique(ref_count)
    ref_count <- factor(ref_count$parents, levels = parent_levels)
    ref_count <- table(ref_count)
    ref_len <- tapply(width(map$ref), map$ref$parent, sum)
    ref_norm_count <- ref_count / ref_len

    alt_count <- findOverlaps(read, map$alt)
    alt_count <- as.data.frame(alt_count)
    alt_count$parents <- map$alt$parent[alt_count$subjectHits]
    alt_count <- subset(alt_count, select = -subjectHits)
    alt_count <- unique(alt_count)
    alt_count <- factor(alt_count$parents, levels = parent_levels)
    alt_count <- table(alt_count)
    alt_len <- tapply(width(map$alt), map$alt$parent, sum)
    alt_norm_count <- alt_count / alt_len

    return(cbind(ref_count, alt_count, ref_norm_count, alt_norm_count))
}


.writeVCF <- function(sing, file_path) {
    options(scipen = 10)
    gt <- apply(sing, c(1, 3), function(x) {
        x <- as.character(sum(which(x[1:2] > 0)))
        x <- switch(
            x,
            "0" = "./.",
            "1" = "0/0",
            "2" = "1/1",
            "3" = "0/1"
        )
        return(x)
    })
    att <- attributes(sing)
    chr <- att$vvs_info.chr
    pos <- as.integer((att$vvs_info.start + att$vvs_info.end) / 2)
    n_mar <- length(chr)

    id <-
        paste0("Chr",
               sprintf(fmt = "%02d", chr),
               "_",
               pos,
               "_",
               att$vvs_info.end)

    ref <- rep("A", n_mar)
    alt <- rep("G", n_mar)

    rad <- sing[, 1, ]
    aad <- sing[, 2, ]
    n_rad <- sing[, 3, ]
    n_aad <- sing[, 4, ]
    ad <- paste(rad, aad, sep = ",")
    n_ad <- paste(n_rad, n_aad, sep = ",")
    dp <- rad + aad
    genodata <-
        matrix(data = paste(gt, ad, dp, n_ad , sep = ":"), nrow = n_mar)
    n_chr <- length(unique(chr))
    chr <-
        paste0("chr", sprintf(fmt = paste0("%0", nchar(n_chr), "d"), chr))

    out_matrix <- cbind(chr,
                        pos,
                        id,
                        ref,
                        alt,
                        ".",
                        "PASS",
                        ".",
                        "GT:AD:DP:NAD",
                        genodata)

    header <- c(
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        att$dimnames$sample
    )

    out_matrix <- rbind(header, out_matrix)

    if (file_path == "") {
        file_path <- tempfile(pattern = "singOutput", fileext = ".vcf")
    }
    meta <- c(
        '##fileformat=VCFv4.0',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the reference and alternate alleles in the order listed">',
        '##FORMAT=<ID=NAD,Number=.,Type=Float,Description="Normalized allelic depths for the reference and alternate alleles in the order listed">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (only filtered reads used for calling)">'
    )
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

    return(file_path)
}

smoothRatio <- function(x,
                        window,
                        nad = TRUE) {
    chr <- getChromosome(x)
    pos <- getPosition(x)
    validscan <- getValidScan(x)
    validmarker <- getValidSnp(x)
    smoothed_ratio_node <- gdsfmt::add.gdsn(
        node = x@data@handler,
        name = "smooth.ratio",
        storage = "float32",
        compress = "LZMA_RA",
        replace = TRUE
    )
    if (nad) {
        ad_data_node <-
            gdsfmt::index.gdsn(node = x@data@handler,
                               path = "/annotation/format/NAD/data")
    } else {
        if (x@data@genotypeVar == "genotype.filt") {
            ad_data_node <-
                gdsfmt::index.gdsn(node = x@data@handler,
                                   path = "/annotation/format/AD/filt.data")
        } else {
            ad_data_node <-
                gdsfmt::index.gdsn(node = x@data@handler,
                                   path = "/annotation/format/AD/data")
        }
    }

    for (chr_i in unique(chr)) {
        message(paste0('Now smoothing chr ', chr_i))
        target <- chr == chr_i
        pos_chr <- pos[target]
        validmarker_chr <- validmarker & target
        sel_marker <- rep(validmarker_chr, each = 2)
        sel <- list(validscan, sel_marker)
        ad <- gdsfmt::readex.gdsn(ad_data_node, sel = sel)
        ratio <-
            ad[, c(TRUE, FALSE)] / (ad[, c(TRUE, FALSE)] + ad[, c(FALSE, TRUE)])
        smoothed <- NULL
        for (i in 1:length(pos_chr)) {
            tmp_pos <- pos_chr[i]
            in_window <-
                pos_chr >= tmp_pos - window & pos_chr <= tmp_pos + window
            if (sum(in_window) < 10) {
                in_window <-
                    which(1:length(pos_chr) >= i - 5 & 1:length(pos_chr) <= i + 4)
                if (length(in_window) < 10) {
                    if (1 %in% in_window) {
                        in_window <- seq(in_window[1],
                                         by = 1,
                                         length.out = 10)
                    } else if (length(pos_chr) %in% in_window) {
                        in_window <- seq(tail(in_window, 1),
                                         by = -1,
                                         length.out = 10)
                    }
                }
                smoothed <-
                    cbind(smoothed, rowMeans(ratio[, in_window], na.rm = T))
            } else {
                smoothed <- cbind(smoothed, rowMeans(ratio[, in_window], na.rm = T))
            }
        }
        gdsn_dim <- gdsfmt::objdesp.gdsn(smoothed_ratio_node)$dim
        if (sum(gdsn_dim) == 0) {
            smoothed_ratio_node <- gdsfmt::add.gdsn(
                node = x@data@handler,
                name = "smooth.ratio",
                storage = "float32",
                compress = "LZMA_RA",
                replace = TRUE,
                val = smoothed,
                valdim = dim(smoothed)
            )
        } else {
            gdsfmt::append.gdsn(node = smoothed_ratio_node,
                                val = smoothed)
        }
    }
    gdsfmt::readmode.gdsn(smoothed_ratio_node)
    return(x)
}


plotReadRatio <- function(x,
                          coord = NULL,
                          chr = NULL,
                          ind = NULL,
                          valid_only = TRUE,
                          dot_fill = "blue",
                          alpha = 0.1,
                          nad = TRUE,
                          smooth = TRUE,
                          lwd = 1,
                          overlap = FALSE,
                          no_legend = TRUE,
                          font_size = 1,
                          ...) {
    validscan <- getValidScan(x)
    validmarker <- getValidSnp(x)
    if (is.null(ind)) {
        ind <- 1:length(validscan)
    }

    if (is.null(chr)) {
        chr <- rep(TRUE, sum(validmarker))
    } else {
        chr <- getChromosome(x) %in% chr
    }
    validmarker <- validmarker & chr
    sel_marker <- rep(validmarker, each = 2)

    parents <- x@scanAnnot$parents[x@scanAnnot$parents != 0]
    if (!is.null(parents)) {
        ind <- c(ind[ind %in% parents], ind[!ind %in% parents])
    }

    if (nad) {
        ad_data_node <- gdsfmt::index.gdsn(node = x@data@handler,
                                           path = "/annotation/format/NAD/data")
    } else {
        if (x@data@genotypeVar == "genotype.filt") {
            ad_data_node <- gdsfmt::index.gdsn(node = x@data@handler,
                                               path = "/annotation/format/AD/filt.data")
        } else {
            ad_data_node <- gdsfmt::index.gdsn(node = x@data@handler,
                                               path = "/annotation/format/AD/data")
        }
    }

    if (smooth) {
        smooth_node <-
            gdsfmt::index.gdsn(node = x@data@handler, path = "smooth.ratio")
    }


    if (overlap) {
        sel <- list(1:nscan(x) %in% ind, validmarker)
        smoothed <- gdsfmt::readex.gdsn(node = smooth_node, sel = sel)
        smoothed <- t(smoothed)
        if (!is.null(x@snpAnnot$flipped)) {
            smoothed[flipped] <- 1 - smoothed[flipped]
        }
        df <- data.frame(chr = getChromosome(x)[validmarker],
                         pos = getPosition(x)[validmarker],
                         smoothed)
        id <- getScanID(x, FALSE)[1:nscan(x) %in% ind]
        names(df) <- c("chr", "pos", id)
        df <-
            tidyr::pivot_longer(
                df,
                cols = -c(chr, pos),
                names_to = "id",
                values_to = "ratio"
            )

        p <- ggplot(
            data = df,
            mapping = aes(
                x = pos * 10 ^ -6,
                y = ratio,
                group = id,
                color = id
            )
        ) +
            geom_line(size = lwd, alpha = alpha) +
            ylim(0, 1) +
            xlab("Position (Mb)") +
            ylab("Reference allele read ratio") +
            facet_wrap(
                facets = ~ chr,
                nrow = coord[1],
                ncol = coord[2],
                scales = "free_x",
                dir = "v",
                strip.position = "right"
            ) +
            theme(
                axis.title = element_text(size = rel(font_size)),
                axis.text = element_text(size = rel(font_size)),
                plot.title = element_text(size = rel(font_size)),
                legend.text = element_text(size = rel(font_size)),
                legend.title = element_text(size = rel(font_size)),
                strip.text = element_text(size = rel(font_size))
            )
        if (no_legend) {
            p <- p + theme(legend.position = "none")
        }
        print(p)

    } else {
        for (i in ind) {
            if (valid_only & !validscan[i]) {
                next()
            }
            id <- getScanID(x, FALSE)[i]
            sel <- list(1:nscan(x) %in% i, sel_marker)
            ad <- gdsfmt::readex.gdsn(ad_data_node, sel = sel)
            ref <- ad[c(TRUE, FALSE)]
            alt <- ad[c(FALSE, TRUE)]
            if (smooth) {
                sel <- list(1:nscan(x) %in% i, validmarker)
                smoothed <-
                    gdsfmt::readex.gdsn(node = smooth_node, sel = sel)
            }

            if (!is.null(x@snpAnnot$flipped)) {
                flipped <- ref[x@snpAnnot$flipped]
                ref[x@snpAnnot$flipped] <- alt[x@snpAnnot$flipped]
                alt[x@snpAnnot$flipped] <- flipped
                if (smooth) {
                    smoothed[flipped] <- 1 - smoothed[flipped]
                }
            }

            if (smooth) {
                df <- data.frame(
                    chr = getChromosome(x)[validmarker],
                    pos = getPosition(x)[validmarker],
                    ad = ref / (ref + alt),
                    smooth = smoothed
                )
            } else {
                df <- data.frame(
                    chr = getChromosome(x)[validmarker],
                    pos = getPosition(x)[validmarker],
                    ad = ref / (ref + alt)
                )
            }

            p <- ggplot(
                data = df,
                mapping = aes(
                    x = pos * 10 ^ -6,
                    y = ad,
                    group = chr,
                    color = chr
                )
            ) +
                geom_point(alpha = alpha, shape = 20) +
                # geom_line() +
                labs(title = id) +
                ylim(0, 1) +
                xlab("Position (Mb)") +
                ylab("Reference allele read ratio") +
                facet_wrap(
                    facets = ~ chr,
                    nrow = coord[1],
                    ncol = coord[2],
                    scales = "free_x",
                    dir = "v",
                    strip.position = "right"
                ) +
                theme(
                    axis.title = element_text(size = rel(font_size)),
                    axis.text = element_text(size = rel(font_size)),
                    plot.title = element_text(size = rel(font_size)),
                    legend.text = element_text(size = rel(font_size)),
                    legend.title = element_text(size = rel(font_size)),
                    strip.text = element_text(size = rel(font_size))
                )

            if (no_legend) {
                p <- p + theme(legend.position = "none")
            }
            if (smooth) {
                p <- p + geom_line(
                    mapping = aes(
                        x = pos * 10 ^ -6,
                        y = smooth,
                        group = chr
                    ),
                    color = "magenta",
                    size = lwd
                )
            }
            print(p)
        }
    }
}
