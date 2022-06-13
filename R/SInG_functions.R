################################################################################
# Make SNP list from mummer results
#'
#' @importFrom data.table fread
#' @importClassesFrom IRanges IRanges
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#'
#' @export
#'
makeSnpList <- function(snps_fn, coord_fn){
    snps <- fread(file = snps_fn,
                  header = FALSE,
                  skip = 4,
                  sep = "\t",
                  colClasses = c("integer", "factor", "factor",
                                 "integer", "integer", "integer",
                                 "integer", "integer", "integer", "integer",
                                 "factor", "factor"))
    message("Number of SNPs: ", nrow(snps))
    ref_chr <- levels(snps$V11)
    ref_chr <- gsub("[^0-9]", "", ref_chr)
    ref_chr <- as.numeric(ref_chr)
    ref_chr <- ref_chr[as.numeric(snps$V11)]
    alt_chr <- levels(snps$V12)
    alt_chr <- gsub("[^0-9]", "", alt_chr)
    alt_chr <- as.numeric(alt_chr)
    alt_chr <- alt_chr[as.numeric(snps$V12)]
    trans_chr <- which(ref_chr != alt_chr)
    snps <- snps[-trans_chr, ]
    message("Number of SNPs after filtering out",
            " inter-chromosomally translocated SNPs: ",
            nrow(snps))
    unknown <- which(snps$V2 == "N" | snps$V3 == "N")
    snps <- snps[-unknown, ]
    message("Number of SNPs after filtering out SNPs at N: ",
            nrow(snps))
    ref_snps <- GRanges(snps[, V11], IRanges(snps[, V1], width = 1))
    ref_snps$allele <- snps$V2
    ref_snps$buff <- snps$V5
    ref_snps$dist <- snps$V6
    alt_snps <- GRanges(snps[, V12], IRanges(snps[, V4], width = 1))
    alt_snps$allele <- snps$V3

    coord <- fread(file = coord_fn,
                   header = FALSE,
                   skip = 4,
                   sep = "\t",
                   colClasses = c("integer", "integer", "integer", "integer",
                                  "integer", "integer", "numeric", "integer",
                                  "integer", "factor", "factor"))
    ref_chr <- as.numeric(gsub("[^0-9]", "", levels(coord$V10)))
    ref_chr <- ref_chr[as.numeric(coord$V10)]
    alt_chr <- as.numeric(gsub("[^0-9]", "", levels(coord$V11)))
    alt_chr <- alt_chr[as.numeric(coord$V11)]
    trans_chr <- which(ref_chr != alt_chr)
    coord <- coord[-trans_chr, ]

    ref_coord <- GRanges(coord$V10, IRanges(coord$V1, coord$V2))

    hits <- unique(queryHits(findOverlaps(ref_snps, ref_coord)))
    ref_snps <- ref_snps[hits]
    alt_snps <- alt_snps[hits]
    ref_snps$id <- seq_along(ref_snps)
    alt_snps$id <- seq_along(alt_snps)

    snps <- list(ref = ref_snps, alt = alt_snps)
    class(snps) <- c(class(snps), "SNPs")
    return(snps)
}

################################################################################
### Make RE site lists, detectable SNP list, and digested genome fragments
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @import BSgenome
#'
#' @export
digestGenome <- function(snps, ref_fn, alt_fn, prefix, read_len, re){
    stopifnot(inherits(snps, "SNPs"))

    ref_genome <- readDNAStringSet(ref_fn)
    ref_seg <- .insilicoRE(ref_genome, re)
    ref_snps <- .validSnps(ref_seg, snps$ref, read_len)
    ref_reads <- .getReads(ref_genome, ref_seg, read_len)
    fastq_fn <- c(paste0(prefix, "ref_R1.fq.gz"),
                  paste0(prefix, "ref_R2.fq.gz"),
                  paste0(prefix, "alt_R1.fq.gz"),
                  paste0(prefix, "alt_R2.fq.gz"))
    writeXStringSet(ref_reads$read1, fastq_fn[1],
                    compress = TRUE, format = "fastq")
    writeXStringSet(ref_reads$read2, fastq_fn[2],
                    compress = TRUE, format = "fastq")

    alt_genome <- readDNAStringSet(alt_fn)
    alt_seg <- .insilicoRE(alt_genome, re)
    alt_snps <- .validSnps(alt_seg, snps$alt, read_len)
    alt_reads <- .getReads(alt_genome, alt_seg, read_len)
    writeXStringSet(alt_reads$read1, fastq_fn[3],
                    compress = TRUE, format = "fastq")
    writeXStringSet(alt_reads$read2, fastq_fn[4],
                    compress = TRUE, format = "fastq")
    names(fastq_fn) <- c("ref_r1", "ref_r2", "alt_r1", "alt_r2")
    names(ref_genome) <- paste0("ref_", names(ref_genome))
    names(alt_genome) <- paste0("alt_", names(alt_genome))
    genome <- c(ref_genome, alt_genome)
    genome_fn <- paste0(prefix, "merge.fa.gz")
    writeXStringSet(genome, genome_fn,
                    compress = TRUE, format = "fasta")
    dg <- list(seg = list(ref = ref_seg, alt = alt_seg),
               snps = list(ref = ref_snps, alt = alt_snps),
               fq = fastq_fn, genome = genome_fn)
    class(dg) <- c(class(dg), "DG")
    return(dg)
}

#' @importFrom XVector subseq
#' @importFrom Biostrings reverseComplement
.getReads <- function(genome, seg, read_len){
    reads <- genome[seg]
    long_reads <- width(reads) > read_len
    read1 <- reads
    read1[long_reads] <- subseq(read1[long_reads],
                                start = 1, width = read_len)
    names(read1) <- paste(seqnames(seg), start(seg), end(seg), "R1", sep = "_")
    read2 <- reverseComplement(reads)
    read2[long_reads] <- subseq(read2[long_reads],
                                start = 1, width = read_len)
    names(read2) <- paste(seqnames(seg), start(seg), end(seg), "R2", sep = "_")
    return(list(read1 = read1, read2 = read2))
}

#' @importFrom GenomicRanges resize findOverlaps width
#' @importFrom S4Vectors queryHits
.validSnps <- function(seg, snps, read_len){
    seg_w <- width(seg)
    long_seg <- seg_w > read_len * 2
    short_seg <- seg[!long_seg]
    long_seg1 <- resize(seg[long_seg], width = 150)
    long_seg2 <- resize(seg[long_seg], width = 150, "end")
    ol <- findOverlaps(snps, c(short_seg, long_seg1, long_seg2))
    return(snps[unique(queryHits(ol))])
}

#' @importFrom DECIPHER DigestDNA
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start<- end<-
.insilicoRE <- function(genome, re){
    re1 <- DigestDNA(re$site[1], genome, "position", "top")
    re1 <- unlist(re1)
    re1 <- data.frame(chr = sub("\\..*", "", names(re1)),
                      pos = unlist(re1),
                      re = "re1")
    re2 <- DigestDNA(re$site[2], genome, "position", "top")
    re2 <- unlist(re2)
    re2 <- data.frame(chr = sub("\\..*", "", names(re2)),
                      pos = unlist(re2),
                      re = "re2")
    re_sites <- rbind(re1, re2)
    re_sites <- re_sites[order(re_sites$chr, re_sites$pos), ]
    valid <- re_sites$re[1:(nrow(re_sites) - 1)] != re_sites$re[2:nrow(re_sites)]
    valid_seg <- cbind(re_sites[c(valid, FALSE), ], re_sites[c(FALSE, valid), ])
    valid_seg <- valid_seg[valid_seg[, 1] == valid_seg[, 4], ]
    valid_seg$strand <- "+"
    valid_seg$strand[valid_seg[, 3] == "re2"] <- "-"
    valid_seg <- GRanges(valid_seg[, 1],
                         IRanges(valid_seg[, 2], valid_seg[, 5]),
                         valid_seg$strand)
    st <- as.character(strand(valid_seg))
    start(valid_seg[st == "+"]) <- start(valid_seg[st == "+"]) - re$f[1]
    end(valid_seg[st == "+"]) <- end(valid_seg[st == "+"]) + re$f[2]
    start(valid_seg[st == "-"]) <- start(valid_seg[st == "-"]) - re$r[2]
    end(valid_seg[st == "-"]) <- end(valid_seg[st == "-"]) + re$r[1]
    return(valid_seg)
}

################################################################################
# Align reads
#' @importFrom Rsubread buildindex align
#' @importFrom Rsamtools sortBam
#' @export

alignRead <- function(dg,
                      n_threads,
                      indexing){
    stopifnot(inherits(dg, "DG"))
    index_fn <- sub("\\.fa.gz", "", dg$genome)
    bam_fn <- c(sub("\\.fa.gz", "_ref.bam", dg$genome),
                sub("\\.fa.gz", "_alt.bam", dg$genome))

    if(indexing){
        buildindex(index_fn, dg$genome)
    }

    ref_aln <- align(index_fn, dg$fq["ref_r1"], dg$fq["ref_r2"],
                     type = "dna", output_file = bam_fn[1],
                     nthreads = n_threads, maxMismatches = 0, unique = TRUE,
                     indels = 0)
    alt_aln <- align(index_fn, dg$fq["alt_r1"], dg$fq["alt_r2"],
                     type = "dna", output_file = bam_fn[2],
                     nthreads = n_threads, maxMismatches = 0, unique = TRUE,
                     indels = 0)
    for(i in seq_along(bam_fn)){
        sortBam(bam_fn[i], sub("\\.bam", "", bam_fn[i]), byQname=TRUE)
    }
    ar <- list(aln_stats = cbind(ref_aln, alt_aln),
               bam_fn = bam_fn, index_fn = index_fn)
    class(ar) <- c(class(ar), "AR")
    return(ar)
}


################################################################################
#'
#' @importFrom Rsamtools scanBam
#' @export
# Get valid TAG markers
getTAG <- function(ar, dg){
    stopifnot(inherits(ar, "AR"))
    stopifnot(inherits(dg, "DG"))
    dg <- .getPairedTag(dg)
    bam1 <- scanBam(ar$bam_fn[1])[[1]]
    seg_ref <- .pickValidTag(bam1, dg$seg$ref)

    bam2 <- scanBam(ar$bam_fn[2])[[1]]
    seg_alt <- .pickValidTag(bam2, dg$seg$alt)

    tags <- .getValidTagPair(seg_ref, seg_alt, dg)
    tags <- .getTagSeq(tags, dg$genome)
    class(tags) <- c(class(tags), "TAG")
    return(tags)
}

#' @importFrom Biostrings readDNAStringSet writeXStringSet writeXStringSet
#' @importFrom GenomeInfoDb seqlevels<- seqlevels
#' @importFrom Rsubread buildindex
.getTagSeq <- function(tags, genome_fn){
    genome <- readDNAStringSet(genome_fn)
    seqlevels(tags$ref) <- paste0("ref_", seqlevels(tags$ref))
    seqlevels(tags$alt) <- paste0("alt_", seqlevels(tags$alt))
    tag_seq <- c(genome[tags$ref], genome[tags$alt])
    names(tag_seq) <- c(paste0("ref_", seq_along(tags$ref)),
                        paste0("alt_", seq_along(tags$alt)))
    tags$ref$id <- paste0("ref_", seq_along(tags$ref))
    tags$alt$id <- paste0("alt_", seq_along(tags$alt))
    stag_fn <- sub("_merge", "_STAG", genome_fn)
    writeXStringSet(tag_seq, stag_fn)
    index_fn <- sub("\\.fa.gz", "", stag_fn)
    buildindex(index_fn, stag_fn)
    tags <- c(tags, list(index_fn = index_fn))
    return(tags)
}

#' @importFrom GenomicRanges findOverlaps reduce
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start
#' @importFrom dplyr left_join
#' @importFrom S4Vectors queryHits
.getValidTagPair <- function(seg_ref, seg_alt, dg){
    ol <- findOverlaps(seg_ref, dg$snps$ref)
    valiable_tag <- seq_along(seg_ref) %in% unique(queryHits(ol))
    seg_ref <- seg_ref[valiable_tag]
    seg_ref <- reduce(seg_ref)
    seg_ref <- seg_ref[order(as.character(seqnames(seg_ref)), start(seg_ref))]

    ol <- findOverlaps(seg_alt, dg$snps$alt)
    valiable_tag <- seq_along(seg_alt) %in% unique(queryHits(ol))
    seg_alt <- seg_alt[valiable_tag]
    seg_alt <- reduce(seg_alt)
    seg_alt <- seg_alt[order(as.character(seqnames(seg_alt)), start(seg_alt))]

    ref_map <- as.data.frame(findOverlaps(seg_ref, dg$snps$ref))
    alt_map <- as.data.frame(findOverlaps(seg_alt, dg$snps$alt))
    map <- left_join(ref_map, alt_map, "subjectHits")
    map <- unique(map[, -2])
    map <- map[!is.na(map$queryHits.x) & !is.na(map$queryHits.y), ]
    seg_ref <- seg_ref[map$queryHits.x]
    seg_alt <- seg_alt[map$queryHits.y]

    return(list(ref = seg_ref, alt = seg_alt))
}

.getPairedTag <- function(dg){
    valid <- (dg$snps$ref$id %in% dg$snps$alt$id)
    valid_id <- dg$snps$ref$id[valid]
    dg$snps$ref <- dg$snps$ref[dg$snps$ref$id %in% valid_id]
    dg$snps$alt <- dg$snps$alt[dg$snps$alt$id %in% valid_id]
    valid_ref_seg <- findOverlaps(dg$seg$ref, dg$snps$ref)
    dg$seg$ref$valid <- FALSE
    dg$seg$ref$valid[unique(queryHits(valid_ref_seg))] <- TRUE
    valid_alt_seg <- findOverlaps(dg$seg$alt, dg$snps$alt)
    dg$seg$alt$valid <- FALSE
    dg$seg$alt$valid[unique(queryHits(valid_alt_seg))] <- TRUE
    return(dg)
}

#' @importFrom GenomicRanges resize
#' @importFrom BiocGenerics width
.pickValidTag <- function(bam, seg){
    max_len <- max(bam$qwidth, na.rm = TRUE)
    long_seg <- width(seg) > max_len
    seg1 <- seg2 <- seg
    seg1[long_seg] <- resize(seg1[long_seg], max_len)
    seg2[long_seg] <- resize(seg2[long_seg], max_len, "end")
    seg <- c(seg1, seg2)
    seg[c(TRUE, FALSE)] <- seg1
    seg[c(FALSE, TRUE)] <- seg2

    valid <- grepl("^[0-9]*M$", bam$cigar)
    valid <- valid & sub("ref_|alt_", "", bam$rname) == as.character(seqnames(seg))
    valid <- valid & bam$pos == start(seg)
    seg$qname <- bam$qname
    seg$seq <- bam$seq
    seg$tag <- rep(seq_len(length(seg)/2), each = 2)
    seg <- seg[valid]
    return(seg)
}

################################################################################
#' @export
#'
#' @importFrom Rsubread align buildindex
#' @importFrom Rsamtools scanBamFlag ScanBamParam scanBam
#' @importFrom BiocGenerics width
#'
mergeGenome <- function(snps, ref_fn, alt_fn, prefix, indexing = TRUE){
    stopifnot(inherits(snps, "SNPs"))
    ref_genome <- readDNAStringSet(ref_fn)
    alt_genome <- readDNAStringSet(alt_fn)
    rgn <- names(ref_genome)
    agn <- names(alt_genome)
    rsn <- seqlevels(snps$ref)
    asn <- seqlevels(snps$alt)
    if(!all(rsn %in% rgn)){
        rsn_i <- as.numeric(sub(".*[^0-9]", "", sub("[^0-9].*", "", rsn)))
        rsn <- paste0("chr", sprintf("%02d", rsn_i))
        rgn_i <- as.numeric(sub(".*[^0-9]", "", sub("[^0-9].*", "", rgn)))
        rgn <- paste0("chr", sprintf("%02d", rgn_i))
    }
    if(!all(asn %in% agn)){
        asn_i <- as.numeric(gsub(".*[^0-9]", "", gsub("[^0-9].*", "", asn)))
        asn <- paste0("chr", sprintf("%02d", asn_i))
        agn_i <- as.numeric(gsub(".*[^0-9]", "", gsub("[^0-9].*", "", agn)))
        agn <- paste0("chr", sprintf("%02d", agn_i))
    }

    names(ref_genome) <- paste0("ref_", rgn)
    names(alt_genome) <- paste0("alt_", agn)
    seqlevels(snps$ref) <- paste0("ref_", rsn)
    seqlevels(snps$alt) <- paste0("alt_", asn)
    genome <- c(ref_genome, alt_genome)
    genome_fn <- paste0(prefix, "merge.fa.gz")
    writeXStringSet(genome, genome_fn,
                    compress = TRUE, format = "fasta")
    index_fn <- paste0(prefix, "merge")
    if(indexing){
        buildindex(index_fn, genome_fn)

    } else {
        check <- file.exists(paste0(index_fn, ".files"))
        if(!check){
            warning(index_fn, " does not exist.\n Need inexing,")
        }
    }
    out <- list(snps = snps, index_fn = index_fn)
    class(out) <- c(class(out), "GSNPs")
    return(out)
}

################################################################################
# Count reads
#'
#'
#' @export
#'
doSInG <- function(object, fq_list, prefix, n_threads, vcf, rnaseq = FALSE){
    if(inherits(object, "TAG")  | inherits(object, "GSNPs")){
        count_list <- .count.reads(object, fq_list, n_threads, rnaseq)
    } else {
        stop("object should be a TAG or GSNPs class object.")
    }
    colnames(count_list) <- fq_list[, ncol(fq_list)]

    message("Creating a GDS file as output...")
    if(inherits(object, "GSNPs")){
        object <- object$snps
    }
    out_fn <- .create.gds(object, count_list, prefix, vcf)
    return(out_fn)
}

.count.reads <- function(object, fq_list, n_threads, rnaseq){
    count_list <- NULL
    for(i in seq_len(nrow(fq_list))){
        count_list <- cbind(count_list,
                            .getCounts(object, fq_list[i, ], n_threads, rnaseq))
    }
    return(count_list)
}

#'
#' @importFrom Rsubread align
#' @importFrom Rsamtools scanBamFlag ScanBamParam scanBam
#' @importFrom BiocGenerics width
#'
.getCounts <- function(object, fq_fn, n_threads, rnaseq){
    dir <- getwd()
    aln_dir <- paste(dir, "aln", sep = "/")
    dir.create(aln_dir, showWarnings = FALSE)

    if(inherits(object, "TAG")){
        count <- .alignToTag(object, fq_fn, n_threads, aln_dir)
    } else {
        count <- .alignToGenome(object, fq_fn, n_threads, aln_dir, rnaseq)
    }

    return(count)
}

.alignToTag <- function(tags, fq_fn, n_threads, aln_dir){
    if(length(fq_fn) == 2){
        bam_fn <- paste0(aln_dir, "/", fq_fn[2], ".bam")
        aln <- align(tags$index_fn, fq_fn[1], type = "dna",
                     output_file = bam_fn,
                     nthreads = n_threads, maxMismatches = 0, unique = FALSE,
                     nBestLocations = 2, indels = 0)

    } else {
        bam_fn <- paste0(aln_dir, "/", fq_fn[3], ".bam")
        aln <- align(tags$index_fn, fq_fn[1], fq_fn[2], type = "dna",
                     output_file = bam_fn,
                     nthreads = n_threads, maxMismatches = 0, unique = FALSE,
                     nBestLocations = 2, indels = 0)
    }

    param <- ScanBamParam(scanBamFlag(isUnmappedQuery = FALSE),
                          simpleCigar = TRUE, what = c("rname", "cigar"))
    bam <- scanBam(bam_fn, param = param)[[1]]
    bam$cigar <- as.numeric(sub("M", "", bam$cigar))
    tag_w <- c(width(tags$ref), width(tags$alt))
    id <- c(tags$ref$id, tags$alt$id)
    map <- match(bam$rname, id)
    full_map <- bam$cigar == 150
    full_map <- full_map | bam$cigar == tag_w[map]
    bam <- bam$rname[full_map]
    return(table(bam))
}

#' @importFrom GenomicAlignments readGAlignments cigar
#' @importFrom Rsubread align
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges grglist findOverlaps
#' @importFrom S4Vectors subjectHits
.alignToGenome <- function(gsnps, fq_fn, n_threads, aln_dir, rnaseq){
    if(length(fq_fn) == 2){
        bam_fn <- paste0(aln_dir, "/", fq_fn[2], ".bam")
        aln <- align(gsnps$index_fn, fq_fn[1], type = "dna",
                     output_file = bam_fn,
                     nthreads = n_threads, maxMismatches = 0, unique = TRUE,
                     indels = 0)

    } else {
        bam_fn <- paste0(aln_dir, "/", fq_fn[3], ".bam")
        aln <- align(gsnps$index_fn, fq_fn[1], fq_fn[2], type = "dna",
                     output_file = bam_fn,
                     nthreads = n_threads, maxMismatches = 0, unique = TRUE,
                     indels = 0)
    }

    ga <- readGAlignments(file = bam_fn)
    check <- all(seqlevels(ga) %in% c(seqlevels(gsnps$snps$ref),
                                      seqlevels(gsnps$snps$alt)))
    if(!check){
        stop('Chromosome names in the given BAM file and the mummer files ',
             'do not match.', '\nChromosome names in the BAM file:\n',
             paste(seqlevels(ga), collapse = "\ "),
             '\nChromosome names in the mummer files:\n',
             paste(c(seqlevels(gsnps$snps$ref), seqlevels(gsnps$snps$alt)),
                   collapse = "\ "))
    }

    if(rnaseq){
        incomple <- grepl("S|I|D|H|P|X", cigar(ga))
        ga <- ga[!incomple]
        ga <- grglist(ga, order.as.in.query=TRUE)
        ga <- unlist(ga)

    } else {
        incomple <- grepl("I|D|H|P|X|N", cigar(ga))
        ga <- ga[!incomple]
    }

    ref_count <- table(factor(subjectHits(findOverlaps(ga, gsnps$snps$ref)),
                              seq_along(gsnps$snps$ref)))
    alt_count <- table(factor(subjectHits(findOverlaps(ga, gsnps$snps$alt)),
                              seq_along(gsnps$snps$alt)))
    count <- c(ref_count, alt_count)
    return(count)
}


#'
#' @importFrom SNPRelate snpgdsCreateGeno
#' @importFrom SeqArray seqSNP2GDS seqGDS2VCF
#' @importFrom BiocGenerics start
.create.gds <- function(object, count_list, prefix, vcf){
    tmp_out <- tempfile(tmpdir = tempdir(), fileext = ".snp")
    on.exit({unlink(tmp_out)})
    genmat <- .make.geno(count_list)
    data <- .make.ADdata(count_list)
    n_mar <- length(object$ref)

    snpgdsCreateGeno(tmp_out,
                     genmat,
                     colnames(count_list),
                     seq_len(n_mar),
                     seq_len(n_mar),
                     sub("ref_", "", as.character(seqnames(object$ref))),
                     start(object$ref),
                     rep("A/G", n_mar),
                     FALSE,
                     "",
                     "",
                     list(alt.pos = start(object$alt)))

    out_fn <- seqSNP2GDS(tmp_out, paste0(prefix, ".gds"), verbose = TRUE)
    .addAD(data, out_fn)

    if(vcf){
        out_fn <- c(out_fn,
                    seqGDS2VCF(out_fn, paste0(prefix, ".vcf"), verbose = FALSE,
                               use_Rsamtools = TRUE))
    }
    return(out_fn)
}

.make.geno <- function(count_list){
    n_row <- nrow(count_list)/2
    n_index <- seq_len(n_row)
    genmat <- matrix(0, n_row, ncol(count_list))
    genmat[count_list[n_index, ] == 0 & count_list[-n_index, ] == 0] <- 3
    genmat[count_list[n_index, ] > 0 & count_list[-n_index, ] > 0] <- 1
    genmat[count_list[n_index, ] > 0 & count_list[-n_index, ] == 0] <- 2
    return(t(genmat))
}

.make.ADdata <- function(count_list){
    n_row <- nrow(count_list)/2
    n_index <- seq_len(n_row)
    data <- matrix(0, ncol(count_list), n_row * 2)
    data[, c(T, F)] <- t(count_list[n_index, ])
    data[, c(F, T)] <- t(count_list[-n_index, ])
    return(data)
}

#' @importFrom methods new
#' @importFrom gdsfmt index.gdsn openfn.gds put.attr.gdsn addfolder.gdsn add.gdsn closefn.gds
.addAD <- function(data, out_fn){
    gds <- openfn.gds(out_fn, FALSE)

    addfolder.gdsn(index.gdsn(gds, "annotation/format"), "AD")
    add.gdsn(index.gdsn(gds, "annotation/format/AD"),
             "data", data, storage = "vl_int", valdim = dim(data),
             compress = "LZMA_RA")
    add.gdsn(index.gdsn(gds, "annotation/format/AD"),
             "@data", rep(2, ncol(data)/2), storage = "int32",
             compress = "LZMA_RA", visible = FALSE)
    put.attr.gdsn(index.gdsn(gds, "annotation/format/AD"),
                  "Number", "R")
    put.attr.gdsn(index.gdsn(gds, "annotation/format/AD"),
                  "Type", "Integer")
    put.attr.gdsn(index.gdsn(gds, "annotation/format/AD"),
                  "Description",
                  "Allelic depths for the reference and alternate alleles in the order listed")
    closefn.gds(gds)
}

#'
#'
#'
#' @importFrom Rsamtools pileup PileupParam
#' @importFrom dplyr left_join
#'
#' @export
#'
doGBS <- function(ref_fn, fq_list, prefix, n_threads, vcf, indexing){
    index_fn <- sub("\\.fas*t*a*\\.*g*z*", "", ref_fn)
    if(indexing){
        buildindex(index_fn, ref_fn)

    } else {
        check <- file.exists(paste0(index_fn, ".files"))
        if(!check){
            warning(index_fn, " does not exist.\n Need inexing,")
        }
    }

    dir <- getwd()
    aln_dir <- paste(dir, "aln", sep = "/")
    dir.create(aln_dir, showWarnings = FALSE)
    aln_dir <- paste(dir, "vcf", sep = "/")
    dir.create(aln_dir, showWarnings = FALSE)

    count_list <- NULL
    for(i in seq_len(nrow(fq_list))){
        if(ncol(fq_list) == 2){
            bam_fn <- paste0(aln_dir, "/", fq_list[i, 2], ".bam")
            aln <- align(index_fn, fq_list[i, 1], type = "dna",
                         output_file = bam_fn, sortReadsByCoordinates = TRUE,
                         readGroup = fq_list[i, 2], readGroupID = fq_list[i, 2],
                         nthreads = n_threads, unique = TRUE)

        } else {
            bam_fn <- paste0(aln_dir, "/", fq_list[i, 3], ".bam")
            aln <- align(index_fn, fq_list[i, 1], fq_list[i, 2], type = "dna",
                         output_file = bam_fn, sortReadsByCoordinates = TRUE,
                         readGroup = fq_list[i, 3], readGroupID = fq_list[i, 3],
                         nthreads = n_threads, unique = TRUE)
        }
        pp <- PileupParam(min_mapq = 20, distinguish_strands = FALSE,
                          include_insertions = TRUE)
        bai_fn <- paste0(bam_fn, ".bai")
        count_list <- c(count_list,
                        list(pileup(bam_fn, bai_fn, pileupParam = pp)))
    }

    by_cols <- c("seqnames", "pos", "nucleotide")
    for(i in seq_along(count_list)){
        if(i == 1){
            df <- count_list[[i]]
        } else {
            df <- left_join(df, count_list[[i]], by_cols)
        }
    }
    df <- df[order(df$seqnames, df$pos), ]
    dup <- duplicated(df[, c("seqnames", "pos")])
    variant <- which(dup)
    variant <- sort(unique(c(variant, variant - 1)))

    df <- df[variant, ]
    chr <- df[c(TRUE, FALSE), "seqnames"]
    pos <- df[c(TRUE, FALSE), "pos"]
    ref_allele <- df[c(TRUE, FALSE), "nucleotide"]
    alt_allele <- df[c(FALSE, TRUE), "nucleotide"]
    ref <- t(df[c(TRUE, FALSE), -(1:3)])
    alt <- t(df[c(FALSE, TRUE), -(1:3)])
    n_mar <- ncol(ref)
    n_sam <- nrow(ref)
    genmat <- matrix(0, n_sam, n_mar)
    genmat[ref == 0 & alt == 0] <- 3
    genmat[ref > 0 & alt > 0] <- 1
    genmat[ref > 0 & alt == 0] <- 2
    ad <- matrix(0, n_sam, n_mar * 2)
    ad[, c(T, F)] <- ref
    ad[, c(F, T)] <- alt

    tmp_out <- tempfile(tmpdir = tempdir(), fileext = ".snp")
    on.exit({unlink(tmp_out)})
    snpgdsCreateGeno(tmp_out,
                     genmat,
                     fq_list[, ncol(fq_list)],
                     seq_len(n_mar),
                     seq_len(n_mar),
                     chr,
                     pos,
                     paste(ref_allele, alt_allele, sep = "/"),
                     FALSE,
                     "",
                     "")

    out_fn <- seqSNP2GDS(tmp_out, paste0(prefix, ".gds"), verbose = FALSE)
    .addAD(ad, out_fn)

    if(vcf){
        out_fn <- c(out_fn,
                    seqGDS2VCF(out_fn, paste0(prefix, ".vcf"), verbose = FALSE,
                               use_Rsamtools = TRUE))
    }
    return(out_fn)
}
