#' Calculate read ratio and smoothed read ratio.
#'　
#' Calculate read ratio at each interval based on normalized read counts
#' obtained via [doSInG()]. Read ratios calculated at each interval
#' usually contain data points with stochastic error which are
#' represented as unexpected changes in read ratio between adjacent
#' intervals. To remove those errorneous intervals, this function
#' also calculate average read ratio within the specified window.
#'
#' @param object A `SInG` object after running [doSInG()].
#' @param window A integer value to specify the window size in base pairs.
#' @param min_marker A integer to indicate the minimum number of markers
#' with in each window. If the number of markers within a window with the
#' size `window` is less than `min_marker`, use `round(min_markers/2)`
#' markers at upstream and downstream from the given marker.
#'
#' @return A `SInG` object with smoothed read ratios.
#'
#' @importFrom GenomicRanges start end
#' @importFrom GenomeInfoDb seqnames
#'
#' @examples
#' # Load a sample SInG object.
#' sing_object <- system.file("extdata", "sample.Rdata", package = "SInG")
#' load(sing_object)
#'
#' object <- calcRatio(object, window = 5000)
#'
#' @export
#'
calcRatio <- function(object, window = 10000, min_marker = 0) {

    ratio_mat <- NULL
    for(i in seq_along(object$counts)){
        counts <- object$counts[[i]]$sin
        ratio <- counts$ref_norm / rowSums(counts[, c("ref_norm", "alt_norm")])
        object$counts[[i]]$ratio <- ratio
        ratio_mat <- rbind(ratio_mat, ratio)
    }

    chr <- as.numeric(seqnames(object$varInt$ref_vs))
    pos <- start(object$varInt$ref_vs) + round((end(object$varInt$ref_vs) - start(object$varInt$ref_vs)) / 2)
    smoothed <- .smoothRatio(ratio_mat, chr, pos, window, min_marker)
    for(i in seq_along(object$counts)){
        object$counts[[i]]$smoothed <- smoothed[i, ]
    }

    return(object)
}

.smoothRatio <- function(ratio, chr, pos, window, min_marker){
    smoothed <- NULL
    for (chr_i in unique(chr)) {
        message(paste0('Now smoothing chr ', chr_i))
        target <- chr == chr_i
        pos_chr <- pos[target]
        ratio_chr <- ratio[, target]
        for (i in 1:length(pos_chr)) {
            tmp_pos <- pos_chr[i]
            in_window <- pos_chr >= tmp_pos - window & pos_chr <= tmp_pos + window
            if(min_marker > 0){
                if (sum(in_window) < min_markers) {
                    in_window <-
                        which(1:length(pos_chr) >= i - round(min_markers/2) &
                                  1:length(pos_chr) <= i + round(min_markers/2))
                    if (length(in_window) < min_markers) {
                        if (1 %in% in_window) {
                            in_window <- seq(in_window[1],
                                             by = 1,
                                             length.out = min_markers)
                        } else if (length(pos_chr) %in% in_window) {
                            in_window <- seq(tail(in_window, 1),
                                             by = -1,
                                             length.out = min_markers)
                        }
                    }
                }
            }
            if(sum(in_window) == 1){
                smoothed <- cbind(smoothed, ratio_chr[, in_window])
            } else {
                smoothed <- cbind(smoothed, rowMeans(ratio_chr[, in_window], na.rm = T))
            }
        }
    }
    return(smoothed)
}
