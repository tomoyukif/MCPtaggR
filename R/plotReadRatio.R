#' Draw line plot of read ratio per interval
#'
#' Draw line plot(s) of read ratio per interval for samples.
#'
#' @param object A SInG object after running [calcRaio()].
#' @param ind An integer vector to specify indices of samples to be
#' visualized. If multiple integers were specified, lines plots for
#' the samples will be overlaid with different colors.
#' @param smooth If TRUE, as default, used smoothed read ratio.
#' Otherwise, raw read ratio.
#' @param genome A string of either "ref" or "alt" to indicate which
#' interval position information should be used.
#' @param smooth_only A logical value to plot only lines of smoothed ratios
#' or both lines of smoothed ratios and dots of raw ratios.
#' @param coord A vector with two integer to specify the number of rows and
#' columns of tiled line plots for each chromosome.
#' @param dot_color A string to specify the color of dots in the plot.
#' @param dot_size A numeric value to specify the size of dots in the plot.
#' @param dot_alpha A numeric value to specify the level of transparency
#' of dots in the plot.
#' @param line_color A string to specify the color of lines in the plot.
#'  @param line_size A numeric value to specify the width of lines in the plot.
#' @param line_alpha A numeric value to specify the level of transparency
#' of lines in the plot.
#' @param font_size A numeric value to specify font size in the plot.
#' @param ggargs A string to modify ggplot.
#'
#' @return Invisibly returns a ggplot object.
#' @import ggplot2
#' @importFrom GenomicRanges start end
#' @importFrom GenomeInfoDb seqnames
#'
#' @examples
#' # Load a sample SInG object.
#' sing_object <- system.file("extdata", "sample.Rdata", package = "SInG")
#' load(sing_object)
#'
#' plotReadRatio(object, ind = 1)
#'
#' @export
#'
#'
plotReadRatio <- function(object,
                          ind = NULL,
                          genome = "ref",
                          smooth_only = FALSE,
                          coord = c(NA, 2),
                          dot_color = "green",
                          dot_size = 0.5,
                          dot_alpha = 0.3,
                          line_color = "magenta",
                          line_size = 0.5,
                          line_alpha = 0.7,
                          font_size = 1,
                          ggards = NULL) {

    if(is.null(ind)){
        ind <- seq_along(object$counts)
    }

    if(genome == "alt"){
        chr <- as.numeric(seqnames(object$varInt$alt_vs))
        start_pos <- start(object$varInt$alt_vs) * 10^-6
        end_pos <- end(object$varInt$alt_vs) * 10^-6
        center_pos <- start_pos + round((end_pos - start_pos) / 2)

    } else {
        chr <- as.numeric(seqnames(object$varInt$ref_vs))
        start_pos <- start(object$varInt$ref_vs) * 10^-6
        end_pos <- end(object$varInt$ref_vs) * 10^-6
        center_pos <- start_pos + round((end_pos - start_pos) / 2)
    }

    df_ratio <- NULL
    df_smooth <- NULL
    for(i in seq_along(object$counts)){
        if(i %in% ind){
            if(is.null(df_ratio)){
                df_ratio <- data.frame(chr = chr,
                                       pos = center_pos,
                                       ratio = object$counts[[i]]$ratio,
                                       Sample = names(object$counts)[i])
                df_smooth <- rbind(data.frame(chr = chr,
                                              pos = start_pos,
                                              ratio = object$counts[[i]]$smoothed,
                                              Sample = names(object$counts)[i]),
                                   data.frame(chr = chr,
                                              pos = end_pos,
                                              ratio = object$counts[[i]]$smoothed,
                                              Sample = names(object$counts)[i]))
            } else {
                df_ratio <- rbind(df_ratio,
                                  data.frame(chr = chr,
                                             pos = center_pos,
                                             ratio = object$counts[[i]]$ratio,
                                             Sample = names(object$counts)[i]))
                df_smooth <- rbind(df_smooth,
                                   rbind(data.frame(chr = chr,
                                                    pos = start_pos,
                                                    ratio = object$counts[[i]]$smoothed,
                                                    Sample = names(object$counts)[i]),
                                         data.frame(chr = chr,
                                                    pos = end_pos,
                                                    ratio = object$counts[[i]]$smoothed,
                                                    Sample = names(object$counts)[i])))
            }
        }
    }

    if(all(is.na(coord))){
        coord[1] <- round(sqrt(length(unique(chr))))
        coord[2] <- ceiling(length(unique(chr)) / coord[1])
    } else if(any(is.na(coord))){
        non_na <- coord[!is.na(coord)]
        coord[is.na(coord)] <- ceiling(length(unique(chr)) / non_na)
    }

    p <- ggplot()

    if(length(ind) == 1){
        if(smooth_only){
            p <- p + geom_line(mapping = aes(x = pos,
                                             y = ratio,
                                             group = Sample),
                               data = df_smooth,
                               size = line_size,
                               alpha = line_alpha,
                               color = line_color)

        } else {
            p <- p +
                geom_point(mapping = aes(x = pos,
                                         y = ratio,
                                         group = Sample),
                           data = df_ratio,
                           size = dot_size,
                           alpha = dot_alpha,
                           color = dot_color,
                           stroke = 0) +
                geom_line(mapping = aes(x = pos,
                                        y = ratio,
                                        group = Sample),
                          data = df_smooth,
                          size = line_size,
                          alpha = line_alpha,
                          color = line_color)
        }
    } else {
        if(smooth_only){
            p <- p + geom_line(mapping = aes(x = pos,
                                             y = ratio,
                                             group = Sample,
                                             color = Sample),
                               data = df_smooth,
                               size = line_size,
                               alpha = line_alpha)

        } else {
            p <- p +
                geom_point(mapping = aes(x = pos,
                                         y = ratio,
                                         group = Sample,
                                         color = Sample),
                           data = df_ratio,
                           size = dot_size,
                           alpha = dot_alpha,
                           stroke = 0) +
                geom_line(mapping = aes(x = pos,
                                        y = ratio,
                                        group = Sample,
                                        color = Sample),
                          data = df_smooth,
                          size = line_size,
                          alpha = line_alpha)
        }
    }

    p <- p + ylim(0, 1) +
        xlab("Position (Mb) in the reference genome") +
        ylab("Reference allele read ratio") +
        facet_wrap(facets = ~ chr,
                   nrow = coord[1],
                   ncol = coord[2],
                   scales = "free_x",
                   dir = "v",
                   strip.position = "right") +
        theme(axis.title = element_text(size = rel(font_size)),
              axis.text = element_text(size = rel(font_size)),
              plot.title = element_text(size = rel(font_size)),
              legend.text = element_text(size = rel(font_size)),
              legend.title = element_text(size = rel(font_size)),
              strip.text = element_text(size = rel(font_size)))
    print(p)
    invisible(p)
}
