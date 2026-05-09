#' plot a genome scan as a circle
#'
#' Plot LOD curves for a genome scan as a circle
#'
#' @param x Output of [qtl2::scan1()]
#' @param map A list of vectors of marker positions, as produced by
#' [qtl2::insert_pseudomarkers()].
#' @param lodcolumn LOD score column to plot (a numeric index, or a
#' character string for a column name). Only one value allowed.
#' @param chr Selected chromosomes to plot; a vector of character
#' strings.
#' @param gap Gap between chromosomes. The default is 1% of the total genome length.
#' @param rlim radius limits
#' @param start_angle angle to start
#' @param clockwise If true, go clockwise
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param chr_labels Whether to add chromosome labels
#' @param lwd line width for LOD curves
#' @param col color for LOD curves
#' @param lty line type for LOD curves
#' @param add Add curves to previous plot (need care about matching `rlim` and `ylim`)
#' @param ... Additional graphics parameters
#'
#' @return No return value
#'
#' @importFrom qtl2 subset_scan1 chr_lengths align_scan1_map
#' @importFrom graphics par text lines
#' @export
#'
#' @examples
#' library(qtl2)
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#' kinship <- calc_kinship(probs, "loco")
#' out_lmm <- scan1(probs, pheno, kinship, covar, Xcovar)
#'
#' max_lod <- max(c(maxlod(out), maxlod(out_lmm)))
#' plot_scan1_circ(out, map, ylim=c(0, max_lod), rlim=c(5, 7))
#' plot_scan1_circ(out_lmm, map, col="violetred", lty=2, add=TRUE,
#'                 ylim=c(0, max_lod), rlim=c(5, 7))

plot_scan1_circ <-
    function(x, map, lodcolumn=1, chr=NULL, gap=NULL, rlim=c(5,6),
             start_angle=pi, clockwise=TRUE, xlim=NULL, ylim=NULL,
             chr_labels=TRUE, lwd=2, col="slateblue", lty=1, add=FALSE, ...)
{

    if(is.null(map)) stop("map is NULL")
    if(!is.list(map)) map <- list(" "=map) # if a vector, treat it as a list with no names

    if(!is.matrix(x) && !is.data.frame(x) && is.numeric(x)) {
        x <- as.matrix(x)
    }

    # subset chromosomes
    if(!is.null(chr)) {
        chri <- match(chr, names(map))
        if(any(is.na(chri)))
            stop("Chromosomes ", paste(chr[is.na(chri)], collapse=", "), " not found")
        x <- qtl2::subset_scan1(x, map, chr)
        map <- map[chri]
    }

    if(is.null(gap)) gap <- sum(qtl2::chr_lengths(map))/100
    if(length(gap) != 1 || gap < 0) stop("gap should be a single non-negative number")

    tmp <- qtl2::align_scan1_map(x, map)
    x <- tmp$scan1
    map <- tmp$map

    if(length(lodcolumn)==0) stop("lodcolumn has length 0")
    if(length(lodcolumn) > 1) { # If length > 1, take first value
        warning("lodcolumn should have length 1; only first element used.")
        lodcolumn <- lodcolumn[1]
    }
    if(is.character(lodcolumn)) { # turn column name into integer
        tmp <- match(lodcolumn, colnames(x))
        if(is.na(tmp)) stop('lodcolumn "', lodcolumn, '" not found')
        lodcolumn <- tmp
    }
    if(lodcolumn < 1 || lodcolumn > ncol(x))
        stop("lodcolumn [", lodcolumn, "] out of range (should be in 1, ..., ", ncol(x), ")")
    lod <- unclass(x)[,lodcolumn]

    xpos <- map_to_xpos(map, gap)
    chrbound <- map_to_boundaries(map, gap)
    indexes <- map_to_index(map)

    if(is.null(xlim)) xlim <- c(min(xpos), max(xpos)+gap)
    if(is.null(ylim)) ylim <- c(0, max(lod, na.rm=TRUE))

    pts <- xy2circ(xpos, lod, xlim=xlim, ylim=ylim,
                   rlim=rlim, start_angle=start_angle, clockwise=clockwise)
    pts0 <- xy2circ(xpos, rep(0, length(lod)), xlim=xlim, ylim=ylim,
                    rlim=rlim, start_angle=start_angle, clockwise=clockwise)

    plot_scan1_circ_internal <-
        function(pts, type="p", xlim=NULL, ylim=NULL, pty="s", bty="n",
                 xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="", lty=1,
                 lwd=2, col="slateblue", col_axis="black", lwd_axis=1, rlim=c(5,6),
                 add=FALSE, ...)
    {
        par(pty=pty, bty=bty)
        xl <- c(-rlim[2], rlim[2])*1.05

#        cat("add:", add)
        if(!add) plot(pts, type="n", xlim=xl, ylim=xl, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab="")

        for(chr in seq_along(map)) {
            if(!add) lines(pts0[indexes[[chr]], ,drop=FALSE], col=col_axis, lwd=lwd_axis)
            lines(pts[indexes[[chr]],,drop=FALSE], col=col, lwd=lwd, lty=lty)

            if(!add & chr_labels) {
                label_pos <- mean(range(xpos[indexes[[chr]]]))
                label_pos <- xy2circ(label_pos, -ylim[2]/5, xlim=xlim, ylim = -ylim,
                                     rlim=c(rlim[1], rlim[1]-(rlim[2]-rlim[1])), start_angle=start_angle,
                                     clockwise=clockwise)
            }

            if(!add) text(label_pos[1], label_pos[2], names(map)[chr])
        }
        invisible()
    }

    plot_scan1_circ_internal(pts, lwd=lwd, col=col, rlim=rlim, add=add, xlim=xlim, ylim=ylim,
                             xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",
                             lty=lty, ...)
}
