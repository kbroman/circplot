# convert map to x-axis positions for plot_scan1
map_to_xpos <-
    function(map, gap=NULL)
{
    if(is.null(gap)) sum(qtl2::chr_lengths(map))/100

    if(length(map)==1) return(map[[1]])

    chr_range <- vapply(map, range, c(0,1), na.rm=TRUE)

    result <- map[[1]]-chr_range[1,1] + gap/2
    for(i in 2:length(map)) {
        result <- c(result,
                    map[[i]] - chr_range[1,i] + gap + max(result, na.rm=TRUE))
    }
    result
}

# boundaries of chromosomes in plot_scan1
# first row: left edges
# second row: right edges
map_to_boundaries <-
    function(map, gap=NULL)
{
    if(is.null(gap)) gap <- sum(qtl2::chr_lengths(map))/100

    if(length(map)==1)
        return(cbind(range(map[[1]], na.rm=TRUE)))

    # range of each chromosome
    chr_range <- lapply(map, range, na.rm=TRUE)

    # corresponding xpos, as matrix with two rows
    startend <- matrix(map_to_xpos(chr_range, gap), nrow=2)

    startend[1,] <- startend[1,] - gap/2
    startend[2,] <- startend[2,] + gap/2

    startend
}

# convert map to list of indexes to LOD vector
map_to_index <-
    function(map)
{
    if(length(map)==1) {
        map[[1]] <- seq(along=map[[1]])
        return(map)
    }

    lengths <- vapply(map, length, 0)
    split(seq_len(sum(lengths)), rep(seq(along=map), lengths))
}
