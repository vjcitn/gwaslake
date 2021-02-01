#' produce a manhattan plot for a given study, based on specific coordinates or gene symbol with radius for flanking region
#' @param studyID character(1) id of study in ieugwasr::gwasinfo() catalog, defaults to "ubm-a-524", which is
#' a study described in Elliott et al Nature October 2018
#' @param symbol a gene symbol (to be used with `genesym_to_string`, implicitly); if absent, chromosome, start and end must be supplied, defaults to 'BCAN'
#' @param start chromosomal address to start region (will have radius subtracted if supplied)
#' @param stop chromosomal address to end region (will have radius added if supplied)
#' @param chromosome name of chromosome in Ensembl nomenclature (e.g., 1, 2, ...)
#' @param radius numeric default to zero, specifies flanking region
#' @examples
#' manhattanPlot(radius=500000)
#' @export
manhattanPlot = function (studyID="ubm-a-524", symbol = "BCAN", start = NULL, stop = NULL, 
    chromosome = NULL, radius=0) 
{
    if (is.character(symbol)) {
        cur = strsplitter(genesym_to_string(symbol))
        start = cur["start"]
        stop = cur["end"]
        chromosome = cur["chr"]
    }
    if (is.null(start)) 
        stop("neither symbol nor start provided")
    location = paste(as.character(chromosome), paste(as.character(as.numeric(start)-radius), 
        as.character(as.numeric(stop)+radius), sep = "-"), sep = ":")
    mplot = ieugwasr::associations(location, studyID)
    ggplot2::ggplot(mplot, ggplot2::aes(x = position, y = -log10(p))) + ggplot2::geom_point()
}
