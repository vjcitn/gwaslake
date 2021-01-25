
#' IEU tables conversion to GRanges
#' @param ieutab data.frame-like entity as returned by e.g., variant_rsid, associations, phewas
#' @param build character(1) defaults to "GRCh37"
#' @note assumes SNP; renames 'pos' to 'position' if necessary
#' @examples
#' if (requireNamespace("ieugwasr")) {
#' ex = ieugwasr::variants_rsid("rs6060535")
#' ieutab_to_GRanges(ex)
#' }
#' @export
ieutab_to_GRanges = function(ieutab, build="GRCh37") {
 stopifnot("chr" %in% names(ieutab))
 if (!requireNamespace("GenomicRanges")) stop("install GenomicRanges to use this function")
 if (!requireNamespace("S4Vectors")) stop("install S4Vectors to use this function")
 if (!requireNamespace("GenomeInfoDb")) stop("install GenomeInfoDb to use this function")
 if ("pos" %in% names(ieutab)) {
    ind = which(names(ieutab)=="pos")
    names(ieutab)[ind] = "position"
    }
 tmp = GenomicRanges::GRanges(ieutab$chr, IRanges::IRanges(ieutab$position, width=1))
 S4Vectors::mcols(tmp) = ieutab
 GenomeInfoDb::genome(tmp) = build
 tmp
}
