
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

#' create a "chr:start-end" string for associations()
#' @param genesym character(1)
#' @param radius numeric(1) flanking range, defaults to zero
#' @param addr_src data.frame defaults to `genes_v75`
#' @examples
#' genesym_to_string("YY1")
#' @export
genesym_to_string = function(genesym, radius=0, addr_src = gwaslake::genes_v75) {
  tmp = addr_src[which(addr_src$gene_name==genesym),]
  if (nrow(tmp)!=1) stop("no unique address available")
  paste0(tmp$seqnames, ":", tmp$start, "-", tmp$end)
}
