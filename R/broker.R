makeUCSC = function(x) {
  haschr = grep("^chr", x)
  if (length(haschr)==length(x)) return(x)
  if (length(haschr)==0) return(paste0("chr", x))
  stop("x is inconsistent NCBI/UCSC annotation for chromosome.")
}
#' format the variant location, content and identity information from MRC IEU API tables for OpenCRAVAT
#' @param ieutab a data.frame-like entity with `chr`, `position`, `ea` (to be regarded
#' @param chain a liftOver chain structure as defined in rtracklayer
#' @param genome character(1) label, if present, attached to output as a column `build`
#' as ALT allele and `nea` to be regarded as REF allele
#' @note If a liftOver chain is provided, some input records may be lost.  We attempt
#' to keep allele content information with each locus.  
#' @examples
#' # Here we assume it is desired to give hg38 variant positions to cravat
#' if (requireNamespace("ieugwasr")) {
#'   pw = try(ieugwasr::phewas("17:37000000-37020000"))
#'   if (!inherits(pw, "try-error")) {
#'     octab = ieutab_to_cravat(pw, chain=gwaslake::ch19to38, genome="hg38" )
#'     print(head(octab))
#'     }
#'   aa = try(ieugwasr::associations("rs6060535", id="ukb-b-10787"))
#'   if (!inherits(aa, "try-error")) {
#'     octab2 = ieutab_to_cravat(aa, chain=gwaslake::ch19to38, genome="hg38")
#'     head(octab2)
#'     }
#' }
#' @export
ieutab_to_cravat = function(ieutab, chain=NULL, genome=NULL) {
  stopifnot(all(c("chr", "position", "ea", "nea") %in% names(ieutab)))
  if (!requireNamespace("GenomicRanges")) stop("install GenomicRanges to use this function")
  if (!requireNamespace("IRanges")) stop("install IRanges to use this function")
  if (!requireNamespace("rtracklayer")) stop("install rtracklayer to use this function")
  ac = as.character
  if (!is.null(chain)) {
    if (is.null(genome)) stop("if chain is supplied, genome cannot be NULL")
    tmp = GenomicRanges::GRanges(makeUCSC(ieutab$chr), IRanges::IRanges(ieutab$position, width=1),
         rsid=ieutab$rsid, trait=ieutab$trait)
    tmp = unlist(rtracklayer::liftOver(tmp, chain))
    w = GenomicRanges::width(tmp)
    if (any(w>1)) tmp = tmp[-which(w>1)]
    tmp = as.data.frame(tmp)
    ieutab = merge(ieutab, tmp, by=c("rsid", "trait"))
    ieutab$position = ieutab$start
    }
  ans = data.frame(chr=ac(ieutab$chr), pos=ieutab$position, ref=ieutab$nea,
    alt=ieutab$ea, var=ieutab$rsid, samp="IEU")
  if (!is.null(genome)) ans$build = genome
  ans
}
