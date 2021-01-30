
#' helper for converting ENSG to symbols, cbinding symbol and chromosome name for the genes to the input data.frame-like entity
#' @param dflike a data.frame-like entity with a column of ENSG gene identifiers; column
#' name must be the value of `bind_sym` parameter `ensgvbl`
#' @param gene_resource a data.frame with columns `gene_name` and `seqnames`,
#' defaulting to `gwaslake::genes_v75`
#' @param ensgvbl character(1) the column of `dflike` holding the ENSG ids
#' @param output_symvbl character(1) name to be used for symbol in output data.frame, defaults to `symbol`
#' @param output_chromvbl character(1) name to be used for chromosome in output data.frame, defaults to `gchrom`
#' @examples
#' head(little_eq)
#' head(bind_sym(little_eq))
#' @export
bind_sym = function( dflike, gene_resource = genes_v75, ensgvbl = "trait",
   output_symvbl = "symbol", output_chromvbl = "gchrom" ) {
 stopifnot(ensgvbl %in% names(dflike))
 stopifnot(all(c("gene_name", "seqnames") %in% names(gene_resource)))
 dflike[[output_symvbl]] = gene_resource[ dflike[[ ensgvbl ]], "gene_name" ] 
 dflike[[output_chromvbl]] = as.character(factor(gene_resource[ dflike[[ ensgvbl ]], "seqnames" ] ))
 dflike
}
