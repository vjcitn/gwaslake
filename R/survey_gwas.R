#' plot sample size and other information about studies in IEU GWAS ecosystem
#' @import ggplot2
#' @param phrases character() vector to grep with ignore.case=TRUE in `trait` field, unique matching rows
#' will be collected
#' @param datasource data.frame with `trait` as a column, defaults to 
#' gwaslake::gwidf_2021_01_30; result of ieugwasr::gwasinfo() is intended input
#' @param title_pref character(1) defaults to "IEU GWAS"
#' @param title_text character(1) defaults to a paste of phrases with collapse=", "
#' @param xlab character(1) x axis label defaults to "N controls"
#' @param ylab character(1) x axis label defaults to "N cases"
#' @param \dots passed to grep
#' @note returns ggplot, intended for ggplotly
#' @examples
#' requireNamespace("ggplot2")
#' sg = survey_gwas(c("asthma", "asthmatic")) + 
#'  ggplot2::theme(text=ggplot2::element_text(size=16))
#' if (interactive()) {
#'  requireNamespace("plotly")
#'  plotly::ggplotly(sg)
#' }
#' @export
survey_gwas = function(phrases, datasource=gwaslake::gwidf_2021_01_30, title_pref="IEU GWAS", 
    title_text=paste(sQuote(phrases), collapse=", "), xlab="N controls", ylab="N cases", ...) {
  hits = lapply(phrases, function(x) grep(x, datasource$trait, ignore.case=TRUE, ...))
  stopifnot(length(hits)>0)
  hits = datasource[unique(unlist(hits)),]
  newdf = data.frame(ncont=hits$ncontrol, 
          ncase=hits$ncase, 
          id=hits$id, year = hits$year, 
          trait=hits$trait, 
          pmid=hits$pmid)
  newdf$tag = paste("id:", newdf$id, "<br>yr:", newdf$year, "<br>trait:", newdf$trait,
          "<br>pmid:", newdf$pmid, sep="")
  if (nchar(title_text)>0) title_text = paste("[", paste(title_text, 
                             ifelse(length(title_text)>1, ", ", ""), sep=""), "]", sep="")
  ggplot(newdf, aes(x=ncont, y=ncase, text=tag)) + 
         geom_point() + ggtitle(paste(title_pref, title_text)) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
}
