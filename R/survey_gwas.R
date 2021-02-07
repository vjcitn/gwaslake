#' plot sample size and other information about studies in IEU GWAS ecosystem
#' @rawNamespace import(ggplot2, except=last_plot)
#' @param phrases character() vector to grep with ignore.case=TRUE in `trait` field, unique matching rows
#' will be collected
#' @param datasource data.frame with `trait` as a column, defaults to 
#' gwaslake::gwidf_2021_01_30; result of ieugwasr::gwasinfo() is intended input
#' @param title_pref character(1) defaults to "IEU GWAS"
#' @param title_text character(1) defaults to a paste of phrases with collapse=", "
#' @param xlab character(1) x axis label defaults to "N controls"
#' @param ylab character(1) x axis label defaults to "N cases"
#' @param \dots passed to grep
#' @note Returns ggplot, intended for ggplotly.  When used with ggplotly, hover over points to get
#' details such as exact definition of trait, study ID, PMID if available.  
#' For some phenotypes, numbers of cases or controls are unavailable.
#' If this occurs for all studies assessing the selected phenotypes, an error is thrown.
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
  ncas = hits$ncase
  if (all(is.na(ncas))) stop("no non-missing case counts present")
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


hitco = function(phrases="asthma", datasource=gwidf_2021_01_30, ...) {
 hits = lapply(phrases, function(x) grep(paste0("^", x, "$"), tolower(datasource$trait), ...))
 stopifnot(length(hits)>0)
 tmp = datasource[unique(unlist(hits)),]
 studs = unique(tmp$id)
 do.call(rbind, lapply(studs, tophits))
}


#' simple app to help survey diverse phenotypes
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @import plotly
#' @export
survey_app = function() {
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText("gwaslake survey app"),
    helpText("Note: until further notice, results of this app are limited to those studies with traits mapping exactly to Disease Ontology terms."),
    helpText("Direct use of ieugwasr can lead to many more relevant studies."),
    selectInput("pheno", "pheno", choices=sort(as.character(mapped_traits_demo)))
   ),
   mainPanel(
    tabsetPanel(
     tabPanel("survey",
      plotlyOutput("surv")
     ),
     tabPanel("top hits",
      dataTableOutput("hittab")
     )
    )
   )
  )
 )
 server = function(input, output) {
  output$surv = renderPlotly({
   ggplotly(survey_gwas(input$pheno))
  })
  output$hittab = renderDataTable({
   as.data.frame(hitco(phrases=input$pheno))
  })
 }
 runApp(list(ui=ui, server=server))
}
