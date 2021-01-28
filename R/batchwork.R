make_batch = function(x) sapply(strsplit(x,"-"), function(z) paste(z[1], z[2], sep="-"))

#' a simple app supporting filtering and exploration of MRC IEU GWAS batches and associated studies
#' @importFrom dplyr filter mutate
#' @importFrom DT renderDataTable dataTableOutput
#' @importFrom ieugwasr gwasinfo batches
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @param dat a data.frame-like thing expected to be like the result of ieugwasr::gwasinfo(); if NULL,
#' the output of gwasinfo() will be used
#' @return just invokes a shiny app at this point, no return
#' @export
browse_ieu = function(dat=NULL) {
 if (is.null(dat)) dat = ieugwasr::gwasinfo()
 dat = dplyr::mutate(dat, batch=make_batch(id))
 ba = ieugwasr::batches()
 desc = ba$description
 bid = ba$id
 names(bid) = desc
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText("browse_ieu app"),
    selectInput("batch", "batch", choices=bid)
    ),
  mainPanel(
   tabsetPanel(
    tabPanel("sizes",
     DT::dataTableOutput("siz")
     ),
    tabPanel("data",
     DT::dataTableOutput("btab")
     )
    )
   )
  )
 )
 server = function(input, output) {
  output$siz = DT::renderDataTable(ba)
  output$btab = DT::renderDataTable({
   dplyr::filter(dat, batch == input$batch)
  })
 }
 runApp(list(ui=ui, server=server))
}
 
