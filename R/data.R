#' a data.frame of genes(EnsDb.Hsapiens.v75) with build appended
#' @docType data
#' @format data.frame
#' @examples
#' head(genes_v75,3)
"genes_v75"

#' a data.frame of eQTL results for testing
#' @docType data
#' @format data.frame
#' @examples
#' head(little_eq)
"little_eq"

#' a UCSC liftOver chain
#' @docType data
#' @format chain
#' @examples
#' ch19to38
"ch19to38"

#' a snapshot of output of ieugwasr::gwasinfo()
#' @docType data
#' @format data.frame
#' @note Includes 'batchcode' computed by code commented out in example
#' @examples
#' # batchcode = strsplit(gwidf_2021_01_30$id, "-")
#' # batchcode = sapply(batchcode,function(x) paste0(x[1], "-", x[2]))
#' head(gwidf_2021_01_30)
"gwidf_2021_01_30"
