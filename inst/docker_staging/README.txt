
The docker image vjcitn/bioc_staged_rel:0.0.2 has, in addition
to base packages

c("abind", "AnnotationDbi", "arrangements", "askpass", "assertthat", 
"backports", "base64enc", "BH", "Biobase", "BiocFileCache", "BiocGenerics", 
"BiocManager", "BiocParallel", "BiocStyle", "BiocVersion", "biomaRt", 
"Biostrings", "bit", "bit64", "bitops", "blob", "BMA", "bookdown", 
"brew", "brio", "broom", "BSgenome", "callr", "car", "carData", 
"cellranger", "cli", "clipr", "colorspace", "commonmark", "CompQuadForm", 
"conquer", "covr", "cowplot", "cpp11", "crayon", "crosstalk", 
"curl", "data.table", "DBI", "dbplyr", "DelayedArray", "DEoptimR", 
"desc", "devtools", "diffobj", "digest", "docopt", "downlit", 
"dplyr", "DT", "ellipsis", "evaluate", "fansi", "farver", "forcats", 
"foreach", "formatR", "fs", "futile.logger", "futile.options", 
"generics", "GenomeInfoDb", "GenomeInfoDbData", "GenomicAlignments", 
"GenomicFeatures", "GenomicRanges", "ggplot2", "ggrepel", "gh", 
"git2r", "glmnet", "glue", "gmp", "gridExtra", "gtable", "haven", 
"highr", "hms", "htmltools", "htmlwidgets", "httr", "ini", "inline", 
"IRanges", "isoband", "iterators", "iterpc", "jsonlite", "knitr", 
"labeling", "lambda.r", "later", "lazyeval", "leaps", "lifecycle", 
"littler", "lme4", "magrittr", "maptools", "markdown", "MatrixGenerics", 
"MatrixModels", "matrixStats", "memoise", "meta", "metafor", 
"mime", "minqa", "mnormt", "munsell", "mvtnorm", "nloptr", "nortest", 
"numDeriv", "openssl", "openxlsx", "pbkrtest", "pcaPP", "pillar", 
"pkgbuild", "pkgconfig", "pkgdown", "pkgload", "plogr", "plotly", 
"plyr", "praise", "prettyunits", "processx", "progress", "promises", 
"ps", "psych", "purrr", "quantreg", "R6", "ragg", "randomForest", 
"rappdirs", "rcmdcheck", "RColorBrewer", "Rcpp", "RcppArmadillo", 
"RcppEigen", "RCurl", "readr", "readxl", "rematch", "rematch2", 
"remotes", "reshape", "reshape2", "rex", "Rhtslib", "rio", "rjson", 
"rlang", "rmarkdown", "robustbase", "roxygen2", "rprojroot", 
"rrcov", "Rsamtools", "RSQLite", "rstudioapi", "rtracklayer", 
"rversions", "S4Vectors", "scales", "sessioninfo", "shape", "snow", 
"snpStats", "sp", "SparseM", "startup", "statmod", "stringi", 
"stringr", "SummarizedExperiment", "sys", "systemfonts", "testthat", 
"textshaping", "tibble", "tidyr", "tidyselect", "tinytex", "tmvnsim", 
"usethis", "utf8", "vctrs", "viridisLite", "waldo", "whisker", 
"withr", "xfun", "XML", "xml2", "xopen", "XVector", "yaml", "zip", 
"zlibbioc")

As of 7 Feb 2021, the vjcitn/gwaslake:latest has, in addition to these,

c("bslib", "cachem", "coloc", "fastmap", "finemapr", "gargle", 
"genetics.binaRies", "googleAuthR", "graph", "gwaslake", "gwasvcf", 
"httpuv", "ieugwasr", "igraph", "jquerylib", "ontologyIndex", 
"ontologyPlot", "ontoProc", "paintmap", "Rgraphviz", "sass", 
"shiny", "shinytoastr", "sourcetools", "VariantAnnotation", "xtable"
)

Now we would want to exclude the packages that are in the remotes

    mrcieu/ieugwasr,
    vjcitn/gwasvcf,
    vjcitn/finemapr,
    cran/ontologyPlot

So the  0.0.3 of staged_rel will add installation of

c("bslib", "cachem", "coloc", "fastmap", "gargle", 
"genetics.binaRies", "googleAuthR", "graph", 
"httpuv", "igraph", "jquerylib", "ontologyIndex", 
"ontoProc", "paintmap", "Rgraphviz", "sass", 
"shiny", "shinytoastr", "sourcetools", "VariantAnnotation", "xtable"
)

in a dockerfile building from 0.0.2
