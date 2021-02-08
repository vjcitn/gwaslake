Development of a Terra-compatible container for gwaslake is complicated
by the fact that the MRC gwasglue package has complex github dependencies.

We used dockerfile1 to introduce most dependencies on top of

us.gcr.io/anvil-gcr-public/anvil-rstudio-bioconductor:0.0.9

dockerfile1 got us to vjcitn/gwlakeanvil:0.0.2

dockerfile2 adds ontologyPlot (which was dropped temporarily from CRAN) and
other packages needed as we built up gwaslake to 0.0.14

dockerfile2 gets us to vjcitn/gwlakeanvil:0.0.3
