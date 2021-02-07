FROM vjcitn/bioc_staged_312:0.0.3

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); devtools::install('.', build_vignettes=TRUE, repos = BiocManager::repositories())"

