FROM vjcitn/bioc_staged_rel:0.0.2

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); devtools::install('.', build_vignettes=TRUE, repos = BiocManager::repositories())"

