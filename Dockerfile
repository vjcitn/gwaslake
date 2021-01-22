FROM vjcitn/bioc_staged_rel:0.0.2

ENV GITHUB_PAT 99b8832b717858e0bed22425d812902ad3f37459

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); devtools::install('.', build_vignettes=TRUE, repos = BiocManager::repositories())"

ENV GITHUB_PAT=
