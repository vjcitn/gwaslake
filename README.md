# gwaslake -- convenience containers and documents for exploring the MRC Integrative Epidemiology Unity "gwasglue" ecosystem for GWAS

We are indebted to the University of Bristol MRC Integrative Epidemiology Unit for
infrastructure underlying this ecosystem:

![dockerhub](https://github.com/vjcitn/gwaslake/raw/master/inst/images/mrcieuglue.png)

This repo includes code for a 'workshop' consisting of a docker container and
Rmarkdown files that explore the use of this infrastructure.  On
any internet-connected machine with a docker client, you can use

```
docker run -ti -e PASSWORD=abc -p 8787:8787 vjcitn/gwaslake:latest
```

to acquire the docker container image and start it up so that it serves
Rstudio at port localhost:8787.  There will be a vignettes folder in
the working directory, and all software necessary to build the vignettes
is already installed in the container.

Note that the MRC IEU `gwasglue` package is not installed -- it has too
many dependencies on remote packages in github.  However, you can
install it and its dependencies manually with remote reference `mrcieu/gwasglue`.
