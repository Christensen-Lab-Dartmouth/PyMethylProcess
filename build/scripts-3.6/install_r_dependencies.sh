#!/bin/bash

python basic_installer.py -p Rcpp simpleCache BiocManager devtools Cairo remotes tidyverse knitr markdown gridExtra multcomp fastICA statmod lme4 base64enc ggpubr forcats
python basic_installer.py -b -p minfi ENmix geneplotter IlluminaHumanMethylation450kanno.ilmn12.hg19
python basic_installer.py -b -p sva S4Vectors DNAcopy gdsfmt illuminaio ggbio
python basic_installer.py -b -p rtracklayer GEOquery LOLA limma missMethyl TCGAbiolinks
python basic_installer.py -g perishky/meffil
