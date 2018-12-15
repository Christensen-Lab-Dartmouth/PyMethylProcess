FROM continuumio/miniconda3:4.5.11

RUN apt-get -y update

RUN apt-get -y install gcc g++ libcairo2-dev libxt-dev xvfb

RUN pip install numpy==1.15.4

RUN conda install -y -c r r=3.5.1 r-cairo=1.5_9 r-devtools=1.13.6 python=3.6.7

RUN conda install -y scikit-learn=0.20.1 rpy2=2.9.1 sqlite=3.25.3 readline=7.0 click=6.7 cairo=1.14.12 python=3.6.7

RUN conda install -y -c bioconda bioconductor-biocinstaller=1.30.0 bioconductor-geneplotter=1.58.0

RUN conda install -y -c conda-forge unzip=6.0 xorg-libx11=1.6.6 tar=1.29 python=3.6.7

RUN pip install fancyimpute==0.4.2 pandas==0.23.4

RUN mkdir /scripts/

COPY installer.py /scripts/

RUN python /scripts/installer.py install_r_packages -p BiocManager -p remotes -p knitr -p markdown -p gridExtra -p multcomp -p fastICA -p statmod -p lme4 -p Cairo

RUN python /scripts/installer.py install_custom -m -p GEOquery=3.8 -p sva -p S4Vectors -p DNAcopy -p gdsfmt -p ENmix -p illuminaio

RUN wget https://github.com/perishky/meffil/archive/master.zip && unzip master.zip && mv meffil-master meffil && R CMD INSTALL meffil

RUN apt-get install -y libopenmpi-dev

RUN conda install -y curl=7.62.0 openssl=1.1.1 python=3.6.7

RUN conda install -y -c conda-forge matplotlib=3.0.2 mpi=1.0 openmpi=3.1.2 pathos=0.2.1 python=3.6.7

RUN pip install pyina==0.2.0

RUN python /scripts/installer.py install_r_packages -p openssl

RUN python /scripts/installer.py install_custom -m -p rtracklayer=3.8 -p GEOquery=3.8

COPY *.py /scripts/

WORKDIR /root

ENTRYPOINT ["/usr/bin/tini","-s","--"]
