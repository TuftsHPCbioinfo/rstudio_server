FROM rocker/tidyverse:4.4.0

ENV PATH="$PATH:/usr/lib/rstudio-server/bin"
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/R/lib"

## Configures to avoid session data filling $HOME
Run apt-get update &&  apt-get install -y vim && \
    echo session-timeout-minutes=0 >> /etc/rstudio/rsession.conf  && \
    echo session-timeout-kill-hours=96 >> /etc/rstudio/rsession.conf && \
    echo session-save-action-default=no >> /etc/rstudio/rsession.conf 

RUN apt-get -y update && apt-get -y install python3 python3-dev python3-pip zlib1g libgmp3-dev libglpk-dev gfortran \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    default-jre \
    default-jdk \
    libmysqlclient-dev \
    texlive-latex-base \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-extra \
    libmagick++-dev fftw3 fftw-dev \
    libboost-all-dev \ 
    python-is-python3 \
    && apt-get -y autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

## Install fftw
RUN cd /opt && wget http://www.fftw.org/fftw-3.3.10.tar.gz \
    && tar -xvf fftw-3.3.10.tar.gz \
    && cd fftw-3.3.10 \
    && ./configure --enable-shared \
    && make CFLAGS=-fPIC \
    && make install \
    && rm /opt/fftw-3.3.10.tar.gz

## Install latest cmake
RUN cd /opt \
    && wget https://github.com/Kitware/CMake/releases/download/v3.28.0-rc3/cmake-3.28.0-rc3.tar.gz \
    && tar -zxvf cmake-3.28.0-rc3.tar.gz \
    && cd cmake-3.28.0-rc3 \
    && ./bootstrap \
    && make \
    && make install

## Install pandoc
RUN cd /opt \
    && wget https://github.com/jgm/pandoc/releases/download/3.1.13/pandoc-3.1.13-1-amd64.deb \
    && dpkg -i pandoc-3.1.13-1-amd64.deb \
    && rm pandoc-3.1.13-1-amd64.deb

RUN && Rscript -e "install.packages('amap')" \
    && Rscript -e "install.packages('apcluster')" \
    && Rscript -e "install.packages('base64')" \
    && Rscript -e "install.packages('basedosdados')" \
    && Rscript -e "install.packages('bdsmatrix')" \
    && Rscript -e "install.packages('beanplot')" \
    && Rscript -e "install.packages('beeswarm')" \
    && Rscript -e "install.packages('biglm')" \
    && Rscript -e "install.packages('bigrquery')" \
    && Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")' \
    && Rscript -e "install.packages('admix')" \
    && Rscript -e "install.packages('BGLR')" \
    && Rscript -e "install.packages('caret')" \
    && Rscript -e "install.packages('chemometrics')" \
    && Rscript -e "install.packages('circlize')" \
    && Rscript -e "install.packages('classInt')" \
    && Rscript -e "install.packages('clue')" \
    && Rscript -e "install.packages('coda')" \
    && Rscript -e "install.packages('confintr')" \
    && Rscript -e "install.packages('corpcor')" \
    && Rscript -e "install.packages('cplm')" \
    && Rscript -e "install.packages('clustree')" \
    && Rscript -e "install.packages('cowplot')" \
    && Rscript -e 'install.packages("data.table")' \ 
    && Rscript -e "install.packages('devtools')" \
    && Rscript -e "install.packages('DT')" \
    && Rscript -e "install.packages('deldir')" \
    && Rscript -e "install.packages('diptest')" \
    && Rscript -e "install.packages('doParallel')" \
    && Rscript -e "install.packages('doSNOW')" \
    && Rscript -e "install.packages('dotCall64')" \
    && Rscript -e "install.packages('dotenv')" \
    && Rscript -e "install.packages('downloader')" \
    && Rscript -e "install.packages('dqrng')" \
    && Rscript -e "install.packages('dynamicTreeCut')" \
    && Rscript -e "install.packages('ellipse')" \
    && Rscript -e "install.packages('e1071')" \
    && Rscript -e "install.packages('FNN')" \
    && Rscript -e "install.packages('fastDummies')" \
    && Rscript -e "install.packages('fastICA')" \
    && Rscript -e "install.packages('fields')" \
    && Rscript -e "install.packages('fitdistrplus')" \
    && Rscript -e "install.packages('flexclust')" \
    && Rscript -e "install.packages('flexmix')" \
    && Rscript -e "install.packages('fpc')" \
    && Rscript -e "install.packages('Formula')" \
    && Rscript -e "install.packages('GlobalOptions')" \
    && Rscript -e "install.packages('gapminder')" \
    && Rscript -e "install.packages('gdata')" \
    && Rscript -e "install.packages('getopt')" \ 
    && Rscript -e "install.packages('ggbeeswarm')" \
    && Rscript -e "install.packages('ggdendro')" \
    && Rscript -e "install.packages('ggm')" \
    && Rscript -e "install.packages('ggpmisc')" \
    && Rscript -e "install.packages('ggpp')" \
    && Rscript -e "install.packages('ggprism')" \
    && Rscript -e "install.packages('glmmTMB')" \
    && Rscript -e "install.packages('glmnet')" \
    && Rscript -e "install.packages('gmodels')" \
    && Rscript -e "install.packages('goftest')" \
    && Rscript -e "install.packages('gridBase')" \
    && Rscript -e "install.packages('ggbiplot')" \
    && Rscript -e "install.packages('ggpubr')" \
    && Rscript -e "install.packages('ggridges')" \
    && Rscript -e "install.packages('ggrepel')" \
    && Rscript -e "install.packages('ggsci')" \ 
    && Rscript -e "install.packages('ggraph')" \
    && Rscript -e "install.packages('ggthemes')" \ 
    && Rscript -e "install.packages('ggvenn')" \ 
    && Rscript -e "install.packages('gridExtra')" \ 
    && Rscript -e "install.packages('Hmisc')" \
    && Rscript -e "install.packages('haven')" \
    && Rscript -e "install.packages('h2o')" \
    && Rscript -e "install.packages('JADE')" \
    && Rscript -e 'install.packages("janitor")' \
    && Rscript -e "install.packages('kableExtra')" \
    && Rscript -e "install.packages('magick')" \

    && Rscript -e "install.packages('markdown')" \ 
    && Rscript -e "install.packages('openxlsx', dependencies = TRUE)" \
    && Rscript -e "install.packages('pals')" \
    && Rscript -e "install.packages('pandoc')" \
    && Rscript -e "install.packages('pals')" \
    && Rscript -e "install.packages('patchwork')" \
    && Rscript -e "install.packages('PerformanceAnalytics')" \
    && Rscript -e "install.packages('plotly')" \
    && Rscript -e "install.packages('poolr')" \
    && Rscript -e "install.packages('randomForest')" \
    && Rscript -e "install.packages('RColorBrewer')" \
    && Rscript -e "install.packages('RANN')" \
    && Rscript -e "install.packages('ROCR')" \
    && Rscript -e "install.packages('RPMM')" \
    && Rscript -e "install.packages('RUnit')" \
    && Rscript -e "install.packages('RcppAnnoy')" \
    && Rscript -e "install.packages('RcppML')" \
    && Rscript -e "install.packages('RcppProgress')" \
    && Rscript -e "install.packages('remotes')" \
    && Rscript -e "install.packages('reshape2')" \
    && Rscript -e "install.packages('rmdformats')" \ 
    && Rscript -e "install.packages('rJava')" \
    && Rscript -e "install.packages('rmarkdown')" \ 
    && Rscript -e "install.packages('reticulate')" \ 
    && Rscript -e "install.packages('Rtsne')" \ 
    && Rscript -e "install.packages('scales')" \ 
    && Rscript -e "install.packages('stringi')" \ 
    && Rscript -e "install.packages('TMB')" \
    && Rscript -e "install.packages('umap')" \
    && Rscript -e "install.packages('vioplot')" \
    && Rscript -e "install.packages('wesanderson')" \
    && Rscript -e "install.packages('xlsx')" 


