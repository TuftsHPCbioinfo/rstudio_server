FROM rocker/tidyverse:4.3.3

ENV PATH="$PATH:/usr/lib/rstudio-server/bin"
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/R/lib"

RUN apt-get -y update && apt-get -y install python3 python3-dev python3-pip zlib1g libgmp3-dev libglpk-dev gfortran \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libmysqlclient-dev \
    texlive-latex-base \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-extra \
    libmagick++-dev fftw3 fftw-dev \
    libboost-all-dev \
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

RUN Rscript -e "install.packages('devtools')" \
    && Rscript -e "install.packages('rJava')" \
    && Rscript -e "install.packages('openxlsx', dependencies = TRUE)" \
    && Rscript -e "install.packages('remotes')" \
    && Rscript -e "install.packages('rmdformats')" \ 
    && Rscript -e "install.packages('vioplot')" \
    && Rscript -e "install.packages('clustree')" \
    && Rscript -e "install.packages('ggraph')" \
    && Rscript -e "install.packages('xlsx')" \
    && Rscript -e "install.packages('ggthemes')" \ 
    && Rscript -e 'install.packages("rmarkdown")' \ 
    && Rscript -e 'install.packages("markdown")' \ 
    && Rscript -e 'install.packages("reticulate")' \ 
    && Rscript -e 'install.packages("data.table")' \ 
    && Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")' \
    && Rscript -e "install.packages('RColorBrewer')" \
    && Rscript -e "install.packages('ggsci')" \ 
    && Rscript -e "BiocManager::install('qusage')" \
    && Rscript -e "BiocManager::install('batchelor')" \ 
    && Rscript -e "BiocManager::install(c('GEOquery'))" \
    && Rscript -e "install.packages('pandoc')" \
    && Rscript -e "install.packages('wesanderson')" \
    && Rscript -e "BiocManager::install('EGSEA')" \
    && Rscript -e "BiocManager::install('EGSEAdata')" \
    && Rscript -e "install.packages('ggbiplot')" \
    && Rscript -e "BiocManager::install('ReactomePA')" 
