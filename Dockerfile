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

