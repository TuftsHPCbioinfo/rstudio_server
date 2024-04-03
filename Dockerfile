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

## Install Seurat5, seurat-wrappers and singleR
RUN Rscript -e 'devtools::install_version("Matrix")' \
    && Rscript -e "install.packages('SeuratObject')" \
    && Rscript -e "install.packages('Seurat')" \
    && Rscript -e "BiocManager::install('SingleR')" 

RUN Rscript -e "BiocManager::install(c('org.Cf.eg.db', 'org.Sc.sgd.db', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db'))" \
    && Rscript -e "BiocManager::install('EnhancedVolcano')" \
    && Rscript -e "install.packages('kableExtra')" \
    && Rscript -e "install.packages('magick')" \
    && Rscript -e "install.packages('pals')" \
    && Rscript -e "install.packages('hgnc')" \
    && Rscript -e "BiocManager::install('orthogene')" \
    && Rscript -e "BiocManager::install('GSVA')" \
    && Rscript -e "BiocManager::install('biomaRt')" \
    && Rscript -e "BiocManager::install('clustifyr')" \
    && Rscript -e "BiocManager::install(c('escape'))" \
    && Rscript -e "BiocManager::install(c('DelayedMatrixStats'))" \
    && Rscript -e "BiocManager::install(c('DelayedArray'))" \
    && Rscript -e "BiocManager::install(c('escape'))" \
    && Rscript -e "BiocManager::install(c('clusterProfiler'))" \
    && Rscript -e "BiocManager::install(c('ReactomePA'))" \
    && Rscript -e "BiocManager::install(c('DOSE'))" \
    && Rscript -e "install.packages('sctransform')" \
    && Rscript -e 'BiocManager::install("slingshot")'

## scDblFinder and DropletUtils
RUN Rscript -e 'BiocManager::install("scDblFinder")' \
    && Rscript -e 'BiocManager::install("DropletUtils")' \
    && Rscript -e 'BiocManager::install("dittoSeq")'

## scCATCH
RUN Rscript -e 'BiocManager::install("scCATCH")'


## CoGAPS, tricycle, celldex, miQC, Nebulosa, schex, and rliger
RUN Rscript -e 'BiocManager::install("CoGAPS")' \
    && Rscript -e 'BiocManager::install("tricycle")' \
    && Rscript -e 'BiocManager::install("celldex")' \ 
    && Rscript -e 'BiocManager::install("miQC")' \    
    && Rscript -e 'BiocManager::install("Nebulosa")' \
    && Rscript -e 'BiocManager::install("schex")' \ 
    && Rscript -e "install.packages('rliger')"

## scMappR
RUN Rscript -e 'BiocManager::install("pcaMethods")' \
    && Rscript -e 'BiocManager::install("preprocessCore")' \
    && Rscript -e 'BiocManager::install("GSVA")' \
    && Rscript -e 'install.packages("ADAPTS")' \
    && Rscript -e "install.packages('scMappR')"


## Harmony, PCAtools, SoupX, scde, sincell
RUN Rscript -e "install.packages('harmony')" \
    && Rscript -e 'BiocManager::install("PCAtools")' \
    && Rscript -e "install.packages('SoupX')" \
    && Rscript -e 'BiocManager::install("scde")' \
    && Rscript -e 'BiocManager::install("sincell")'

## BioSingular, glmGamPoi, GSVA, splatter, muscat, BayesSpace, pagoda2
RUN Rscript -e "BiocManager::install('BiocSingular')" \ 
    && Rscript -e "BiocManager::install('glmGamPoi')" \
    && Rscript -e 'BiocManager::install("GSVA")' \
    && Rscript -e 'BiocManager::install("splatter")' \
    && Rscript -e 'BiocManager::install("muscat")' \ 
    && Rscript -e 'BiocManager::install("BayesSpace")' \
    && Rscript -e "install.packages('pagoda2')"

## signac, scHot, M3Drop, cellkonverter, iCellR
RUN Rscript -e "install.packages('Signac')" \
    && Rscript -e 'BiocManager::install("scHOT")' \
    && Rscript -e 'BiocManager::install("M3Drop")' \
    && Rscript -e 'BiocManager::install("zellkonverter")'\
    && Rscript -e "install.packages('scatterplot3d')" \
    && Rscript -e "install.packages('NbClust')" \
    && Rscript -e "install.packages('iCellR')"

# monocle
RUN Rscript -e "devtools::install_github('cysouw/qlcMatrix')" \
    && Rscript -e 'BiocManager::install("monocle")'

## cellxgene.census
RUN pip install tiledb \
    && Rscript -e "install.packages('tiledb')" \
    && Rscript -e "install.packages('cellxgene.census',repos=c('https://chanzuckerberg.r-universe.dev', 'https://cloud.r-project.org'))"

# Seuratwrappers
RUN Rscript -e "remotes::install_github('satijalab/seurat-wrappers')" \
    && Rscript -e "devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')" \
    && Rscript -e "remotes::install_github('omnideconv/immunedeconv')" \
    && Rscript -e "devtools::install_github('arc85/singleseqgset')" \
    && Rscript -e "remotes::install_github('mojaveazure/seurat-disk')"

## Install  monocle3
RUN Rscript -e "BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', \
                       'limma', 'lme4', 'S4Vectors', 'batchelor', 'SingleCellExperiment', \
                       'SummarizedExperiment', 'HDF5Array', \
                       'terra', 'ggrastr'))" \
    && Rscript -e "devtools::install_github('cole-trapnell-lab/monocle3')"

# SCOPfunctions
RUN Rscript -e 'BiocManager::install("MAST")' \
    && Rscript -e 'devtools::install_github("CBMR-Single-Cell-Omics-Platform/SCOPfunctions")'

## Install ProjecTILs
RUN Rscript -e "install.packages('scGate')" \
    && Rscript -e "remotes::install_github('carmonalab/STACAS')" \
    && Rscript -e "remotes::install_github('carmonalab/ProjecTILs')" \
    && Rscript -e "remotes::install_github('10XGenomics/loupeR')"


## CellAnnotatoR
RUN Rscript -e "install.packages('conos')" \
    && Rscript -e "devtools::install_github('khodosevichlab/CellAnnotatoR')"


## Garnett
RUN Rscript -e "BiocManager::install(c('DelayedArray', 'DelayedMatrixStats'))" \
    && Rscript -e "devtools::install_github('cole-trapnell-lab/garnett')"


## velocyto.R
RUN Rscript -e "devtools::install_github('velocyto-team/velocyto.R')"

## CellChat
RUN Rscript -e "install.packages('NMF')" \
    && Rscript -e "devtools::install_github('jokergoo/circlize')" \
    && Rscript -e "BiocManager::install('ComplexHeatmap')" \
    && Rscript -e "devtools::install_github('sqjin/CellChat')"

## Dependencies for SCHNAPPs
RUN Rscript -e 'devtools::install_github("nghiavtr/BPSC")' \
    && Rscript -e 'devtools::install_github("statOmics/zingeR")' \
    && Rscript -e 'devtools::install_github("bnprks/BPCells")'

## SCHNAPPs
RUN Rscript -e "devtools::install_github('C3BI-pasteur-fr/UTechSCB-SCHNAPPs',dependencies = TRUE)"

