FROM rocker/tidyverse:4.5.2

# Pass in GitHub PAT via build argument and set environment variable
ARG GPAT
ENV GITHUB_PAT=${GPAT}

# Extend PATH and set library path for R
ENV PATH="$PATH:/usr/lib/rstudio-server/bin"
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/R/lib"

# Configure RStudio session settings
RUN echo "session-timeout-minutes=0" >> /etc/rstudio/rsession.conf && \
    echo "session-save-action-default=no" >> /etc/rstudio/rsession.conf && \
    echo "copilot-enabled=1" >> /etc/rstudio/rsession.conf

# Install OS-level dependencies in one layer (using --no-install-recommends to avoid extra packages)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        vim \
        python3 python3-dev python3-pip \
        zlib1g libgmp3-dev libglpk-dev gfortran \
        jags \
        libudunits2-dev \
        libgdal-dev libgeos-dev libproj-dev \
        default-jre default-jdk \
        libmysqlclient-dev \
        texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra \
        libmagick++-dev fftw-dev \
        libgsl-dev \
        libboost-all-dev \
        libzmq5 \
        python-is-python3 && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# --------------------------------------------------
# Build FFTW from source
RUN cd /opt && \
    wget http://www.fftw.org/fftw-3.3.10.tar.gz && \
    tar -xvf fftw-3.3.10.tar.gz && \
    cd fftw-3.3.10 && \
    ./configure --enable-shared && \
    make CFLAGS=-fPIC && \
    make install && \
    rm /opt/fftw-3.3.10.tar.gz && \
    rm -rf /opt/fftw-3.3.10

# --------------------------------------------------
# Install latest CMake from source
RUN cd /opt && \
    wget https://github.com/Kitware/CMake/releases/download/v4.0.3/cmake-4.0.3.tar.gz && \
    tar -zxvf cmake-4.0.3.tar.gz && \
    cd cmake-4.0.3 && \
    ./bootstrap && \
    make && \
    make install && \
    rm /opt/cmake-4.0.3.tar.gz && \
    rm -rf /opt/cmake-4.0.3

# --------------------------------------------------
# Install Pandoc (DEB package)
RUN cd /opt && \
    wget https://github.com/jgm/pandoc/releases/download/3.7.0.2/pandoc-3.7.0.2-1-amd64.deb && \
    dpkg -i pandoc-3.7.0.2-1-amd64.deb && \
    rm pandoc-3.7.0.2-1-amd64.deb

# --------------------------------------------------
# Install Miniforge (Conda)
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh -O Miniforge3.sh && \
    bash Miniforge3.sh -b -p /opt/conda && \
    rm Miniforge3.sh
ENV PATH="/opt/conda/bin:$PATH"


# --------------------------------------------------
# Bioconductor packages
# Copy the Bioconductor package list into the container
COPY bioc_package_list.txt /tmp/

# Install Bioconductor packages using the list
RUN Rscript -e "install.packages('BiocManager')" && Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' \
    && Rscript -e "BiocManager::install(scan('/tmp/bioc_package_list.txt', what = '', sep = '\n'))"

## Install CRAN packages that require Bioconductor dependencies
RUN Rscript -e "install.packages(c( \
    'ADAPTS', 'codebook', 'conos', 'ggm', 'isva', 'metap', 'mutoss', 'NMF', \
    'pcalg', 'scGate', 'scMappR', 'Signac', 'restfulr', 'rliger', 'WGCNA' \
))"


# --------------------------------------------------
# CRAN packages
# Copy the package list into the container
COPY cran_package_list.txt /tmp/


# Install CRAN packages using the list
RUN Rscript -e "install.packages(scan('/tmp/cran_package_list.txt', what = '', sep = '\n'))"

# --------------------------------------------------
# Github packages
# Copy GitHub package list into the container
COPY github_packages.txt /tmp/

# Install all GitHub packages in a single Rscript call
RUN Rscript -e "pkgs <- readLines('/tmp/github_packages.txt'); \
                for (pkg in pkgs) { \
                    if (grepl('=', pkg)) { \
                        parts <- unlist(strsplit(pkg, '=')); \
                        remotes::install_github(parts[1], ref = parts[2]); \
                    } else { \
                        remotes::install_github(pkg); \
                    } \
                }"

## Install cellxgene.census and dependencies
RUN pip install --no-cache-dir tiledb && \
    Rscript -e "install.packages('tiledb')" && \
    Rscript -e "install.packages('cellxgene.census', repos = c(\
        'https://chanzuckerberg.r-universe.dev', \
        'https://cloud.r-project.org' \
    ))"

## Install MCPcounter
RUN Rscript -e "devtools::install_github('ebecht/MCPcounter', ,ref='master', subdir='Source')"

## Install ProjecTILs from GitHub
RUN Rscript -e "remotes::install_github('carmonalab/ProjecTILs')"

## Install velocyto.R
RUN wget -q https://github.com/velocyto-team/velocyto.R/archive/refs/tags/0.6.tar.gz && \
    tar -xf 0.6.tar.gz && \
    sed -i '48s/^/\/\//' velocyto.R-0.6/src/routines.cpp && \
    Rscript -e "install.packages('./velocyto.R-0.6', repos = NULL, type = 'source')" && \
    rm -rf 0.6.tar.gz velocyto.R-0.6

# --------------------------------------------------
# Clean up caches and unnecessary files to reduce image size
RUN rm -rf /tmp/* ~/.cache/R /root/.R/ && \
    rm -rf /usr/local/lib/R/site-library/*/help && \
    find /usr/local/lib/R/site-library/ -type f -name "*.o" -delete && \
    find /usr/local/lib/R/site-library/ -type f -name "*.pdf" -delete && \
    find /usr/local/lib/R/site-library/ -type f -name "*.so" -exec strip --strip-unneeded {} \; && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /root/.local /root/.cache /root/.npm /root/.pip /root/.conda
