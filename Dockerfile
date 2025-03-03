FROM rocker/tidyverse:4.4.3

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
    wget https://github.com/Kitware/CMake/releases/download/v3.31.6/cmake-3.31.6.tar.gz && \
    tar -zxvf cmake-3.31.6.tar.gz && \
    cd cmake-3.31.6 && \
    ./bootstrap && \
    make && \
    make install && \
    rm /opt/cmake-3.31.6.tar.gz && \
    rm -rf /opt/cmake-3.31.6

# --------------------------------------------------
# Install Pandoc (DEB package)
RUN cd /opt && \
    wget https://github.com/jgm/pandoc/releases/download/3.6.3/pandoc-3.6.3-1-amd64.deb && \
    dpkg -i pandoc-3.6.3-1-amd64.deb && \
    rm pandoc-3.6.3-1-amd64.deb

# --------------------------------------------------
# Install Miniforge (Conda)
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh -O Miniforge3.sh && \
    bash Miniforge3.sh -b -p /opt/conda && \
    rm Miniforge3.sh
ENV PATH="/opt/conda/bin:$PATH"

## Install CRAN packages
RUN Rscript -e "install.packages(c(
    'abind', 'ade4', 'admix', 'afex', 'alarmdata', 'amap', 'animation', 'aod', 'apcluster', 'ape', 
    'aplot', 'arrangements', 'arrayhelpers', 'arrow', 'arsenal', 'ash', 'askpass', 'assertthat', 'AssocTests', 'attempt',
    'aws.s3', 'aws.signature', 'babelgene', 'backports', 'base64', 'base64enc', 'basedosdados', 'BayesLogit', 'bayesm', 'bayesplot', 
    'bayestestR', 'bbmle', 'bdsmatrix', 'beanplot', 'BEDMatrix', 'beeswarm', 'betareg', 'bezier', 'BGLR', 'BH', 
    'BiasedUrn', 'biglm', 'bigmemory', 'bigmemory.sri', 'bigrquery', 'binom', 'BiocManager', 'bipartite', 'bit', 'bit64', 
    'bitops', 'biwt', 'blme', 'blob', 'bmp', 'bookdown', 'BoolNet', 'brew', 'bridgesampling', 'brio', 
    'brms', 'Brobdingnag', 'broom', 'broom.helpers', 'broom.mixed', 'bslib', 'BuyseTest', 'C50', 'ca', 'cachem',
    'Cairo', 'calibrate', 'callr', 'car', 'carData', 'caret', 'caTools', 'cccd', 'cellranger', 'censable',
    'censusapi', 'checkmate', 'chemometrics', 'chirps', 'chron', 'circlize', 'CircStats', 'circular', 'classInt', 'cli', 
    'clipr', 'clock', 'clue', 'ClusterR', 'clustree', 'cmdfun', 'CMplot', 'cmprsk', 'coda', 'codetools',
    'coin', 'colorBlindness', 'colorRamps', 'colorspace', 'colourpicker', 'combinat', 'ComICS', 'commonmark', 'compositions', 'CompQuadForm',
    'concaveman', 'config', 'confintr', 'conflicted', 'conquer', 'corpcor', 'corrplot', 'covr', 'cowplot', 'coxme',
    'cplm', 'cpp11', 'crawl', 'crayon', 'credentials', 'crochet', 'crosstalk', 'crul', 'cubature', 'Cubist',
    'curl', 'cvar', 'dagitty', 'dashboardthemes', 'DataExplorer', 'data.table', 'data.tree', 'dataverse', 'datawizard', 'DBI',
    'dbplyr', 'dbscan', 'DDRTree', 'deldir', 'dendextend', 'dendsort', 'densEstBayes', 'DEoptimR', 'desc', 'DescTools',
    'deSolve', 'devEMF', 'devtools', 'diagram', 'dials', 'DiceDesign', 'dichromat', 'diffobj', 'digest', 'diptest',
    'DirichletReg', 'distillery', 'distr', 'distributional', 'docopt', 'doFuture', 'doMC', 'doParallel', 'doRNG', 'doSNOW',
    'dotCall64', 'dotenv', 'downlit', 'downloader', 'dplyr', 'dqrng', 'drat', 'DT', 'dtplyr', 'dtt',
    'dtw', 'duckdb', 'dygraphs', 'dynamicTreeCut', 'e1071', 'ecmwfr', 'editData', 'effectsize', 'egg', 'elevatr',
    'ellipse', 'ellipsis', 'emdbook', 'emmeans', 'energy', 'english', 'entropy', 'EnvStats', 'estimability', 'evaluate',
    'evd', 'Exact', 'exactextractr', 'exactRankTests', 'expm', 'expss', 'extrafont', 'extrafontdb', 'extRemes', 'factoextra',
    'FactoMineR', 'fANCOVA', 'fansi', 'farver', 'fastcluster', 'fastDummies', 'fastICA', 'fastmap', 'fastmatch', 'fastqcr',
    'fBasics', 'fdrtool', 'fExtremes', 'fftw', 'fGarch', 'fields', 'filelock', 'filematrix', 'fitdistrplus', 'flashClust',
    'flexclust', 'flexmix', 'flextable', 'fMultivar', 'FNN', 'fontawesome', 'forcats', 'foreach', 'forestplot', 'formatR',
    'Formula', 'formula.tools', 'fpc', 'fresh', 'fs', 'fst', 'fstcore', 'furrr', 'futile.logger', 'futile.options',
    'future', 'future.apply', 'GA', 'gamlss', 'gamlss.data', 'gamlss.dist', 'gamlss.tr', 'gapminder', 'gargle', 'gbutils',
    'gclus', 'gdata', 'gdtools', 'generics', 'geomander', 'geometries', 'geometry', 'geonames', 'geos', 'geosphere',
    'gert', 'getopt', 'GetoptLong', 'getPass', 'ggalluvial', 'GGally', 'ggalt', 'ggbeeswarm', 'ggbiplot', 'ggdendro',
    'ggdist', 'ggeffects', 'ggforce', 'ggfortify', 'ggfun', 'ggiraph', 'ggmap', 'ggnetwork', 'ggnewscale', 'ggplot2',
    'ggplotify', 'ggpmisc', 'ggpointdensity', 'ggpp', 'ggprism', 'ggpubr', 'ggraph', 'ggrastr', 'ggrepel', 'ggridges',
    'ggseqlogo', 'ggsci', 'ggsignif', 'ggstats', 'ggtext', 'ggthemes', 'ggvenn', 'gh', 'git2r', 'gitcreds',
    'glasso', 'gld', 'glmmTMB', 'glmnet', 'glmpca', 'GlobalOptions', 'globals', 'glue', 'GMMAT', 'gmodels',
    'gmp', 'gnm', 'goftest', 'golem', 'googleAuthR', 'googledrive', 'googlesheets4', 'gower', 'GPfit', 'gplots',
    'gProfileR', 'gprofiler2', 'graphlayouts', 'gridBase', 'gridExtra', 'gridGraphics', 'gridtext', 'grr', 'GSA', 'gsl',
    'gson', 'gss', 'gsubfn', 'gsw', 'gt', 'gtable', 'gtools', 'gtsummary', 'GWASExactHW', 'h2o',
    'hablar', 'hardhat', 'harmony', 'hash', 'haven', 'hdf5r', 'HDInterval', 'heatmaply', 'here', 'hexbin'
), Ncpus = parallel::detectCores())"


## Install Bioconductor packages
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' \
    && Rscript -e "BiocManager::install(c(
        'affy', 'affyio', 'alabaster.base', 'alabaster.matrix', 'alabaster.ranges', 'alabaster.schemas', 'alabaster.se', 'ALL', 'annotate', 'AnnotationDbi',
        'AnnotationFilter', 'AnnotationForge', 'AnnotationHub', 'AnnotationHubData', 'APAlyzer', 'apeglm', 'AUCell', 'bamsignals', 'basilisk', 'basilisk.utils',
        'batchelor', 'BayesSpace', 'beachmat', 'Biobase', 'BiocBaseUtils', 'BiocCheck', 'BiocFileCache', 'BiocGenerics', 'BiocIO', 'BiocNeighbors',
        'BiocParallel', 'BiocSingular', 'BiocStyle', 'BiocVersion', 'biocViews', 'bioDist', 'biomaRt', 'biomformat', 'Biostrings', 'biovizBase',
        'BiSeq', 'bluster', 'BSgenome', 'BSgenome.Sscrofa.UCSC.susScr3', 'BSgenome.Sscrofa.UCSC.susScr3.masked', 'bsseq', 'bumphunter', 'Category', 'ccrepe', 'celldex',
        'CHETAH', 'clusterProfiler', 'clustifyr', 'CoGAPS', 'ComplexHeatmap', 'ConsensusClusterPlus', 'dada2', 'DECIPHER', 'decontam', 'DEGreport',
        'DelayedArray', 'DelayedMatrixStats', 'DeMixT', 'densvis', 'depmap', 'derfinder', 'derfinderHelper', 'DESeq2', 'DESpace', 'dir.expiry',
        'DirichletMultinomial', 'discordant', 'dittoSeq', 'DMRcate', 'DMRcatedata', 'DNAcopy', 'DO.db', 'DOSE', 'DropletUtils', 'DSS',
        'edgeR', 'EGSEA', 'EGSEAdata', 'EnhancedVolcano', 'ENmix', 'enrichplot', 'EnsDb.Hsapiens.v79', 'ensembldb', 'epiNEM', 'escape',
        'ExperimentHub', 'ExperimentHubData', 'FastqCleaner', 'fastseg', 'FDb.InfiniumMethylation.hg19', 'fgsea', 'fishpond', 'FlowSorted.Blood.EPIC', 'gage', 'gdsfmt',
        'genefilter', 'geneplotter', 'GENESIS', 'GENIE3', 'genomation', 'GenomeInfoDb', 'GenomeInfoDbData', 'GenomicAlignments', 'GenomicFeatures', 'GenomicFiles',
        'GenomicRanges', 'GeomxTools', 'GeoMxWorkflows', 'GEOquery', 'ggtree', 'glmGamPoi', 'globaltest', 'GO.db', 'GOSemSim', 'GOstats',
        'graph', 'graphite', 'GSEABase', 'GSVA', 'Gviz', 'GWASTools', 'gypsum', 'HDF5Array', 'HDO.db', 'hgu133a.db',
        'hgu133plus2.db', 'HSMMSingleCell', 'HybridMTest', 'hypeR', 'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylation450kmanifest', 
        'IlluminaHumanMethylationEPICanno.ilm10b4.hg19', 'IlluminaHumanMethylationEPICmanifest', 'illuminaio', 'impute',
        'interactiveDisplayBase', 'IRanges', 'karyoploteR', 'KEGGdzPathwaysGEO', 'KEGGgraph', 'KEGGREST', 'LEA', 'lfa', 'limma', 'Linnorm',
        'lionessR', 'LOLA', 'LoomExperiment', 'lpsymphony', 'lumi', 'M3Drop', 'Maaslin2', 'MAGeCKFlute', 'MAST', 'MatrixGenerics',
        'megadepth', 'metagenomeSeq', 'metapod', 'MethylAid', 'methylclock', 'methylclockData', 'methylKit', 'methylumi', 'microbiome', 'minet',
        'minfi', 'minfiData', 'miQC', 'missMethyl', 'mixOmics', 'mnem', 'monocle', 'msa', 'MultiAssayExperiment', 'multiMiR',
        'multtest', 'muscat', 'NanoStringNCTools', 'Nebulosa', 'nempi', 'OmnipathR', 'OrganismDbi', 'org.Cf.eg.db', 'org.Hs.eg.db', 'org.Mm.eg.db',
        'org.Rn.eg.db', 'org.Sc.sgd.db', 'org.Ss.eg.db', 'orthogene', 'PADOG', 'PAIRADISE', 'pandaR', 'pathview', 'pcaMethods', 'PCAtools',
        'phyloseq', 'planet', 'plyranges', 'preprocessCore', 'projectR', 'ProtGenerics', 'quantiseqr', 'quantsmooth', 'qusage', 'qvalue',
        'RBGL', 'RCy3', 'reactome.db', 'ReactomePA', 'recount', 'regioneR', 'ResidualMatrix', 'Rgraphviz', 'rhdf5', 'rhdf5filters',
        'Rhdf5lib', 'Rhtslib', 'rnaEditr', 'ROC', 'Rsamtools', 'Rsubread', 'rtracklayer', 'S4Arrays', 'S4Vectors', 'safe',
        'ScaledMatrix', 'scater', 'scDblFinder', 'scde', 'schex', 'scHOT', 'scran', 'scRNAseq', 'scry', 'scuttle', 'SeqArray',
        'seqPattern', 'SeqVarTools', 'sesame', 'sesameData', 'ShortRead', 'siggenes', 'sincell', 'SingleCellExperiment', 'SingleR', 'singscore',
        'slingshot', 'SNPRelate', 'snpStats', 'SparseArray', 'sparseMatrixStats', 'SpatialDecon', 'SpatialExperiment', 'splatter', 'STRINGdb', 'SummarizedExperiment',
        'sva', 'TFBSTools', 'topGO', 'TrajectoryUtils', 'treeio', 'TreeSummarizedExperiment', 'tricycle', 'TxDb.Hsapiens.UCSC.hg19.knownGene', 'txdbmaker', 
        'TxDb.Sscrofa.UCSC.susScr3.refGene', 'tximeta', 'tximport', 'UCell', 'UCSC.utils', 'variancePartition', 'VariantAnnotation', 'vsn', 'wateRmelon', 'Wrench',
        'XVector', 'zellkonverter', 'zlibbioc'
    ), Ncpus = parallel::detectCores())"

## Install CRAN packages that require Bioconductor dependencies
RUN Rscript -e "install.packages(c(
    'ADAPTS', 'codebook', 'conos', 'ggm', 'isva', 'metap', 'mutoss', 'NMF', 
    'pcalg', 'scGate', 'scMappR', 'Signac', 'restfulr', 'rliger', 'WGCNA'
))"

## Install GitHub packages
RUN Rscript -e "remotes::install_github(c(
    'satijalab/azimuth',
    'stefpeschel/NetCoMi',
    'cnfoley/hyprcoloc',
    'bnprks/BPCells/r',
    'bnprks/BPCells',
    'nghiavtr/BPSC',
    'itsrainingdata/sparsebnUtils',
    'itsrainingdata/ccdrAlgorithm',
    'khodosevichlab/CellAnnotatoR',
    'sqjin/CellChat',
    'jokergoo/circlize',
    'stan-dev/cmdstanr',
    'cansysbio/ConsensusTME',
    'CBMR-Single-Cell-Omics-Platform/SCOPfunctions',
    'CostaLab/CrossTalkeR',
    'DARTH-git/dampack',
    'chris-mcginnis-ucsf/DoubletFinder',
    'GfellerLab/EPIC',
    'cole-trapnell-lab/garnett',
    'JGCRI/gcamextractor',
    'RubD/Giotto@cless',
    'ghcarlalan/graveler',
    'omnideconv/immunedeconv',
    '10XGenomics/loupeR',
    'dklinges9/mcera5',
    'ebecht/MCPcounter',
    'ilyamaclean/microclima',
    'cit-bioinfo/mMCP-counter',
    'cole-trapnell-lab/monocle3',
    'MRCIEU/MRInstruments',
    'gqi/MRMix',
    'rondolab/MR-PRESSO',
    'immunogenomics/presto',
    'carmonalab/ProjecTILs',
    'jbisanz/qiime2R',
    'cysouw/qlcMatrix',
    'WSpiller/RadialMR',
    'rmcelreath/rethinking',
    'JGCRI/rgcam',
    'JGCRI/rpackageutils',
    'cellgeni/sceasy',
    'CBMR-Single-Cell-Omics-Platform/SCOPfunctions',
    'igrabski/sc-SHC',
    'satijalab/seurat-data',
    'pinin4fjords/shinyngs',
    'satijalab/seurat-wrappers',
    'arc85/singleseqgset',
    'EliGurarie/smoove',
    'zdk123/SpiecEasi',
    'GraceYoon/SPRING',
    'carmonalab/STACAS',
    'poisonalien/trackplot',
    'MRCIEU/TwoSampleMR',
    'timokelder/UNSEEN',
    'velocyto-team/velocyto.R',
    'dviraran/xCell',
    'statOmics/zingeR'
))"


## Install cellxgene.census and dependencies
RUN pip install --no-cache-dir tiledb && \
    Rscript -e "install.packages('tiledb')" && \
    Rscript -e "install.packages('cellxgene.census', repos = c(
        'https://chanzuckerberg.r-universe.dev', 
        'https://cloud.r-project.org'
    ))"

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