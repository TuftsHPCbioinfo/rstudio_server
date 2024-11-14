FROM rocker/tidyverse:4.4.2

ARG GPAT
ENV GITHUB_PAT=ENV GITHUB_PAT=${GPAT}

ENV PATH="$PATH:/usr/lib/rstudio-server/bin"
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/R/lib"

## Configures to avoid session data filling $HOME
RUN apt-get update &&  apt-get install -y vim && \
    echo session-timeout-minutes=0 >> /etc/rstudio/rsession.conf  && \
    echo session-save-action-default=no >> /etc/rstudio/rsession.conf && \
    echo copilot-enabled=1 >> /etc/rstudio/rsession.conf 

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
    libgsl-dev \
    libboost-all-dev \ 
    libzmq5 \
    python-is-python3 \
    && apt-get -y autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Install fftw
RUN cd /opt && wget http://www.fftw.org/fftw-3.3.10.tar.gz \
    && tar -xvf fftw-3.3.10.tar.gz \
    && cd fftw-3.3.10 \
    && ./configure --enable-shared \
    && make CFLAGS=-fPIC \
    && make install \
    && rm /opt/fftw-3.3.10.tar.gz \
    && rm -rf /opt/fftw-3.3.10

## Install latest cmake
RUN cd /opt \
    && wget https://github.com/Kitware/CMake/releases/download/v3.31.0/cmake-3.31.0.tar.gz\
    && tar -zxvf cmake-3.31.0.tar.gz\
    && cd cmake-3.31.0 \
    && ./bootstrap \
    && make \
    && make install \
    && rm /opt/cmake-3.31.0.tar.gz \
    && rm -rf /opt/cmake-3.31.0

## Install pandoc
RUN cd /opt \
    && wget https://github.com/jgm/pandoc/releases/download/3.5/pandoc-3.5-1-amd64.deb\
    && dpkg -i pandoc-3.5-1-amd64.deb \
    && rm pandoc-3.5-1-amd64.deb

## Install Miniforge
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh -O Miniforge3.sh \
    && bash Miniforge3.sh -b -p /opt/conda \
    && rm Miniforge3.sh
ENV PATH="/opt/conda/bin:$PATH"

## Install CRAN packages
RUN Rscript -e "install.packages(c('abind','ade4','admix','afex','alarmdata','amap','animation','aod','apcluster','ape','aplot','arrangements','arrayhelpers','arrow','arsenal','ash','askpass','assertthat','AssocTests','attempt','aws.s3','aws.signature','babelgene','backports','base64','base64enc','basedosdados','BayesLogit','bayesm','bayesplot','bayestestR','bbmle','bdsmatrix','beanplot','BEDMatrix','beeswarm','betareg','bezier','BGLR','BH','BiasedUrn','biglm','bigmemory','bigmemory.sri','bigrquery','binom','BiocManager','bipartite','bit','bit64','bitops','biwt','blme','blob','bmp','bookdown','BoolNet','brew','bridgesampling','brio','brms','Brobdingnag','broom','broom.helpers','broom.mixed','bslib','BuyseTest','C50','ca','cachem','Cairo','calibrate','callr','car','carData','caret','caTools','cccd','cellranger','censable','censusapi','checkmate','chemometrics','chirps','chron','circlize','CircStats','circular','classInt','cli','clipr','clock','clue','ClusterR','clustree','cmdfun','CMplot','cmprsk','coda','codetools','coin','colorBlindness','colorRamps','colorspace','colourpicker','combinat','ComICS','commonmark','compositions','CompQuadForm','concaveman','config','confintr','conflicted','conquer','corpcor','corrplot','covr','cowplot','coxme','cplm','cpp11','crawl','crayon','credentials','crochet','crosstalk','crul','cubature','Cubist','curl','cvar','dagitty','dashboardthemes','DataExplorer','data.table','data.tree','dataverse','datawizard','DBI','dbplyr','dbscan','DDRTree','deldir','dendextend','dendsort','densEstBayes','DEoptimR','desc','DescTools','deSolve','devEMF','devtools','diagram','dials','DiceDesign','dichromat','diffobj','digest','diptest','DirichletReg','distillery','distr','distributional','docopt','doFuture','doMC','doParallel','doRNG','doSNOW','dotCall64','dotenv','downlit','downloader','dplyr','dqrng','drat','DT','dtplyr','dtt','dtw','duckdb','dygraphs','dynamicTreeCut','e1071','ecmwfr','editData','effectsize','egg','elevatr','ellipse','ellipsis','emdbook','emmeans','energy','english','entropy','EnvStats','estimability','evaluate','evd','Exact','exactextractr','exactRankTests','expm','expss','extrafont','extrafontdb','extRemes','factoextra','FactoMineR','fANCOVA','fansi','farver','fastcluster','fastDummies','fastICA','fastmap','fastmatch','fastqcr','fBasics','fdrtool','fExtremes','fftw','fGarch','fields','filelock','filematrix','fitdistrplus','flashClust','flexclust','flexmix','flextable','fMultivar','FNN','fontawesome','forcats','foreach','forestplot','formatR','Formula','formula.tools','fpc','fresh','fs','fst','fstcore','furrr','futile.logger','futile.options','future','future.apply','GA','gamlss','gamlss.data','gamlss.dist','gamlss.tr','gapminder','gargle','gbutils','gclus','gdata','gdtools','generics','geomander','geometries','geometry','geonames','geos','geosphere','gert','getopt','GetoptLong','getPass','ggalluvial','GGally','ggalt','ggbeeswarm','ggbiplot','ggdendro','ggdist','ggeffects','ggforce','ggfortify','ggfun','ggiraph','ggmap','ggnetwork','ggnewscale','ggplot2','ggplotify','ggpmisc','ggpointdensity','ggpp','ggprism','ggpubr','ggraph','ggrastr','ggrepel','ggridges','ggseqlogo','ggsci','ggsignif','ggstats','ggtext','ggthemes','ggvenn','gh','git2r','gitcreds','glasso','gld','glmmTMB','glmnet','glmpca','GlobalOptions','globals','glue','GMMAT','gmodels','gmp','gnm','goftest','golem','googleAuthR','googledrive','googlesheets4','gower','GPfit','gplots','gProfileR','gprofiler2','graphlayouts','gridBase','gridExtra','gridGraphics','gridtext','grr','GSA','gsl','gson','gss','gsubfn','gsw','gt','gtable','gtools','gtsummary','GWASExactHW','h2o','hablar','hardhat','harmony','hash','haven','hdf5r','HDInterval','heatmaply','here','hexbin','hgnc','highr','Hmisc','hms','Hmsc','hoardr','homologene','htmlTable','htmltools','HTMLUtils','htmlwidgets','httpcode','httpuv','httr','httr2','huge','hunspell','hwriter','HyperG','ica','iCellR','ids','ieugwasr','igraph','imager','infer','infotheo','ini','inline','insight','interp','intervals','inum','ipred','irlba','irr','isdparser','Iso','isoband','iterators','iterpc','itertools','JADE','janeaustenr','janitor','jomo','jpeg','jquerylib','jsonify','jsonlite','jsonvalidate','kableExtra','kernlab','keyring','kinship2','km.ci','KMsurv','knitr','ks','labeling','labelled','Lahman','lambda.r','LambertW','lamW','LaplacesDemon','lars','latentcor','later','latex2exp','lattice','latticeExtra','lava','lavaan','lazyeval','leaps','LearnBayes','leiden','leidenAlg','leidenbase','lexicon','lhs','libcoin','libgeos','lidR','lifecycle','likert','limSolve','linprog','listenv','littler','lme4','lmerTest','lmodel2','lmom','Lmoments','lmtest','locfit','logger','logging','logistf','logNormReg','lokern','loo','lpSolve','lubridate','MAGEE','magic','magick','magrittr','mapdata','mapproj','maps','maptree','markdown','MASS','MatchIt','mathjaxr','Matrix','matrixcalc','MatrixModels','matrixStats','matrixTests','maxLik','maxstat','mboost','mclust','mcmc','MCMCpack','mdatools','memoise','MendelianRandomization','meta','metadat','metafor','mets','mgsub','mice','miceadds','microbenchmark','mime','miniUI','minqa','miscTools','missForest','missMDA','mitml','mitools','mixedCCA','mixtools','mnormt','modeldata','modelenv','ModelMetrics','modelr','modeltools','moments','moonBook','move','MplusAutomation','mr.raps','msigdbr','multcomp','multcompView','multicool','munsell','mvtnorm','N2R','NADA','naniar','nanoarrow','nanotime','naturalsort','NbClust','ncdf4','ncmeta','ndtv','nestedcv','netdiffuseR','network','networkD3','networkDynamic','neurobase','nleqslv','nlme','nloptr','NLP','nnet','nnls','nor1mix','norm','nortest','Nozzle.R1','numDeriv','nycflights13','oce','officer','openssl','openxlsx','operator.tools','optparse','orca','ordinal','oro.nifti','orthopolynom','outliers','pack','packrat','pagoda2','paletteer','pals','pan','pander','pandoc','parallelly','parameters','parsnip','partykit','patchwork','pbapply','pbivnorm','pbkrtest','pbmcapply','PBSmodelling','pcaPP','penalized','performance','PerformanceAnalytics','permute','phangorn','pheatmap','pillar','pingr','pixmap','pkgbuild','pkgconfig','pkgdown','pkgload','pkgmaker','plogr','plotly','plotrix','pls','plsRglm','plyr','png','PoissonMultinomial','polspline','polyclip','polynom','poolr','posterior','prabclus','pracma','praise','PredictABEL','prettyunits','princurve','prismatic','pROC','processx','prodlim','profvis','progress','progressr','proj4','promises','proto','protolite','proxy','ps','pscl','psych','Publish','pulsar','purrr','pwr','qap','qdapRegex','qgraph','qlcMatrix','qpdf','qqconf','qqman','quadprog','quantmod','quantreg','QuickJSR','qvcalc','R2HTML','R6','ragg','randomForest','ranger','RANN','rapidjsonr','rappdirs','rARPACK','raster','rasterVis','rbibutils','R.cache','rcmdcheck','RColorBrewer','rcompanion','Rcpp','RcppAnnoy','RcppArmadillo','RcppCCTZ','RcppDate','RcppDist','RcppEigen','RcppGSL','RcppHNSW','RcppInt64','RcppML','RcppNumerical','RcppParallel','RcppProgress','RcppRoll','RcppSpdlog','RcppThread','RcppTOML','RcppZiggurat','RCurl','Rdpack','reactable','reactR','readbitmap','readr','readxl','recipes','redist','redistmetrics','refseqR','registry','reldist','relimp','remaCor','rematch','rematch2','remote','remotes','rentrez','repmis','repr','reprex','reshape','reshape2','reticulate','rex','Rfast','rgl','RhpcBLASctl','rio','riskRegression','rjags','rJava','rjson','RJSONIO','rlang','rlas','rle','rly','RMariaDB','rmarkdown','rmdformats','rmdpartials','R.methodsS3','Rmpfr','rms','RMTstat','RNCEP','RNetCDF','rngtools','RNifti','robustbase','ROCR','R.oo','Rook','rootSolve','ropenblas','roxygen2','RPMM','RPostgres','rprojroot','rrBLUP','rrtable','rsample','rsconnect','RSpectra','RSQLite','rstan','rstanarm','rstantools','rstatix','rstudioapi','rsvd','Rtsne','Rttf2pt1','RUnit','R.utils','ruv','rvcheck','rversions','rvest','rvg','s2','sandwich','sass','scales','scattermore','scatterpie','scatterplot3d','scCATCH','sccore','scCustomize','scrime','sctransform','segmented','selectr','seriation','servr','sessioninfo','Seurat','SeuratObject','sf','sfheaders','sfsmisc','shadowtext','shape','shiny','shinyBS','shinydashboard','shinydashboardPlus','shinyjs','shinystan','shinythemes','shinyWidgets','sitmo','sjlabelled','sjmisc','sjPlot','sjstats','skimr','slam','slider','slippymath','sm','sn','sna','snakecase','snow','SnowballC','snowfall','sodium','som','SoupX','sourcetools','sp','spam','SparseM','sparsesvd','spatstat','spatstat.data','spatstat.explore','spatstat.geom','spatstat.linnet','spatstat.random','spatstat.sparse','spatstat.utils','spData','spdep','spdl','speedglm','splus2R','sqldf','SQUAREM','stabledist','stabs','StanHeaders','stars','startupmsg','statmod','statnet.common','stringdist','stringi','stringr','survey','survminer','survMisc','svglite','svUnit','synchronicity','sys','systemfonts','syuzhet','table1','tableone','taxa','taxonomizr','TeachingDemos','tensor','tensorA','terra','testit','testthat','texreg','textclean','textshape','textshaping','TFisher','tgp','TH.data','threejs','tibble','tictoc','tidybayes','tidygraph','tidymodels','tidync','tidyr','tidyselect','tidytext','tidytree','tidyverse','tiff','tigris','tiledb','timechange','timeDate','timereg','timeSeries','tinytex','tinytiger','TMB','tmvnsim','tokenizers','toOrdinal','triangle','triebeard','truncdist','truncnorm','tsne','TSP','TTR','tune','tweedie','tweenr','typed','tzdb','ucminf','umap','units','UpSetR','urlchecker','urltools','useful','usethis','utf8','uuid','uwot','V8','vcd','vcdExtra','vctrs','vegan','VennDiagram','verification','VGAM','vioplot','vipor','viridis','viridisLite','visdat','visNetwork','vroom','vtable','waiter','waldo','warp','webshot','wesanderson','wheatmap','whisker','withr','wk','workflows','workflowsets','writexl','xfun','xgboost','xlsx','xlsxjars','XML','xml2','xopen','xtable','xts','yaml','yardstick','yulab.utils','zCompositions','zeallot','zip','zoo','ztable'))"

## Install Bioconductor packages
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")' \
    && Rscript -e "BiocManager::install(c('affy','affyio','alabaster.base','alabaster.matrix','alabaster.ranges','alabaster.schemas','alabaster.se','ALL','annotate','AnnotationDbi','AnnotationFilter','AnnotationForge','AnnotationHub','AnnotationHubData','APAlyzer','apeglm','AUCell','bamsignals','basilisk','basilisk.utils','batchelor','BayesSpace','beachmat','Biobase','BiocBaseUtils','BiocCheck','BiocFileCache','BiocGenerics','BiocIO','BiocNeighbors','BiocParallel','BiocSingular','BiocStyle','BiocVersion','biocViews','bioDist','biomaRt','biomformat','Biostrings','biovizBase','BiSeq','bluster','BSgenome','BSgenome.Sscrofa.UCSC.susScr3','BSgenome.Sscrofa.UCSC.susScr3.masked','bsseq','bumphunter','Category','ccrepe','celldex','CHETAH','clusterProfiler','clustifyr','CoGAPS','ComplexHeatmap','ConsensusClusterPlus','dada2','DECIPHER','decontam','DEGreport','DelayedArray','DelayedMatrixStats','DeMixT','densvis','depmap','derfinder','derfinderHelper','DESeq2','DESpace','dir.expiry','DirichletMultinomial','discordant','dittoSeq','DMRcate','DMRcatedata','DNAcopy','DO.db','DOSE','DropletUtils','DSS','edgeR','EGSEA','EGSEAdata','EnhancedVolcano','ENmix','enrichplot','EnsDb.Hsapiens.v79','ensembldb','epiNEM','escape','ExperimentHub','ExperimentHubData','FastqCleaner','fastseg','FDb.InfiniumMethylation.hg19','fgsea','fishpond','FlowSorted.Blood.EPIC','gage','gdsfmt','genefilter','geneplotter','GENESIS','GENIE3','genomation','GenomeInfoDb','GenomeInfoDbData','GenomicAlignments','GenomicFeatures','GenomicFiles','GenomicRanges','GeomxTools','GeoMxWorkflows','GEOquery','ggtree','glmGamPoi','globaltest','GO.db','GOSemSim','GOstats','graph','graphite','GSEABase','GSVA','Gviz','GWASTools','gypsum','HDF5Array','HDO.db','hgu133a.db','hgu133plus2.db','HSMMSingleCell','HybridMTest','hypeR','IlluminaHumanMethylation450kanno.ilmn12.hg19','IlluminaHumanMethylation450kmanifest','IlluminaHumanMethylationEPICanno.ilm10b4.hg19','IlluminaHumanMethylationEPICmanifest','illuminaio','impute','interactiveDisplayBase','IRanges','karyoploteR','KEGGdzPathwaysGEO','KEGGgraph','KEGGREST','LEA','lfa','limma','Linnorm','lionessR','LOLA','LoomExperiment','lpsymphony','lumi','M3Drop','Maaslin2','MAGeCKFlute','MAST','MatrixGenerics','megadepth','metagenomeSeq','metapod','MethylAid','methylclock','methylclockData','methylKit','methylumi','microbiome','minet','minfi','minfiData','miQC','missMethyl','mixOmics','mnem','monocle','msa','MultiAssayExperiment','multiMiR','multtest','muscat','NanoStringNCTools','Nebulosa','nempi','OmnipathR','OrganismDbi','org.Cf.eg.db','org.Hs.eg.db','org.Mm.eg.db','org.Rn.eg.db','org.Sc.sgd.db','org.Ss.eg.db','orthogene','PADOG','PAIRADISE','pandaR','pathview','pcaMethods','PCAtools','phyloseq','planet','plyranges','preprocessCore','projectR','ProtGenerics','quantiseqr','quantsmooth','qusage','qvalue','RBGL','RCy3','reactome.db','ReactomePA','recount','regioneR','ResidualMatrix','Rgraphviz','rhdf5','rhdf5filters','Rhdf5lib','Rhtslib','rnaEditr','ROC','Rsamtools','Rsubread','rtracklayer','S4Arrays','S4Vectors','safe','ScaledMatrix','scater','scDblFinder','scde','schex','scHOT','scran','scry','scuttle','SeqArray','seqPattern','SeqVarTools','sesame','sesameData','ShortRead','siggenes','sincell','SingleCellExperiment','SingleR','singscore','slingshot','SNPRelate','snpStats','SparseArray','sparseMatrixStats','SpatialDecon','SpatialExperiment','splatter','STRINGdb','SummarizedExperiment','sva','TFBSTools','topGO','TrajectoryUtils','treeio','TreeSummarizedExperiment','tricycle','TxDb.Hsapiens.UCSC.hg19.knownGene','txdbmaker','TxDb.Sscrofa.UCSC.susScr3.refGene','tximeta','tximport','UCell','UCSC.utils','variancePartition','VariantAnnotation','vsn','wateRmelon','Wrench','XVector','zellkonverter','zlibbioc'))"

## Install CRAN packages requiring bioconductor pakcages as dependencies
RUN Rscript -e "install.packages(c('ADAPTS','codebook','conos','ggm','isva','metap','mutoss','NMF','pcalg','scGate','scMappR','Signac','restfulr','rliger','WGCNA'))"

## Install GitHub packages
RUN Rscript -e "remotes::install_github('satijalab/azimuth')" \
   && Rscript -e "devtools::install_github('stefpeschel/NetCoMi', dependencies = c('Depends', 'Imports', 'LinkingTo'))" \
    && Rscript -e "devtools::install_github('cnfoley/hyprcoloc')" \
    && Rscript -e "remotes::install_github('bnprks/BPCells/r')" \
    && Rscript -e "devtools::install_github('bnprks/BPCells')" \
    && Rscript -e "devtools::install_github('nghiavtr/BPSC')" \
    && Rscript -e "devtools::install_github('itsrainingdata/sparsebnUtils')" \
    && Rscript -e "devtools::install_github('itsrainingdata/ccdrAlgorithm')" \
    && Rscript -e "devtools::install_github('khodosevichlab/CellAnnotatoR')" \
    && Rscript -e "devtools::install_github('sqjin/CellChat')" \
    && Rscript -e "devtools::install_github('jokergoo/circlize')" \
    && Rscript -e "remotes::install_github('stan-dev/cmdstanr')" \
    && Rscript -e "devtools::install_github('cansysbio/ConsensusTME')" \
    && Rscript -e 'devtools::install_github("CBMR-Single-Cell-Omics-Platform/SCOPfunctions")' \
    && Rscript -e "devtools::install_github('https://github.com/CostaLab/CrossTalkeR', build_vignettes = TRUE)" \
    && Rscript -e "devtools::install_github('DARTH-git/dampack')" \
    && Rscript -e "devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')" \
    && Rscript -e "devtools::install_github('GfellerLab/EPIC', build_vignettes=TRUE)" \
    && Rscript -e "devtools::install_github('cole-trapnell-lab/garnett')" \
    && Rscript -e "devtools::install_github('JGCRI/gcamextractor')" \
    && Rscript -e "remotes::install_github('RubD/Giotto@cless')" \    
    && Rscript -e "devtools::install_github('ghcarlalan/graveler')" \
    && Rscript -e "remotes::install_github('omnideconv/immunedeconv')" \
    && Rscript -e "remotes::install_github('10XGenomics/loupeR')" \
    && Rscript -e "devtools::install_version('Matrix')" \
    && Rscript -e "remotes::install_github('dklinges9/mcera5')" \  
    && Rscript -e "devtools::install_github('ebecht/MCPcounter', ,ref='master', subdir='Source')" \
    && Rscript -e "devtools::install_github('ilyamaclean/microclima')" \
    && Rscript -e "devtools::install_github('cit-bioinfo/mMCP-counter')" \
    && Rscript -e "devtools::install_github('cole-trapnell-lab/monocle3')" \
    && Rscript -e "devtools::install_github('MRCIEU/MRInstruments')" \
    && Rscript -e "devtools::install_github('gqi/MRMix')" \
    && Rscript -e "devtools::install_github('rondolab/MR-PRESSO')" \
    && Rscript -e "devtools::install_github('immunogenomics/presto')" \
    && Rscript -e "remotes::install_github('carmonalab/ProjecTILs')" \
    && Rscript -e "devtools::install_github('jbisanz/qiime2R')" \
    && Rscript -e "devtools::install_github('cysouw/qlcMatrix')" \
    && Rscript -e "remotes::install_github('WSpiller/RadialMR')" \
    && Rscript -e "remotes::install_github('rmcelreath/rethinking')" \
    && Rscript -e "devtools::install_github('JGCRI/rgcam')" \
    && Rscript -e "devtools::install_github('JGCRI/rpackageutils')" \
    && Rscript -e "devtools::install_github('cellgeni/sceasy')" \
    && Rscript -e 'devtools::install_github("CBMR-Single-Cell-Omics-Platform/SCOPfunctions")' \
    && Rscript -e "devtools::install_github('igrabski/sc-SHC')" \
    && Rscript -e "devtools::install_github('satijalab/seurat-data')" \
    && Rscript -e "devtools::install_github('pinin4fjords/shinyngs')" \
    && Rscript -e "remotes::install_github('satijalab/seurat-wrappers')" \
    && Rscript -e "remotes::install_github('satijalab/seurat-wrappers')" \
    && Rscript -e "devtools::install_github('arc85/singleseqgset')" \
    && Rscript -e "devtools::install_github('EliGurarie/smoove')" \
    && Rscript -e "devtools::install_github('zdk123/SpiecEasi')" \ 
    && Rscript -e "devtools::install_github('GraceYoon/SPRING')" \
    && Rscript -e "remotes::install_github('carmonalab/STACAS')" \
    && Rscript -e "remotes::install_github('poisonalien/trackplot')" \
    && Rscript -e "remotes::install_github('MRCIEU/TwoSampleMR')" \
    && Rscript -e "devtools::install_github('timokelder/UNSEEN')" \
    && Rscript -e "devtools::install_github('velocyto-team/velocyto.R')" \
    && Rscript -e "devtools::install_github('dviraran/xCell')" \
    && Rscript -e "devtools::install_github('statOmics/zingeR')" 

## cellxgene.census
RUN pip install tiledb \
    && Rscript -e "install.packages('tiledb')" \
    && Rscript -e "install.packages('cellxgene.census',repos=c('https://chanzuckerberg.r-universe.dev', 'https://cloud.r-project.org'))"
 
## ProjecTILs
RUN Rscript -e "remotes::install_github('carmonalab/ProjecTILs')" 

## velocyto.R
RUN wget https://github.com/velocyto-team/velocyto.R/archive/refs/tags/0.6.tar.gz \
    && tar -xvf 0.6.tar.gz \
    && sed -i '48s/^/\/\//' velocyto.R-0.6/src/routines.cpp \
    && Rscript -e "install.packages('velocyto.R-0.6',repos=NULL,type='source')" \
    && rm 0.6.tar.gz  \
    && rm -rf velocyto.R-0.6


# Clean up package installation cache to reduce image size
RUN rm -rf /tmp/* \
    && rm -rf ~/.cache/R \
    && rm -rf /root/.R/ \
    && rm -rf /usr/local/lib/R/site-library/*/help \
    && find /usr/local/lib/R/site-library/ -type f -name "*.o" -delete \
    && find /usr/local/lib/R/site-library/ -type f -name "*.so" -exec strip --strip-unneeded {} \;
