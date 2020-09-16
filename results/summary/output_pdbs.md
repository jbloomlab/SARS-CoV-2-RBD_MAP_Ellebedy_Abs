Output pdbs with escape scores as b-factor
================
8/7/2020

-   [Setup](#setup)

This notebook annotates pdb structures by escape metrics by filling in b
factor, enabling visualization in PyMol.

    #list of packages to install/load
    packages = c("yaml","data.table","tidyverse","bio3d","knitr")
    #install any packages not already installed
    installed_packages <- packages %in% rownames(installed.packages())
    if(any(installed_packages == F)){
      install.packages(packages[!installed_packages])
    }
    #load packages
    lapply(packages, library, character.only=T)

    ## [[1]]
    ## [1] "yaml"      "stats"     "graphics"  "grDevices" "utils"     "datasets" 
    ## [7] "methods"   "base"     
    ## 
    ## [[2]]
    ## [1] "data.table" "yaml"       "stats"      "graphics"   "grDevices" 
    ## [6] "utils"      "datasets"   "methods"    "base"      
    ## 
    ## [[3]]
    ##  [1] "forcats"    "stringr"    "dplyr"      "purrr"      "readr"     
    ##  [6] "tidyr"      "tibble"     "ggplot2"    "tidyverse"  "data.table"
    ## [11] "yaml"       "stats"      "graphics"   "grDevices"  "utils"     
    ## [16] "datasets"   "methods"    "base"      
    ## 
    ## [[4]]
    ##  [1] "bio3d"      "forcats"    "stringr"    "dplyr"      "purrr"     
    ##  [6] "readr"      "tidyr"      "tibble"     "ggplot2"    "tidyverse" 
    ## [11] "data.table" "yaml"       "stats"      "graphics"   "grDevices" 
    ## [16] "utils"      "datasets"   "methods"    "base"      
    ## 
    ## [[5]]
    ##  [1] "knitr"      "bio3d"      "forcats"    "stringr"    "dplyr"     
    ##  [6] "purrr"      "readr"      "tidyr"      "tibble"     "ggplot2"   
    ## [11] "tidyverse"  "data.table" "yaml"       "stats"      "graphics"  
    ## [16] "grDevices"  "utils"      "datasets"   "methods"    "base"

    knitr::opts_chunk$set(echo = T)
    knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

    #read in config file
    config <- read_yaml("config.yaml")

    #read in file giving concordance between RBD numbering and SARS-CoV-2 Spike numbering
    RBD_sites <- data.table(read.csv(file=config$RBD_sites,stringsAsFactors=F))

    #make output directories
    if(!file.exists(config$pdb_outputs_dir)){
      dir.create(file.path(config$pdb_outputs_dir))
    }

Session info for reproducing environment:

    sessionInfo()

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /app/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_haswellp-r0.3.7.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] knitr_1.26        bio3d_2.4-0       forcats_0.4.0     stringr_1.4.0    
    ##  [5] dplyr_0.8.3       purrr_0.3.3       readr_1.3.1       tidyr_1.0.0      
    ##  [9] tibble_3.0.1      ggplot2_3.3.0     tidyverse_1.3.0   data.table_1.12.8
    ## [13] yaml_2.2.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.5 xfun_0.11        haven_2.2.0      lattice_0.20-38 
    ##  [5] colorspace_1.4-1 vctrs_0.2.4      generics_0.0.2   htmltools_0.4.0 
    ##  [9] rlang_0.4.5      pillar_1.4.3     glue_1.3.1       withr_2.1.2     
    ## [13] DBI_1.1.0        dbplyr_1.4.2     modelr_0.1.5     readxl_1.3.1    
    ## [17] lifecycle_0.2.0  munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0
    ## [21] rvest_0.3.5      evaluate_0.14    parallel_3.6.2   fansi_0.4.0     
    ## [25] broom_0.5.6      Rcpp_1.0.3       scales_1.1.0     backports_1.1.5 
    ## [29] jsonlite_1.6     fs_1.3.1         hms_0.5.2        digest_0.6.23   
    ## [33] stringi_1.4.3    grid_3.6.2       cli_2.0.0        tools_3.6.2     
    ## [37] magrittr_1.5     crayon_1.3.4     pkgconfig_2.0.3  ellipsis_0.3.0  
    ## [41] xml2_1.2.2       reprex_0.3.0     lubridate_1.7.4  assertthat_0.2.1
    ## [45] rmarkdown_2.0    httr_1.4.1       rstudioapi_0.10  R6_2.4.1        
    ## [49] nlme_3.1-143     compiler_3.6.2

Setup
-----

Read in tables of escape scores, the list of samples, and the structure
that we want to map these scores to.

    pdb <- read.pdb(file=config$pdb)

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

    escape <- data.table(read.csv(file=config$escape_fracs, stringsAsFactors=F))
    #keep only the average scores
    escape <- escape[library=="average",]

To visualize site-wise mutational sensitivity on the 3D structure, we
output `.pdb` files for the ACE2-bound RBD structure in which we replace
the B factor column with escape metrics for each antibody. In PyMol, we
can then visualize the structure in various ways, colored by these
metrics of escape. In particular, we will color by the sum of per-site
mutation escape fractions, as well as the maximum escape fraction of
mutations at each site. For sites where no mutations have filtered
escape scores, because all mutations were filtered out due to being too
deleterious, we fill the b-factor value as -1 to enable collapsing to 0
or callout as a separate class, depending how we choose to color sites
for different visualizations.

    #get list of the experimental scores we'll output
    samples <- unique(escape$selection)

    #iterate through samples, replacing b-factor with the total escape fraction of muts at a site, or the max escape at a site
    #for sites where all mutations are deleterious and were filtered before escape expt, ascribe -1 for total escape (b/c no mut is tolerated)
    for(s in samples){
      # pdb_tot <- pdb #for replacement with total escape fraction at a site
      # pdb_max <- pdb #for replacement with maximum escape fraction at a site
      # pdb_tot$atom$b <- NA
      # pdb_max$atom$b <- NA
      b_tot <- rep(0, length(pdb$atom$b))
      b_max <- rep(0, length(pdb$atom$b))
      for(i in 1:nrow(pdb$atom)){
        if(pdb$atom$chain[i]=="E"){
          res <- pdb$atom$resno[i]
          escape_scores <- escape[selection==s & protein_site==res,mut_escape_frac_epistasis_model]
          if(length(escape_scores)>0){
            b_tot[i] <- sum(escape_scores,na.rm=T)
            b_max[i] <- max(escape_scores,na.rm=T)
          }else{
            b_tot[i] <- -1
            b_max[i] <- -1
          }
        }
      }
      write.pdb(pdb=pdb,file=paste(config$pdb_outputs_dir,"/",s,"_6m0j_total_escape.pdb",sep=""), b=b_tot)
      write.pdb(pdb=pdb,file=paste(config$pdb_outputs_dir,"/",s,"_6m0j_max_escape.pdb",sep=""), b=b_max)
    }

For loading in Chimera, it is easier if all structures annotated by
max-escape are already normalized 0 to 1. Therefore, we also output a
set of pdbs colored by max escape that are pre-scaled from 0 to 1.

    #create subdirectory
    if(!file.exists(paste(config$pdb_outputs_dir,"/rescale_0_1/",sep=""))){
      dir.create(file.path(paste(config$pdb_outputs_dir,"/rescale_0_1/",sep="")))
    }

    #iterate through samples, replacing b-factor with the total escape fraction of muts at a site, or the max escape at a site
    #for sites where all mutations are deleterious and were filtered before escape expt, ascribe -1 for total escape (b/c no mut is tolerated)
    for(s in samples){
      b_max <- rep(0, length(pdb$atom$b))
      for(i in 1:nrow(pdb$atom)){
        if(pdb$atom$chain[i]=="E"){
          res <- pdb$atom$resno[i]
          escape_scores <- escape[selection==s & protein_site==res,mut_escape_frac_epistasis_model]
          if(length(escape_scores)>0){
            b_max[i] <- max(escape_scores,na.rm=T)
          }else{
            b_max[i] <- 0
          }
        }
      }
      b_max <- b_max/(max(b_max))
      write.pdb(pdb=pdb,file=paste(config$pdb_outputs_dir,"/rescale_0_1/",s,"_6m0j_max_escape.pdb",sep=""), b=b_max)
    }

Though we will want more elaborate series of commands to codify our
visualization of these RBD structures colored by escape, the series of
commands below, when executed in a PyMol session with one of these PDBs
open, will color the RBD surface according to escape scores.

For example, to normalize each structure colored by the max mut effect,
we might want to have a white to red scale from 0 to 1:

    create RBD, chain E
    hide all; show cartoon, chain A; color gray20, chain A
    show surface, RBD; spectrum b, white red, RBD, minimum=0, maximum=1

For something like total escape, maybe we want each structure normalized
to the maximum total escape in that structure, in which case we can just
leave the maximum argument empty.

    create RBD, chain E
    hide all; show cartoon, chain A; color gray20, chain A
    show surface, RBD; spectrum b, white red, RBD, minimum=0