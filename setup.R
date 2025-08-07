# the purpose of this script is to load all libraries in R and to setup the appropriate 
# directories for storing the results

library(rtracklayer)
library(Seurat)
library(devtools)
library(presto)
library(tidyverse)
library(ggrepel)
library(Matrix)
library(cowplot)
library(grid)
library(patchwork)
library(gridExtra)
library(gprofiler2)
library(EnhancedVolcano)
library(glmGamPoi)


rootdir="/media/inserm-root/P4/Single_cell_rna_seq_analyses_eline/sc3samples_whitelist"
setwd(rootdir)  # Set this to correct location

dir.create("./flames", recursive = TRUE, showWarnings = FALSE)
dir.create("./QC", recursive = TRUE, showWarnings = FALSE)
dir.create("./counts", recursive = TRUE, showWarnings = FALSE)
dir.create("./merge_integration", recursive = TRUE, showWarnings = FALSE)
dir.create("./analyses", recursive = TRUE, showWarnings = FALSE)
