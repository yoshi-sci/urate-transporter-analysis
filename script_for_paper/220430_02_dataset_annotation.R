setwd("")

#set the package
pacman::p_load(rgl, knitr, Matrix, scater, DropletUtils, scran, limma, 
               GSVA, GSEABase, genefilter, AnnotationDbi, org.Mm.eg.db, 
               org.Hs.eg.db, gplots, gdata, RColorBrewer, scales, reshape2,
               sctransform, Seurat, magrittr, dplyr, reticulate, ggplot2, 
               SeuratData, biomaRt, gProfileR, patchwork, cowplot, future)

plan()
# change the current plan to access parallelization
plan("multisession", workers = 4)
plan()
#fall back to 
plan("sequential")

options(future.globals.maxSize = 3000 * 1024^2)
#(Default: 500 * 1024 ^ 2 = 500 MiB)

#load("/Users/yoshi/R/scRNAseq/kidney/URAT1/Expression_profile/Physiological urate handling in human male/rawdata/220430_rawdata_afterbiomaRt.RData")
#
####################################################################################
###Make the list of human healthy kidney datasets
Healthy_62_GSE118184_Human_Male$ID <- "Ident.1"
Healthy_54_GSE131882_Human_Male$ID <- "Ident.2"
Healthy_62_GSE131882_Human_Male$ID <- "Ident.3"

Healthy_62_GSE118184_Human_Male$species <- "Human"
Healthy_54_GSE131882_Human_Male$species <- "Human"
Healthy_62_GSE131882_Human_Male$species <- "Human"

Healthy_62_GSE118184_Human_Male$sex <- "Male"
Healthy_54_GSE131882_Human_Male$sex <- "Male"
Healthy_62_GSE131882_Human_Male$sex <- "Male"

Healthy_62_GSE118184_Human_Male$GSE <- "GSE118184"
Healthy_54_GSE131882_Human_Male$GSE <- "GSE131882"
Healthy_62_GSE131882_Human_Male$GSE <- "GSE131882"

Healthy_62_GSE118184_Human_Male$disease <- "Healthy"
Healthy_54_GSE131882_Human_Male$disease <- "Healthy"
Healthy_62_GSE131882_Human_Male$disease <- "Healthy"

Healthy_Human_kidney_list <- list(Healthy_62_GSE118184_Human_Male, Healthy_54_GSE131882_Human_Male, 
                                  Healthy_62_GSE131882_Human_Male)

names(Healthy_Human_kidney_list) <- c("Healthy_62_GSE118184_Human_Male", "Healthy_54_GSE131882_Human_Male", 
                                      "Healthy_62_GSE131882_Human_Male")

