#set working directinary
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


#############################################################################################################
library('biomaRt')
Ctrl54M <- readRDS("./rawdata/GSE131882_RAW/control.s1.dgecounts")
mart <- useEnsembl(biomart = "ensembl", 
                           dataset = "hsapiens_gene_ensembl", 
                           mirror = "www")
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"), mirror = "asia")
umicounts <- Ctrl54M[["umicount"]][["exon"]][["all"]]
ensembl_gene_id <- rownames(umicounts)

umicounts_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                        values = ensembl_gene_id, mart= mart)
head(umicounts_list)
#write.csv(umicounts_list,"Ctrl54Mumicounts_list.csv",quote=FALSE, row.names=FALSE)
#write.csv(umicounts,"Ctrl54Mumicounts.csv",quote=FALSE, row.names=T)

for (i in 1:nrow(umicounts_list)){
  if (umicounts_list[i,"hgnc_symbol"] == "") umicounts_list[i,"hgnc_symbol"]  <- umicounts_list[i,"ensembl_gene_id"]
}
dim(umicounts_list)
ensembl_gene_id <- as.matrix(ensembl_gene_id)
colnames(ensembl_gene_id) <- c("ensembl_gene_id")
umicounts_list <- merge(umicounts_list , ensembl_gene_id, all.y = T)
dim(umicounts_list)

for (i in 1:nrow(umicounts_list)){
  if (is.na(umicounts_list[i,"hgnc_symbol"]) == T) umicounts_list[i,"hgnc_symbol"]  <- umicounts_list[i,"ensembl_gene_id"]
}

for (i in 1:nrow(umicounts_list)){
  if (umicounts_list[i,"ensembl_gene_id"] != ensembl_gene_id[i,"ensembl_gene_id"]) umicounts_list[i,] <- NA 
  if (umicounts_list[i,"ensembl_gene_id"] != ensembl_gene_id[i,"ensembl_gene_id"]) break
}
umicounts_list <- na.omit(umicounts_list)
dim(umicounts_list)

rownames(umicounts_list) <- 1:nrow(ensembl_gene_id)

View(umicounts@Dimnames[[1]])
umicounts@Dimnames[[1]] <- umicounts_list$hgnc_symbol

Ctrl54M <- CreateSeuratObject(counts = umicounts, min.cells = 3, min.features = 200, project = "Human_adult_kidney")
Ctrl54M
##28375 features across 3740 samples

rm(ensembl_gene_id)
rm(umicounts)
rm(umicounts_list)
rm(i)
rm(mart)

#############################22222222222222222222222222222222222########################################
library('biomaRt')
Ctrl62M <- readRDS("./rawdata/GSE131882_RAW/control.s2.dgecounts")
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "www")
umicounts <- Ctrl62M[["umicount"]][["exon"]][["all"]]
ensembl_gene_id <- rownames(umicounts)

umicounts_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values = ensembl_gene_id, mart= mart)

#write.csv(umicounts_list,"Ctrl62Mumicounts_list.csv",quote=FALSE, row.names=FALSE)
#write.csv(umicounts,"Ctrl62Mumicounts.csv",quote=FALSE, row.names=T)

for (i in 1:nrow(umicounts_list)){
  if (umicounts_list[i,"hgnc_symbol"] == "") umicounts_list[i,"hgnc_symbol"]  <- umicounts_list[i,"ensembl_gene_id"]
}

ensembl_gene_id <- as.matrix(ensembl_gene_id)
colnames(ensembl_gene_id) <- c("ensembl_gene_id")
umicounts_list <- merge(umicounts_list , ensembl_gene_id, all.y = T)
dim(umicounts_list)

for (i in 1:nrow(umicounts_list)){
  if (is.na(umicounts_list[i,"hgnc_symbol"]) == T) umicounts_list[i,"hgnc_symbol"]  <- umicounts_list[i,"ensembl_gene_id"]
}

for (i in 1:nrow(umicounts_list)){
  if (umicounts_list[i,"ensembl_gene_id"] != ensembl_gene_id[i,"ensembl_gene_id"]) umicounts_list[i,] <- NA 
  if (umicounts_list[i,"ensembl_gene_id"] != ensembl_gene_id[i,"ensembl_gene_id"]) break
}
umicounts_list <- na.omit(umicounts_list)
dim(umicounts_list)
for (i in 1:nrow(umicounts_list)){
  if (umicounts_list[i,"ensembl_gene_id"] != ensembl_gene_id[i,"ensembl_gene_id"]) umicounts_list[i,] <- NA 
  if (umicounts_list[i,"ensembl_gene_id"] != ensembl_gene_id[i,"ensembl_gene_id"]) break
}
umicounts_list <- na.omit(umicounts_list)
dim(umicounts_list)

rownames(umicounts_list) <- 1:nrow(ensembl_gene_id)

View(umicounts@Dimnames[[1]])
umicounts@Dimnames[[1]] <- umicounts_list$hgnc_symbol

Ctrl62M <- CreateSeuratObject(counts = umicounts, min.cells = 3,  min.features = 200, project = "Human_adult_kidney")

rm(ensembl_gene_id)
rm(umicounts)
rm(umicounts_list)
rm(i)
rm(mart)
Ctrl62M
##22205 features across 4043 samples 

#############################################################################################################
Healthy_54_GSE131882_Human_Male <- Ctrl54M
rm(Ctrl54M)
Healthy_62_GSE131882_Human_Male <- Ctrl62M
rm(Ctrl62M)

Human_kidney1 <- read.table(
  "./rawdata/GSE118184/GSE118184human.dge.txt", row.names = 1, header = T)%>%
  as.matrix() %>%
  CreateSeuratObject(min.cells = 3, min.features = 200)
Human_kidney1
#20456 features across 4524 samples
Healthy_62_GSE118184_Human_Male　<- Human_kidney1
#
rm(Human_kidney1)

save.image("./Object/rawdata_afterbiomaRt.RData")
######################################################################################3333


