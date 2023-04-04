setwd("")
#set the package
pacman::p_load(rgl, knitr, Matrix, scater, DropletUtils, scran, limma, 
               GSVA, GSEABase, genefilter, AnnotationDbi, org.Mm.eg.db, 
               org.Hs.eg.db, gplots, gdata, RColorBrewer, scales, reshape2,
               sctransform, Seurat, magrittr, dplyr, reticulate, ggplot2, 
               SeuratData, biomaRt, gProfileR, patchwork, cowplot, future,ggsci)

plan()
# change the current plan to access parallelization
plan("multisession", workers = 4)
plan()
#fall back to 
plan("sequential")

options(future.globals.maxSize = 3000 * 1024^2)
#(Default: 500 * 1024 ^ 2 = 500 MiB)

######################################################################################
#Check the QC
#calculate percent.mt
for (i in 1:length(Healthy_Human_kidney_list)) {
  Healthy_Human_kidney_list[[i]] <- PercentageFeatureSet(Healthy_Human_kidney_list[[i]],　pattern = "^MT-", col.name = "percent.mt")
}

for (i in 1:length(Healthy_Human_kidney_list)) {
  ploti <- VlnPlot(Healthy_Human_kidney_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(ploti)
}
for (i in 1:length(Healthy_Human_kidney_list)) {
  ploti <-FeatureScatter(Healthy_Human_kidney_list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  print(ploti)
}
for (i in 1:length(Healthy_Human_kidney_list)) {
  ploti <-FeatureScatter(Healthy_Human_kidney_list[[i]], feature1 = "nFeature_RNA", feature2 = "percent.mt")
  print(ploti)
}
for (i in 1:length(Healthy_Human_kidney_list)) {
  ploti <-FeatureScatter(Healthy_Human_kidney_list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(ploti)
}

rm(ploti)
#subset the cells having good quality 
Healthy_Human_kidney_list[[1]] <- subset(Healthy_Human_kidney_list[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Healthy_Human_kidney_list[[2]] <- subset(Healthy_Human_kidney_list[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Healthy_Human_kidney_list[[3]] <- subset(Healthy_Human_kidney_list[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)

# normalize data with SCTransform()
for (i in 1:length(Healthy_Human_kidney_list)) {
  Healthy_Human_kidney_list[[i]] <- SCTransform(Healthy_Human_kidney_list[[i]],
                                                assay = 'RNA', new.assay.name = 'SCT',
                                                vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'))
}

#calculate cell cycle score
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cc.genes <- c(s.genes, g2m.genes)

#check the genes which  are used in cell cycle score analysis
cc.genes[cc.genes %in% rownames(Healthy_54_GSE131882_Human_Male)] #the cell cycle gene which match the gene in my dataset
cc.genes[!cc.genes %in% rownames(Healthy_54_GSE131882_Human_Male)] #the cell cycle gene which do not match the gene in my dataset

for (i in 1:length(Healthy_Human_kidney_list)) {
  Healthy_Human_kidney_list[[i]] <- CellCycleScoring(Healthy_Human_kidney_list[[i]], 
                                                     s.features = s.genes, g2m.features = g2m.genes,
                                                     assay = 'SCT',set.ident = TRUE) 
}
for (i in 1:length(Healthy_Human_kidney_list)) {
  Healthy_Human_kidney_list[[i]]$CC.Difference <- Healthy_Human_kidney_list[[i]]$S.Score - Healthy_Human_kidney_list[[i]]$G2M.Score
}
rm(s.genes)
rm(g2m.genes)
rm(cc.genes)

# normalise again but this time including also the cell cycle scores
for (i in 1:length(Healthy_Human_kidney_list)) {
  Healthy_Human_kidney_list[[i]] <- SCTransform(Healthy_Human_kidney_list[[i]],
                                                assay = 'RNA', new.assay.name = 'SCT',
                                                vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
}

########################################################################################################################
#Perform an integrated analysis
Healthy_Human_kidney_features <- SelectIntegrationFeatures(object.list = Healthy_Human_kidney_list, nfeatures = 3000)
Healthy_Human_kidney_list <- PrepSCTIntegration(object.list = Healthy_Human_kidney_list, anchor.features = Healthy_Human_kidney_features, 
                                                verbose = T)
Healthy_Human_kidney_anchors <- FindIntegrationAnchors(object.list = Healthy_Human_kidney_list, normalization.method = "SCT", 
                                                       anchor.features = Healthy_Human_kidney_features, verbose = FALSE)

Healthy_Human_kidney_integrated <- IntegrateData(anchorset = Healthy_Human_kidney_anchors, normalization.method = "SCT", 
                                                 verbose = T) 

Healthy_Human_kidney_integrated <- RunPCA(object = Healthy_Human_kidney_integrated, npcs = 200, verbose = F)

Healthy_Human_kidney_integrated <- JackStraw(object = Healthy_Human_kidney_integrated, num.replicate = 100, dims = 50) %>%
  ScoreJackStraw(dims = 1:50)
JackStrawPlot(object = Healthy_Human_kidney_integrated, dims = 1:30)
ElbowPlot(object = Healthy_Human_kidney_integrated, ndims = 30)

Healthy_Human_kidney_integrated <- RunUMAP(object = Healthy_Human_kidney_integrated, reduction = "pca", dims = 1:200)
Healthy_Human_kidney_integrated <- FindNeighbors(Healthy_Human_kidney_integrated, dims = 1:200, verbose = FALSE)
Healthy_Human_kidney_integrated <- FindClusters(Healthy_Human_kidney_integrated, verbose = FALSE)

DimPlot(Healthy_Human_kidney_integrated, reduction = "umap", label = T, split.by = "GSE")
DimPlot(Healthy_Human_kidney_integrated, reduction = "umap", label = T)
VlnPlot(Healthy_Human_kidney_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Healthy_Human_kidney_integrated, features = c("S.Score", "G2M.Score", "CC.Difference"), ncol = 3)

rm(Healthy_Human_kidney_features)
rm(Healthy_Human_kidney_anchors)

#DefaultAssay(object = Healthy_Human_kidney_integrated) <- "integrated"

# Create dummy new assay to demo switching default assays
new.assay <- Healthy_Human_kidney_integrated[["RNA"]]

Key(object = new.assay) <- "DefaultRNA_"
Healthy_Human_kidney_integrated[["DefaultRNA"]] <- new.assay
# switch default assay to DefaultRNA
DefaultAssay(object = Healthy_Human_kidney_integrated) <- "DefaultRNA"
DefaultAssay(object = Healthy_Human_kidney_integrated)
rm(new.assay)

#return the assay to RNA and logNormalization
DefaultAssay(Healthy_Human_kidney_integrated) <- "RNA"
Healthy_Human_kidney_integrated <- NormalizeData(Healthy_Human_kidney_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
   #stored in Healthy_Human_kidney_integrated[["RNA"]]@data -->  Vlnplot use the lognormalization data as default.  
DefaultAssay(object = Healthy_Human_kidney_integrated)

########################################################################################################
DefaultAssay(object = Healthy_Human_kidney_integrated) <- "DefaultRNA"
markers.to.plot <- c("NPHS1","NPHS2", "PODXL",
                     "CUBN","LRP2","CRYAB",
                     "SLC12A1","UMOD" ,
                     "SLC12A3","CALB1","ATP1B3","AQP2", 
                     "SLC4A1", "AQP6","SLC26A7","SLC26A4",
                     "EMCN","PLAT","ITGA8","PDGFRB","CFH", "CLDN1")

DotPlot(Healthy_Human_kidney_integrated, features = rev(markers.to.plot), cols = c("blue", "red"),
        col.min = 0,col.max = 4.0, dot.scale = 6, dot.min = 0.05
) + RotatedAxis() + labs(title = "Cluster annotation") +ylab("Cluster") +xlab("Marker gene")
rm(markers.to.plot)

#Rename the cluster name
DefaultAssay(Healthy_Human_kidney_integrated) <- "RNA"
Healthy_Human_kidney_integrated <- NormalizeData(Healthy_Human_kidney_integrated, normalization.method = "LogNormalize", scale.factor = 10000)

mother_population <- Healthy_Human_kidney_integrated
mother_population <- RenameIdents(mother_population, 
                                                `0` = "PTb", `1` = "LOH (AL)", `2` = "DCT", 
                                                `3` = "PTa", `4` = "PC", `5` = "PTc", 
                                                `6` = "LOH (AL)", `7` = "PC", `8` = "ICA",
                                                `9` = "EDC", `10`= "ICB", `11`= "Podo",
                                                `12`= "PEC", `13` = "LOH (DL)", `14` = "CNT",
                                                `15` = "MGC")
Idents(mother_population) <- factor(Idents(mother_population), 
                                    levels = unique(c("Podo","PTa", "PTb", "PTc",
                                                      "LOH (DL)","LOH (AL)","DCT","CNT", 
                                                      "PC","ICA", "ICB", "EDC","MGC","PEC")))

mother_population <- Healthy_Human_kidney_integrated
markers.to.plot <- c("NPHS1","NPHS2", "PODXL",
                     "CUBN","LRP2","CRYAB",
                     "SLC12A1","UMOD" ,
                     "SLC12A3","CALB1","ATP1B3","AQP2", 
                     "SLC4A1", "AQP6","SLC26A7","SLC26A4",
                     "EMCN","PLAT","ITGA8","PDGFRB","CFH", "CLDN1")
DotPlot(mother_population, features = rev(markers.to.plot), cols = c("blue", "red"),
        col.min = 0,col.max = 4.0, dot.scale = 6, dot.min = 0.05
) + RotatedAxis() +ylab("Region") +xlab("Marker gene")
rm(markers.to.plot)

#Figure 1B & Supplementary Figure 1A 
Col_cluster <- c("#FF4B00","#F6AA00","#F6AA00","#F6AA00",
                 "#03AF7A","#03AF7A","#4DC4FF","#4DC4FF","#005AFF",
                 "#005AFF","#005AFF","#84919E","#84919E","#84919E")

DimPlot(mother_population, reduction = "umap", label = F, cols = Col_cluster)  + 
  FontSize(x.title = 10, y.title = 10)
DimPlot(mother_population, reduction = "umap", label = T, cols = Col_cluster, label.size = 3) + labs(title = "Human adult kidney") + 
  FontSize(x.title = 10, y.title = 10) + ggmin::theme_powerpoint()
rm(Col_cluster)

#Supplementary Figure 1C
Col_cluster <- c("#FF4B00","#F6AA00","#F6AA00","#F6AA00",
                 "#03AF7A","#03AF7A","#4DC4FF","#4DC4FF","#005AFF",
                 "#005AFF","#005AFF","#84919E","#84919E","#84919E")
DimPlot(mother_population, reduction = "umap", label = F, split.by = "ID", cols = Col_cluster) + 
  FontSize(x.title = 10, y.title = 10) + ggmin::theme_powerpoint()
rm(Col_cluster)
FeatureScatter(Healthy_Human_kidney_integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
               group.by = "GSE", pt.size = 0.1 ,cols = c("#00B0F6","#00B0F6")) + labs(title = "Adult human male kidney")  + ggmin::theme_powerpoint()


#"Podo"                   283
#"PT_S1"/"PT_S2"/"PT_S3"  1000/1686/767 = 3453
#"LOH (DL)"               257
#"LOHa (AL)"/"LOHb (AL)"  1430/563 = 1993
#"DCT"                    1204
#"CNT"                    206
#"PCa"/"PCb"              486/818 =1304
#"ICA"/"ICB"              467/287 =754
#"EDC"/"MGC"/"PEC"        291/75/260 = 626

#Whole_kidney Supplementary Figure 1B
Cluster_name <- c("Podocyte","PT","LOH (DL)","LOH (AL)","DCT","CNT","PC","IC","Others")
Cluster_ncells <- c(283, 3453, 257, 1993, 1204, 206, 1304, 754, 626)
Cluster_sex <- c("Male","Male","Male","Male","Male","Male","Male","Male","Male")

Cluster_proportion <- data.frame(Cluster = Cluster_name, ncell = Cluster_ncells, sex = Cluster_sex)
Cluster_proportion$Cluster <- factor(Cluster_proportion$Cluster, levels=c(Cluster_name))

g <- ggplot(Cluster_proportion, aes(x = sex, y = ncell, fill = Cluster)) + 
  geom_bar(stat = "identity",position = "fill")+ 
  scale_y_continuous(labels = percent)+ 
  scale_fill_manual(values = c("#FF82B2", "#FFFF88", "light green","#00CC00", "light blue","light blue", "#8EB8FF","#B384FF", "black"))
plot(g)+ theme_bw()
rm(g)
rm(Cluster_proportion)
rm(Cluster_name)
rm(Cluster_sex)
rm(Cluster_ncells)

#PT-Annotation
#PT_S1-malker:"SLC2A2", "SLC5A2"
#PT_S1-malker:"SLC25A1", "CA4"
#PT_S3-maker:"SLC3A1", "SLC7A13"
Genename_focused_Chr <- list("SLC5A1","SLC5A2", "SLC36A2", "SLC6A13","PTH1R","SLC2A1","SLC2A2","CASR")
Titlename_Chr <- list("SLC5A1/SGLT1","SLC5A2/SGLT2","SLC36A2/PAT2","SLC6A13/GAT2","PTH1R/PTH1R","SLC2A1/GLUT1","SLC2A2/GLUT2","CASR/CaSR")
PT_color <- c("#FFCC66","#FFFF00", "#CCCC33","#008800")
for (i in 1:length(Genename_focused_Chr)) {
  print(VlnPlot(PT_Healthy_Human_kidney_integrated, features = Genename_focused_Chr[[i]],
                cols = PT_color, pt.size = 0.1)+  
          coord_cartesian(ylim = c(0, 5.0)) +labs(title = Titlename_Chr[[i]]))+ 
    ggmin::theme_powerpoint()
}
for (i in 1:length(Genename_focused_Chr)) {
  print(VlnPlot(Healthy_Human_kidney_integrated, features = Genename_focused_Chr[[i]],
                pt.size = 0.1)+  
          coord_cartesian(ylim = c(0, 4.5)) +labs(title = Titlename_Chr[[i]]))
}
rm(Genename_focused_Chr)
rm(Titlename_Chr)

#Annotate the cluster
Healthy_Human_kidney_integrated <- RenameIdents(Healthy_Human_kidney_integrated, 
                                                `0` = "PT_S2", `1` = "LOH (AL)", `2` = "DCT", 
                                                `3` = "PT_S1", `4` = "PC", `5` = "PT_S3", 
                                                `6` = "LOH (AL)", `7` = "PC", `8` = "ICA",
                                                `9` = "EDC", `10`= "ICB", `11`= "Podo",
                                                `12`= "PEC", `13` = "LOH (DL)", `14` = "CNT",
                                                `15` = "MGC")
Idents(Healthy_Human_kidney_integrated) <- factor(Idents(Healthy_Human_kidney_integrated), 
                                                  levels = unique(c("Podo","PT_S1", "PT_S2", "PT_S3",
                                                                    "LOH (DL)","LOH (AL)","DCT","CNT", "PC",
                                                                    "ICA", "ICB", "EDC","MGC","PEC")))
Col_cluster <- c("#FF4B00","#F6AA00","#F6AA00","#F6AA00",
                 "#03AF7A","#03AF7A","#4DC4FF","#4DC4FF","#005AFF",
                 "#005AFF","#005AFF","#84919E","#84919E","#84919E")
DimPlot(Healthy_Human_kidney_integrated, reduction = "umap", label = F, 
        cols = Col_cluster)
DimPlot(Healthy_Human_kidney_integrated, reduction = "umap", label = T, 
        cols = Col_cluster)
rm(Col_cluster)

#Marker gene plot Figure 1A
markers.to.plot <- c("NPHS1","NPHS2", "PODXL",
                     "CUBN","LRP2","CRYAB",
                     "SLC12A1","UMOD" ,
                     "SLC12A3","CALB1","ATP1B3","AQP2", 
                     "SLC4A1", "AQP6","SLC26A7","SLC26A4",
                     "EMCN","PLAT","ITGA8","PDGFRB","CFH", "CLDN1")

DotPlot(Healthy_Human_kidney_integrated, features = rev(markers.to.plot), cols = c("blue", "red"),
        col.min = 0,col.max = 4.0, dot.scale = 6, dot.min = 0.05
) + RotatedAxis() + labs(title = "Cluster annotation") +ylab("Region") +xlab("Marker gene")
rm(markers.to.plot)

#Supplementary Figure 1B & Supplementary Table2
#Cluster proportion
Cluster_name <- list("Podo","PT_S1", "PT_S2", "PT_S3",
                     "LOH (DL)","LOH (AL)","DCT","CNT", "PC", 
                     "ICA", "ICB", "EDC","MGC","PEC")
dim(mother_population) 
for (i in 1:length(Cluster_name)) {
  print(dim(subset(mother_population, idents = Cluster_name[[i]])))
}
rm(Cluster_name)
#34145   283
#34145  1000
#34145  1686
#34145   767
#34145   257
#34145  1993
#34145  1204
#34145   206
#34145  1304
#34145   467
#34145   287
#34145   291
#34145    75
#34145   260

#set the object 
PT_Healthy_Human_kidney_integrated <- subset(Healthy_Human_kidney_integrated,  idents = c("PT_S1", "PT_S2", "PT_S3","LOH (DL)"))
dim(Healthy_Human_kidney_integrated )
Idents(PT_Healthy_Human_kidney_integrated) <- factor(Idents(PT_Healthy_Human_kidney_integrated), 
                                                  levels = unique(c("PT_S1", "PT_S2", "PT_S3","LOH (DL)")))
PT_Healthy_Human_kidney_integrated <- RenameIdents(PT_Healthy_Human_kidney_integrated, 
                                    "PT_S1" = "S1", "PT_S2" = "S2", "PT_S3" = "S3", 
                                    "LOH (DL)" = "DL")
Idents(PT_Healthy_Human_kidney_integrated) <- factor(Idents(PT_Healthy_Human_kidney_integrated), 
                                    levels = unique(c("S1","S2","S3","DL")))

#34145 10080
dim(PT_Healthy_Human_kidney_integrated)
#34145  3710
# PT_color:("#E58700","#C99800","#A3A500","#6BB100")

#Figure 1C
markers.to.plot <-　c("SLC5A2","SLC2A2","SLC3A1","SLC7A13")
mother_population <- subset(PT_Healthy_Human_kidney_integrated,  idents = c("S1","S2","S3"))
DotPlot(mother_population, features = rev(markers.to.plot), cols = c("blue", "red"),
        col.min = 0,col.max = 1.0, dot.scale = 6,  scale.min = 5 , scale.max = 80,
) + RotatedAxis()  +ylab("Region") +xlab("Marker gene")
rm(markers.to.plot)
rm(mother_population)
############################################################################################################
DefaultAssay(Healthy_Human_kidney_integrated) <- "RNA"

#Figure 1D
markers.to.plot <- c("ABCG2","ABCC4","SLC2A9", "SLC17A1","SLC17A3","SLC22A6","SLC22A7","SLC22A8","SLC22A11","SLC22A12","SLC22A13")
DotPlot(Healthy_Human_kidney_integrated, features = rev(markers.to.plot), cols = c("blue", "red"),
        col.min = 0,col.max = 2.0, dot.scale = 6, dot.min = 0.05
) + RotatedAxis()+ylab("Region") +xlab("Urate transporter") + FontSize(x.title = 12, y.title = 12)
rm(markers.to.plot)

#Supplementary Figure 5B & Supplementary Table 3
#PT_LOH(DL)
Cluster_name <- c("S1","S2","S3","DL")
Cluster_ncells <- c(1000, 1686, 767, 257)
Cluster_sex <- c("Male","Male","Male","Male")
PT_color <- c("#FFCC66","#FFFF00", "#CCCC33","#008800")

Cluster_proportion <- data.frame(Cluster = Cluster_name, ncell = Cluster_ncells, sex = Cluster_sex)
Cluster_proportion$Cluster <- factor(Cluster_proportion$Cluster, levels=c(Cluster_name))

g <- ggplot(Cluster_proportion, aes(x = sex, y = ncell, fill = Cluster)) + 
  geom_bar(stat = "identity",position = "fill")+ 
  scale_y_continuous(labels = percent)+ 
  scale_fill_manual(values = PT_color )
plot(g)+ theme_bw()
rm(g)
rm(Cluster_proportion)
rm(Cluster_ncells)
rm(Cluster_sex)
rm(Cluster_name)

FeatureScatter(PT_Healthy_Human_kidney_integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.5 ,
               cols = PT_color) + labs(title = "Adult human male kidney")  + ggmin::theme_powerpoint()

#show the color
show_col(hue_pal()(15))
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
ggplotColours(n=15)
#"#F8766D" "#E58700" "#C99800" "#A3A500" "#6BB100" "#00BA38" "#00BF7D" "#00C0AF" "#00BCD8"
#"#00B0F6" "#619CFF" "#B983FF" "#E76BF3" "#FD61D1" "#FF67A4"
rm(ggplotColours)

#Cell_population_definition
mother_population <-　PT_Healthy_Human_kidney_integrated
AIP_cell_population <- subset(mother_population, subset = SLC22A11 >= 0.5 | SLC22A12 >= 0.5) %>%
  subset(subset = SLC22A6 < 0.5 & SLC22A7 < 0.5 & SLC22A8 < 0.5)
DIP_cell_population <- subset(mother_population, subset = SLC22A11 >= 0.5 | SLC22A12 >= 0.5) %>%
  subset(subset = SLC22A6 >= 0.5 |SLC22A7 >= 0.5 |SLC22A8 >= 0.5)
BIP_cell_population <- subset(mother_population, subset = SLC22A6 >= 0.5 |SLC22A7 >= 0.5 |SLC22A8 >= 0.5) %>%
  subset(subset = SLC22A11 < 0.5 & SLC22A12 < 0.5)
DIN_cell_population <- subset(mother_population, subset = SLC22A11 < 0.5 & SLC22A12 < 0.5) %>%
  subset(subset = SLC22A6 < 0.5 & SLC22A7 < 0.5 & SLC22A8 < 0.5)
AIP_cell_population$model <- "AIP" #190
DIP_cell_population$model <- "DIP" #1174
BIP_cell_population$model <- "BIP" #1517
DIN_cell_population$model <- "DIN" #829

Cute_model <- merge(AIP_cell_population, c(DIP_cell_population,BIP_cell_population,DIN_cell_population))
Idents(Cute_model) <- factor(Idents(Cute_model), levels = unique(c("S1","S2","S3","DL")))


Cell_population_list <- SplitObject(Cute_model, split.by = "model")
Cell_population_list <- list(Cell_population_list[["BIP"]], Cell_population_list[["DIP"]],Cell_population_list[["AIP"]],Cell_population_list[["DIN"]])
names(Cell_population_list) <- c("BIP", "DIP", "AIP","DIN")
for (i in 1:length(Cell_population_list)) {
  print(dim(Cell_population_list[[i]]))
}
#1517, 1174, 190, 829
dim(PT_Healthy_Human_kidney_integrated)
#34145  3710

rm(BIP_cell_population)
rm(AIP_cell_population)
rm(DIP_cell_population)
rm(DIN_cell_population)
rm(mother_population)
