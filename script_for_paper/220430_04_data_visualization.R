setwd("")

#set the package
pacman::p_load(rgl, knitr, Matrix, scater, DropletUtils, scran, limma, 
               GSVA, GSEABase, genefilter, AnnotationDbi, org.Mm.eg.db, 
               org.Hs.eg.db, gplots, gdata, RColorBrewer, scales, reshape2,
               sctransform, Seurat, magrittr, dplyr, reticulate, ggplot2, 
               SeuratData, biomaRt, gProfileR, patchwork, cowplot, future, ggsci)

plan()
# change the current plan to access parallelization
plan("multisession", workers = 4)
plan()
#fall back to 
plan("sequential")

options(future.globals.maxSize = 3000 * 1024^2)
#(Default: 500 * 1024 ^ 2 = 500 MiB)

#############################################################################################################
#Urate transporters
Genename_focused_Chr <- list("SLC22A11","SLC22A12","SLC22A13","SLC2A9",
                             "SLC22A6","SLC22A7","SLC22A8",
                             "ABCG2","ABCC2","ABCC4","SLC17A1","SLC17A3")
Titlename_Chr <- list("SLC22A11/OAT4","SLC22A12/URAT1","SLC22A13/OAT10","SLC2A9/GLUT9",
                      "SLC22A6/OAT1","SLC22A7/OAT2","SLC22A8/OAT3",
                      "ABCG2/BCRP","ABCC2/MRP2","ABCC4/MRP4","SLC17A1/NPT1","SLC17A3/NPT4")
PT_color <- c("#FFCC66","#FFFF00", "#CCCC33","#008800")

for (i in 1:length(Genename_focused_Chr)) {
  print(VlnPlot(Cute_model, features = Genename_focused_Chr[[i]],
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

#Average the expression 
Functional_PT_Male_Healthy_Human_integrated <- merge(Cell_population_list[["BIP"]], c(Cell_population_list[["DIP"]],Cell_population_list[["AIP"]]))
Idents(Functional_PT_Male_Healthy_Human_integrated) <- factor(Idents(Functional_PT_Male_Healthy_Human_integrated), levels = unique(c("S1","S2","S3","DL")))

Urate_transporter_Chr <- c("SLC22A11","SLC22A12","SLC2A9","SLC22A6","SLC22A7","SLC22A8","SLC17A1","SLC17A3")

#1. select the population
Male_average <- AverageExpression(Functional_PT_Male_Healthy_Human_integrated , use.scale = F, use.counts = F,
                                  features = Urate_transporter_Chr)
Male_average <- AverageExpression(Cell_population_list[["AIP"]] , use.scale = F, use.counts = F,
                                  features = Urate_transporter_Chr)
Male_average <- AverageExpression(Cell_population_list[["DIP"]] , use.scale = F, use.counts = F,
                                  features = Urate_transporter_Chr)
Male_average <- AverageExpression(Cell_population_list[["BIP"]] , use.scale = F, use.counts = F, 
                                  features = Urate_transporter_Chr)
Male_average <- AverageExpression(Cell_population_list[["DIN"]] , use.scale = F, use.counts = F, 
                                  features = Urate_transporter_Chr)
#2. export the result
Male_average[["RNA"]]
#BIP, AIP, DIP cell population Table1
#S1         S2        S3        DL
#SLC22A11  2.9917552  3.9675884  4.520269 1.9352763
#SLC22A12  0.3630302  0.9212517  1.418967 0.6965603
#SLC2A9    4.6814368  4.5289187  7.352984 5.4708192
#SLC22A6  10.2461186 12.9317463 12.496015 1.5039022
#SLC22A7   1.2889020  1.9107321  4.737629 0.7940940
#SLC22A8   3.7525337  2.9620344  3.130045 1.2540011
#SLC17A1  10.5992864  8.5419028  7.658221 5.8940376
#SLC17A3   6.0233234  3.8219564  3.563649 2.1811390

rm(Functional_PT_Male_Healthy_Human_integrated)
rm(Urate_transporter_Chr)
rm(Male_average)

###############################################################################################################
#Transporter_location Supplementary Figure 2
Cute_model$model <- factor(Cute_model$model, levels=rev(c("BIP", "DIP", "AIP", "DIN")))
#1. Select the cluster
Cluster_name <- "S1"
Cluster_name <- "S2"
Cluster_name <- "S3"
Cluster_name <- "DL"

#2. Select the tranporter type
#AI transporter Supplementary Figure 2A,B 
mother_population <- merge(Cell_population_list[["DIP"]], Cell_population_list[["AIP"]]) %>% subset(idents = Cluster_name)
mother_population$model <- factor(mother_population$model, levels=rev(c("AIP", "DIP")))

AI_transporter <- c("SLC22A12")
Transporter_name <- AI_transporter

AI_transporter <- c("SLC22A11")
Transporter_name <- AI_transporter

#BI transporter Supplementary Figure 2C, D, E
mother_population <- merge(Cell_population_list[["DIP"]], Cell_population_list[["BIP"]]) %>% subset(idents = Cluster_name)
mother_population$model <- factor(mother_population$model, levels=rev(c("BIP", "DIP")))
BI_transporter <- c("SLC22A8")
Transporter_name <- BI_transporter

BI_transporter <- c("SLC22A6")
Transporter_name <- BI_transporter

BI_transporter <- c("SLC22A7")
Transporter_name <- BI_transporter

#Efflux transporter Figure 3C & Supplementary Figure 3B, C
mother_population <-　subset(Cute_model, idents = Cluster_name)
mother_population$model <- factor(mother_population$model, levels=rev(c("BIP", "DIP", "AIP", "DIN")))
Efflux_transporter <- c("SLC17A3")
Transporter_name <- Efflux_transporter

Efflux_transporter <- c("SLC17A1")
Transporter_name <- Efflux_transporter

Efflux_transporter <- c("SLC2A9")
Transporter_name <- Efflux_transporter

#3. export the figure
DotPlot(mother_population, features = rev(Transporter_name), cols = c("blue", "red"),
        col.min = -0.8 ,col.max = 0.8, dot.scale = 6, dot.min = 0.05, group.by = "model",
        scale.min = 5 , scale.max = 80
) + RotatedAxis() +ylab("Cell population") +xlab(Cluster_name) 
dim(mother_population)

##Influx&Efflux transporter
mother_population <- subset(Cell_population_list[["DIP"]],idents = "DL" )
Influx_transporter <- c("SLC22A11","SLC22A12","SLC22A6","SLC22A7","SLC22A8")
DotPlot(Cute_model, features = rev(Influx_transporter), cols = c("blue", "purple","red","grey"),
        split.by = "model"
) + RotatedAxis() + labs(title = ) +ylab("Cell population") 

DotPlot(Cute_model, features = rev(Efflux_transporter), cols = c("blue", "purple","red","grey"),
        col.min = -0.8 ,col.max = 0.8, dot.scale = 6, dot.min = 0.05, 
        scale.min = 5 , scale.max = 80, split.by = "model"
) + RotatedAxis() + labs(title = ) +ylab("Cell population")

rm(mother_population)
rm(Cluster_name)
rm(Transporter_name)
rm(AI_transporter)
rm(BI_transporter)
rm(Influx_transporter)
rm(Efflux_transporter)

###########################################################################################################3333
#Cell_population_proportion Figure 3B
#Whole kidney (1000, 1686, 767, 257)
Cluster_name <- c("S1","S2","S3","DL")
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cute_model, idents = Cluster_name[[i]])))
}
rm(Cluster_name)

#1. Select the cluster
Cluster_name <- "S1" #442, 473, 29, 56
Cluster_name <- "S2" #702, 374, 102, 508
Cluster_name <- "S3" #338, 313, 40, 76
Cluster_name <- "DL" #35, 14, 19, 189

#2. calculate the nubmers  of the cells
for (i in 1:length(Cell_population_list)) {
  print(dim(subset(Cell_population_list[[i]],  idents = Cluster_name)))
}
rm(Cluster_name)

#3 create the dataset
Cluster_proportion <- data.frame(
  model   = c("BIP_cell_population","DIP_cell_population","AIP_cell_population","DIN_cell_population","BIP_cell_population","DIP_cell_population","AIP_cell_population","DIN_cell_population",
              "BIP_cell_population","DIP_cell_population","AIP_cell_population","DIN_cell_population","BIP_cell_population","DIP_cell_population","AIP_cell_population","DIN_cell_population"),
  Cluster = c("S1","S1","S1","S1", "S2","S2","S2","S2",
              "S3","S3","S3","S3","DL","DL","DL","DL"),
  ncell = c(442, 473, 29, 56, 702, 374, 102, 508,
            338, 313, 40, 76, 35, 14, 19, 189)
)
Cluster_proportion$Cluster <- factor(Cluster_proportion$Cluster, levels=c("S1","S2","S3","DL"))
Cluster_proportion$model <- factor(Cluster_proportion$model, levels=c("BIP_cell_population","DIP_cell_population","AIP_cell_population","DIN_cell_population"))

#4 export the figure
g <- ggplot(Cluster_proportion, aes(x = Cluster, y = ncell, fill = model)) + 
  geom_bar(stat = "identity",position = "fill")+ 
  scale_y_continuous(labels = percent)+ 
  scale_fill_manual(values = c("#5D99FF", "#A16EFF", "#FF69A3", "#444444"))+labs(title = "Proportion of each cell population")
plot(g)+ theme_bw()
rm(Cluster_proportion)
rm(g)

CellPopulationCount <- matrix(c(442, 473, 29, 56, 702, 374, 102, 508,338, 313, 40, 76,35, 14, 19, 189), nrow = 4, ncol = 4) 
CellPopulationCount <- t(CellPopulationCount)
rownames(CellPopulationCount) <- c("S1",	"S2",	"S3",	"DL")
colnames(CellPopulationCount) <- c("BIP","DIP","AIP","DIN")
mosaicplot(CellPopulationCount, main = "Model distribution", color = TRUE, cex.axis = 1.3)
rm(CellPopulationCount)

#################################################################33333
#DIP->AIP->BIP->DIN
Cell_population_list <- SplitObject(Cute_model, split.by = "model")
Cell_population_list <- list(Cell_population_list[["DIP"]], Cell_population_list[["AIP"]],Cell_population_list[["BIP"]],Cell_population_list[["DIN"]])
names(Cell_population_list) <- c("DIP", "AIP", "BIP","DIN")

#AI_Overlapping analyses_DIP cell population Figure 2C
#1. Calculate the number of cells
Cluster_name <- c("S1","S2","S3","DL")
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["DIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A11 >= 0.5 & SLC22A12 < 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["DIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A11 < 0.5 & SLC22A12 >= 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["DIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A11 >= 0.5 & SLC22A12 >= 0.5)))
}
dim(subset(Cell_population_list[["DIP"]],  idents = "S1")%>%
      subset(subset = SLC22A11 >= 0.5 & SLC22A12 < 0.5))

S1_ncells <- c(381, 42, 50)
S2_ncells <- c(279, 67, 28)
S3_ncells <- c(205, 57, 51)
LOH_ncells <- c(9,4,1)

#2. Select the cluster
Cluster_name <- "S1"
data_ncells <- S1_ncells

Cluster_name <- "S2"
data_ncells <- S2_ncells

Cluster_name <- "S3"
data_ncells <- S3_ncells

Cluster_name <- "DL"
data_ncells <- LOH_ncells

#3. export the figure
selected_cell_population <- c("DIP_cell_population","DIP_cell_population","DIP_cell_population")
DIP_AI_detail_proportion <- data.frame(model   = selected_cell_population, ncell = data_ncells,
  AI_transporter = c("SLC22A11/OAT4"	,"SLC22A12/URAT1"	,"SLC22A11/OAT4 & SLC22A12/URAT1"))
DIP_AI_detail_proportion$AI_transporter <- factor(DIP_AI_detail_proportion$AI_transporter, levels=c("SLC22A11/OAT4"	,"SLC22A12/URAT1"	,"SLC22A11/OAT4 & SLC22A12/URAT1"))

g <- ggplot(DIP_AI_detail_proportion, aes(x = model, y = ncell, fill = AI_transporter)) + 
  geom_bar(stat = "identity",position = "fill")+ 
  scale_y_continuous(labels = percent)+ 
  scale_fill_manual(values = c("#FFDDFF", "#FFABCE", "#DD0000"))+
  labs(title = Cluster_name)
plot(g)

rm(Cluster_name)
rm(g)
rm(DIP_AI_detail_proportion)
rm(selected_cell_population)
rm(data_ncells)

rm(S1_ncells)
rm(S2_ncells)
rm(S3_ncells)
rm(LOH_ncells)

#AI_Overlapping analyses_AIP cell population Supplementary Figure 2F
#1. Calculate the number of cells
Cluster_name <- c("S1","S2","S3","DL")
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["AIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A11 >= 0.5 & SLC22A12 < 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["AIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A11 < 0.5 & SLC22A12 >= 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["AIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A11 >= 0.5 & SLC22A12 >= 0.5)))
}
S1_ncells <- c(25, 3, 1)
S2_ncells <- c(85, 15, 2)
S3_ncells <- c(31, 7, 2)
LOH_ncells <- c(15, 4, 0)

#2. Select the cluster
Cluster_name <- "S1"
data_ncells <- S1_ncells

Cluster_name <- "S2"
data_ncells <- S2_ncells

Cluster_name <- "S3"
data_ncells <- S3_ncells

Cluster_name <- "DL"
data_ncells <- LOH_ncells

#3. export the figure
selected_cell_population <- c("AIP_cell_population","AIP_cell_population","AIP_cell_population")
AIP_AI_detail_proportion <- data.frame(model   = selected_cell_population, ncell = data_ncells,
                                      AI_transporter = c("SLC22A11/OAT4"	,"SLC22A12/URAT1"	,"SLC22A11/OAT4 & SLC22A12/URAT1"))
AIP_AI_detail_proportion$AI_transporter <- factor(AIP_AI_detail_proportion$AI_transporter, levels=c("SLC22A11/OAT4"	,"SLC22A12/URAT1"	,"SLC22A11/OAT4 & SLC22A12/URAT1"))

g <- ggplot(AIP_AI_detail_proportion, aes(x = model, y = ncell, fill = AI_transporter)) + 
  geom_bar(stat = "identity",position = "fill")+ 
  scale_y_continuous(labels = percent)+ 
  scale_fill_manual(values = c("#FFDDFF", "#FFABCE", "#DD0000"))+
  labs(title = Cluster_name)
plot(g)

rm(g)
rm(AIP_AI_detail_proportion)
rm(Cell_population_list)
rm(selected_cell_population)
rm(data_ncells)
rm(Cluster_name)

rm(S1_ncells)
rm(S2_ncells)
rm(S3_ncells)
rm(LOH_ncells)
##################################################
#BIP->DIP>AIP->DIN
Cell_population_list <- SplitObject(Cute_model, split.by = "model")
Cell_population_list <- list(Cell_population_list[["BIP"]], Cell_population_list[["DIP"]],Cell_population_list[["AIP"]],Cell_population_list[["DIN"]])
names(Cell_population_list) <- c("BIP", "DIP", "AIP","DIN")

#BI_Overlapping analyses_DIP cell population Figure 2D
#1. Calculate the number of cells
Cluster_name <- c("S1","S2","S3","DL")
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["DIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 >= 0.5 & SLC22A7 < 0.5 & SLC22A8 < 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["DIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 < 0.5 & SLC22A7 >= 0.5 & SLC22A8 < 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["DIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 < 0.5 & SLC22A7 < 0.5 & SLC22A8 >= 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["DIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 >= 0.5 & SLC22A7 >= 0.5 & SLC22A8 < 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["DIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 >= 0.5 & SLC22A7 < 0.5 & SLC22A8 >= 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["DIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 < 0.5 & SLC22A7 >= 0.5 & SLC22A8 >= 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["DIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 >= 0.5 & SLC22A7 >= 0.5 & SLC22A8 >= 0.5)))
}
S1_ncells <- c(99,3,15,38,182,6,130)
S2_ncells <- c(168,14,33,40,77,7,35)
S3_ncells <- c(68,29,6,56,49,9,96)
LOH_ncells <- c(3,3,3,2,2,0,1)

#2. Select the cluster
Cluster_name <- "S1"
data_ncells <- S1_ncells

Cluster_name <- "S2"
data_ncells <- S2_ncells

Cluster_name <- "S3"
data_ncells <- S3_ncells

Cluster_name <- "DL"
data_ncells <- LOH_ncells

#select the model
selected_cell_population <- c("DIP_cell_population","DIP_cell_population","DIP_cell_population","DIP_cell_population","DIP_cell_population","DIP_cell_population","DIP_cell_population")
DIP_BI_detail_proportion <- data.frame(
  model   = selected_cell_population,
  BI_transporter = c("SLC22A6/OAT1",	"SLC22A7/OAT2",	"SLC22A8/OAT3", 
          "SLC22A6/OAT1 & SLC22A7/OAT2","SLC22A6/OAT1 & SLC22A8/OAT3",	"SLC22A7/OAT2 & SLC22A8/OAT3"	,
          "SLC22A6/OAT1 & SLC22A7/OAT2 & SLC22A8/OAT3"),
  ncell = data_ncells
)
DIP_BI_detail_proportion$BI_transporter <- factor(DIP_BI_detail_proportion$BI_transporter, levels=c("SLC22A6/OAT1",	"SLC22A7/OAT2",	"SLC22A8/OAT3", 
                                                                        "SLC22A6/OAT1 & SLC22A7/OAT2","SLC22A6/OAT1 & SLC22A8/OAT3",	"SLC22A7/OAT2 & SLC22A8/OAT3"	,
                                                                        "SLC22A6/OAT1 & SLC22A7/OAT2 & SLC22A8/OAT3"))

g <- ggplot(DIP_BI_detail_proportion, aes(x = model, y = ncell, fill = BI_transporter)) + 
  geom_bar(stat = "identity",position = "fill")+ 
  scale_y_continuous(labels = percent)+ 
  scale_fill_manual(values = c("#EEFFFF", "#CEF9DC", "#BAD3FF","#8EF1FF",  "#4689FF","#30F9B2","#000088"))+
  labs(title = Cluster_name)
plot(g)

remove(g)
remove(DIP_BI_detail_proportion)
rm(data_ncells)
rm(Cluster_name)
rm(selected_cell_population)

rm(S1_ncells)
rm(S2_ncells)
rm(S3_ncells)
rm(LOH_ncells)

#BI_Overlapping analyses_BIP cell population Supplementary Figure 2G
#1. Calculate the number of cells
Cluster_name <- c("S1","S2","S3","DL")
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["BIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 >= 0.5 & SLC22A7 < 0.5 & SLC22A8 < 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["BIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 < 0.5 & SLC22A7 >= 0.5 & SLC22A8 < 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["BIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 < 0.5 & SLC22A7 < 0.5 & SLC22A8 >= 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["BIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 >= 0.5 & SLC22A7 >= 0.5 & SLC22A8 < 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["BIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 >= 0.5 & SLC22A7 < 0.5 & SLC22A8 >= 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["BIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 < 0.5 & SLC22A7 >= 0.5 & SLC22A8 >= 0.5)))
}
for (i in 1:length(Cluster_name)) {
  print(dim(subset(Cell_population_list[["BIP"]],  idents = Cluster_name[i]) %>%
              subset(subset = SLC22A6 >= 0.5 & SLC22A7 >= 0.5 & SLC22A8 >= 0.5)))
}

S1_ncells <- c(127,4,25,29,159,7,91)
S2_ncells <- c(364,36,71,73,107,9,42)
S3_ncells <- c(102,38,15,80,32,9,62)
LOH_ncells <- c(12,4,13,2,1,2,1)

#2. Select the cluster
Cluster_name <- "S1"
data_ncells <- S1_ncells

Cluster_name <- "S2"
data_ncells <- S2_ncells

Cluster_name <- "S3"
data_ncells <- S3_ncells

Cluster_name <- "DL"
data_ncells <- LOH_ncells

#select the model
selected_cell_population <- c("BIP_cell_population","BIP_cell_population","BIP_cell_population","BIP_cell_population","BIP_cell_population","BIP_cell_population","BIP_cell_population")
BIP_BI_detail_proportion <- data.frame(
  model   = selected_cell_population,
  BI_transporter = c("SLC22A6/OAT1",	"SLC22A7/OAT2",	"SLC22A8/OAT3", 
          "SLC22A6/OAT1 & SLC22A7/OAT2","SLC22A6/OAT1 & SLC22A8/OAT3",	"SLC22A7/OAT2 & SLC22A8/OAT3"	,
          "SLC22A6/OAT1 & SLC22A7/OAT2 & SLC22A8/OAT3"),
  ncell = data_ncells
)
BIP_BI_detail_proportion$BI_transporter <- factor(BIP_BI_detail_proportion$BI_transporter, levels=c("SLC22A6/OAT1",	"SLC22A7/OAT2",	"SLC22A8/OAT3", 
                                                                        "SLC22A6/OAT1 & SLC22A7/OAT2","SLC22A6/OAT1 & SLC22A8/OAT3",	"SLC22A7/OAT2 & SLC22A8/OAT3"	,
                                                                        "SLC22A6/OAT1 & SLC22A7/OAT2 & SLC22A8/OAT3"))

g <- ggplot(BIP_BI_detail_proportion, aes(x = model, y = ncell, fill = BI_transporter)) + 
  geom_bar(stat = "identity",position = "fill")+ 
  scale_y_continuous(labels = percent)+ 
  scale_fill_manual(values = c("#EEFFFF", "#CEF9DC", "#BAD3FF","#8EF1FF",  "#4689FF","#30F9B2","#000088"))+
  labs(title = Cluster_name)
plot(g)

remove(g)
remove(BIP_BI_detail_proportion)
rm(data_ncells)
rm(Cluster_name)
rm(selected_cell_population)

rm(S1_ncells)
rm(S2_ncells)
rm(S3_ncells)
rm(LOH_ncells)
#####################################################################################################################
#1.select the region
Region <- "S1"
Region <- "S2"
Region <- "S3"
Region <- "DL"

#2. select the turn of the cell poplation
#AIP->DIP>BIP->DIN
Cell_population_list <- SplitObject(Cute_model, split.by = "model")
Cell_population_list <- list(Cell_population_list[["AIP"]], Cell_population_list[["DIP"]],Cell_population_list[["BIP"]],Cell_population_list[["DIN"]])
names(Cell_population_list) <- c("AIP", "DIP", "BIP","DIN")

#DIP->AIP->BIP->DIN
Cell_population_list <- SplitObject(Cute_model, split.by = "model")
Cell_population_list <- list(Cell_population_list[["DIP"]], Cell_population_list[["AIP"]],Cell_population_list[["BIP"]],Cell_population_list[["DIN"]])
names(Cell_population_list) <- c("DIP", "AIP", "BIP","DIN")

#BIP->DIP>AIP->DIN
Cell_population_list <- SplitObject(Cute_model, split.by = "model")
Cell_population_list <- list(Cell_population_list[["BIP"]], Cell_population_list[["DIP"]],Cell_population_list[["AIP"]],Cell_population_list[["DIN"]])
names(Cell_population_list) <- c("BIP", "DIP", "AIP","DIN")

#3.calculate the numbers of the cells
for (i in 1:length(Cell_population_list)) {
  print(dim(subset(Cell_population_list[[i]]) %>%
              subset(idents = Region)
  ))}
#BIP DIP AIP DIN
#S1 442, 473, 29, 56
#S2 702, 374, 102, 508
#S3 338, 313, 40, 76
#LOH 35, 14, 19, 189

#4. select the population
Single_basolateral <- subset(Cell_population_list[["DIP"]], 
                             subset = (SLC22A6 >= 0.5 & SLC22A7 < 0.5 & SLC22A8 < 0.5) |
                               (SLC22A6 < 0.5 & SLC22A7 >= 0.5 & SLC22A8 < 0.5) |
                               (SLC22A6 < 0.5 & SLC22A7 < 0.5 & SLC22A8 >= 0.5) 
)
Single_basolateral$basolateral <- "Single"
dim(Single_basolateral) #DIP 444, BIP 811
Double_basolateral <- subset(Cell_population_list[["DIP"]], 
                             subset = (SLC22A6 >= 0.5 & SLC22A7 >= 0.5 & SLC22A8 < 0.5) |
                               (SLC22A6 < 0.5 & SLC22A7 >= 0.5 & SLC22A8 >= 0.5) |
                               (SLC22A6 >= 0.5 & SLC22A7 < 0.5 & SLC22A8 >= 0.5)
)
Double_basolateral$basolateral <- "Double"
dim(Double_basolateral) #DIP 468, BIP 510
Triple_basolateral <- subset(Cell_population_list[["DIP"]], 
                             subset = (SLC22A6 >= 0.5 & SLC22A7 >= 0.5 & SLC22A8 >= 0.5)
)
Triple_basolateral$basolateral <- "Triple"
dim(Triple_basolateral) #DIP 262, BIP 196
Multiple_basolateral <-  merge(Double_basolateral, c(Triple_basolateral))
Multiple_basolateral$basolateral <- "Multiple"
dim(Multiple_basolateral) #730
rm(Double_basolateral)
rm(Triple_basolateral)

mother_population <- list(Multiple_basolateral, Single_basolateral)
rm(Single_basolateral)
rm(Multiple_basolateral)


#############################################################
##Supplementary Table 4 & Figure 5C
mother_population <- Cell_population_list
#No efflux population (Influx) Figure 
for (i in 1:length(mother_population)) {
  print(dim(subset(mother_population[[i]], subset = SLC2A9 < 0.5) %>%
              subset(subset = SLC17A1 < 0.5 & SLC17A3 < 0.5)%>%
              subset(idents = Region)
  ))}
print(dim(subset(mother_population[["AIP"]], subset = SLC2A9 < 0.5) %>%
            subset(subset = SLC17A1 < 0.5 & SLC17A3 < 0.5)%>%
            subset(idents = Region)
))
print(dim(subset(mother_population[["DIN"]], subset = SLC2A9 < 0.5) %>%
            subset(subset = SLC17A1 < 0.5 & SLC17A3 < 0.5)%>%
            subset(idents = Region)
))
#BIP DIP AIP DIN/ DIP PDZK1+ -
#S1 26, 25, 8, 16 / 12, 13
#S2 164, 67, 40, 225 / 19, 48
#S3 67, 37, 22, 32 / 16, 21
#LOH 4, 0, 4, 71 / 0, 0

#AE-postive population (Influx+AE) Figure 3D
for (i in 1:length(mother_population)) {
  print(dim(subset(mother_population[[i]], subset = SLC2A9 < 0.5) %>%
              subset(subset = SLC17A1 >= 0.5 | SLC17A3 >= 0.5)%>%
              subset(idents = Region)
  ))}
#S1 130, 126, 11, 31 / 78, 48
#S2 226, 149, 36, 166 / 60, 89
#S3 69, 63, 10, 24 / 30, 33
#LOH 4, 3, 8, 35 / 1, 2

#BE-postive population (Influx+BE)
for (i in 1:length(mother_population)) {
  print(dim(subset(mother_population[[i]], subset = SLC2A9 >= 0.5) %>%
              subset(subset = SLC17A1 < 0.5 & SLC17A3 < 0.5)%>%
              subset(idents = Region)
  ))}
#S1 8, 11, 5, 2/ 6, 5
#S2 48, 29, 5, 52 / 12, 17
#S3 35, 27, 3, 7 / 15 ,12
#LOH 7, 2, 1, 27 / 1, 1

#AE&BE-postive population (Influx+AE+BE)
for (i in 1:length(mother_population)) {
  print(dim(subset(mother_population[[i]], subset = SLC2A9 >= 0.5) %>%
              subset(subset = SLC17A1 >= 0.5 | SLC17A3 >= 0.5)%>%
              subset(idents = Region)
  ))}
#S1 278, 311, 5, 7 / 260, 51
#S2 264, 129, 21, 65 / 68, 61
#S3 167, 186, 5, 13 / 149, 37
#LOH 20, 9, 6, 56 / 3, 6

rm(Region)
rm(mother_population)

Cell_population_list <- SplitObject(Cute_model, split.by = "model")
Cell_population_list <- list(Cell_population_list[["BIP"]], Cell_population_list[["DIP"]],Cell_population_list[["AIP"]],Cell_population_list[["DIN"]])
names(Cell_population_list) <- c("BIP", "DIP", "AIP","DIN")

#####################################################################################################################
for (i in 1:length(Cell_population_list)) {
  print(FeatureScatter(Cell_population_list[[i]] , feature1 = "SLC22A11", feature2 = "SLC22A12", slot = "data", pt.size = 1.0)) + NoLegend() + 
    scale_y_continuous(breaks=seq(0,5,1/2))+   scale_x_continuous(breaks=seq(0,5,1/2))+  
    coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 4.5))
}

for (i in 1:length(Cell_population_list)) {
  ploti <- VlnPlot(Cell_population_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(ploti)
}