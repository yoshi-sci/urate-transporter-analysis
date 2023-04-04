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
#1. set the population
 #Figure 4C & Supplementary Figure 4C
population <- "DIP cell population"
mother_population <- Cell_population_list[["DIP"]]
focused_gene <- c("SLC22A11","SLC22A12","SLC2A9","SLC22A6","SLC22A7","SLC22A8","SLC17A1","SLC17A3")

 #
population <- "BIP cell population"
mother_population <- Cell_population_list[["BIP"]]
focused_gene <- c("SLC2A9","SLC22A6","SLC22A7","SLC22A8","SLC17A1","SLC17A3")

#2. set the group kind
Group <- c("Single", "Double", "Triple") 
Group <- c("Single", "Multiple") 

#3. group the cell
membrane <- "basolateral"
Single_basolateral <- subset(mother_population, 
                             subset = (SLC22A6 >= 0.5 & SLC22A7 < 0.5 & SLC22A8 < 0.5) |
                               (SLC22A6 < 0.5 & SLC22A7 >= 0.5 & SLC22A8 < 0.5) |
                               (SLC22A6 < 0.5 & SLC22A7 < 0.5 & SLC22A8 >= 0.5) 
)
Single_basolateral$basolateral <- "Single"
dim(Single_basolateral) #DIP 444, BIP 811

Double_basolateral <- subset(mother_population, 
                             subset = (SLC22A6 >= 0.5 & SLC22A7 >= 0.5 & SLC22A8 < 0.5) |
                               (SLC22A6 < 0.5 & SLC22A7 >= 0.5 & SLC22A8 >= 0.5) |
                               (SLC22A6 >= 0.5 & SLC22A7 < 0.5 & SLC22A8 >= 0.5)
)
Double_basolateral$basolateral <- "Double" # or "Multiple" 
dim(Double_basolateral) #DIP 468, BIP 510

Triple_basolateral <- subset(mother_population, 
                             subset = (SLC22A6 >= 0.5 & SLC22A7 >= 0.5 & SLC22A8 >= 0.5)
)
Triple_basolateral$basolateral <- "Triple" # or "Multiple" 
dim(Triple_basolateral) #DIP 262, BIP 196

mother_population_DGE <-  merge(Single_basolateral, c(Double_basolateral,Triple_basolateral))
dim(mother_population)
dim(mother_population_DGE)
Idents(mother_population_DGE) <- factor(Idents(mother_population_DGE), 
                                        levels = unique(Group))
mother_population_DGE$basolateral <- factor(mother_population_DGE$basolateral, levels=Group)
rm(Single_basolateral)
rm(Double_basolateral)
rm(Triple_basolateral)

DotPlot(mother_population_DGE, features = rev(focused_gene), group.by = "basolateral", 
        cols = c("blue", "red"),col.min = -0.8 ,col.max = 0.8, scale.min = 5 , scale.max = 80,
        dot.scale = 6) + 
  RotatedAxis() +ylab("BI transporter") +xlab(population)
rm(mother_population_DGE)
rm(mother_population)
rm(population)
rm(focused_gene)
rm(Group)
rm(membrane)

#Figure 4A
Cute_model$model <- factor(Cute_model$model, levels=rev(c("BIP", "DIP", "AIP", "DIN")))
mother_population <-　Cute_model
Gene_name <- "PDZK1"
DotPlot(Cute_model, features = rev(Gene_name), cols = c("blue", "red"),
        col.min = -0.8 ,col.max = 0.8, dot.scale = 6, dot.min = 0.05, group.by = "model",
        scale.min = 5 , scale.max = 80
) + RotatedAxis() 
rm(Gene_name)
rm(mother_population)

#Figure 4B
population <- "DIP cell population"
membrane <- "apical"
Group <- c("Single", "Double", "Triple", "Fourth")
focused_gene <- c("PDZK1")

mother_population <- Cell_population_list[["DIP"]]
population <- "DIP cell population"
dim(mother_population) #1174
Single_apical <- subset(mother_population, 
                        subset = (SLC22A11 >= 0.5 & SLC22A12 < 0.5 & SLC17A1 < 0.5 & SLC17A3 < 0.5) | 
                          (SLC22A11 < 0.5 & SLC22A12 >= 0.5 & SLC17A1 < 0.5 & SLC17A3 < 0.5)) #479
Single_apical$apical <- "Single"
dim(Single_apical) #173
Double_apical <- subset(mother_population, (SLC22A11 >= 0.5 & SLC22A12 < 0.5 & SLC17A1 >= 0.5 & SLC17A3 < 0.5) |
                          (SLC22A11 >= 0.5 & SLC22A12 < 0.5 & SLC17A1 < 0.5 & SLC17A3 >= 0.5) |
                          (SLC22A11 < 0.5 & SLC22A12 >= 0.5 & SLC17A1 >= 0.5 & SLC17A3 < 0.5) |
                          (SLC22A11 < 0.5 & SLC22A12 >= 0.5 & SLC17A1 < 0.5 & SLC17A3 >= 0.5) |
                          (SLC22A11 >= 0.5 & SLC22A12 >= 0.5 & SLC17A1 < 0.5 & SLC17A3 < 0.5)) #593
Double_apical$apical <- "Double"
dim(Double_apical)#393

Triple_apical <- subset(mother_population, (SLC22A11 >= 0.5 & SLC22A12 < 0.5 & SLC17A1 >= 0.5 & SLC17A3 >= 0.5) |
                          (SLC22A11 < 0.5 & SLC22A12 >= 0.5 & SLC17A1 >= 0.5 & SLC17A3 >= 0.5) |
                          (SLC22A11 >= 0.5 & SLC22A12 >= 0.5 & SLC17A1 >= 0.5 & SLC17A3 < 0.5)  |
                          (SLC22A11 >= 0.5 & SLC22A12 >= 0.5 & SLC17A1 < 0.5 & SLC17A3 >= 0.5)) #593
Triple_apical$apical <- "Triple"
dim(Triple_apical) #546

Fourth_apical <- subset(mother_population, 
                        subset = (SLC22A11 >= 0.5 & SLC22A12 >= 0.5 & SLC17A1 >= 0.5 & SLC17A3 >= 0.5)) #479
Fourth_apical$apical <- "Fourth"
dim(Fourth_apical) #62

mother_population_DGE <-  merge(Single_apical, c(Double_apical, Triple_apical, Fourth_apical))
dim(mother_population_DGE) #1174
Idents(mother_population_DGE) <- factor(Idents(mother_population_DGE), 
                                        levels = unique(Group))
mother_population_DGE$apical <- factor(mother_population_DGE$apical, levels= c(Group))
rm(Single_apical)
rm(Double_apical)
rm(Triple_apical)
rm(Fourth_apical)

DotPlot(mother_population_DGE, features = rev(focused_gene), group.by = membrane, 
        cols = c("blue", "red"),col.min = -0.8 ,col.max = 0.8, scale.min = 5 , scale.max = 80, 
        dot.scale = 6) + 
  RotatedAxis()  +ylab("Apical transporter") +xlab(population) 
rm(mother_population)
rm(mother_population_DGE)
rm(population)
rm(focused_gene)
rm(membrane)
rm(Group)

#
membrane <- "apical"
Group <- c("AI+No AE", "AI+SLC17A3", "AI+SLC17A1", "AI+Double AE")
focused_gene <- c("SLC22A11","SLC22A12","SLC2A9","SLC22A6","SLC22A7","SLC22A8")
mother_population <- Cell_population_list[["DIP"]]
population <- "DIP cell population"

No_apical <- subset(mother_population, subset = SLC17A1 < 0.5 & SLC17A3 < 0.5) #445
No_apical$apical <- "AI+No AE"
Single_apical_A <- subset(mother_population, subset = SLC17A1 >= 0.5 & SLC17A3 < 0.5) #445
Single_apical_A$apical <- "AI+SLC17A1"
Single_apical_B <- subset(mother_population, subset = SLC17A1 < 0.5 & SLC17A3 >= 0.5) #138
Single_apical_B$apical <- "AI+SLC17A3"
Single_apical <- subset(mother_population, 
                        subset = (SLC17A1 >= 0.5 & SLC17A3 < 0.5) | (SLC17A1 < 0.5 & SLC17A3 >= 0.5)) #593
Single_apical$apical <- "AI+Single"
Double_apical <- subset(mother_population, subset = SLC17A1 >= 0.5 & SLC17A3 >= 0.5) #557
Double_apical$apical <- "AI+Double AE"

mother_population_DGE <-  merge(No_apical, c(Single_apical_A, Single_apical_B, Double_apical))
dim(mother_population)
dim(mother_population_DGE)
Idents(mother_population_DGE) <- factor(Idents(mother_population_DGE), 
                                        levels = unique(Group))
mother_population_DGE$apical <- factor(mother_population_DGE$apical, levels= c(Group))
rm(No_apical)
rm(Single_apical)
rm(Single_apical_A)
rm(Single_apical_B)
rm(Double_apical)

DotPlot(mother_population_DGE, features = rev(focused_gene), group.by = membrane, 
        cols = c("blue", "red"),col.min = -0.8 ,col.max = 0.8, scale.min = 5 , scale.max = 80, 
        dot.scale = 6) + 
  RotatedAxis()  +ylab("Apical transporter") +xlab("Urate transporter") #+ labs(title = population)
rm(mother_population)
rm(mother_population_DGE)
rm(focused_gene)
rm(membrane)
rm(Group)
rm(population)

#Figure 4C & Supplementary Figure 4A
Group <- c("Neg", "Pos")
mother_population <- Cell_population_list[["DIP"]]
population <- "DIP cell population"

focused_gene <- c("SLC22A11","SLC22A12","SLC17A1","SLC17A3") #apical
focused_gene <- c("SLC2A9","SLC22A6","SLC22A7","SLC22A8") #basolateral

Scaffold_positive_model <- subset(mother_population, subset = PDZK1 >= 0.5) 
Scaffold_negative_model <- subset(mother_population, subset = PDZK1 < 0.5) 
dim(Scaffold_positive_model) #DIP 744
dim(Scaffold_negative_model) #DIP 430
Scaffold_positive_model$scaffold <- "Pos"
Scaffold_negative_model$scaffold <- "Neg"
mother_population_DGE <- merge(Scaffold_positive_model, c(Scaffold_negative_model))
Idents(mother_population_DGE) <- factor(Idents(mother_population_DGE), 
                                        levels = unique(Group))
mother_population_DGE$scaffold <- factor(mother_population_DGE$scaffold, levels= c(Group))

DotPlot(mother_population_DGE, features = rev(focused_gene), group.by = "scaffold", 
        cols = c("blue", "red"),col.min = -0.8 ,col.max = 0.8, scale.min = 5 , scale.max = 80, 
        dot.scale = 6) + 
  RotatedAxis()  +ylab("PDZK1") +xlab("Urate transporter") #+ labs(title = population)

#Figure 4D
mother_population <- Scaffold_positive_model
mother_population <- Scaffold_negative_model
#Count the single BI Pos: 187, Neg:257
dim(subset(mother_population, 
           subset = (SLC22A6 >= 0.5 & SLC22A7 < 0.5 & SLC22A8 < 0.5) |
             (SLC22A6 < 0.5 & SLC22A7 >= 0.5 & SLC22A8 < 0.5) |
             (SLC22A6 < 0.5 & SLC22A7 < 0.5 & SLC22A8 >= 0.5)))
#Count the double BI Pos: 346, Neg:122
dim(subset(mother_population,
           subset = (SLC22A6 >= 0.5 & SLC22A7 >= 0.5 & SLC22A8 < 0.5) |
             (SLC22A6 < 0.5 & SLC22A7 >= 0.5 & SLC22A8 >= 0.5) |
             (SLC22A6 >= 0.5 & SLC22A7 < 0.5 & SLC22A8 >= 0.5)))
#Count the Triple BI Pos: 211, Neg:51
dim(subset(mother_population, 
           subset = (SLC22A6 >= 0.5 & SLC22A7 >= 0.5 & SLC22A8 >= 0.5)))

Cluster_proportion <- data.frame(
  PDZK1   = c("Pos", "Neg","Pos", "Neg","Pos", "Neg"),
  BI_kinds = c("Single", "Single", "Double", "Double", "Triple", "Triple"),
  BI_transporter = c("Single", "Single", "Multiple", "Multiple", "Multiple", "Multiple"),
  ncell = c(187, 257, 346, 122, 211, 51)
)
Cluster_proportion$PDZK1 <- factor(Cluster_proportion$PDZK1 , levels=c("Neg","Pos"))
Cluster_proportion$BI_transporter <- factor(Cluster_proportion$BI_transporter, levels=c("Multiple","Single"))
Cluster_proportion$BI_kinds <- factor(Cluster_proportion$BI_kinds, levels=c("Triple","Double", "Single"))

#4 export the figure
g <- ggplot(Cluster_proportion, aes(x = PDZK1, y = ncell, fill = BI_transporter)) + 
  geom_bar(stat = "identity",position = "fill")+ 
  scale_y_continuous(labels = percent)+ 
  scale_fill_manual(values = c("#D3DEF1", "#000088"))
plot(g)+ theme_bw()
rm(Cluster_proportion)
rm(g)


rm(mother_population)
rm(mother_population_DGE)
rm(focused_gene)

########################################################################
mother_population <- Scaffold_positive_model
mother_population <- Scaffold_negative_model
#Supplementary Table 4 & Figure 4E & Figure 5C
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

#AE-postive population (Influx+AE) 
for (i in 1:length(mother_population)) {
  print(dim(subset(mother_population[[i]], subset = SLC2A9 < 0.5) %>%
              subset(subset = SLC17A1 >= 0.5 | SLC17A3 >= 0.5)%>%
              subset(idents = Region)
  ))}

#BE-postive population (Influx+BE)
for (i in 1:length(mother_population)) {
  print(dim(subset(mother_population[[i]], subset = SLC2A9 >= 0.5) %>%
              subset(subset = SLC17A1 < 0.5 & SLC17A3 < 0.5)%>%
              subset(idents = Region)
  ))}

#AE&BE-postive population (Influx+AE+BE)
for (i in 1:length(mother_population)) {
  print(dim(subset(mother_population[[i]], subset = SLC2A9 >= 0.5) %>%
              subset(subset = SLC17A1 >= 0.5 | SLC17A3 >= 0.5)%>%
              subset(idents = Region)
  ))}

rm(Region)
rm(mother_population)
rm(Group)
rm(population)
rm(Scaffold_positive_model)
rm(Scaffold_negative_model)

Cell_population_list <- SplitObject(Cute_model, split.by = "model")
Cell_population_list <- list(Cell_population_list[["BIP"]], Cell_population_list[["DIP"]],Cell_population_list[["AIP"]],Cell_population_list[["DIN"]])
names(Cell_population_list) <- c("BIP", "DIP", "AIP","DIN")

###########################################################################
Scaffoldgroup_markers <- FindMarkers(mother_population, ident.1 = "Pos",ident.2 = "Neg",group.by = "scaffold", only.pos = F,
                                     min.pct = 0.1, logfc.threshold = 0.25,min.diff.pct = -Inf, test.use = "wilcox")
Sigfini_Scaffoldgroup_markers <- subset(Scaffoldgroup_markers, subset =  Scaffoldgroup_markers[,"p_val_adj"] < 0.05)

write.table(Sigfini_Scaffoldgroup_markers, file = "/Users/yoshi/Downloads/Significant_maker2.tsv")
rm(Scaffoldgroup_markers)
rm(Sigfini_Scaffoldgroup_markers)
rm(mother_population)

