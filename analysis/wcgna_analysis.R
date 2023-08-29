osama_main <- readRDS("bulk_rna_mat.rds")

library(readr)
library(tidyverse)
expression_matrix <- osama_main@assays@data$counts
dim(expression_matrix) #20276 genes #76 samples
phenotype_data <-  as.data.frame(osama_main@colData)
View(phenotype_data)
all(colnames(expression_matrix) %in% phenotype_data$Sample) # FALSE
sum(colnames(expression_matrix) %in% phenotype_data$Sample) # 75
missing_sample <-  colnames(expression_matrix) %in% phenotype_data$Sample
View(missing_sample) # 25 is false
colnames(expression_matrix)[25] # H5608_Cing

missing_sample2 <- phenotype_data$Sample %in%  colnames(expression_matrix) 
View(missing_sample2)    
phenotype_data$Sample[25]   # Saved as "H5608" thus the difference
phenotype_data$Sample[25] = colnames(expression_matrix)[25]
    
all(colnames(expression_matrix) %in% phenotype_data$Sample) ## TRUE
all(colnames(expression_matrix) == phenotype_data$Sample)  # TRUE
rownames(phenotype_data) = phenotype_data$Sample
phenotype_data <- phenotype_data %>% mutate(Status = case_when(
  Juvenile =="Y" ~ "Juvenile_HD",
  Condition == "HD" ~ "HD",
  TRUE ~ as.character(Condition)  
))
library(tidyr)
library(stringr)
phenotype_data <- phenotype_data %>% mutate(Status_Location = str_c(Status,Location, sep = "_"))
View(phenotype_data$Status)


#Binarize Categorical Variables
library(WGCNA)
options(stringsAsFactors = F)
# Binarize Location
var_location <- binarizeCategoricalVariable(phenotype_data$Location, includePairwise = FALSE, includeLevelVsAll = TRUE)
View(var_location)

# Binarize gender 
var_gender <- binarizeCategoricalVariable(phenotype_data$Gender,includePairwise = T,includeLevelVsAll = F)
View(binarize_gender)

#Binarize Condition
var_condition <- binarizeCategoricalVariable(phenotype_data$Condition, includePairwise = T, includeLevelVsAll = F)
View(var_condition)

#Binarize Status
var_status <- binarizeCategoricalVariable(phenotype_data$Status, includePairwise = F, includeLevelVsAll = T)
View(var_status)


#B

colnames(var_location) <- c("Accumbens", "Caudate", "Cingulate")
colnames(var_gender) <- "Male" # Gender, 1 is male, 0 is female
colnames(var_status) <- c("Control", "HD", "Juvenile HD")

var_combined<- cbind(var_status, var_location, var_gender, phenotype_data$Age)
rownames(var_combined) <- phenotype_data$Sample
colnames(var_combined)[8] <- "Age"
View(var_combined)

BiocManager::install("edgeR")
library(edgeR)

keeps <- filterByExpr(expression_matrix, design = var_combined)
filtered_expression_matrix<- expression_matrix[keeps,]
dim(filtered_expression_matrix) # 15413 genes remaining

normalized_expression_matix <- limma::voom(filtered_expression_matrix, design = var_combined, plot = T, normalize.method = "quantile")$E   # Normalize counts with voom

# Visualize Batch Effects
library(factoextra)
library(FactoMineR)
??PCA
pca.plot.norm <- PCA(t(normalized_expression_matix), scale= TRUE, graph = F)
pdf(file = "pca_plot_precorrection.pdf", height = 9, width = 6)
fviz_pca_ind(pca.plot.norm, geom.ind = "point", pointshape = 21, pointsize = 2, fill.ind = as.character(phenotype_data$batch), col.ind = "black", palette = "jco", addEllipses = TRUE, label = "var", col.var = "black",repel = TRUE, legend.title = "Batch") + ggtitle("2D PCA-plot to Visualize Batch Effects") +theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "pca_plot_precorrection_by_staus.pdf", height = 9, width = 6)
fviz_pca_ind(pca.plot.norm, geom.ind = "point", pointshape = 21, pointsize = 2, fill.ind = as.character(phenotype_data$Status), col.ind = "black", palette = "jco", addEllipses = TRUE, label = "var", col.var = "black",repel = TRUE, legend.title = "Status") + ggtitle("2D PCA-plot to Visualize Batch Effects") +theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "pca_plot_precorrection_by_location.pdf", height = 9, width = 6)
fviz_pca_ind(pca.plot.norm, geom.ind = "point", pointshape = 21, pointsize = 2, fill.ind = as.character(phenotype_data$Location), col.ind = "black", palette = "jco", addEllipses = TRUE, label = "var", col.var = "black",repel = TRUE, legend.title = "Location") + ggtitle("2D PCA-plot to Visualize Batch Effects") +theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Perform Batch correction 
batch_corrected_matrix <- limma::removeBatchEffect(normalized_expression_matix, batch = phenotype_data$batch, design = var_combined)

#Post batch correction PCA

pca.plot.corrected <- PCA(t(batch_corrected_matrix), scale= TRUE, graph = F)
pdf(file = "pca_plot_corrected.pdf", height = 9, width = 6)
fviz_pca_ind(pca.plot.corrected, geom.ind = "point", pointshape = 21, pointsize = 2, fill.ind = as.character(phenotype_data$batch), col.ind = "black", palette = "jco", addEllipses = TRUE, label = "var", col.var = "black",repel = TRUE, legend.title = "Batch") + ggtitle("2D PCA Plot Corrected for Batch") +theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "pca_plot_corrected_status.pdf", height = 9, width = 6)
fviz_pca_ind(pca.plot.corrected, geom.ind = "point", pointshape = 21, pointsize = 2, fill.ind = as.character(phenotype_data$Status), col.ind = "black", palette = "jco", addEllipses = TRUE, label = "var", col.var = "black",repel = TRUE, legend.title = "Status") + ggtitle("2D PCA_Plot Corrected for Batch") +theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "pca_plot_corrected_location.pdf", height = 9, width = 6)
fviz_pca_ind(pca.plot.corrected, geom.ind = "point", pointshape = 21, pointsize = 2, fill.ind = as.character(phenotype_data$Location), col.ind = "black", palette = "jco", addEllipses = TRUE, label = "var", col.var = "black",repel = TRUE, legend.title = "Location") + ggtitle("2D PCA_Plot Corrected for Batch") +theme(plot.title = element_text(hjust = 0.5))
dev.off()




#1 Visualization of global changes

#PCA and Correlation Matrix 
library(factoextra)
library(ggpubr)

summarised_phenodata <- phenotype_data %>% group_by(Status) %>% summarise(count=n(), Accumbens = sum(Location=="Accumbens"), Cingulate = sum(Location=="Cingulate"), Caudate = sum(Location=="Caudate"))
View(summarised_phenodata)


pdf(file = "summary_of_samples.pdf", height = 3,width = 6)
ggtexttable(summarised_phenodata, row=NULL, theme = ttheme("mBlue"))
dev.off()

# Sample Correlation Heatmaps
pdf(file="correlation_heatmap_status_n_location.pdf", height = 12, width = 8)
heatmap(cor(batch_corrected_matrix),cexCol=0.7,cexRow=0.7, labRow =phenotype_data$Status, labCol = phenotype_data$Location)
dev.off()


# Analysis of Genes Resulting in Greatest Variation
# Scree plot for corrected data
pdf(file = "scree_plot_of_corrected_pca.pdf", height = 6, width = 6) #will you 7 principal components for further analysis, based on the scree plot
fviz_eig (pca.plot.corrected, addlabels = TRUE , ylim = c ( 0 , 30 )) 
dev.off()


# Contributions of variables to PC1-PC5
pdf(file = "top_20_contributing_variables_to_PC1-PC5.pdf", height = 9, width = 9)
fviz_contrib (pca.plot.corrected, choice = "var" , axes = 1:5, top = 30 )
dev.off()

#Biplot for PCA for variables with top 10 contribution
pdf(file = "biplot_pca_corrected.pdf", height = 12, width = 12)
fviz_pca_biplot (pca.plot.corrected, select.var = list(contrib =10),
                 # Individuals 
                 geom.ind = "point" ,
                 fill.ind = phenotype_data$Condition, col.ind = "black" ,
                 pointshape = 21 , pointsize = 2 ,
                 palette = "jama" ,
                 addEllipses = TRUE ,
                 # Variables 
                  col.var = "contrib" ,
                 gradient.cols = "RdYlBu" ,
                 repel = T,
                 legend.title = list ( fill = "Condition" , color = "Contrib" ,
                                       repel = T)
) 
dev.off()

pdf(file = "biplot_pca_corrected_location.pdf", height = 12, width = 12)
fviz_pca_biplot (pca.plot.corrected, select.var = list(contrib =10),
                 # Individuals 
                 geom.ind = "point" ,
                 fill.ind = phenotype_data$Location, col.ind = "black" ,
                 pointshape = 21 , pointsize = 2 ,
                 palette = "jco" ,
                 addEllipses = TRUE ,
                 # Variables 
                 col.var = "contrib" ,
                 gradient.cols = "RdYlBu" ,
                 repel = T,
                 legend.title = list ( fill = "Location" , color = "Contrib" ,
                                       repel = T)
) 
dev.off()

pdf(file = "biplot_pca_corrected_status.pdf", height = 12, width = 12)
fviz_pca_biplot (pca.plot.corrected, select.var = list(contrib =10),
                 # Individuals 
                 geom.ind = "point" ,
                 fill.ind = phenotype_data$Status, col.ind = "black" ,
                 pointshape = 21 , pointsize = 2 ,
                 palette = "jama" ,
                 addEllipses = TRUE ,
                 # Variables 
                 col.var = "contrib" ,
                 gradient.cols = "RdYlBu" ,
                 repel = T,
                 legend.title = list ( fill = "Status" , color = "Contrib" ,
                                       repel = T)
) 
dev.off()

# Hierarchical Clusteering of 5 Principal Components

try_hpcc <- PCA(t(batch_corrected_matrix), ncp=5, graph=F)
res.hcpc <- HCPC(try_hpcc, graph=F)
pdf(file = "hierarchical_clustering_pca.pdf", height = 20, width=12)
fviz_dend (res.hcpc, cex = 0.6 , # Label size 
          palette = "jco" ,
          rect = TRUE , rect_fill = TRUE , rect_border = "jco" , labels_track_height = 0.8) 
dev.off()

View(head(res.hcpc$data.clust, 50))

phenotype_data$HCPC_Cluster <- res.hcpc$data.clust$clust
View(phenotype_data)

summarised_clusters <- phenotype_data %>% group_by(Status) %>% summarise(Count=n(), Cluster_1=sum(HCPC_Cluster==1), Cluster_2=sum(HCPC_Cluster==2), Cluster_3=sum(HCPC_Cluster==3))
View(summarised_clusters)
colnames(summarised_clusters) <- c("Status","Count", "Cluster 1", "Cluster 2", "Cluster 3")

pdf(file = "summary_of_hcpc_clusters.pdf", height = 4,width = 6)
ggtexttable(summarised_clusters, row=NULL, theme = ttheme("mBlue"))
dev.off()


pdf(file = "summary_of_samples.pdf", height = 3,width = 6)
ggtexttable(summarised_phenodata, row=NULL, theme = ttheme("mBlue"))
dev.off()

# Test of association between clusters and status

summarised_clusters <- data.frame(summarised_clusters, row.names = summarised_clusters$Status)
View(summarised_clusters)
summarised_clusters <- summarised_clusters[,-c(1,2)]

dt <- as.table(as.matrix(summarised_clusters))
library(gplots)
pdf("balloon_plot.pdf", height = 5, width = 5)
balloonplot(t(dt), xlab = "", ylab = "", label= False, show.margins = F)
dev.off()

chisq_clusters <- chisq.test(summarised_clusters)
chisq_clusters # pvalue of chi square test is 1.701 e^-07
View(res.hcpc$desc.var$quanti$`1`)
top_10_cluster1 <- rownames((res.hcpc$desc.var$quanti$`1`))[1:10]
top_10_cluster1
top_10_cluster2 <- rownames((res.hcpc$desc.var$quanti$`2`))[1:10]
top_10_cluster3 <- rownames((res.hcpc$desc.var$quanti$`3`))[1:10]
heatmap_df_genes <- c(top_10_cluster1, top_10_cluster2, top_10_cluster3)
View(heatmap_df_genes)
heat_df_genes <- data.frame(heatmap_df_genes, Cluster=rep(c(1,2,3),each=10), row.names = heatmap_df_genes)
write.csv(heat_df_genes,file = "hcpc_cluster_genes.csv")
heatmap_df <- batch_corrected_matrix[rownames(heat_df_genes),]
View(heatmap_df)
annotation_df <- data.frame(Status=phenotype_data$Status, row.names = rownames(phenotype_data))
annotation_df_row <- data.frame(Cluster=as.character(heat_df_genes$Cluster), row.names = rownames(heat_df_genes))
View(annotation_df)
library(RColorBrewer)
heat_colors <- brewer.pal(6, "YlOrRd")
library(pheatmap)
### Run pheatmap
pdf(file = "heatmap_of_cluster_flagship_genes.pdf", height = 8, width = 6)
pheatmap(heatmap_df, color=heat_colors, cluster_rows = F, show_rownames=T,show_colnames = F,cluster_cols = T, annotation_row = annotation_df_row,
         annotation_col = annotation_df, annotation_border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)
dev.off()

save(batch_corrected_matrix, pca.plot.corrected, var_combined, pca.plot.norm, phenotype_data, res.hcpc, file = "New_Preprocessing_and_Global_Visualization.RData")




# Weighted Gene Correlation Network Analyss

library(WGCNA)
library(matrixStats)
library(DGCA)
options(stringsAsFactors = F)
lnames = load("New_Preprocessing_and_Global_Visualization.RData")



# Filtergenes base on level of expression and variance
nrow(batch_corrected_matrix) #15413 genes
library(matrixStats, quietly = TRUE)
?filterGenes
datExpr0= filterGenes(batch_corrected_matrix, 
                                             filterTypes = c("central", "dispersion"), filterCentralType = "median", filterDispersionType = "cv", sequential = T,
                                             filterDispersionPercentile = 0.3, filterCentralPercentile = 0.3)
nrow(datExpr0)  #7552 for network
View(phenotype_data)
all(rownames(var_combined)==rownames(phenotype_data))  # TRUE
all(rownames(var_combined)==colnames(datExpr0)) # TRUE

dat_traits_hd_0 <- var_combined[which(phenotype_data$Condition=="HD"),]  # HD only clinical traits
View(dat_traits_hd_0)
CAG_repeat <- data.frame(CAG = phenotype_data$CAG, row.names = phenotype_data$Sample)
View(CAG_repeat)
CAG_repeat <- CAG_repeat[which(CAG_repeat$CAG !="NA"),, drop=F]
View(CAG_repeat)
dat_traits_hd <- dat_traits_hd_0[rownames(CAG_repeat),]
View(dat_traits_hd)
all(rownames(dat_traits_hd)==rownames(CAG_repeat))  # TRUE
CAG_repeat <- as.numeric(substring(CAG_repeat$CAG, 1,2)) # Takes care of CAGs with range. We use only the higher number
dat_traits_hd<- cbind(dat_traits_hd, CAG_repeat)
View(dat_traits_hd)
dat_traits_hd <- dat_traits_hd[,-1] #removes the control covariate
dat_traits_hd <- dat_traits_hd[, -c(1,2)] # also takes away HD vs Juvenile
dat_traits_hd <- data.frame(dat_traits_hd, row.names = rownames(dat_traits_hd))
datExpr <- datExpr0[,rownames(dat_traits_hd)]
all(rownames(dat_traits_hd)==colnames(datExpr)) # TRUE


# We first check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK # TRUE

# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
datExpr <- t(datExpr)
sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches

pdf(file = "sampleClustering.pdf", width = 12, height = 16);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# No outlier

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(dat_traits_hd, signed = T);
# Plot the sample dendrogram and the colors underneath.
pdf(file = "Sample dendrogram and trait heatmap.pdf", height = 12, width = 9)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(dat_traits_hd),
                    main = "Sample dendrogram and trait heatmap")

dev.off()


save(datExpr, dat_traits_hd, dat_traits_hd, file = "HDOnly-01-dataInput.RData")


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr,networkType = "signed", powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=50,col="red")

# we will use a softthreshold of 16
?blockwiseModules

net = blockwiseModules(datExpr, power = 16, networkType = "signed",
                       TOMType = "signed", minModuleSize = 30, maxBlockSize = 8000,
                       deepSplit = 2,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "HDonlyTOM",
                       verbose = 3)

table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
pdf(file = "clusterdendrogram.pdf", height = 12, width = 9)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "HDOnly-02-networkConstruction-auto.RData")

range(table(moduleColors))
module_and_sizes <- data.frame(table(moduleColors))
pdf(file = "module_sizes.pdf", height = 8, width = 6)
ggtexttable(module_and_sizes,rows = NULL, cols = c("Module", "Size"), theme = ttheme("mOrange"))
dev.off()


#Module trait association


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, dat_traits_hd, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


pdf(file = "module_trait_relaitonships.pdf", height = 10, width = 6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(dat_traits_hd),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()

View(dat_traits_hd)

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable CAG containing the weight column of datTrait
CAG = as.data.frame(dat_traits_hd$CAG_repeat)
names(CAG) = "CAG"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, dat_traits_hd$CAG_repeat, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(CAG), sep="");
names(GSPvalue) = paste("p.GS.", names(CAG), sep="");

# Intramodular analysis: identifying genes with high GS and MM
module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
pdf(file = "gene_significance_vs_module_membership.pdf", height = 7, width = 7)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for CAG Length",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

module = "pink"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
pdf(file = "gene_significance_vs_pink_module_membership.pdf", height = 7, width = 7)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for CAG Length",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()



geneInfo0 = data.frame(geneSymbol= colnames(datExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, CAG, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
  {
  oldNames = colnames(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.CAG));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")
View(geneInfo)



#Visualizing the gene network
library(WGCNA)
dissTOM = 1-TOMsimilarityFromExpr(datExpr,TOMType = "signed", power = 16);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^17;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotTOM, geneTree, moduleColors,col=myheatcol,main = "Network heatmap plot")



# Module Enrichment Analysis using ModuleGO in DGCA

View(geneInfo)
modules <- geneInfo[,c("geneSymbol","moduleColor")]
View(modules)
colnames(modules) <- c("values", "ind")
rownames(modules)= c()
str(modules)
View(modules)

library(GOstats, quietly = TRUE)
library(HGNChelper, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
library(DGCA)
moduleGO_res <- moduleGO(genes = modules$values, labels = modules$ind, universe = colnames(datExpr), pval_GO_cutoff = 1)
moduleGO_df = extractModuleGO(moduleGO_res)

library(ggplot2, quietly = T)
pdf(file = "module_enrichment.pdf", height = 12, width = 9)
plotModuleGO(moduleGO_df, nTerms = 8, text_size = 8,coord_flip = F)
dev.off()


### Differential Correlation Analysis
# We intend to determine if there is a difference in the correlation of genes in datExpr with HTT in controls vs hd

# Determine if HTT  is in the list of 5069 filtered genes in datExpr

sum(str_detect(colnames(datExpr0), "HTT", negate=F) ) # 0 HTT not in the list
sum(str_detect(rownames(batch_corrected_matrix), "HTT", negate=F) ) # Found
# Retrieve HTT from the expression matrix and add to the 5069 genes
library(dplyr)
library(tidyverse)
HTT  <-  as.data.frame(batch_corrected_matrix["HTT",], rownames = colnames(normalized_expression_matix))
colnames(HTT)= "HTT"
View(HTT)
nrow(HTT)
nrow(datExpr0)
library(DGCA)
dgca_matrix= filterGenes(batch_corrected_matrix, 
                      filterTypes = c("central", "dispersion"), filterCentralType = "median", filterDispersionType = "cv", sequential = T,
                      filterDispersionPercentile = 0.5, filterCentralPercentile = 0.5)
dgca_matrix <- cbind(HTT, t(dgca_matrix))
View(dgca_matrix) #3854 genes
View(dgca_matrix) 

## Identify genes whose expression correlate differently between controls vs hd## Identify genes whose expression correlate differently between controls vs hd

disease_status <- as.data.frame(phenotype_data[, "Condition", drop=F], row.names = phenotype_data$Sample)
View(disease_status)
all(rownames(disease_status)==rownames(dgca_matrix))  #TRUE
design_mat = model.matrix(~ disease_status$Condition + 0)
View(design_mat)
colnames(design_mat) = c("control", "hd")

library(GOstats, quietly = TRUE)
library(HGNChelper, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
ddcor_res_HTT = ddcorAll(inputMat = t(dgca_matrix), design = design_mat,
                         compare = c("control", "hd"),
                         splitSet = "HTT",
                         classify = T,
                         adjust = "perm", nPerms = 1000,
                         heatmapPlot = F)

View(ddcor_res_HTT)
write.csv(ddcor_res_HTT, file = "differential_correlation_with_HTT.csv")

signicant_change_in_corrleation <-ddcor_res_HTT[ddcor_res_HTT$pValDiff_adj<0.01,] #489
View(signicant_change_in_corrleation)
table(signicant_change_in_corrleation$Classes)
final_gain_of_correlation <- signicant_change_in_corrleation[signicant_change_in_corrleation$zScoreDiff>0,] #312 # this list was sent to gprofiler for enrichment map 
final_loss_of_correlation <- signicant_change_in_corrleation[signicant_change_in_corrleation$zScoreDiff<0,] #177
write.csv(signicant_change_in_corrleation, file = "significant_change_in_correlation.csv")
write.csv(final_gain_of_correlation, file = "final_gain_of_correlation.csv")
write.csv(final_loss_of_correlation, file = "final_loss_of_correlation.csv")

top_12_gain <- arrange(signicant_change_in_corrleation, desc(signicant_change_in_corrleation$zScoreDiff))[1:12,]
View(top_12_gain)
top_12_loss <-arrange(signicant_change_in_corrleation, zScoreDiff)[1:12,]
View(top_12_loss)

# T0p 10 Altered Correlation by adj_pvalue
library(ggplot2)
library(ggpubr)
top_12_cor <- top_12_gain$Gene1
down_12_cor <- top_12_loss$Gene1
View(dgca_matrix)
View(design_mat)

# Plot with top 12 lesions with a significant positive gain with their correlation with HTT
p1 <- ggscatter(dgca_matrix[which(disease_status$Condition=="Control" ),], x = "HTT", y = top_12_cor, 
          add = "reg.line", conf.int = TRUE, combine = T,
          cor.coef = TRUE, cor.method = "pearson",title ="Correlation of Genes with HTT in Control",
          color= "blue",
          xlab = "Normalized HTT Expression", ylab = "Normalized Gene Expression")

p2 <- ggscatter(dgca_matrix[which(disease_status$Condition=="HD" ),], x = "HTT", y = top_12_cor, 
          add = "reg.line", conf.int = TRUE, combine = T,
          cor.coef = TRUE, cor.method = "pearson",title ="Correlation of Genes with HTT in HD",
          color= "red",
          xlab = "Normalized HTT Expression", ylab = "Normalized Gene Expression")



pdf(file = "gain_in_correlation.pdf", height = 9, width = 16)
ggarrange(p1, p2, nrow = 1, ncol = 2)
dev.off()

p3 <- ggscatter(dgca_matrix[which(disease_status$Condition=="Control" ),], x = "HTT", y = down_12_cor, 
                add = "reg.line", conf.int = TRUE, combine = T,
                cor.coef = TRUE, cor.method = "pearson",title ="Correlation of Genes with HTT in Control",
                color= "blue",
                xlab = "Normalized HTT Expression", ylab = "Normalized Gene Expression")

p4 <- ggscatter(dgca_matrix[which(disease_status$Condition=="HD" ),], x = "HTT", y = down_12_cor, 
                add = "reg.line", conf.int = TRUE, combine = T,
                cor.coef = TRUE, cor.method = "pearson",title ="Correlation of Genes with HTT in HD",
                color= "red",
                xlab = "Normalized HTT Expression", ylab = "Normalized Gene Expression")



pdf(file = "loss_in_correlation.pdf", height = 9, width = 16)
ggarrange(p3, p4, nrow = 1, ncol = 2)
dev.off()

# Classes of gained correlation visualized
types_of_gained_correlation <- data.frame(table(final_gain_of_correlation$Classes))
View(types_of_gained_correlation)
types_of_gained_correlation <- types_of_gained_correlation[which(types_of_gained_correlation$Freq>0),]
colnames(types_of_gained_correlation) <- c("Class", "Frequency")

View(types_of_gained_correlation)

labs <- paste0(types_of_gained_correlation$Class, " (", types_of_gained_correlation$Frequency, ")")
View(labs)
types_of_gained_correlation$labs <- labs

pdf(file = "classes_of_gained_correlation.pdf", height = 9, width = 9)
ggdonutchart(data = types_of_gained_correlation, x = "Frequency", label = "labs", color = "white",
                     lab.font = "white",fill = "Class",lab.pos = "in", 
                     orientation = "horizontal")+
  ggsci::scale_fill_jama()+
  theme(legend.position = "left")

dev.off()

# Classes of lost correlation visualized

types_of_lost_correlation <- data.frame(table(final_loss_of_correlation$Classes))
types_of_lost_correlation <- types_of_lost_correlation[which(types_of_lost_correlation$Freq>0),]
colnames(types_of_lost_correlation) <- c("Class", "Frequency")

View(types_of_lost_correlation)

labs <- paste0(types_of_lost_correlation$Class, " (", types_of_lost_correlation$Frequency, ")")
View(labs)
types_of_lost_correlation$labs <- labs

pdf(file = "classes_of_lost_correlation.pdf", height = 9, width = 9)
ggdonutchart(data = types_of_lost_correlation, x = "Frequency", label = "labs", color = "white",
             lab.font = "white",fill = "Class",lab.pos = "in", 
             orientation = "horizontal")+
  ggsci::scale_fill_jco()+
  theme(legend.position = "left")

dev.off()

gprofiler_enrichment <- read.csv(file = "/Users/nanaaffoh/Documents/Documents - Kenneth’s MacBook Pro/Research/Osama_HD_New/osama_hd_new/New_DGCA/gain_of_correlation_gProfiler_hsapiens_12-18-2021_5-46-28 PM__intersections.csv", header=T)
View(gprofiler_enrichment)



tf_barblot <- gprofiler_enrichment[which(gprofiler_enrichment$source=="TF"),]
tf_barblot <- as_tibble(tf_barblot)
View(tf_barblot)
library(ggpubr)

pdf(file = "transcription_factor_enriched_gain_of_correlation.pdf", height = 3, width = 8)
ggbarplot(tf_barblot[1:15,], x = "term_name", y = "negative_log10_of_adjusted_p_value",
          fill="turquoise",
          palette="jama",
          width= 0.2,
          color= "white", # Set bar border colors to white
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = F,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Negative log (adjusted_p_value)",
          legend.title = "Transcription Factor Binding Motifs Enriched in Genes with gain of Correlation",
          rotate = T,
          ggtheme = theme_minimal()
)
          
dev.off()


ddcorGO_res = ddcorGO(ddcor_res_HTT, universe = colnames(dgca_matrix), adjusted = T,pval_gene_thresh = 0.05,
                      gene_ontology = "all", HGNC_clean = TRUE, HGNC_switch = TRUE, annotation = "org.Hs.eg.db", calculateVariance = TRUE)
names(ddcorGO_res)


significant_gain_of_correlation_genes <- ddcorGO_res[[1]]
significant_loss_of_correlation_genes <- ddcorGO_res[[2]]
View(significant_gain_of_correlation_genes)
View(ddcor_res_HTT[ddcor_res_HTT$Gene1==significant_gain_of_correlation_genes,])
all(significant_gain_of_correlation_genes %in% ddcor_res_HTT$Gene1)
sig_gain <- ddcor_res_HTT$Gene1[significant_gain_of_correlation_genes]
View(sig_gain)
View(ddcor_res_HTT[sig_gain,])
enrichment_gain_of_correlation_BP <- ddcorGO_res[[3]][[1]]
 
enrichment_loss_of_correlation_BP <- ddcorGO_res[[4]][[1]]

View(ddcor_res_HTT)

View(ddcorGO_res)

write_csv(as.data.frame(significant_gain_of_correlation_genes), file = "significant_gain_of Correlation_with_HTT.csv")
write_csv(as.data.frame(significant_loss_of_correlation_genes), file = "significant_loss_of Correlation_with_HTT.csv")
write_csv(as.data.frame(enrichment_gain_of_correlation_BP), file = "BP_enrichment_significant_gain_of correlation_with_HTT.csv")
write_csv(as.data.frame(enrichment_loss_of_correlation_BP), file = "BP_enrichment_significant_loss_of correlation_with_HTT.csv")

rm(enrichment_gain_of_correlation_BP, enrichment_loss_of_correlation_BP)
datalist_gain_of_correlation_enrichment = list()
for(i in 1:3){
  dat <- ddcorGO_res[[3]][[i]]  # extracts type of G). i=1 is BP, i=2 is MF, i=3 is CC
  dat$i <- i # keeps track of which iteration produced whatt
  datalist_gain_of_correlation_enrichment[[i]] <- dat # adds to list
}
View(datalist_gain_of_correlation_enrichment)
enrichment_of_gain_of_correlation <- bind_rows(datalist_gain_of_correlation_enrichment)
View(enrichment_of_gain_of_correlation)
enrichment_of_gain_of_correlation <-  enrichment_of_gain_of_correlation %>% arrange(Pvalue)

write_csv(enrichment_of_gain_of_correlation, file = "enrichment_gain_of_correlation_with_HTT.csv")


datalist_loss_of_correlation_enrichment = list()
for(i in 1:3){
  dat <- ddcorGO_res[[4]][[i]]  # extracts type of G). i=1 is BP, i=2 is MF, i=3 is CC
  dat$i <- i # keeps track of which iteration produced whatt
  datalist_loss_of_correlation_enrichment[[i]] <- dat # adds to list
}
View(datalist_loss_of_correlation_enrichment)
enrichment_of_loss_of_correlation <- bind_rows(datalist_gain_of_correlation_enrichment)
View(enrichment_of_loss_of_correlation)
enrichment_of_loss_of_correlation <-  enrichment_of_loss_of_correlation %>% arrange(Pvalue)

write_csv(enrichment_of_loss_of_correlation, file = "enrichment_loss_of_correlation_with_HTT.csv")

save(dgca_matrix, design, geneInfo, phenotype_data,
     file = "hd_allsamples-01-dgca.RData")
write.csv(t(dgca_matrix), file = "background_dgca.csv")





# Regression Analysis

View(phenotype_data)
all(colnames(batch_corrected_matrix)==phenotype_data$Sample) ##TRUE
regression_hd_exp <- batch_corrected_matrix[, which(phenotype_data$Condition=="HD")]
View(regression_hd_exp)
regression_hd_traits <- phenotype_data[phenotype_data$Condition=="HD",]
View(regression_hd_traits)
regression_hd_exp <- regression_hd_exp[c("HTT", rownames(regression_hd_exp!="HTT")),]
View(regression_hd_exp)
regression_hd_exp <- t(regression_hd_exp)

# Add age, sex, length of repeat, and locaton to the expression dataset
all(rownames(regression_hd_exp)==regression_hd_traits$Sample) #TRUE
View(regression_hd_traits)
library(WGCNA)
options(stringsAsFactors = F)
hd_trait_CAGs <- regression_hd_traits[which(regression_hd_traits$CAG!="NA"),]
View(hd_trait_CAGs)
hd_trait_CAGs["CAG"] <- as.numeric(substring(hd_trait_CAGs$CAG, 0,2)) # since some CAG lengths are given with 2 peaks, we will choose only the higer number
View(hd_trait_CAGs)

#Binarize Categorical Variables
# Binarize Location
binarize_location <- binarizeCategoricalVariable(hd_trait_CAGs$Location, includePairwise = FALSE, includeLevelVsAll = TRUE)
View(binarize_location)

# Binarize gender 
binarize_gender <- binarizeCategoricalVariable(hd_trait_CAGs$Gender,includePairwise = T,includeLevelVsAll = F)
View(binarize_gender)
colnames(binarize_location) <- c("Accumbens", "Caudate", "Cingulate")
colnames(binarize_gender) <- "Male" # Gender, 1 is male, 0 is female
View(hd_trait_CAGs)
combined_regression_traits <- cbind(hd_trait_CAGs$CAG, hd_trait_CAGs$Age, binarize_location, binarize_gender)
rownames(combined_regression_traits) <- hd_trait_CAGs$Sample
View(combined_regression_traits)
colnames(combined_regression_traits)[1] <- "CAG"
colnames(combined_regression_traits)[2] <- "Age"
View(regression_hd_exp)
new_hd_regression_exp <- regression_hd_exp[rownames(combined_regression_traits),]
View(new_hd_regression_exp)
regression_combined <- cbind(combined_regression_traits, new_hd_regression_exp)
View(regression_combined)

library(broom)
library(dplyr)
library(tidyr)
library(purrr)


regression_combined <- as_tibble(regression_combined, rownames=NA)
View(regression_combined)
models_regression <- regression_combined %>% 
  pivot_longer(
    cols = colnames(new_hd_regression_exp),
    names_to = "gene_name",
    values_to = "gene_value"
  )
View(models_regression)

regression_table_hd <- models_regression %>% split(.$gene_name) %>%
  map(~lm(gene_value ~ CAG + Age + Male +Accumbens + Caudate + Cingulate, data = .)) %>%
  tibble(
    dvsub = names(.),
    untidied = .
    ) %>%
  mutate(tidy = map(untidied, broom::tidy)) %>%
  unnest(tidy)
View(regression_table_hd)
regression_table_hd$p_adjust <- p.adjust(regression_table_hd$p.value, method = "fdr")
View(regression_table_hd)
write.csv(as.data.frame(regression_table_hd[,-2]), file = "regression_hd.multivariable.csv")

glance_table_hd <- models_regression %>% split(.$gene_name) %>%
  map(~lm(gene_value ~ CAG + Age + Male +Accumbens + Caudate + Cingulate, data = .)) %>%
  tibble(
    dvsub = names(.),
    untidied = .
  ) %>%
  mutate(tidy = map(untidied, broom::glance)) %>%
  unnest(tidy)
View(glance_table_hd)
glance_table_hd$p_adjust <- p.adjust(glance_table_hd$p.value, method = "fdr")

save(regression_table_hd, models_regression, glance_table_hd, regression_combined, file = "Regression_HD.RData")

write.csv(as.data.frame(glance_table_hd[,-2]), file = "glance_hd.multivariable.csv")
significant_models <- glance_table_hd[,-2]
View(significant_models)
significant_models <- significant_models[which(significant_models$p_adjust <0.01),] #3122
View(significant_models)
significant_genes <- regression_table_hd[regression_table_hd$term=="CAG",-2] #20276
View(significant_genes)
significant_genes <- significant_genes[which(significant_genes$p_adjust <0.01),] #685
View(significant_genes)
genes <- significant_genes$dvsub %in% significant_models$dvsub 
sig_genes_sig_models <- significant_genes[genes,] #548
View(sig_genes_sig_models)
range(sig_genes_sig_models$estimate)  #ranges from -0.3334747 to 0.446870
write.csv(sig_genes_sig_models, file = "significant_genes_in_significant_models.csv") #model is gene = int+cag+age+gender+location
models <- significant_models$dvsub %in% significant_genes$dvsub
sig_models_for_sig_genes <- significant_models[models,] 
View(sig_models_for_sig_genes) #548
range(sig_models_for_sig_genes$adj.r.squared) #0.2507877 to 0.5749975
write.csv(sig_models_for_sig_genes, file = "significant_models_for_significant_genes.csv")
increase_with_CAG <- sig_genes_sig_models[which(sig_genes_sig_models$estimate>0),]
View(increase_with_CAG) #376
write.csv(increase_with_CAG, file = "increase_with_CAG.csv")
decrease_with_CAG <- sig_genes_sig_models[which(sig_genes_sig_models$estimate<0),]
write.csv(decrease_with_CAG, file = "decrease_with_CAG.csv")
View(decrease_with_CAG) #172



#Regreession for normal

View(normalized_expression_matix)
View(phenotype_data)
all(colnames(normalized_expression_matix)==phenotype_data$Sample) ##TRUE
regression_ctrl_exp <- normalized_expression_matix[, which(phenotype_data$Condition=="Control")]
View(regression_hd_exp)
regression_ctrl_traits <- phenotype_data[phenotype_data$Condition=="Control",]
View(regression_ctrl_traits)
regression_ctrl_exp <- regression_ctrl_exp[c("HTT", rownames(regression_ctrl_exp!="HTT")),]
View(regression_ctrl_exp)
regression_ctrl_exp <- t(regression_ctrl_exp)

# Add age, sex, and locaton to the expression dataset
all(rownames(regression_ctrl_exp)==regression_ctrl_traits$Sample) #TRUE
View(regression_ctrl_traits)
library(WGCNA)
options(stringsAsFactors = F)

#Binarize Categorical Variables
# Binarize Location
binarize_location_ctrl <- binarizeCategoricalVariable(regression_ctrl_traits$Location, includePairwise = FALSE, includeLevelVsAll = TRUE)
View(binarize_location_ctrl)

# Binarize gender 
binarize_gender_ctrl <- binarizeCategoricalVariable(regression_ctrl_traits$Gender,includePairwise = FALSE,includeLevelVsAll = TRUE)
View(binarize_gender_ctrl)
colnames(binarize_location_ctrl) <- c("Accumbens", "Caudate", "Cingulate")
colnames(binarize_gender_ctrl) <- c("Female", "Male")
combined_regression_traits_ctrl <- cbind(regression_ctrl_traits$Age, binarize_location_ctrl, binarize_gender_ctrl)
rownames(combined_regression_traits_ctrl) <- regression_ctrl_traits$Sample
View(combined_regression_traits_ctrl)
colnames(combined_regression_traits_ctrl)[1] <- "Age"
View(regression_ctrl_exp)
new_ctrl_regression_exp <- regression_ctrl_exp[rownames(combined_regression_traits_ctrl),]
View(new_ctrl_regression_exp)
regression_combined_ctrl <- cbind(combined_regression_traits_ctrl, new_ctrl_regression_exp)
View(regression_combined)

library(broom)
library(dplyr)
library(tidyr)
library(purrr)


regression_combined_ctrl <- as_tibble(regression_combined_ctrl, rownames=NA)
View(regression_combined_ctrl)
models_regression_ctrl <- regression_combined_ctrl %>% 
  pivot_longer(
    cols = colnames(new_ctrl_regression_exp),
    names_to = "gene_name",
    values_to = "gene_value"
  )
View(models_regression_ctrl)

regression_table_ctrl <- models_regression_ctrl %>% split(.$gene_name) %>%
  map(~lm(gene_value ~ Age+ Female + Male +Accumbens + Caudate + Cingulate, data = .)) %>%
  tibble(
    dvsub = names(.),
    untidied = .
  ) %>%
  mutate(tidy = map(untidied, broom::tidy)) %>%
  unnest(tidy)
View(regression_table_ctrl)
write.csv(as.data.frame(regression_table_ctrl[,-2]), file = "regression_ctrl.multivariable.csv")

glance_table_ctrl <- models_regression_ctrl %>% split(.$gene_name) %>%
  map(~lm(gene_value ~ Age + Female + Male +Accumbens + Caudate + Cingulate, data = .)) %>%
  tibble(
    dvsub = names(.),
    untidied = .
  ) %>%
  mutate(tidy = map(untidied, broom::glance)) %>%
  unnest(tidy)
View(glance_table_ctrl)
write.csv(as.data.frame(glance_table_ctrl[,-2]), file = "glance_ctrl.multivariable.csv")
significant_models_ctrl <- glance_table_ctrl[,-2]
View(significant_models_ctrl)
significant_models_ctrl <- significant_models_ctrl[which(significant_models_ctrl$p.value<0.05),]  #11,384
write.csv(as.data.frame(significant_models_ctrl), file = "ctrl_genes_varying_with_ctrl_model.csv")


library(ggplot2)
ggplot(regression_combined, aes(x=CAG, y=(HTT/ACTB))) + 
  geom_point()+
  geom_smooth(method=lm, color="black")+
  labs(title="HTT Expression with CAG Expansion",
       x="CAG Repeat Length", y = "HTT Expression")+
  theme_classic() 

# Visualizing the regression

# Top 25 increasing with CAG and Top 25 decreasing with CAG
top_25_incr <- arrange(increase_with_CAG, desc(estimate))[1:25,1]
top_25_decr <- arrange(decrease_with_CAG, estimate)[1:25,1]
View(decrease_with_CAG)
dot_plot_genes <- rbind(top_25_incr, top_25_decr)
dot_plot_select <- match(dot_plot_genes$dvsub, sig_genes_sig_models$dvsub)
View(dot_plot_genes)
df_dot_plot <- sig_genes_sig_models[dot_plot_select,]
View(df_dot_plot)
df_dot_plot$group <- if_else(df_dot_plot$estimate>0, "increases with CAG repeat expansion", "decreases with CAG repeat expansion")

library(ggpubr)


pdf(file="plot_of_rgeression_weights_50_genes.pdf", width = 8, height = 8)
ggbarplot(df_dot_plot, x = "dvsub", y = "estimate",
          fill = "group",           # change fill color by group            
          color= "white", # Set bar border colors to white
          palette = "jama",            
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Regression Weight of CAG Repeat Length",
          xlab="Genes",
          legend.title = "Effect of CAG Expansion on Transcript Levels",
          rotate = TRUE,
          ggtheme = theme_minimal()
)
dev.off()

# Plot of amount of variation explained by significant models, where CAG had a significant regression weight
View(sig_models_for_sig_genes)
top_50_var <- arrange(sig_models_for_sig_genes, desc(adj.r.squared))[1:50,1]
var_genes <- match(top_50_var$dvsub, sig_models_for_sig_genes$dvsub)
df_var <- sig_models_for_sig_genes[var_genes,]
df_var$adj.r.squared <- df_var$adj.r.squared*100
View(df_var)

var_groups <- data_frame(if_else(sig_genes_sig_models$estimate>0, "Increases with CAG Expansion", "Decreases with CAG Expansion"),dvsub = sig_genes_sig_models$dvsub)
View(var_groups)
colnames(var_groups)[1] <- "Effect of CAG Expansion on Transcript Levels"
df_var <- left_join(df_var,var_groups)
View(df_var)
pdf(file="plot_of_rsquaresof_50_genes.pdf", width = 8, height = 8)
ggdotchart(df_var, x = "dvsub", y = "adj.r.squared",
           color = "Effect of CAG Expansion on Transcript Levels",                              # Color by groups
           palette = "jama",            # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           rotate = TRUE,                                # Rotate vertically
           dot.size = 5,                                 # Large dot size
           y.text.col = TRUE,                            # Color y text by groups
           ylab = "Adjusted R Square (%)",
           xlab = "Genes",
           ggtheme = theme_pubr()                        # ggplot2 theme
)+
  theme_cleveland() 
dev.off()

enrichment_with_CAG <- read.csv("/Users/nanaaffoh/Documents/Documents - Kenneth’s MacBook Pro/Research/Osama_HD_New/osama_hd_new/regression_analysis/HD/increase_with_CAG/gProfiler_hsapiens_10-13-2021_10-08-05 PM__intersections.csv", header = T)
View(enrichment_with_CAG)
biological_processes <- enrichment_with_CAG[which(enrichment_with_CAG$source=="GO:BP"),] %>% arrange(desc(negative_log10_of_adjusted_p_value))
View(biological_processes)
biological_processes <- biological_processes[1:30,]
reactome_pathways<- enrichment_with_CAG[which(enrichment_with_CAG$source=="REAC"),] %>% arrange(desc(negative_log10_of_adjusted_p_value))
reactome_pathways <- reactome_pathways[1:20,]
enrichment_plot <- rbind(biological_processes, reactome_pathways)
View(enrichment_plot)

pdf(file="plot_of_enrichment.pdf", width = 8, height = 8)
ggbarplot(enrichment_plot, x = "term_name", y = "negative_log10_of_adjusted_p_value",
          fill = "source",           # change fill color by group            
          color= "white", # Set bar border colors to white
          palette = "jco",            
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = T,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Negative log(adjusted p_value)",
          legend.title = "Enrichment Term",
          rotate = TRUE,
          ggtheme = theme_minimal()
)
dev.off()

transcription_factor<- enrichment_with_CAG[which(enrichment_with_CAG$source=="TF"),] %>% arrange(desc(negative_log10_of_adjusted_p_value))
View(transcription_factor)
transcription_factor <- transcription_factor[c(1:17,19:30),]

pdf(file="plot_of_transcription_factor.pdf", width = 8, height = 8)
ggbarplot(transcription_factor, x = "term_name", y = "negative_log10_of_adjusted_p_value",
          fill = "turquoise",           # change fill color by group            
          color= "white", # Set bar border colors to white
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = F,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Negative log(adjusted p_value)",
          xlab="Transcription Factor",
          legend.title = "Enrichment Term",
          rotate = TRUE,
          ggtheme = theme_minimal()
)
dev.off()
transcritio_decreased_with_CAG <- read.csv("/Users/nanaaffoh/Documents/Documents - Kenneth’s MacBook Pro/Research/Osama_HD_New/osama_hd_new/regression_analysis/HD/decrease_with_CAG/gProfiler_hsapiens_10-14-2021_5-48-03 AM__intersections.csv", header = T)
View(transcritio_decreased_with_CAG)


transcritio_decreased_with_CAG <- transcritio_decreased_with_CAG[which(enrichment_with_CAG$source=="TF"),] %>% arrange(desc(negative_log10_of_adjusted_p_value))
transcritio_decreased_with_CAG <- transcritio_decreased_with_CAG[1:30,]

pdf(file="plot_of_transcription_factor_decreased_with_CAG.pdf", width = 8, height = 8)
ggbarplot(transcritio_decreased_with_CAG, x = "term_name", y = "negative_log10_of_adjusted_p_value",
          fill = "#005CAB",           # change fill color by group            
          color= "white", # Set bar border colors to white
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = F,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Negative log(adjusted p_value)",
          xlab="Transcription Factor",
          legend.title = "Enrichment Term",
          rotate = TRUE,
          ggtheme = theme_minimal()
)
dev.off()


library(Rtsne)
?Rtsne
tsne_out <- Rtsne(t(batch_corrected_expression_matrix),perplexity = 10)
plot(tsne_out$Y,col=as.factor(phenotype_data$Location),
     pch=19)
legend("bottomleft",
       legend=unique(phenotype_data$Location),
       fill =palette("default"),
       border=NA,box.col=NA)


























