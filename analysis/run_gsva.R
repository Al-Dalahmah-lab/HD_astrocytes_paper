
library(GSVA)
library(Seurat)
library(magrittr)
library(stringr)
library(ggplot2)
source("/functions//pseudo_gsva_general.R")

#Loading object 
comb_hd_cd44_pos <- readRDS("comb_hd_cd44_positive.rds")
comb_hd_cd44_neg <- readRDS("comb_hd_cd44_negative.rds")

#Gene set loading
astro_genemark <- read.csv("/data/astro_genesets.csv")
astro_mark_list <- list(astro_genemark$Protoplasmic_genes[1:47], astro_genemark$MT_score[1:60], astro_genemark$Reactive_score)

astro_mark_list <- list(astro_mark_list[[1]][!is.na(match(astro_mark_list[[1]], rownames(comb_hd_cd44_pos)))],
                        astro_mark_list[[2]][!is.na(match(astro_mark_list[[2]], rownames(comb_hd_cd44_pos)))],
                        astro_mark_list[[3]][!is.na(match(astro_mark_list[[3]], rownames(comb_hd_cd44_pos)))])

increased_cag <- read.csv("/data/increased_cag.csv")
decreased_cag <- read.csv("/mnt/mfs/hgrcgrid/homes/fp2409/hd_snRNA/decreased_cag.csv")
cag_list <- list("increased_cag" = increased_cag[[1]],
                 "decreased_cag" = decreased_cag[[1]])
cag_list <- list("increased_cag" = cag_list$increased_cag[!is.na(match(cag_list$increased_cag, rownames(comb_hd_cd44_pos)))],
                 "decreased_cag" = cag_list$decreased_cag[!is.na(match(cag_list$decreased_cag, rownames(comb_hd_cd44_pos)))])



lipidomic_gene <- c("HSPA1A", "MKNK2", "DNAJB1", "MRPL18", "HSPB1", "PCBD1",
                    "HSPA1B", "DEDD2", "SIAH2", "SNHG7", "ATF4", "RP9P", 
                    "FKBP4", "ZDHHC18", "PPID", "GATAD2A", "HSF1", "FTL",
                    "NAA60", "SLC25A39", "UBB", "STK40", "VPS37B", "ASB6",
                    "RPS19BP1", "SRSF9", "JTB", "TAF10", "TDP52L2")


all_gene_set <- list("Quiescent Geneset" = astro_mark_list[[1]],
                     "Neuroprotective Geneset" = astro_mark_list[[2]],
                     "CAG Correlated Geneset" = cag_list$increased_cag,
                     "RNA Correlated Lipid Geneset" = lipidomic_gene)


#Running the pseudo_gsva algorithm on astrocyte data with location splitting
gsva_res_cd44_pos_loc <- pseudo_gsva_general(comb_hd_cd44_pos,all_gene_set, split_loc = T)
gsva_res_cd44_neg_loc <- pseudo_gsva_general(comb_hd_cd44_neg,all_gene_set, split_loc = T)


#Plotting
gsva_res_cd44_pos_loc$gsva_heat_scaled$mean_var <- factor(gsva_res_cd44_pos_loc$gsva_heat_scaled$mean_var, levels = c("Quiescent Geneset", "Neuroprotective Geneset",
                                                                                                                      "CAG Correlated Geneset", "RNA Correlated Lipid Geneset")
                                                          
gsva_res_cd44_neg_loc$gsva_heat_scaled$mean_var <- factor(gsva_res_cd44_neg_loc$gsva_heat_scaled$mean_var, levels = c("Quiescent Geneset", "Neuroprotective Geneset",
                                                                                                                      "CAG Correlated Geneset", "RNA Correlated Lipid Geneset"))

gsva_res_cd44_pos_loc$gsva_heat_scaled$cluster <- paste0("F",gsva_res_cd44_pos_loc$gsva_heat_scaled$cluster)
ggplot(gsva_res_cd44_pos_loc$gsva_heat_scaled, aes(cluster, Location)) +
  geom_tile(aes(fill = score)) +
  scale_fill_gradientn(colours = c("blue", "yellow", "red")) +
  facet_wrap(~mean_var)+
  labs(title = "GSVA Score Fibrous Astrocytes") +
  theme(axis.text.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 12),
        axis.title = element_blank(),
        title = element_text(size = 15)) 

gsva_res_cd44_neg_loc$gsva_heat_scaled$cluster <- paste0("P",gsva_res_cd44_neg_loc$gsva_heat_scaled$cluster)
ggplot(gsva_res_cd44_neg_loc$gsva_heat_scaled, aes(cluster, Location)) +
  geom_tile(aes(fill = score)) +
  scale_fill_gradientn(colours = c("blue", "yellow", "red")) +
  facet_wrap(~mean_var)+
  labs(title = "GSVA Score Protoplasmic Astrocytes") +
  theme(axis.text.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 12),
        axis.title = element_blank(),
        title = element_text(size = 15)) 



#Running the pseudo_gsva algorithm on astrocyte data with grade splitting
gsva_res_loc <- pseudo_gsva_general(comb_hd,all_gene_set, split_loc = T)
gsva_res_cd44_pos_grade <- pseudo_gsva_general(comb_hd_cd44_pos,all_gene_set, split_loc = T, split_var_loc = "grade")
gsva_res_cd44_neg_grade <- pseudo_gsva_general(comb_hd_cd44_neg,all_gene_set, split_loc = T, split_var_loc = "grade")

ggplot(gsva_res_cd44_pos_grade$gsva_heat_scaled, aes(cluster, grade)) +
  geom_tile(aes(fill = score)) +
  scale_fill_gradientn(colours = c("blue", "yellow", "red")) +
  facet_wrap(~mean_var)+
  labs(title = "GSVA Score Fibrous Astrocytes")

ggplot(gsva_res_cd44_neg_grade$gsva_heat_scaled, aes(cluster, grade)) +
  geom_tile(aes(fill = score)) +
  scale_fill_gradientn(colours = c("blue", "yellow", "red")) +
  facet_wrap(~mean_var)+
  labs(title = "GSVA Score Protoplasmic Astrocytes")

