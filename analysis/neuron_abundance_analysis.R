library(ANCOMBC)
library(phyloseq)
library(reshape2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(ggpubr)


#obj loading 
neuron_hd_obj_acc_caud <- readRDS("neuron_hd_obj_acc_caud_sub_3_28.rds")
neuron_hd_obj_cing <- readRDS("neuron_hd_obj_cing.rds")
comb_hd_cd44_pos <- readRDS("comb_hd_cd44_positive.rds")
comb_hd_cd44_neg <- readRDS("comb_hd_cd44_negative.rds")


#Prep metadata by donor
meta_variable <- c("Batch", "Age", "Gender", "Condition", "Donor", "Grade")
split_meta <- split(neuron_hd_obj_acc_caud@meta.data, neuron_hd_obj_acc_caud@meta.data[["Donor"]])

final_meta <- do.call(rbind,lapply(split_meta, function(df){
  df[1,meta_variable]
}) )

final_meta$Condition <- as.factor(final_meta$Condition)
final_meta$Age <- as.numeric(final_meta$Age)

cts <- table(neuron_hd_obj_acc_caud$sub_type_4, neuron_hd_obj_acc_caud$Donor, neuron_hd_obj_acc_caud$Region)

acc_cts <- cts[,,1]
caud_cts <- cts[,,2]


input_cts_acc <- otu_table(as.matrix.data.frame(acc_cts), taxa_are_rows = T)
rownames(input_cts_acc) <- neuron_hd_obj_acc_caud$sub_type_4 %>% table() %>% names()
colnames(input_cts_acc) <- rownames(final_meta)

input_cts_caud <- otu_table(as.matrix.data.frame(caud_cts), taxa_are_rows = T)
rownames(input_cts_caud) <- neuron_hd_obj_acc_caud$sub_type_4 %>% table() %>% names()
colnames(input_cts_caud) <- rownames(final_meta)


ancom_in_acc <- phyloseq(input_cts_acc, sample_data(final_meta))
ancom_in_caud <- phyloseq(input_cts_caud, sample_data(final_meta))

#Abundance analysis
out_acc <- ancombc(phyloseq = ancom_in_acc, formula = "Age + Gender + Batch + Condition", tax_level = NULL,
                   p_adj_method = "BH", prv_cut = 0, lib_cut = 1,
                   group = "Condition", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                   max_iter = 100, conserve = FALSE, alpha = 0.05, global = TRUE,
                   n_cl = 1, verbose = TRUE)

out_caud <- ancombc(phyloseq = ancom_in_caud, formula = "Age + Gender + Batch + Condition", tax_level = NULL,
                    p_adj_method = "BH", prv_cut = 0, lib_cut = 20,
                    group = "Condition", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = FALSE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)

res_acc <- out_acc$res
res_caud <- out_caud$res



prep_df_ancom <- function(res_ancom, region_in = NULL){
  
  df_cond_hd <- lapply(res_ancom, function(df){
    as.vector(df[,"ConditionHD"])
  }) 
  df_cond_hd <- do.call(cbind,df_cond_hd) %>% as.data.frame()
  df_cond_hd$cell_type <- res_ancom$lfc$taxon
  
  df_fig <- df_cond_hd %>% arrange(desc(lfc)) %>%
    mutate(direct = ifelse(lfc > 0, "Positive LFC", "Negative LFC"))
  df_fig$cell_type = factor(df_fig$cell_type, levels = df_fig$cell_type)
  df_fig$direct = factor(df_fig$direct, 
                         levels = c("Positive LFC", "Negative LFC"))
  
  if(is.null(region_in) == F){
    df_fig$region <- region_in
  }
  return(df_fig)
}


df_fig_acc <- prep_df_ancom(res_acc, region_in = "Accumbens")
df_fig_caud <- prep_df_ancom(res_caud, region_in = "Caudate")


############## Plotting individual results #######################
ggplot(df_fig_acc,aes(x = cell_type, y = lfc, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  # geom_text(aes(label = ifelse(diff_abn, "*", "")), 
  #           position = position_dodge(width = .9), vjust = 1, size = 20 / .pt) +
  labs(x = NULL, y = "Log fold change", 
       title = "Abundance analysis of neurons in Accumbens") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

ggplot(df_fig_caud,aes(x = cell_type, y = lfc, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Abundance analysis of neurons in Caudate") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

#############################################################

df_fig <- rbind(df_fig_acc, df_fig_caud)

ggplot(df_fig,aes(x = cell_type, y = lfc, fill = region)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), 
                width = 0.2, position = position_dodge(0.7), color = "black") +
  labs(x = NULL, y = "Log fold change", 
       title = "Abundance analysis of neurons in Accumbens/Caudate") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))




cd44_pos_clus_prop <- table(comb_hd_cd44_pos$seurat_clusters, comb_hd_cd44_pos$Donor) %>% apply(2,function(col){col/sum(col)})
cd44_neg_clus_prop <- table(comb_hd_cd44_neg$seurat_clusters, comb_hd_cd44_neg$Donor) %>% apply(2,function(col){col/sum(col)})
colnames(cd44_pos_clus_prop) <- stringr::str_replace(colnames(cd44_pos_clus_prop), "T","T-")
colnames(cd44_neg_clus_prop) <- stringr::str_replace(colnames(cd44_neg_clus_prop), "T","T-")

#Proportion correlation results 

prep_cor_prop <- function(neuron_cts, cd44_pos_prop, cd44_neg_prop, res_ancom, region){
  
  # res_caud$diff_abn$taxon[which(res_caud$diff_abn$ConditionHD == T)]
  
  sig_clus <- res_ancom$diff_abn$taxon[which(res_ancom$diff_abn$ConditionHD == T)]
  
  neuron_clus_prop <- neuron_cts[,colSums(neuron_cts) != 0] %>% apply(2,function(col){col/sum(col)})
  ord_astro_donor <- colnames(cd44_pos_prop) %in% colnames(neuron_clus_prop)
  
  cd44_pos_clus_prop <- cd44_pos_clus_prop[,ord_astro_donor]
  cd44_neg_clus_prop <- cd44_neg_clus_prop[,ord_astro_donor]
  
  check_ord <- colnames(cd44_pos_clus_prop) == colnames(neuron_clus_prop)
  
  if( sum(check_ord) != length(check_ord) ){
    stop("Donor order is not correct")
  } else{
    print("Donor order pass")
  }
  
  if(length(sig_clus) < 2){
    
    df_neuron_prop <- data.frame("Donor" = colnames(neuron_clus_prop),
                                 "neuron_type" = paste0(region,".",sig_clus),
                                 "neuron_proportion" = neuron_clus_prop[rownames(neuron_clus_prop) %in% sig_clus,])
    
  } else {
    
    df_neuron_prop <- data.frame("Donor" = colnames(neuron_clus_prop),
                                 t(neuron_clus_prop[rownames(neuron_clus_prop) %in% sig_clus,]))
    colnames(df_neuron_prop)[2:dim(df_neuron_prop)[2]] <- paste0(region,".",colnames(df_neuron_prop)[2:dim(df_neuron_prop)[2]] )
    
    
    df_neuron_prop <- df_neuron_prop %>%  reshape2::melt(variable.name = "neuron_type",
                                                         value.name = "neuron_proportion")
    
  }
  
  
  
  # print(df_neuron_prop)
  
  
  df_astro_prop <-  data.frame("Donor" = colnames(neuron_clus_prop),
                               "fibrous"= t(cd44_pos_clus_prop),
                               "protoplasmic" = t(cd44_neg_clus_prop)) %>%  reshape2::melt(variable.name = "astro_type",
                                                                                           value.name = "astro_proportion", 
                                                                                           id = c("Donor"))
  
  df_prop <- dplyr::inner_join(df_neuron_prop, df_astro_prop)
  
  df_prop_split <- split(df_prop, df_prop$neuron_type) %>% lapply(function(df){split(df,df$astro_type)}) %>% unlist(recursive = F)
  
  
  all_fit <- lapply(df_prop_split, function(df){
    
    cor_res <- cor.test(x = df$astro_proportion, y = df$neuron_proportion)
    
    return_df <- data.frame("neuron_type" = df$neuron_type[1],
                            "astro_type" = df$astro_type[1],
                            "cor" = round(cor_res$estimate, digits = 3),
                            "cor_p_val" = round(cor_res$p.value,digits = 3))
    
    
    return(return_df)
    
  }) 
  
  all_fit <- do.call(rbind, all_fit) %>% as.data.frame()
  all_fit <- all_fit[!is.na(all_fit$cor),]
  
  p1 <- ggplot(all_fit, aes(astro_type,neuron_type)) +
    geom_tile(aes(fill = cor)) +
    scale_fill_gradientn(colours = c("green", "white", "yellow")) + 
    geom_text(aes(label = cor)) +
    geom_text(aes(label = paste0("(",cor_p_val,")")), vjust = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size = 12),
          axis.title = element_blank())
  
  p2 <- ggplot(df_prop, aes(x = astro_proportion, y = neuron_proportion)) +
    geom_point() +
    geom_smooth(method = "lm", formula = "neuron_proportion ~ Grade +  astro_proportion") +
    stat_cor(method="pearson") +
    facet_grid(cols = vars(astro_type), rows = vars(neuron_type)) + 
    theme(strip.text.x = element_text(size = 12)) +
    theme(strip.text.y = element_text(size = 12)) +
    labs(title = paste0( region," Neuron - Astro Proportion"))
  
  sum_res <- list("raw_prop_df" = df_prop,
                  "cor_df" = all_fit,
                  "cor_heatmap" = p1,
                  "cor_diagram" = p2)
  
  return(sum_res)
  
  
}

acc_cor_res <- prep_cor_prop(neuron_cts = acc_cts, cd44_pos_clus_prop, cd44_neg_clus_prop, res_acc, "Accumbens")
caud_cor_res <- prep_cor_prop(neuron_cts = caud_cts, cd44_pos_clus_prop, cd44_neg_clus_prop, res_caud, "Caudate")



acc_cor_res$cor_heatmap
caud_cor_res$cor_heatmap




##################################################
############ Cingulate ###########################
##################################################

meta_variable <- c("Batch", "Age", "Gender", "Condition", "Donor", "Grade")

split_meta <- split(neuron_hd_obj_cing@meta.data, neuron_hd_obj_cing@meta.data[["Donor"]])
final_meta <- do.call(rbind,lapply(split_meta, function(df){
  df[1,meta_variable]
})) 


final_meta$Condition <- as.factor(final_meta$Condition)
final_meta$Age <- as.numeric(final_meta$Age)

input_cts_cing <- table(neuron_hd_obj_cing$sub_type_3, neuron_hd_obj_cing$Donor)

input_cts_cing <- otu_table(as.matrix.data.frame(input_cts_cing), taxa_are_rows = T)
rownames(input_cts_cing) <- neuron_hd_obj_cing$sub_type_3 %>% table() %>% names()
colnames(input_cts_cing) <- rownames(final_meta)

ancom_in_cing <- phyloseq(input_cts_cing, sample_data(final_meta))


out_cing <- ancombc(phyloseq = ancom_in_cing, formula = "Age + Gender + Batch + Condition", tax_level = NULL,
                    p_adj_method = "holm", prv_cut = 0, lib_cut = 0,
                    group = "Condition", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                    max_iter = 100, conserve = FALSE, alpha = 0.05, global = TRUE,
                    n_cl = 1, verbose = TRUE)

res_cing <- out_cing$res



df_fig_cing <- prep_df_ancom(res_cing, region_in = "Cingulate")


############## Plotting individual results #######################
ggplot(df_fig_cing,aes(x = cell_type, y = lfc)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", fill = "purple", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Abundance analysis of neurons in Cingulate") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

#############################################################



cd44_pos_clus_prop <- table(comb_hd_cd44_pos$seurat_clusters, comb_hd_cd44_pos$Donor) %>% apply(2,function(col){col/sum(col)})
cd44_neg_clus_prop <- table(comb_hd_cd44_neg$seurat_clusters, comb_hd_cd44_neg$Donor) %>% apply(2,function(col){col/sum(col)})
colnames(cd44_pos_clus_prop) <- stringr::str_replace(colnames(cd44_pos_clus_prop), "T","T-")
colnames(cd44_neg_clus_prop) <- stringr::str_replace(colnames(cd44_neg_clus_prop), "T","T-")


colnames(input_cts_cing)[!(colnames(input_cts_cing) %in% colnames(cd44_pos_clus_prop))]

cts_cing <- table(neuron_hd_obj_cing$sub_type_3, neuron_hd_obj_cing$Donor)
cts_cing <- cts_cing[,colnames(input_cts_cing) != "T-4110"]


cing_cor_res <- prep_cor_prop(neuron_cts = cts_cing, cd44_pos_clus_prop, cd44_neg_clus_prop, res_cing, "Cingulate")


cing_cor_res$cor_heatmap


#pulls out significant cluster 
sig_clus <- res_cing$diff_abn$taxon[which(res_cing$diff_abn$ConditionHD == T)]

neuron_clus_prop <- cts_cing[,colSums(cts_cing) != 0] %>% apply(2,function(col){col/sum(col)})
ord_astro_donor <- colnames(cd44_pos_clus_prop) %in% colnames(neuron_clus_prop)

cd44_pos_clus_prop <- cd44_pos_clus_prop[,ord_astro_donor]
cd44_neg_clus_prop <- cd44_neg_clus_prop[,ord_astro_donor]

check_ord <- colnames(cd44_pos_clus_prop) == colnames(neuron_clus_prop)

if( sum(check_ord) != length(check_ord) ){
  stop("Donor order is not correct")
} else{
  print("Donor order pass")
}


df_neuron_prop <- data.frame("Donor" = colnames(neuron_clus_prop),
                             "neuron_type" = paste0("Cingulate.",sig_clus),
                             "neuron_proportion" = neuron_clus_prop[rownames(neuron_clus_prop) %in% sig_clus,])


print(df_neuron_prop)


df_astro_prop <-  data.frame("Donor" = colnames(neuron_clus_prop),
                             "fibrous"= t(cd44_pos_clus_prop),
                             "protoplasmic" = t(cd44_neg_clus_prop)) %>%  reshape2::melt(variable.name = "astro_type",
                                                                                         value.name = "astro_proportion", 
                                                                                         id = c("Donor"))

df_prop <- dplyr::inner_join(df_neuron_prop, df_astro_prop)
df_prop_split <- split(df_prop, df_prop$neuron_type) %>% lapply(function(df){split(df,df$astro_type)}) %>% unlist(recursive = F)

#apply correlation function
all_fit <- lapply(df_prop_split, function(df){
  
  cor_res <- cor.test(x = df$astro_proportion, y = df$neuron_proportion)
  
  return_df <- data.frame("neuron_type" = df$neuron_type[1],
                          "astro_type" = df$astro_type[1],
                          "cor" = round(cor_res$estimate, digits = 3),
                          "cor_p_val" = round(cor_res$p.value,digits = 3))
  
  
  return(return_df)
  
}) 

all_fit <- do.call(rbind, all_fit) %>% as.data.frame()
all_fit <- all_fit[!is.na(all_fit$cor),]

p1 <- ggplot(all_fit, aes(astro_type,neuron_type)) +
  geom_tile(aes(fill = cor)) +
  scale_fill_gradientn(colours = c("green", "white", "yellow")) + 
  geom_text(aes(label = cor)) +
  geom_text(aes(label = paste0("(",cor_p_val,")")), vjust = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size = 12),
        axis.title = element_blank())

p2 <- ggplot(df_prop, aes(x = astro_proportion, y = neuron_proportion)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "neuron_proportion ~ Grade +  astro_proportion") +
  stat_cor(method="pearson") +
  facet_grid(cols = vars(astro_type), rows = vars(neuron_type)) + 
  theme(strip.text.x = element_text(size = 12)) +
  theme(strip.text.y = element_text(size = 12)) +
  labs(title = paste0( region," Neuron - Astro Proportion"))



