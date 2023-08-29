
########################
######### Input ########
########################
# *obj: Refers to Seurat object
# *gene_list: A list with multiple gene sets to run GSVA algorithm on
# *split_loc: Boolean which refers whether the algorithm will run on 2 levels of pseudo-bulking
# ie (pseudo-bulking cells from the same donor and running GSVA is 1 level of aggregation and 2 level
# is pseudo-bulking cells from the same donor and region then running GSVA is 2 level)
# *split_var: variable to pseudo-bulk (must be in the meta-data)
# *split_var_loc: The 2nd variable to pseudo-bulk 


pseudo_gsva_general <- function(obj, gene_list,split_loc = F, split_var = "Donor", split_var_loc = "Location"){
  
  
  
  if(split_loc == F){
    
    #preparing first level of aggregation or grouping
    obj$index <- 1:dim(obj)[2]
    cluster_meta <- split(obj@meta.data, obj$seurat_clusters)
    meta_clus_don <- lapply(cluster_meta,function(x){split(x,x[,split_var])})
    
    
    total_clus <- length(table(obj$seurat_clusters))
    clus_name <- names(meta_clus_don)
    
    meta_don_list <- list()
    pseudo_don_clus <- list()
    for(i in 1:total_clus){
      
      donor_sum <- lapply(meta_clus_don[[i]], function(x){
        don_ind <- x$index
        
        #Pseudobulking the samples
        #Take account into cases when only 1 sample of a donor is in a cluster
        if(length(don_ind) != 1){
          apply(obj@assays$RNA@counts[,don_ind],1,sum)
        } else {
          obj@assays$RNA@counts[,don_ind]
        }
        
      })
      print(paste("Cluster",clus_name[i],"donor pseudo bulk completed"))
      expr_mat <- do.call(cbind, donor_sum)
      
      pseudo_don_clus[[i]] <- expr_mat
      meta_don_list[[i]] <- paste0(names(meta_clus_don[[i]]),"_clus",clus_name[i])
      
    }
    
    pseudo_don_clus_counts <- do.call(cbind, pseudo_don_clus)
    colnames(pseudo_don_clus_counts) <- unlist(meta_don_list)
    
    #This is to remove column with no count
    pseudo_don_clus_counts <- pseudo_don_clus_counts[,colSums(pseudo_don_clus_counts) != 0]
    
    #creates Seurat object for normalization 
    pseudo_don_clus_seurat <- Seurat::CreateSeuratObject(counts = pseudo_don_clus_counts)
    pseudo_don_clus_seurat <- Seurat::NormalizeData(pseudo_don_clus_seurat)
    pseudo_don_clus_seurat <- Seurat::SCTransform(pseudo_don_clus_seurat)
    
    #GSVA package takes in contonuous normalized values
    gsva_res <- gsva(as.matrix(pseudo_don_clus_seurat@assays$SCT@data), gene_list, mx.diff=F, method = "gsva")
    
    #extracting metadata name
    meta_don_pseudo <- colnames(pseudo_don_clus_counts)
    pseudo_clus <- str_extract(meta_don_pseudo,'(?<=clus?)\\d+') %>% as.numeric()
    pseudo_don <- str_extract(meta_don_pseudo, "[^_]+")
    
    #scaling procedure 
    scaled_pseudo_gsva <- apply(gsva_res,2,scale) 
    rownames(scaled_pseudo_gsva) <- rownames(gsva_res)
    
    scaled_pseudo_gsva_df <- scaled_pseudo_gsva %>% t() %>% as.data.frame()
    scaled_pseudo_gsva_df$Donor <- pseudo_don
    scaled_pseudo_gsva_df$cluster <- pseudo_clus
    meta_pseudo <- split(scaled_pseudo_gsva_df, scaled_pseudo_gsva_df$cluster)
    
    #Takes average of each pseudo-bulk samples
    mean_meta_pseudo <- lapply(meta_pseudo, function(df){
      gsva_avg <- lapply(names(gene_list), function(gene_set){
        mean(df[,gene_set])
      }) %>% unlist() %>% as.vector()
      names(gsva_avg) <- names(gene_list)
      gsva_avg
    })
    
    mean_meta_pseudo <- do.call(rbind, mean_meta_pseudo) %>% reshape2::melt()
    names(mean_meta_pseudo) <- c("cluster", "mean_var", "score")
    
    
    gsva_res_list <- list("gsva_raw" = gsva_res,
                          "gsva_raw_scaled" = scaled_pseudo_gsva,
                          "df_gsva" = scaled_pseudo_gsva_df,
                          "df_gsva_mean" = mean_meta_pseudo)
    return(gsva_res_list)
  }
  
  else{
    
    obj$index <- 1:dim(obj)[2]
    source("/home/fp2409/useful_func/uneven_list_to_df.R")
    meta_clus_don_loc <- split(obj@meta.data, list(obj@meta.data$seurat_clusters,
                                                   obj@meta.data[[split_var]],
                                                   obj@meta.data[[split_var_loc]]), drop = TRUE)
    #Removes CC region from analysis
    remove_cc <- which(stringr::str_detect( names(meta_clus_don_loc),"CC"))
    if(length(remove_cc) != 0){
      
      meta_clus_don_loc <- meta_clus_don_loc[-c(remove_cc)]
      
    }
    meta_clus_don_loc_df <- uneven_list_to_df(strsplit(names(meta_clus_don_loc),"\\.")) %>% t() %>% as.data.frame()
    colnames(meta_clus_don_loc_df) <- c("Cluster",split_var, split_var_loc)
    # names(meta_clus_don_loc) <- paste0("clus_",names(meta_clus_don_loc))
    
    print("Pseudo bulking starting")
    pseudo_sum <- lapply(meta_clus_don_loc, function(x){
      don_ind <- x$index
      if(length(don_ind) != 1){
        apply(obj@assays$RNA@counts[,don_ind],1,sum)
      } else {
        obj@assays$RNA@counts[,don_ind]
      }
      
    })
    
    print("Pseudo bulk completed")
    
    pseudo_don_clus_counts <- do.call(cbind, pseudo_sum)
    pseudo_don_clus_seurat <- Seurat::CreateSeuratObject(counts = pseudo_don_clus_counts)
    pseudo_don_clus_seurat <- Seurat::NormalizeData(pseudo_don_clus_seurat)
    pseudo_don_clus_seurat <- Seurat::SCTransform(pseudo_don_clus_seurat)
    gsva_res <- gsva(as.matrix(pseudo_don_clus_seurat@assays$SCT@data), gene_list, mx.diff=F, method = "gsva")
    
    scaled_pseudo_gsva <- apply(gsva_res,2,scale) 
    rownames(scaled_pseudo_gsva) <- rownames(gsva_res)
    
    
    
    #Plug in GSVA output
    prep_gsva_df <- function(gsva_output){
      
      pseudo_gsva_df <- gsva_output %>% t() %>% as.data.frame()
      pseudo_gsva_df[split_var] <- meta_clus_don_loc_df[[split_var]]
      pseudo_gsva_df$cluster <- meta_clus_don_loc_df$Cluster
      pseudo_gsva_df[split_var_loc] <- meta_clus_don_loc_df[[split_var_loc]] 
      
      return(pseudo_gsva_df)
    }
    
    
    prep_heat <- function(df){
      gsva_split <- split(df, list(df$cluster, df[[split_var_loc]]))
      var_clus_split <- strsplit(names(gsva_split),"\\.")
      
      df_store <- lapply(names(gsva_split), function(group){
        
        mean_score <- lapply(names(gene_list), function(gene_set){
          
          if(is.na(gsva_split[[group]][,1][1] == T)){
            return(0)
          } else {
            mean(gsva_split[[group]][,gene_set])
          }
        }) %>% unlist()
        
        var_clus_split <- strsplit(group,"\\.")
        
        names(mean_score) <- names(all_gene_set)
        mean_score <- data.frame(mean_score,  "cluster" = as.factor(var_clus_split[[1]][1]))
        
        mean_score[split_var_loc] <- as.factor(var_clus_split[[1]][2])
        mean_score$mean_var <- rownames(mean_score)
        mean_score <- mean_score %>% reshape2::melt()
        
        #Removes redundant column
        return(mean_score[,-4])
        
        
      })
      gsva_loc_clus <- do.call(rbind, df_store) 
      names(gsva_loc_clus)[dim(gsva_loc_clus)[2]] <- "score"
      
      return(gsva_loc_clus)
    }
    
    df_raw_gsva <- prep_gsva_df(gsva_res)
    df_scale_gsva <- prep_gsva_df(scaled_pseudo_gsva)
    
    gsva_heat_raw <- prep_heat(df_raw_gsva)  
    gsva_heat_scaled <- prep_heat(df_scale_gsva)
    
    gsva_res_list <- list("gsva_raw" = gsva_res,
                          "gsva_raw_scaled" = scaled_pseudo_gsva,
                          "df_gsva_mean_raw" = df_raw_gsva,
                          "df_gsva_mean_scale" = df_scale_gsva,
                          "gsva_heat_raw" = gsva_heat_raw,
                          "gsva_heat_scaled" = gsva_heat_scaled)
    
    return(gsva_res_list)
    
    
    
  }
}
