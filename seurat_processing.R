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


setwd(rootdir)  # Set this to correct location

# QC with Seurat:
perform_Seurat_QC<- function(sample_name,
                             count_matrix,
                             type="isoform", #specify if you are working with an isoform- or a gene-level count matrix
                             max.features =2000,
                             min.features = 200,
                             min.counts = 400,
                             max.counts = 5000,
                             MT = 10,
                             npc = 15,
                             output_path=paste0(rootdir,"/QC")
){
  
  #path_to_counts=file.path(paste0(rootdir,"/counts"), paste0("gene_symbol_", sample_name, "transcript_level_counts.csv"))
  
  
  
  df<-count_matrix
  #head(rownames(df)<-df$X)
  #df$X<-NULL
  print("Creation objet Seurat")
  seurat_obj <- CreateSeuratObject(counts = df , min.features = 0)
  
  if ("X" %in% colnames(seurat_obj)) {
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[colnames(seurat_obj) != "X"])
  }
  
  if (type=="isoform") {
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "--MT") #calculate the percentage of all the counts belonging to a subset of the possible features for each cell. 
    print("calcul du pourcentage de MT")
  } else if (type =="gene") {
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") #calculate the percentage of all the counts belonging to a subset of the possible features for each cell. 
    print("calcul du pourcentage de MT")
  } else {
    print("type argument must be isoform or gene")
  }
  plot_scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm") + NoLegend() +
    labs(title =paste0("Reads vs Unique Genes per Cell before filtering, ",sample_name))
  
  p1 <- VlnPlot(seurat_obj, features = c('nFeature_RNA', 'nCount_RNA'))+
    labs(caption = paste0("Features before filtering, ",sample_name))
  
  
  data <- as.data.frame(seurat_obj[["percent.mt"]])
  p2<-ggplot(data = data, aes(x = percent.mt)) +
    geom_histogram(fill = "steelblue", color = "black") +
    geom_vline(xintercept = 10, color = "red", linetype = "dashed", size = 1)+
    labs(
      title = paste0("% MT before filtering, ",sample_name),
      x = "Percentage of mitochondrial DNA",
      y = "Number of cell"
    ) +
    theme_minimal()
  
  #Filtering
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mt < MT & nCount_RNA <max.counts & nCount_RNA >min.counts)
  #Plots after filtering
  p3<- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"))+
    labs(caption = paste0("% Features after filtering, ",sample_name))
  
  data <- as.data.frame(seurat_obj[["percent.mt"]])
  p4<-ggplot(data = data, aes(x = percent.mt)) +
    geom_histogram(fill = "steelblue", color = "black") +
    labs(
      title = paste0("% MT after filtering, ",sample_name),
      x = "Percentage of mitochondrial DNA",
      y = "Number of cell"
    ) +
    theme_minimal()
  #save object and plots
  seurat_obj$sample<-sample_name
  saveRDS(seurat_obj,paste0(rootdir,"/QC/seu_obj_",sample_name,"_",type,".rds"))
  
  
  
  ggsave(filename = file.path(output_path, paste0("meta_features_",type,"_",sample_name, "_plots_before_qc.png")), 
         plot = p1
  ) 
  p1
  ggsave(filename = file.path(output_path, paste0("meta_features_",type,"_",sample_name, "_plots_after_qc.png")), 
         plot = p3
  ) 
  p3
  ggsave(filename = file.path(output_path, paste0("perc_mt_before_qc_",type,"_",sample_name, "_plots.png")), 
         plot = p2
  ) 
  p2
  ggsave(filename = file.path(output_path, paste0("perc_mt_after_qc_",type,"_",sample_name, "_plots.png")), 
         plot = p4
  ) 
  p4
  ggsave(filename = file.path(output_path, paste0("reads_vs_unique_gene",type,"_",sample_name, "_plots_before_qc.png")), 
         plot = plot_scatter1
  ) 
  res=list()
  res$seu_obj<-seurat_obj
  res$features_before<-p1
  res$features_after<-p2
  res$mt_before<-p3
  res$mt_after<-p4
  res$scatter<-plot_scatter1
  print(names(res))
  return(res)
}


# some general functions to have more information about our isoforms
# these can be used on any seurat object since they run on the count layer of an RNA assay

number_of_isoforms_per_gene <- function(seu_obj,sample_name) { 
  
  
  counts<-LayerData(seu_obj, assay = "RNA", layer = "counts")
  as.data.frame(counts) -> df
  row.names(df) -> df$feature
  
  #1 split transcript ids into gene and transcript id
  pseudobulk_data <- df %>% 
    mutate(
      gene_name = sub(".*--", "",feature))|> 
    mutate(transcript_id=sub("--.*", "",feature))
  
  
  
  # 2. Count the number of isoforms per gene
  isoform_count_per_gene <- pseudobulk_data %>%
    group_by(gene_name) %>%
    summarise(n_isoforms = n_distinct(transcript_id))
  
  # 3. count isforms per category 
  isoform_count_per_gene <- isoform_count_per_gene %>%
    mutate(isoform_category = case_when(
      n_isoforms == 1 ~ "1",
      n_isoforms == 2 ~ "2",
      n_isoforms == 3 ~ "3",
      n_isoforms == 4 ~ "4",
      n_isoforms == 5 ~ "5",
      n_isoforms == 6 ~ "6",
      n_isoforms == 7 ~ "7",
      n_isoforms >= 8 ~ "8+"
    ))
  
  
  # 4. Calculate the percentage of genes in each bin
  isoform_count_summary <- isoform_count_per_gene %>%
    dplyr::count(isoform_category) %>%
    mutate(percent = (n / sum(n)) * 100)
  
  plot<-ggplot(isoform_count_summary, aes(x = isoform_category, y = percent)) +
    geom_bar(stat = "identity", fill = "lightblue", color = "black") +  
    labs(title = paste0("Number of Isoforms per Gene ",sample_name), 
         x = "Isoforms per Gene", 
         y = "Genes, %") +
    theme_minimal() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5))
  ggsave(filename=paste0(rootdir,"/analyses/",sample_name,"_number_iso_gene.png"),plot=plot)
  print(plot)
  return(plot)
}





top10_gene_with_most_iso<- function(seu_obj,sample_name) {
  
  
  counts<-LayerData(seu_obj, assay = "RNA", layer = "counts")
  as.data.frame(counts) -> df
  row.names(df) -> df$feature
  
  #split transcript ids into gene and transcript id
  
  pseudobulk_data <- df %>% 
    mutate(
      gene_name = sub(".*--", "",feature))|> 
    mutate(transcript_id=sub("--.*", "",feature))
  
  
  # Genes ranked by the number of transcript isoforms detected across all samples 
  gene_transcript_counts <- pseudobulk_data %>%
    group_by(gene_name) %>%
    summarise(unique_transcripts = n_distinct(transcript_id)) %>%
    arrange(desc(unique_transcripts))
  write.csv(gene_transcript_counts, paste0(rootdir,"/analyses/",sample_name,"_gene_most_iso.csv"))
  # Select the top 10 genes based on unique transcript counts
  top10 <- gene_transcript_counts %>% slice_max(n=10, order_by=unique_transcripts)
  tail10 <- gene_transcript_counts %>% slice_min(n=10, order_by=unique_transcripts)
  top10
  tail10
  
  res=list()
  res$counts<-gene_transcript_counts
  
  # Plot ranked genes by unique "BambuTx" transcript count
  plot<- ggplot(gene_transcript_counts, aes(x = rank(desc(unique_transcripts)), y = unique_transcripts)) +
    geom_point(color = "darkblue", size = 1) +  # Points for each gene
    
    # Log scale for both axes
    scale_x_log10() +
    scale_y_log10() +
    
    # Title and labels
    labs(
      title = paste0("Unique Transcripts per Gene ",sample_name),
      x = "Rank (log scale)",
      y = "# Transcripts (log scale)"
    ) +
    
    # Highlight and label the top 10 genes with gray background and black border around the text
    geom_label_repel(
      data = gene_transcript_counts %>% filter(gene_name %in% top10$gene_name),
      aes(label = gene_name),
      fill = "gray",          # Gray background for the label
      color = "black",         # Black text color
      label.size = 0.25,       # Border thickness around the label
      label.r = unit(0.15, "lines"),  # Border radius (rounded corners)
      size = 3,
      box.padding = 0.2,
      max.overlaps = 25 
    ) +
    
    # Minimal theme and additional styling
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # Centered title
      axis.text = element_text(size = 10, color = "black"),             # Black axis tick labels
      axis.title = element_text(color = "black"),                       # Black axis titles
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Black border around the graph
    )
  ggsave(filename=paste0(rootdir,"/analyses/",sample_name,"_gene_most_iso.png"),plot=plot)
  
  res$plot<-plot
  print(plot)
  print(names(res))
  return(res)
}


# to run a regression on cell cycle genes, ribosomal genes .
# with a scale data version or a SCtransform version 
#should be used instead of NormalizeData+FindvariableFEatures+Scaledata or of the classical SCTransform
# pour le moment on run SCT et RNA classique
#dir.create(paste0(output_path,"/cell_cycle_regression"))
run_regression<-function(seu_obj,
                         type="isoform",
                        # SCT=TRUE,
                         output_path,
                        integrateRPCA=FALSE) {
  
  dir.create(paste0(output_path,"/cell_cycle_regression"))
  if (type=="isoform"){
    
    # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
    # segregate this list into markers of G2/M phase and markers of S phase
    s.genes <- as.data.frame(cc.genes$s.genes)
    colnames(s.genes)<-c("gene_name")
    g2m.genes <- as.data.frame(cc.genes$g2m.genes)
    colnames(g2m.genes)<-c("gene_name")
    
    g2m.features<-as.data.frame(Features(seu_obj))
    colnames(g2m.features)<-c("feature")
    g2m.features<-g2m.features%>%
      mutate(
        gene_name = sub(".*--", "",feature))|> 
      mutate(transcript_id=sub("--.*", "",feature))%>%
      right_join(y=g2m.genes)
    g2m.features <- g2m.features$feature
    
    s.features<-as.data.frame(Features(seu_obj))
    colnames(s.features)<-c("feature")
    s.features<-s.features%>%
      mutate(
        gene_name = sub(".*--", "",feature))|> 
      mutate(transcript_id=sub("--.*", "",feature))%>%
      right_join(y=s.genes) 

    
    s.features <- s.features$feature
    
    ribosomal_genes <- c(
      "RPSA", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y", "RPS5", "RPS6", "RPS7", "RPS8", 
      "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", 
      "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", 
      "RPS27", "RPS27A", "RPS28", "RPS29", 
      "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPL10", "RPL10A", 
      "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", 
      "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A", 
      "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", 
      "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39", "RPL40", "RPL41"
    )
    ribosomal_genes <- as.data.frame(ribosomal_genes)
    colnames(ribosomal_genes)<-c("gene_name")
    
    ribosomal_features<-as.data.frame(Features(seu_obj))
    colnames(ribosomal_features)<-c("feature")
    ribosomal_features<-ribosomal_features%>%
      mutate(
        gene_name = sub(".*--", "",feature))|> 
      mutate(transcript_id=sub("--.*", "",feature))%>%
      right_join(y=ribosomal_genes)#|> select(feature)
    #g2m.features<-as.vector(g2m.features)
    ribosomal_features <- ribosomal_features$feature
    
    head(ribosomal_features)
    
    
    DefaultAssay(seu_obj)<-"RNA"
    
  
    seu_obj<- NormalizeData(seu_obj)
    seu_obj<- FindVariableFeatures(seu_obj, selection.method = "vst")
    seu_obj <- ScaleData(seu_obj, features = rownames(seu_obj))
    seu_obj <- RunPCA(seu_obj, features = VariableFeatures(seu_obj))
    
    seu_obj <- CellCycleScoring(seu_obj, s.features = s.features, g2m.features = g2m.features, 
                                set.ident = TRUE)
    seu_obj <- AddModuleScore(seu_obj, features = list(ribosomal_features), name = "ribo_score")
    

    # view cell cycle scores and phase assignments
    head(seu_obj[[]])
    
    # and now regress out
    #seu_obj <- ScaleData(seu_obj, vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"), features = rownames(seu_obj))
    seu_obj<-  SCTransform(object =seu_obj,vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"))
    seu_obj<- RunPCA(seu_obj)
    integ=""
    if (integrateRPCA==TRUE) {
      seu_obj <- IntegrateLayers(object = seu_obj, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                                    verbose = FALSE, normalization.method = "SCT")
      
      # now that integration is complete, rejoin layers
      seu_obj <- JoinLayers(seu_obj)
      merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.rpca", dims = 1:15)
      merged_obj <- FindClusters(merged_obj, resolution = 1, cluster.name = "rpca_clusters")
      merged_obj <- RunUMAP(merged_obj, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")
      
      merged_obj$cluster<-Idents(merged_obj) 
      
      dimplot<-DimPlot(merged_obj, reduction = "umap.rpca", group.by = c("sample","rpca_clusters"),
                       combine = TRUE, label.size = 2) + labs(caption="Integrated dimplot")
      ggsave(filename = paste0(output_basis,"_merged_integrated_RCPA_dimplot.png"), 
             plot = dimplot
      ) 
      print(dimplot)
      #enregistrer l’objet merged et integrated
      saveRDS(merged_obj, file = paste0(output_basis,"_merged_integrated_RCPA_seurat.rds"))
      integ="integrated"
      
    }
 
    seu_obj <- FindNeighbors(object = seu_obj, dims = 1:30)
    seu_obj <- FindClusters(object = seu_obj)
    seu_obj <- RunUMAP(object = seu_obj, dims = 1:30)
    
    plot <- DimPlot(seu_obj,group.by = "Phase") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
      guides(colour = guide_legend(override.aes = list(size = 10)))
    ggsave(filename = paste0(output_path,"/cell_cycle_regression/dimplot_phase_after_reg",integ,".jpg"), height = 7, width = 12, plot = plot,
           quality = 50)
    ggsave(filename = paste0(output_path,"/cell_cycle_regression/featureplot_ribo_after_reg",integ,".jpg"), height = 7, width = 12, plot = FeaturePlot(seu_obj,features = "ribo_score1"),
           quality = 50)
    
    print(plot)
    print(FeaturePlot(seu_obj,features = "ribo_score1"))
    saveRDS(seu_obj, file = paste0(output_path,"/cell_cycle_regression/regressed_out_iso_seu_obj",integ,".rds"))
    
  } else if (type=="gene") {
    # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
    # segregate this list into markers of G2/M phase and markers of S phase
    s.genes <- cc.genes$s.genes
    #colnames(s.genes)<-c("gene_name")
    g2m.genes <- cc.genes$g2m.genes
    #colnames(g2m.genes)<-c("gene_name")
    ribosomal_genes <- c(
      "RPSA", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y", "RPS5", "RPS6", "RPS7", "RPS8", 
      "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", 
      "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", 
      "RPS27", "RPS27A", "RPS28", "RPS29", 
      "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPL10", "RPL10A", 
      "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", 
      "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A", 
      "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", 
      "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39", "RPL40", "RPL41"
    )
    ribosomal_features_present <- ribosomal_genes[ribosomal_genes %in% Features(seu_obj)]
    
    DefaultAssay(seu_obj)<-"RNA"
    seu_obj<- NormalizeData(seu_obj)
    seu_obj<- FindVariableFeatures(seu_obj, selection.method = "vst")
    seu_obj <- ScaleData(seu_obj, features = rownames(seu_obj))
    seu_obj <- RunPCA(seu_obj, features = VariableFeatures(seu_obj))
    
    seu_obj <- CellCycleScoring(seu_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    seu_obj <- AddModuleScore(seu_obj, features = list(ribosomal_features_present), name = "ribo_score")
    
    # view cell cycle scores and phase assignments
    head(seu_obj[[]])
    
    # and now regress out
    #seu_obj <- ScaleData(seu_obj, vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"), features = rownames(seu_obj))
    seu_obj<-  SCTransform(object =seu_obj,vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"))
    seu_obj<- RunPCA(seu_obj)
    integ=""
    if (integrateRPCA==TRUE) {
      seu_obj <- IntegrateLayers(object = seu_obj, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                                 verbose = TRUE, normalization.method = "SCT")
      
      # now that integration is complete, rejoin layers
      #seu_obj <- JoinLayers(seu_obj)
      seu_obj <- FindNeighbors(seu_obj, reduction = "integrated.rpca", dims = 1:15)
      seu_obj <- FindClusters(seu_obj, resolution = 1, cluster.name = "rpca_clusters")
      seu_obj <- RunUMAP(seu_obj, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")
      
      seu_obj$cluster<-Idents(seu_obj) 
      
      dimplot<-DimPlot(seu_obj, reduction = "umap.rpca", group.by = c("sample","rpca_clusters"),
                       combine = TRUE, label.size = 2) + labs(caption="Integrated dimplot")
      ggsave(filename = paste0(output_basis,"_merged_integrated_RCPA_dimplot.png"), 
             plot = dimplot
      ) 
      print(dimplot)
      #enregistrer l’objet merged et integrated
      saveRDS(seu_obj, file = paste0(output_basis,"_merged_integrated_RCPA_seurat.rds"))
      integ="_integrated"
    }
    seu_obj <- FindNeighbors(object = seu_obj, dims = 1:30)
    seu_obj <- FindClusters(object = seu_obj)
    seu_obj <- RunUMAP(object = seu_obj, dims = 1:30)
    
    plot <- DimPlot(seu_obj,group.by = "Phase") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
      guides(colour = guide_legend(override.aes = list(size = 10)))
    ggsave(filename = paste0(output_path,"/cell_cycle_regression/dimplot_phase_after_reg_gene",integ,".jpg"), height = 7, width = 12, plot = plot,
           quality = 50)
    
    ggsave(filename = paste0(output_path,"/cell_cycle_regression/featureplot_ribo_after_reg_gene",integ,".jpg"), height = 7, width = 12, plot = FeaturePlot(seu_obj,features = "ribo_score1"),
           quality = 50)
    
    print(plot)
    print(FeaturePlot(seu_obj,features = "ribo_score1"))
    saveRDS(seu_obj, file = paste0(output_path,"/cell_cycle_regression/regressed_out_obj_gene",integ,".rds"))
    
  } else {print("type must be isoform or gene")}
  
  return(seu_obj)
}

# Merge a list of seurat objects: merge, preprocess choose to integrate and plot dimplot



perform_merge_seurat<-function(sample_names, # a list of the sample's names, in the right order
                               seu_objects, # a list of seurat objects after QC (same order as their names)
                               type="isoform",
                               output_path=paste0(rootdir,"/merge_integration"),
                               integrateRPCA=FALSE){
  
  
  if (type=="isoform") { 
    
  }
  
  merged_obj <-merge(
    seu_objects[[1]],
    y=seu_objects[2:length(seu_objects)],
    merge.data = TRUE,
    add.cell.ids = sample_names
  )
  output_basis<-paste0(output_path,"/",type,"_",sample_names[1],"_",sample_names[2],"_",sample_names[3])
  if (integrateRPCA==FALSE) {
    merged_obj<- JoinLayers(merged_obj)
    
    #the whole standard preprocess is replace by the "run_regression" function
    
    #merged_obj <- NormalizeData(merged_obj)
    #merged_obj <- FindVariableFeatures(merged_obj)
    #merged_obj <- ScaleData(merged_obj)
    #replacing with the SC transform version:
    #merged_obj <- SCTransform(object = merged_obj)
    #merged_obj <- RunPCA(merged_obj)
    
    #merged_obj<- FindNeighbors(merged_obj, dims = 1:15, reduction = "pca")
    #merged_obj<- FindClusters(merged_obj, resolution = 1)
    
    #merged_obj<- RunUMAP(merged_obj, dims = 1:15, reduction = "pca")
    
    merged_obj<-run_regression0(seu_obj=merged_obj, type, output_path, integrateRPCA=FALSE)
    
    
    output_basis<-paste0(output_path,"/",type,"_",sample_names[1],"_",sample_names[2],"_",sample_names[3])
    
    dimplot<-DimPlot(merged_obj, reduction = "umap", group.by = c("sample"),
                     combine = TRUE, label.size = 2) + plot_annotation(title = paste0("Umap unintegrated at ",type," level",sample_names[[1]]," ",sample_names[[2]]," ",sample_names[[3]]))
    
    ggsave(filename = paste0(output_basis,"_merged_dimplot.png"), 
           plot = dimplot
    )
    print(dimplot)
    saveRDS(merged_obj, file = paste0(output_basis,"_merged_unintegrated_seurat.rds"))
    
    
  } else { # if we want to integrate
    merged_obj<-run_regression(seu_obj=merged_obj, type, output_path, integrateRPCA=TRUE)
    
    
    
    
    dimplot<-DimPlot(merged_obj, reduction = "umap", group.by = c("sample"),
                     combine = TRUE, label.size = 2) + plot_annotation(title = paste0("Umap integrated at ",type," level",sample_names[[1]]," ",sample_names[[2]]," ",sample_names[[3]]))
    
    ggsave(filename = paste0(output_basis,"_merged_integrated_dimplot.png"), 
           plot = dimplot
    )
    print(dimplot)
    saveRDS(merged_obj, file = paste0(output_basis,"_merged_integrated_seurat.rds"))
    
    
  }
  
  res=list()
  res$merged_obj<-merged_obj
  res$dimplot<-dimplot
  print(names(res))
  
  return(res)}


run_regression_and_merge0<-function(sample_names, # a list of the sample's names, in the right order
                                   seu_objects, # a list of seurat objects after QC (same order as their names)
                                   type="isoform",
                                   output_path=paste0(rootdir,"/merge_integration"),
                                   integrateRPCA=FALSE) {
  
  dir.create(paste0(output_path,"/cell_cycle_regression"))
  
  if (integrateRPCA==FALSE) { 
    if (type=="isoform"){
      
      # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
      # segregate this list into markers of G2/M phase and markers of S phase
     integ=""
      
      seu_objects<-lapply(seu_objects,function(x) {
        
        s.genes <- as.data.frame(cc.genes$s.genes)
        colnames(s.genes)<-c("gene_name")
        g2m.genes <- as.data.frame(cc.genes$g2m.genes)
        colnames(g2m.genes)<-c("gene_name")
        
        g2m.features<-as.data.frame(Features(x))
        colnames(g2m.features)<-c("feature")
        g2m.features<-g2m.features%>%
          mutate(
            gene_name = sub(".*--", "",feature))|> 
          mutate(transcript_id=sub("--.*", "",feature))%>%
          right_join(y=g2m.genes)
        g2m.features <- g2m.features$feature
        
        s.features<-as.data.frame(Features(x))
        colnames(s.features)<-c("feature")
        s.features<-s.features%>%
          mutate(
            gene_name = sub(".*--", "",feature))|> 
          mutate(transcript_id=sub("--.*", "",feature))%>%
          right_join(y=s.genes) 
        
        
        s.features <- s.features$feature
        
        ribosomal_genes <- c(
          "RPSA", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y", "RPS5", "RPS6", "RPS7", "RPS8", 
          "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", 
          "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", 
          "RPS27", "RPS27A", "RPS28", "RPS29", 
          "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPL10", "RPL10A", 
          "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", 
          "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A", 
          "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", 
          "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39", "RPL40", "RPL41"
        )
        ribosomal_genes <- as.data.frame(ribosomal_genes)
        colnames(ribosomal_genes)<-c("gene_name")
        
        ribosomal_features<-as.data.frame(Features(x))
        colnames(ribosomal_features)<-c("feature")
        ribosomal_features<-ribosomal_features%>%
          mutate(
            gene_name = sub(".*--", "",feature))|> 
          mutate(transcript_id=sub("--.*", "",feature))%>%
          right_join(y=ribosomal_genes)#|> select(feature)
        #g2m.features<-as.vector(g2m.features)
        ribosomal_features <- ribosomal_features$feature
        
        head(ribosomal_features)
        integ=""
        DefaultAssay(x)<-"RNA"
        
        
        x<- NormalizeData(x)
        x<- FindVariableFeatures(x, selection.method = "vst")
        x <- ScaleData(x, features = rownames(x))
        x <- RunPCA(x, features = VariableFeatures(x))
        
        x <- CellCycleScoring(x, s.features = s.features, g2m.features = g2m.features, 
                                    set.ident = TRUE)
        x <- AddModuleScore(x, features = list(ribosomal_features), name = "ribo_score")
        
        
        # view cell cycle scores and phase assignments
        head(x[[]])
        
        # and now regress out
        #seu_obj <- ScaleData(seu_obj, vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"), features = rownames(seu_obj))
        x<-  SCTransform(object =x,vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"))
        
        return(x)})
        
      merged_obj <-merge(
        seu_objects[[1]],
        y=seu_objects[2:length(seu_objects)],
        merge.data = TRUE,
        add.cell.ids = sample_names
      )
      output_basis<-paste0(output_path,"/",type,"_",sample_names[1],"_",sample_names[2],"_",sample_names[3])
      #merged_obj <- SeuratObject::JoinLayers(merged_obj)
      merged_obj<-SCTransform(object =merged_obj,vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"))
      merged_obj<- RunPCA(merged_obj, assay="SCT")
      
      merged_obj <- FindNeighbors(object = merged_obj, dims = 1:30)
      merged_obj <- FindClusters(object = merged_obj)
      merged_obj <- RunUMAP(object = merged_obj, dims = 1:30)
      
      plot <- DimPlot(merged_obj,group.by = "Phase") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
        guides(colour = guide_legend(override.aes = list(size = 10)))
      ggsave(filename = paste0(output_path,"/cell_cycle_regression/dimplot_phase_after_reg",integ,".jpg"), height = 7, width = 12, plot = plot,
             quality = 50)
      ggsave(filename = paste0(output_path,"/cell_cycle_regression/featureplot_ribo_after_reg",integ,".jpg"), height = 7, width = 12, plot = FeaturePlot(merged_obj,features = "ribo_score1"),
             quality = 50)
      
      print(plot)
      print(FeaturePlot(merged_obj,features = "ribo_score1"))
      saveRDS(merged_obj, file = paste0(output_path,"/cell_cycle_regression/regressed_out_iso_seu_obj",integ,".rds"))
      
    } else if (type=="gene") {
      # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
      # segregate this list into markers of G2/M phase and markers of S phase
      s.genes <- cc.genes$s.genes
      #colnames(s.genes)<-c("gene_name")
      g2m.genes <- cc.genes$g2m.genes
      #colnames(g2m.genes)<-c("gene_name")
      ribosomal_genes <- c(
        "RPSA", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y", "RPS5", "RPS6", "RPS7", "RPS8", 
        "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", 
        "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", 
        "RPS27", "RPS27A", "RPS28", "RPS29", 
        "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPL10", "RPL10A", 
        "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", 
        "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A", 
        "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", 
        "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39", "RPL40", "RPL41"
      )
      
      
      integ=""
      
      seu_objects<-lapply(seu_objects,function(x) {
        ribosomal_features_present <- ribosomal_genes[ribosomal_genes %in% Features(x)]
        DefaultAssay(x)<-"RNA"
        x<- NormalizeData(x)
        x<- FindVariableFeatures(x, selection.method = "vst")
        x <- ScaleData(x, features = rownames(x))
        x <- RunPCA(x, features = VariableFeatures(x))
        
        x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
        x <- AddModuleScore(x, features = list(ribosomal_features_present), name = "ribo_score")
        
        # view cell cycle scores and phase assignments
        head(x[[]])
        
        # and now regress out
        #seu_obj <- ScaleData(seu_obj, vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"), features = rownames(seu_obj))
        x<-  SCTransform(object =x,vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"))
        x
        return(x)})
        
      merged_obj <-merge(
        seu_objects[[1]],
        y=seu_objects[2:length(seu_objects)],
        merge.data = TRUE,
        add.cell.ids = sample_names
      )
      
      output_basis<-paste0(output_path,"/",type,"_",sample_names[1],"_",sample_names[2],"_",sample_names[3])
      
      #merged_obj <- JoinLayers(merged_obj)
      merged_obj<-SCTransform(object =merged_obj,vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"))
      merged_obj<- RunPCA(merged_obj, assay="SCT")
      merged_obj <- FindNeighbors(object = merged_obj, dims = 1:30)
      merged_obj <- FindClusters(object = merged_obj)
      merged_obj <- RunUMAP(object = merged_obj, dims = 1:30)
      
      plot <- DimPlot(merged_obj,group.by = "Phase") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
        guides(colour = guide_legend(override.aes = list(size = 10)))
      ggsave(filename = paste0(output_path,"/cell_cycle_regression/dimplot_phase_after_reg_gene",integ,".jpg"), height = 7, width = 12, plot = plot,
             quality = 50)
      
      ggsave(filename = paste0(output_path,"/cell_cycle_regression/featureplot_ribo_after_reg_gene",integ,".jpg"), height = 7, width = 12, plot = FeaturePlot(merged_obj,features = "ribo_score1"),
             quality = 50)
      
      print(plot)
      print(FeaturePlot(merged_obj,features = "ribo_score1"))
      saveRDS(merged_obj, file = paste0(output_path,"/cell_cycle_regression/regressed_out_obj_gene",integ,".rds"))
      
    } else {print("type must be isoform or gene")}
    
  } else { #integrateRCPA==TRUE, we integrate#####################################################################"
    if (type=="isoform"){
      
     
      integ="_integrated"
      
      seu_objects<-lapply(seu_objects,function(x) {
        s.genes <- as.data.frame(cc.genes$s.genes)
        colnames(s.genes)<-c("gene_name")
        g2m.genes <- as.data.frame(cc.genes$g2m.genes)
        colnames(g2m.genes)<-c("gene_name")
        
        g2m.features<-as.data.frame(Features(x))
        colnames(g2m.features)<-c("feature")
        g2m.features<-g2m.features%>%
          mutate(
            gene_name = sub(".*--", "",feature))|> 
          mutate(transcript_id=sub("--.*", "",feature))%>%
          right_join(y=g2m.genes)
        g2m.features <- g2m.features$feature
        
        s.features<-as.data.frame(Features(x))
        colnames(s.features)<-c("feature")
        s.features<-s.features%>%
          mutate(
            gene_name = sub(".*--", "",feature))|> 
          mutate(transcript_id=sub("--.*", "",feature))%>%
          right_join(y=s.genes) 
        
        
        s.features <- s.features$feature
        
        ribosomal_genes <- c(
          "RPSA", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y", "RPS5", "RPS6", "RPS7", "RPS8", 
          "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", 
          "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", 
          "RPS27", "RPS27A", "RPS28", "RPS29", 
          "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPL10", "RPL10A", 
          "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", 
          "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A", 
          "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", 
          "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39", "RPL40", "RPL41"
        )
        ribosomal_genes <- as.data.frame(ribosomal_genes)
        colnames(ribosomal_genes)<-c("gene_name")
        
        ribosomal_features<-as.data.frame(Features(x))
        colnames(ribosomal_features)<-c("feature")
        ribosomal_features<-ribosomal_features%>%
          mutate(
            gene_name = sub(".*--", "",feature))|> 
          mutate(transcript_id=sub("--.*", "",feature))%>%
          right_join(y=ribosomal_genes)#|> select(feature)
        #g2m.features<-as.vector(g2m.features)
        ribosomal_features <- ribosomal_features$feature
        
        head(ribosomal_features)
        DefaultAssay(x)<-"RNA"
        
        
        x<- NormalizeData(x)
        x<- FindVariableFeatures(x, selection.method = "vst")
        x <- ScaleData(x, features = rownames(x))
        x <- RunPCA(x, features = VariableFeatures(x))
        
        x <- CellCycleScoring(x, s.features = s.features, g2m.features = g2m.features, 
                              set.ident = TRUE)
        x <- AddModuleScore(x, features = list(ribosomal_features), name = "ribo_score")
        
        
        # view cell cycle scores and phase assignments
        head(x[[]])
        
        # and now regress out
        #seu_obj <- ScaleData(seu_obj, vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"), features = rownames(seu_obj))
        x<-  SCTransform(object =x,vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"))
        #x<- RunPCA(x)
        return(x)})
      
      merged_obj <-merge(
        seu_objects[[1]],
        y=seu_objects[2:length(seu_objects)],
        merge.data = TRUE,
        add.cell.ids = sample_names
      )
      #output_basis<-paste0(output_path,"/",type,"_",sample_names[1],"_",sample_names[2],"_",sample_names[3])
      
      merged_obj<-  SCTransform(object =merged_obj,vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"))
      merged_obj<- RunPCA(merged_obj)
      integ="_integrated"
      
      library(future)
      oopts <- options()  # <- sauvegarde toutes les options actuelles (y compris future.globals.maxSize)
      options(future.globals.maxSize = 15 * 1024^3)  # 15 GiB
      on.exit(options(oopts))  # <- on remet tout à la fin de l'exécution
      
      merged_obj <- IntegrateLayers(object = merged_obj, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                                   verbose = TRUE, normalization.method = "SCT")
        
        # now that integration is complete, rejoin layers
     # merged_obj <- JoinLayers(merged_obj)
        merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.rpca", dims = 1:15)
        merged_obj <- FindClusters(merged_obj, resolution = 1, cluster.name = "rpca_clusters")
        merged_obj <- RunUMAP(merged_obj, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")
        
        merged_obj$cluster<-Idents(merged_obj) 
      
        saveRDS(merged_obj, file = paste0(output_path,"/cell_cycle_regression/regressed_out_iso_seu_obj",integ,".rds"))
      plot <- DimPlot(merged_obj, group.by = "Phase") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
        guides(colour = guide_legend(override.aes = list(size = 10)))
      ggsave(filename = paste0(output_path,"/cell_cycle_regression/dimplot_phase_after_reg",integ,".jpg"), height = 7, width = 12, plot = plot,
             quality = 50)
      ggsave(filename = paste0(output_path,"/cell_cycle_regression/featureplot_ribo_after_reg",integ,".jpg"), height = 7, width = 12, plot = FeaturePlot(merged_obj,features = "ribo_score1"),
             quality = 50)
      
      print(plot)
      print(FeaturePlot(merged_obj,features = "ribo_score1"))
     
      
    } else if (type=="gene") {
      # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
      # segregate this list into markers of G2/M phase and markers of S phase
      s.genes <- cc.genes$s.genes
      #colnames(s.genes)<-c("gene_name")
      g2m.genes <- cc.genes$g2m.genes
      #colnames(g2m.genes)<-c("gene_name")
      ribosomal_genes <- c(
        "RPSA", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y", "RPS5", "RPS6", "RPS7", "RPS8", 
        "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", 
        "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", 
        "RPS27", "RPS27A", "RPS28", "RPS29", 
        "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPL10", "RPL10A", 
        "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", 
        "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A", 
        "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", 
        "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39", "RPL40", "RPL41"
      )
      
      
      integ="_integrated"
      
      seu_objects<-lapply(seu_objects,function(x) {
        ribosomal_features_present <- ribosomal_genes[ribosomal_genes %in% Features(x)]
        DefaultAssay(x)<-"RNA"
        x<- NormalizeData(x)
        x<- FindVariableFeatures(x, selection.method = "vst")
        x <- ScaleData(x, features = rownames(x))
        x <- RunPCA(x, features = VariableFeatures(x))
        
        x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
        x <- AddModuleScore(x, features = list(ribosomal_features_present), name = "ribo_score")
        
        # view cell cycle scores and phase assignments
        head(x[[]])
        
        # and now regress out
        #seu_obj <- ScaleData(seu_obj, vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"), features = rownames(seu_obj))
        x<-  SCTransform(object =x,vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"))
        x<- RunPCA(x)
        return(x)})
      
      merged_obj <-merge(
        seu_objects[[1]],
        y=seu_objects[2:length(seu_objects)],
        merge.data = TRUE,
        add.cell.ids = sample_names
      )
      output_basis<-paste0(output_path,"/",type,"_",sample_names[1],"_",sample_names[2],"_",sample_names[3])
      
      merged_obj<-  SCTransform(object =merged_obj,vars.to.regress = c("S.Score", "G2M.Score","ribo_score1"))
      merged_obj<- RunPCA(merged_obj)
      integ="_integrated"
      
      library(future)
      oopts <- options()  # <- sauvegarde toutes les options actuelles (y compris future.globals.maxSize)
      options(future.globals.maxSize = 15 * 1024^3)  # 15 GiB
      on.exit(options(oopts))  # <- on remet tout à la fin de l'exécution
      merged_obj <- IntegrateLayers(object = merged_obj, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                                    verbose = TRUE, normalization.method = "SCT")
      
      # now that integration is complete, rejoin layers
      #merged_obj <- JoinLayers(merged_obj)
      merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.rpca", dims = 1:15)
      merged_obj <- FindClusters(merged_obj, resolution = 1, cluster.name = "rpca_clusters")
      merged_obj <- RunUMAP(merged_obj, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")
      
      merged_obj$cluster<-Idents(merged_obj) 
      
      plot <- DimPlot(merged_obj,group.by = "Phase") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
        guides(colour = guide_legend(override.aes = list(size = 10)))
      ggsave(filename = paste0(output_path,"/cell_cycle_regression/dimplot_phase_after_reg_gene",integ,".jpg"), height = 7, width = 12, plot = plot,
             quality = 50)
      
      ggsave(filename = paste0(output_path,"/cell_cycle_regression/featureplot_ribo_after_reg_gene",integ,".jpg"), height = 7, width = 12, plot = FeaturePlot(merged_obj,features = "ribo_score1"),
             quality = 50)
      
      print(plot)
      print(FeaturePlot(merged_obj,features = "ribo_score1"))
      saveRDS(merged_obj, file = paste0(output_path,"/cell_cycle_regression/regressed_out_obj_gene",integ,".rds"))
      
    } else {print("type must be isoform or gene")}
    
    dimplot<-DimPlot(merged_obj, reduction = "umap.rpca", group.by = c("sample","rpca_clusters"),
                     combine = TRUE, label.size = 2) + labs(caption="Integrated dimplot")
    ggsave(filename = paste0(output_path,"/cell_cycle_regression/merged_integrated_RCPA_dimplot.png"), 
           plot = dimplot
    ) 
    
      }
  
  dplot<-DimPlot(merged_obj, group.by = "sample")
  print(dplot)
  ggsave(filename = paste0(output_path,"/cell_cycle_regression/dimplot_sample.png"), 
         plot = dplot
  ) 
  return(merged_obj)
}


#########################################


run_regression_and_merge <- function(sample_names, # a list of sample-names in the right order
                                     seu_objects, # a list of seurat objects, one for each sample
                                     type = "isoform", #set to 'gene' or 'isoform'
                                     output_path = paste0(rootdir, "/merge_integration"),
                                     integrateRPCA = FALSE) {
  
  dir.create(paste0(output_path, "/cell_cycle_regression"), recursive = TRUE, showWarnings = FALSE)
  
  ribosomal_genes <-  c(
    "RPSA", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y", "RPS5", "RPS6", "RPS7", "RPS8", 
    "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", 
    "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", 
    "RPS27", "RPS27A", "RPS28", "RPS29", 
    "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPL10", "RPL10A", 
    "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", 
    "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A", 
    "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", 
    "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39", "RPL40", "RPL41"
  )
  integ_suffix <- if (integrateRPCA) "_integrated" else ""
  
  message("Starting preprocessing of individual samples")
  
  seu_objects <- lapply(seu_objects, function(x) {
    DefaultAssay(x) <- "RNA"
    
    if (type == "isoform") {
      features_df <- data.frame(feature = Features(x)) %>%
        mutate(
          gene_name = sub(".*--", "", feature),
          transcript_id = sub("--.*", "", feature)
        )
      s.features <- features_df %>%
        right_join(data.frame(gene_name = cc.genes$s.genes), by = "gene_name") %>% pull(feature)
      g2m.features <- features_df %>%
        right_join(data.frame(gene_name = cc.genes$g2m.genes), by = "gene_name") %>% pull(feature)
      ribo.features <- features_df %>%
        right_join(data.frame(gene_name = ribosomal_genes), by = "gene_name") %>% pull(feature)
    } else if (type == "gene") {
      s.features <- cc.genes$s.genes
      g2m.features <- cc.genes$g2m.genes
      ribo.features <- ribosomal_genes[ribosomal_genes %in% Features(x)]
    } else {
      stop("type must be 'isoform' or 'gene'")
    }
    
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst")
    x <- ScaleData(x, features = rownames(x))
    x <- RunPCA(x, features = VariableFeatures(x))
    x <- CellCycleScoring(x, s.features = s.features, g2m.features = g2m.features, set.ident = TRUE)
    x <- AddModuleScore(x, features = list(ribo.features), name = "ribo_score")
    x <- SCTransform(x, vars.to.regress = c("S.Score", "G2M.Score", "ribo_score1"))
    
    return(x)
  })
  
  message("Merging Seurat objects")
  merged_obj <- merge(seu_objects[[1]],
                      y = seu_objects[-1],
                      add.cell.ids = sample_names,
                      merge.data = TRUE)
  
  if (integrateRPCA) {
    message("Running RPCA integration")
    merged_obj <- SCTransform(merged_obj, vars.to.regress = c("S.Score", "G2M.Score", "ribo_score1"))
    merged_obj <- RunPCA(merged_obj)
    
    library(future)
    oopts <- options()
    options(future.globals.maxSize = 15 * 1024^3)
    on.exit(options(oopts))
    
    merged_obj <- IntegrateLayers(merged_obj,
                                  method = RPCAIntegration,
                                  orig.reduction = "pca",
                                  new.reduction = "integrated.rpca",
                                  normalization.method = "SCT",
                                  verbose = TRUE)
    
    merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.rpca", dims = 1:15)
    merged_obj <- FindClusters(merged_obj, resolution = 0.5, cluster.name = "rpca_clusters")
    merged_obj <- RunUMAP(merged_obj, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")
  } else {
    message("Running merged object's processing")
    merged_obj <- SCTransform(merged_obj, vars.to.regress = c("S.Score", "G2M.Score", "ribo_score1"))
    merged_obj <- RunPCA(merged_obj, assay = "SCT")
    merged_obj <- FindNeighbors(merged_obj, dims = 1:15)
    merged_obj <- FindClusters(merged_obj,resolution = 0.5, cluster.name = "clusters")
    merged_obj <- RunUMAP(merged_obj, reduction = "pca", dims = 1:15, reduction.name = "umap")
  }
  merged_obj<-PrepSCTFindMarkers(merged_obj, assay = "SCT", verbose = TRUE)
  saveRDS(merged_obj, file = paste0(output_path, "/cell_cycle_regression/regressed_obj_", type, integ_suffix, ".rds"))
  message(integ_suffix," object saved in:",output_path, "/cell_cycle_regression/regressed_obj_", type, integ_suffix, ".rds")
  
  # ----- Visualization -----
  message("Plotting results")
  plot_phase <- DimPlot(merged_obj, group.by = "Phase") +
    labs(title = paste("UMAP - Cell Cycle Phase after regression", type, integ_suffix)) +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
    guides(colour = guide_legend(override.aes = list(size = 10)))
  
  ggsave(paste0(output_path, "/cell_cycle_regression/dimplot_phase_after_reg", type, integ_suffix, ".jpg"),
         plot = plot_phase, height = 7, width = 12)
  
  plot_ribo <- FeaturePlot(merged_obj, features = "ribo_score1") +
    labs(title = paste("Ribosomal Score after regression", type, integ_suffix))
  
  ggsave(paste0(output_path, "/cell_cycle_regression/featureplot_ribo_after_reg", type, integ_suffix, ".jpg"),
         plot = plot_ribo, height = 7, width = 12)
  
  print(plot_phase)
  print(plot_ribo)
  

  # Optional: integrated dimplot
  if (integrateRPCA) {
    dimplot <- DimPlot(merged_obj, reduction = "umap.rpca", group.by = c("sample", "rpca_clusters"),
                       combine = TRUE, label.size = 2) +
      labs(title = paste("Integrated UMAP - RPCA", type) )+
      theme_minimal()
    
    ggsave(paste0(output_path, "/cell_cycle_regression/merged_integrated_", type,"_RCPA_dimplot.png"), plot = dimplot)
    print(dimplot)
    } else {
      dplot <- DimPlot(merged_obj, group.by = "clusters") +
        labs(title = paste("UMAP colored by clusters",integ_suffix, type))
      ggsave(paste0(output_path, "/cell_cycle_regression/dimplot_clusters_", type,integ_suffix,".png"), plot = dplot)
      print(dplot)
      
    }
 
  # Sample-level plot
  dplot <- DimPlot(merged_obj, group.by = "sample") +
    labs(title = paste("UMAP colored by sample",integ_suffix,type))
  ggsave(paste0(output_path, "/cell_cycle_regression/dimplot_sample_", type,integ_suffix,".png"), plot = dplot)
  print(dplot)
  
  
  return(merged_obj)
}