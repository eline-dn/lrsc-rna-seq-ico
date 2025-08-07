
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
#library(glmGamPoi)


setwd(rootdir)  # Set this to correct location

DE_samples<-function(seu_obj,
                     output_dir=paste0(rootdir,"/analyses"),
                     sample_name="", #to specify what we are working with in the output file name 
                     topn=10,
                     type="isoform",
                     PrepSCT=FALSE){
  
  DefaultAssay(object = seu_obj) <- "SCT"
  
  if (PrepSCT==TRUE) {  seu_obj<-PrepSCTFindMarkers(seu_obj, assay = "SCT", verbose = TRUE)
}
  
  res=list()
  
  Idents(seu_obj)<-"sample"
  
  all_markers <- FindAllMarkers(seu_obj, do.print = FALSE,
                                assay="SCT", slot="data",
                                # logfc.threshold = 0.5, min.pct = 0.20,
                                only.pos = FALSE) %>% dplyr::filter(p_val_adj < 0.05)
  res$markers<-all_markers
  #write.csv(all_markers, file=paste0(output_dir,"/all_markers_",sample_name,type,".csv"))
  
  head(all_markers)
  
  
  
  #then heatmap 
  all_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = topn) %>%
    ungroup() -> top_markers
  
  heatmap<-DoHeatmap(seu_obj, slot="data", 
                     assay="SCT",
                     features = top_markers$gene,group.by = "sample") + labs(caption = paste0("Top markers for each sample at ",type," level"))
  
 ggsave(filename=paste0(output_dir,"/heatmap_top_markers_",sample_name,type,".png"),
         heatmap)
  #print(heatmap) 
  res$heatmap<-heatmap
  return(res)
}



DE_clusters<-function(seu_obj, # an integrated seurat object
                      #sample_names, # a list of the sample names
                      outdir=paste0(rootdir,"/analyses/DE_clusters"),
                      topn=10,
                      type="isoform"
){
  
  # split our object into a list of 3 objects, one for each sample
  so_list = SplitObject(seu_obj,split.by = "sample")
  names(so_list) <- sapply(so_list, function(x) {return(unique(x$sample))})
  
  #run DE on each sample to find specific markers for each cluster in the three samples
  markers_list <- lapply(so_list, function(x) {
    x$cluster<-Idents(x)
    all_markers<- FindAllMarkers(x, do.print = FALSE,
                                 group.by="rpca_clusters", 
                                 assay="SCT", slot="data",
                                 logfc.threshold = 0.5, min.pct = 0.20, only.pos = FALSE) %>% 
      dplyr::filter(p_val_adj < 0.05)
    
    head(all_markers)
    return(all_markers)
  })
  
  names(markers_list) = names(so_list)
  
  
  # one heatmap will be generated for each sample, featuring the top markers for each cluster in the given sample
  heatmap_list<-lapply(names(markers_list), function(x) { # x is the sample's name
    
    #plot heatmaps
    top_markers<-markers_list[[x]] %>%
      group_by(cluster) %>%
      dplyr::filter(avg_log2FC > 1) %>%
      slice_head(n = topn) %>%
      ungroup() 
    
    heatmap<-DoHeatmap(so_list[[x]],
                       features = top_markers$gene,
                       assay = "SCT",
                       slot="data") + ggtitle(paste0("Top markers for each cluster ",x," at ",type," level") )
    
    ggsave(filename=paste0(outdir,"/heatmap_",x,"_",type,".png"),heatmap)
    print(heatmap)
    return(heatmap)
  })
  names(heatmap_list)=names(markers_list)
  res=list()
  res$markers=markers_list
  res$heatmaps=heatmap_list
  
  return(res)
}




extract_cluster_feature<-function(num_cluster,
                                  marker_table_list, # the $markers output from DE_cluster()
                                  type="isoform",
                                  cluster_names="rpca_clusters", # or "clusters" if the object is not integrated
                                  integrated=TRUE,
                                  seu_obj # used in DE_cluster()
) {
  integ_suffix <- if (integrated) "_integrated" else ""
  
  listx<-lapply(marker_table_list,function(x) {
    markers<-x%>% dplyr::filter(cluster == num_cluster) %>%
      dplyr::filter(avg_log2FC > 1 & p_val_adj<0.05) %>%
      slice_head(n = 20)|>pull(gene)
    return(markers) })
  list<- reduce(listx, union)
  
  
  dir.create(paste0(rootdir,"/analyses/cluster",num_cluster,"_",type))
  setwd(paste0(rootdir,"/analyses/cluster",num_cluster,"_",type))
  
  for (i in seq(1, length(list), by = 6)) {
    vplot<-VlnPlot(seu_obj, features = list[i:min(i+5, length(list))], group.by = cluster_names)+labs(caption = paste0("Cluster ",num_cluster, " ",type,"s"))
    print(vplot)
    ggsave(filename = paste0("vlnplot_",i,"_",min(i+5, length(list)),type,integ_suffix,".png"),vplot)
  } 
  return(list)
  
}
