# here are a few visualisation tolls based mostyl on seurat's, as well as pathway analyses functions 
#running with gprofiler2 or EnrichR

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
library(enrichR)
library(ggVennDiagram)




plot_venn<-function(seu_obj,
                    type="isoform",
                    output=rootdir,
                    sample="") {
  so_list = SplitObject(seu_obj, split.by = "sample")
  names(so_list) <- sapply(so_list, function(x) {return(unique(x$sample))})
  
  feat_list<-lapply(so_list, function(x){
    
    return(rownames(GetAssayData(x, assay = "SCT", layer = "counts"))[
      Matrix::rowSums(GetAssayData(x, assay = "SCT", layer = "counts") > 0) > 0 ])
  })
  
  vd<-ggVennDiagram(feat_list) +labs(title = paste0(sample," : Overlapping ", type))
  ggsave(filename = paste0(output,"/",sample,type,".png"),vd)
  print(vd)
  return (vd)
}
#Here is a useful function for plotting feature plots of every isoform associated with one specific gene:

plot_feature_iso<-function(seu_obj,
                           gene,
                           sample_name,
                 
                           output_dir="./analyses") {
  
  features <- Features(seu_obj)
  
  plot_features_list <- grep(paste0("--",gene,"$"), features, value = TRUE)
  
  seu_obj$sample <- factor(seu_obj$sample , levels = c("ctrl", "S24", "R"))
  plot1<-FeaturePlot(seu_obj, features = plot_features_list,combine=TRUE,split.by="sample",label=TRUE) +theme(strip.text = element_text(size = 1)) #+labs(title = paste0("Differential transcript expression per sample for ",gene))
  
  ggsave(filename = paste0(output_dir,"/featureplot",paste0(sample_name,"_",gene,".png")), 
         plot = plot1, width=8,height=5*length(plot_features_list),limitsize = FALSE)
  
  
  
  print(paste0("plotted in ",paste0(output_dir,"/featureplot",paste0(sample_name,"_",gene,".png"))))
  return(plot1)
}

#The second function allows us to plot only the transcripts that are in the marker list specified in input, for example a list of markers that are significantly differentially expressed. 
#NB: the markers in this list need to be written like this: `transcript_ensembl_id--gene_symbol`

plot_de_feature_iso<-function(seu_obj,
                              gene,
                              sample_name,
                              output_dir="./analyses",
                              marker_list) {
  features <- Features(seu_obj)
  
  plot_features_list <- grep(paste0("--",gene,"$"),( features %in% marker_list), value = TRUE) 
  
  seu_obj$sample <- factor(seu_obj$sample , levels = c("ctrl", "S24", "R"))
  plot1<-FeaturePlot(seu_obj, features = plot_features_list,combine=TRUE,split.by="sample") +theme(strip.text = element_text(size = 1)) #+labs(title = paste0("Differential transcript expression per sample for ",gene))
  
  ggsave(filename = paste0(output_dir,"/featureplot_select",paste0(sample_name,"_",gene,".png")), 
         plot = plot1, width=8,height=5*length(plot_features_list),limitsize = FALSE)
  
  
  
  print(paste0("plotted in ",paste0(output_dir,"/featureplot",paste0(sample_name,"_",gene,".png"))))
  return(plot1)
}



# pathway enrichment analysis using gprofiler2'sgost

perform_pathway_enrichment<-function(seu_obj,
                                     sample_name, # the sample you want to focus on
                                     background=FALSE,
                                     type="isoform" #set this to "isoform" or "gene"
){
  
  DefaultAssay(object = seu_obj) <- "SCT"
  #seu_obj<-NormalizeData(seu_obj)
  #seu_obj <- FindVariableFeatures(object = seu_obj)
  #seu_obj <- ScaleData(object = seu_obj)
  
  Idents(seu_obj) = "sample"
  DE_all = FindAllMarkers(seu_obj,assay="SCT"#,slot="scale.data"
  )
  all_markers <- DE_all %>% dplyr::filter(p_val_adj < 0.05)
  
  if (type=="isoform") {
    #Filter for significant genes in the current cluster
    sig_transcripts <- all_markers %>%
      dplyr::filter(cluster==sample_name & p_val_adj < 0.05) %>%
      pull(gene)  %>% # Extract transcript names
      str_extract("--[[:alnum:].-]+$") #extract gene names
    sig_transcripts <- sub("^--", "", sig_transcripts)
    
    if (background==TRUE){
      background_genes <- rownames(GetAssayData(seu_obj, assay = "RNA", layer = "counts"))[
        Matrix::rowSums(GetAssayData(seu_obj, assay = "RNA", layer = "counts") > 0) > 0
      ]  %>% str_extract("--[[:alnum:].-]+$")
      background_genes <- sub("^--", "", background_genes) }
    
  } else if (type=="gene") {
    #Filter for significant genes in the current cluster
    sig_transcripts <- all_markers %>%
      dplyr::filter(cluster==sample_name & p_val_adj < 0.05) %>%
      pull(gene)  
    
    if (background==TRUE){
      background_genes <- rownames(GetAssayData(seu_obj, assay = "RNA", layer = "counts"))[
        Matrix::rowSums(GetAssayData(seu_obj, assay = "RNA", layer = "counts") > 0) > 0
      ] }
  } else { 
    print("type must be isoform ar gene")}
  
  # Run g:Profiler for pathway enrichment analysis
  
  if (background==TRUE) {
    pathway_results <- gprofiler2::gost(
      query = sig_transcripts,
      ordered_query = TRUE,
      correction_method = 'fdr',
      custom_bg = background_genes,
      sources = c("GO", "KEGG", "REACTOME"),
      evcodes = TRUE,
      #significant=FALSE,
      organism="hsapiens")
  }
  
  else {
    pathway_results <- gprofiler2::gost(
      query = sig_transcripts,
      ordered_query = TRUE,
      correction_method = 'fdr',
      #custom_bg = background_genes,
      sources = c("GO", "KEGG", "REACTOME"),
      evcodes = TRUE,
      #significant=FALSE,
      organism="hsapiens")
  }
  
  # Prepare the data for plotting
  df_path <- as_tibble(pathway_results$result) %>%
    filter(term_size < 3000, term_size > 5) %>%
    filter(!term_id %in% unlist(pathway_results$parents))
  
  # Plot top 5 results per database
  plot1<-df_path %>%
    group_by(source) %>%
    slice_min(p_value, n = 5, with_ties = TRUE) %>%
    ungroup() %>%
    ggplot(aes(x = reorder(term_name, -p_value), y = -log10(p_value), fill = source)) +
    geom_bar(stat = 'identity', position = position_identity()) +
    coord_flip() +
    theme_bw() +
    labs(x = "") +
    facet_grid(source ~ ., space = 'free', scales = 'free') +
    theme(legend.position = 'none',
          axis.text.y = element_text(angle = 0, size = 8)) +  # Rotate and adjust y-axis text size
    ggtitle(paste0("Pathway Enrichment for ",sample_name," ",type)) 
  
  ggsave(filename=paste0(rootdir,"/analyses/pathway_enrichment_analysis_",type,"_",sample_name,".png"),
         plot=plot1
  )
  res<-list()
  res$pathway_results<-pathway_results
  res$plot<-plot1
  
  plot1 
  return(res)
}




# function to look up isoforms of genes associated with a specific GO term
select_GO_term<-function(pathway_results,
                         GO_term){
  gene_GO_associations<-pathway_results$result |>
    filter(p_value<0.05) |>
    select(intersection,term_name)
  
  
  gene_term_associations<-gene_GO_associations |>
    filter(term_name==GO_term) |>
    select(intersection)
  
  
  features <- all_markers_isoforms_sample |>group_by(cluster) |> filter(avg_log2FC>2 & p_val_adj<0.05) |>pull(gene)
  #slice_max(n=15, order_by=avg_log2FC) 
  
  # Liste des gènes d’intérêt
  gene_list <- str_split(gene_term_associations$intersection,pattern=",")
  
  # Créer un tibble avec les features et extraire le nom du gène
  features_tbl <- tibble(feature = features) %>%
    mutate(
      gene_name = str_extract(feature, "[^--]+$")  # extrait ce qui suit le dernier tiret
    )
  
  # Filtrer les lignes dont le gène est dans la liste
  filtered_features <- features_tbl %>%
    filter(gene_name %in% unlist(gene_list))
  return(filtered_features$gene_name)
}


# pathway analysis with enrichR
# input: a list of genes to test

run_enrichR<-function(input,
                      sample,
                      output_path=paste0(rootdir,"/analyses")) {

websiteLive <- getOption("enrichR.live")

if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
}
#select databases
dbs <- c("GO_Molecular_Function_2025", "GO_Cellular_Component_2025", 
         "GO_Biological_Process_2025","Reactome_Pathways_2024")




#query

if (websiteLive) {
    enriched <- enrichr(input, dbs)
}

if (websiteLive) {
  head(enriched[["GO_Biological_Process_2025"]])
}
# GO_Biological_Process_2025
if (websiteLive) {
  title <- paste0("GO_Biological_Process, ", sample)
  plot <- plotEnrich(enriched[["GO_Biological_Process_2025"]], showTerms = 20, numChar = 60,
                     y = "Count", orderBy = "P.value", title = title)
  print(plot)
  filename <- file.path(output_path, paste0(gsub("[ ,]", "_", title), ".png"))
  ggsave(filename, plot, width = 10, height = 8)
}

# GO_Molecular_Function_2025
if (websiteLive) {
  title <- paste0("GO_Molecular_Function, ", sample)
  plot <- plotEnrich(enriched[["GO_Molecular_Function_2025"]], showTerms = 20, numChar = 60,
                     y = "Count", orderBy = "P.value", title = title)
  print(plot)
  filename <- file.path(output_path, paste0(gsub("[ ,]", "_", title), ".png"))
  ggsave(filename, plot, width = 10, height = 8)
}

# GO_Cellular_Component_2025
if (websiteLive) {
  title <- paste0("GO_Cellular_Component, ", sample)
  plot <- plotEnrich(enriched[["GO_Cellular_Component_2025"]], showTerms = 20, numChar = 60,
                     y = "Count", orderBy = "P.value", title = title)
  print(plot)
  filename <- file.path(output_path, paste0(gsub("[ ,]", "_", title), ".png"))
  ggsave(filename, plot, width = 10, height = 8)
}

# Reactome_Pathways_2024
if (websiteLive) {
  title <- paste0("Reactome_Pathways, ", sample)
  plot <- plotEnrich(enriched[["Reactome_Pathways_2024"]], showTerms = 20, numChar = 60,
                     y = "Count", orderBy = "P.value", title = title)
  print(plot)
  filename <- file.path(output_path, paste0(gsub("[ ,]", "_", title), ".png"))
  ggsave(filename, plot, width = 10, height = 8)
}
return(enriched)




}





# pathway annotation heatmap

enrich_heatmap<-function(enrich_df,
                         sample_name,
                         databases=c("Reactome_Pathways_2024","GO_Biological_Process_2025"),
                         rank=10,
                         show_labels=TRUE) {
  
  
  #- création d'une matrice avec en ligne le numéro de cluster, en colonne le go term et dans les cellules le c-score
  enrich_mat<-enrich_df%>%
    filter((Adjusted.P.value <0.05) & (dbs %in% databases) & (Rank %in% 1:rank)) %>%
    select(Term, cluster, Combined.Score) %>%
    pivot_wider(names_from = Term, values_from = Combined.Score)
  
  rownames(enrich_mat)<-enrich_mat$cluster
  # pour gérer les clusters qui n'existent pas et empiler les heatmaps après
  
  missing_clust<-setdiff(0:max(unique(enrich_df$cluster)),enrich_mat$cluster)
  
  # Créer un data.frame vide avec mêmes colonnes et les bons rownames
  na_rows <- as.data.frame(matrix(NA, nrow = length(missing_clust), ncol = ncol(enrich_mat)))
  colnames(na_rows) <- colnames(enrich_mat)
  rownames(na_rows) <- missing_clust
  
  # Fusion
  enrich_mat <- rbind(enrich_mat, na_rows)
  na_idx <- is.na(enrich_mat$cluster)
  enrich_mat$cluster[na_idx] <- as.numeric(rownames(enrich_mat)[na_idx])

  
  row_labels<-enrich_mat$cluster
  
  
  
  enrich_mat <- as.matrix(enrich_mat[,-1])  # enlève la première colonne "cluster"
  enrich_mat <- apply(enrich_mat, 2, as.numeric)  # convertit chaque colonne en numérique
  rownames(enrich_mat) <- row_labels
  enrich_mat<-enrich_mat[order(as.numeric(rownames(enrich_mat))), ]
  
  row_labels<-rownames(enrich_mat)
  #- heatmap avec le combined score
  
  min=min(enrich_df$Combined.Score)
  max=max(enrich_df$Combined.Score)
  median=median(enrich_df$Combined.Score)
  median=20
  min=0
  max=2000
  col_fun = colorRamp2(c(min, median, max), c("green", "white", "red"))
  #col_fun(seq(-3, 3))
  

  column_labels <- if (show_labels) colnames(enrich_mat) else rep(".", length(colnames(enrich_mat)))
  #the transposition of the previous matrix is used for plotting the heatmap for more readability
  ht<-Heatmap(t(enrich_mat), 
              name = sample_name, 
              col = col_fun, 
              na_col = "lightblue",
              column_labels = row_labels,
              row_labels = column_labels,
              cluster_rows = F,
              cluster_columns = F,
              column_names_rot = 0,
              column_names_side = "top",
              column_title =sample_name,
              left_annotation = rowAnnotation(sample = anno_block(gp = gpar(fill = sample(x = 2:5, size = 1)),
                                                                  labels = sample_name, 
                                                                  labels_gp = gpar(col = "black", fontsize = 10))),
              row_km = 1)
  #draw(ht)
  return(ht)
  
}


apopt_heatmap<-function(enrich_df,
                        sample_name,
                        show_labels=TRUE) {

#- création d'une matrice avec en ligne le numéro de cluster, en colonne le go term et dans les cellules le c-score
  enrich_mat<-enrich_df%>%
    filter((Adjusted.P.value <0.05) & (dbs %in% c("Reactome_Pathways_2024","GO_Biological_Process_2025")) ) %>%
    select(Term, cluster, Combined.Score) %>%
    pivot_wider(names_from = Term, values_from = Combined.Score)
  
  rownames(enrich_mat)<-enrich_mat$cluster
  # pour gérer les clusters qui n'existent pas et empiler les heatmaps après
  
  missing_clust<-setdiff(0:max(unique(enrich_df$cluster)),enrich_mat$cluster)
  
  # Créer un data.frame vide avec mêmes colonnes et les bons rownames
  na_rows <- as.data.frame(matrix(NA, nrow = length(missing_clust), ncol = ncol(enrich_mat)))
  colnames(na_rows) <- colnames(enrich_mat)
  rownames(na_rows) <- missing_clust
  
  # Fusion
  enrich_mat <- rbind(enrich_mat, na_rows)
  na_idx <- is.na(enrich_mat$cluster)
  enrich_mat$cluster[na_idx] <- as.numeric(rownames(enrich_mat)[na_idx])
  
  
  row_labels<-enrich_mat$cluster
  
  
  enrich_mat <- as.matrix(enrich_mat[,-1])  # enlève la première colonne "cluster"
  enrich_mat <- apply(enrich_mat, 2, as.numeric)  
  pattern <- "mitocho|death|apopt|stress|ferrop|p53|bcl|caspas|tnf|trai|fas|bak|bid|escape|evas|senescen|autophag|pyropt|necros|cytotox|granzym|perforin|survival|prolif|oxidative"
  enrich_mat <- as.matrix(enrich_mat[,grep(pattern, colnames(enrich_mat), value = TRUE, ignore.case = TRUE)])  
  enrich_mat[is.na(enrich_mat)] <- -1000
  rownames(enrich_mat) <- row_labels
  
  enrich_mat<-enrich_mat[order(as.numeric(rownames(enrich_mat))), ]
  row_labels<-rownames(enrich_mat)
 
  
#- heatmap avec le combined score

min=min(enrich_df$Combined.Score)
max=max(enrich_df$Combined.Score)
median=median(enrich_df$Combined.Score)
median=20
min=0
max=2000
column_labels <- if (show_labels) colnames(enrich_mat) else rep(".", length(colnames(enrich_mat)))

col_fun = colorRamp2(c(min, median, max), c("lightblue", "white", "red"))

#the transposition of the previous matrix is used for plotting the heatmap for more readability

ht<-Heatmap(t(enrich_mat), 
            name = sample_name, 
            col = col_fun, 
            na_col = "lightblue",
            column_labels = row_labels,
            row_labels = column_labels,
            cluster_rows = F,
            cluster_columns = F,
            column_names_rot = 0,
            column_names_side = "top",
            column_title =sample_name,
            left_annotation = rowAnnotation(sample = anno_block(gp = gpar(fill = sample(x = 2:5, size = 1)),
                                                             labels = sample_name, 
                                                             labels_gp = gpar(col = "black", fontsize = 10))),
            row_km = 1)
#draw(ht)

return(ht)

}




# pathway annotation heatmap with dendrogramm and clustering

enrich_heatmap_dendro<-function(enrich_df,
                         sample_name,
                         databases=c("Reactome_Pathways_2024","GO_Biological_Process_2025"),
                         rank=10,
                         show_labels=TRUE,
                         dendro_col=F,
                         dendro_rows=F) {
  
  
  #- création d'une matrice avec en ligne le numéro de cluster, en colonne le go term et dans les cellules le c-score
  enrich_mat<-enrich_df%>%
    filter((Adjusted.P.value <0.05) & (dbs %in% databases) & (Rank %in% 1:rank)) %>%
    select(Term, Sample_name, Combined.Score) %>%
    pivot_wider(names_from = Term, values_from = Combined.Score)
  
  rownames(enrich_mat)<-enrich_mat$Sample_name
  row_labels<-enrich_mat$Sample_name
  enrich_mat <- as.matrix(enrich_mat[,-1])  # enlève la première colonne "Sample_name"
  enrich_mat <- apply(enrich_mat, 2, as.numeric)  # convertit chaque colonne en numérique
  enrich_mat[is.na(enrich_mat)] <- -1000
  enrich_mat[is.infinite(enrich_mat)] <- 500000
  rownames(enrich_mat) <- row_labels
  enrich_mat<-enrich_mat[order(as.numeric(rownames(enrich_mat))), ]
  row_labels<-rownames(enrich_mat)
  
  #- heatmap avec le combined score
  
  min=min(enrich_df$Combined.Score)
  max=max(enrich_df$Combined.Score)
  median=median(enrich_df$Combined.Score)
  median=20
  min=0
  max=2000
  col_fun = colorRamp2(c(min, median, max), c("lightblue", "white", "red"))
  #col_fun(seq(-3, 3))
  
  
  column_labels <- if (show_labels) colnames(enrich_mat) else rep(".", length(colnames(enrich_mat)))
  
  # order the Sample_name labels:
  vec <- colnames(t(enrich_mat))
  x <- as.integer(sub("_.*", "", vec))
  group <- sub(".*_", "", vec)
  group_levels <- c("ctrl", "S24", "R")
  ord <- order(x, match(group, group_levels))
  column_order <- vec[ord]
  
  # split per cluster

  column_split = as.numeric(sub("_.*", "", column_order))
  
  #the transposition of the previous matrix is used for plotting the heatmap for more readability
  ht<-Heatmap(t(enrich_mat), 
              name = sample_name, 
              col = col_fun, 
              na_col = "lightblue",
              column_labels = row_labels,
              row_labels = column_labels,
              cluster_rows = dendro_rows,
              cluster_columns = dendro_col,
              column_names_rot = 0,
              column_names_side = "top",
              column_title =sample_name,
              column_order = column_order,
              #column_split = factor(column_split, levels = unique(column_split)),
              column_dend_height = unit(4, "cm"),
              left_annotation = rowAnnotation(sample = anno_block(gp = gpar(fill = sample(x = 2:5, size = 1)),
                                                                  labels = sample_name, 
                                                                  labels_gp = gpar(col = "black", fontsize = 10))),
              row_km = 1)
  #draw(ht)
  return(ht)
  
}




apopt_heatmap_dendro<-function(enrich_df,
                        sample_name,
                        show_labels=TRUE) {
  
  #- création d'une matrice avec en ligne le numéro de cluster, en colonne le go term et dans les cellules le c-score
  enrich_mat<-enrich_df%>%
    filter((Adjusted.P.value <0.05) & (dbs %in% c("Reactome_Pathways_2024","GO_Biological_Process_2025")) ) %>%
    select(Term, Sample_name, Combined.Score) %>%
    pivot_wider(names_from = Term, values_from = Combined.Score)
  
  rownames(enrich_mat)<-enrich_mat$Sample_name
  row_labels<-enrich_mat$Sample_name
  enrich_mat <- as.matrix(enrich_mat[,-1])  # enlève la première colonne "Sample_name"
  enrich_mat <- apply(enrich_mat, 2, as.numeric)  # convertit chaque colonne en numérique
  enrich_mat[is.na(enrich_mat)] <- -1000
  enrich_mat[is.infinite(enrich_mat)] <- 500000
  rownames(enrich_mat) <- row_labels
  enrich_mat<-enrich_mat[order(as.numeric(rownames(enrich_mat))), ]
  row_labels<-rownames(enrich_mat) 
  pattern <- "mitocho|death|apopt|stress|ferrop|p53|bcl|caspas|tnf|trai|fas|bak|bid|escape|evas|senescen|autophag|pyropt|necros|cytotox|granzym|perforin|survival|prolif|oxidative"
  enrich_mat <- as.matrix(enrich_mat[,grep(pattern, colnames(enrich_mat), value = TRUE, ignore.case = TRUE)])  
  rownames(enrich_mat) <- row_labels
  
  
  #- heatmap avec le combined score
  
  min=min(enrich_df$Combined.Score)
  max=max(enrich_df$Combined.Score)
  median=median(enrich_df$Combined.Score)
  median=20
  min=0
  max=2000
  column_labels <- if (show_labels) colnames(enrich_mat) else rep(".", length(colnames(enrich_mat)))
  
  col_fun = colorRamp2(c(min, median, max), c("lightblue", "white", "red"))
  
  
  # order the Sample_name labels:
  vec <- colnames(t(enrich_mat))
  x <- as.integer(sub("_.*", "", vec))
  group <- sub(".*_", "", vec)
  group_levels <- c("ctrl", "S24", "R")
  ord <- order(x, match(group, group_levels))
  column_order <- vec[ord]
  
  
  #the transposition of the previous matrix is used for plotting the heatmap for more readability
  
  ht<-Heatmap(t(enrich_mat), 
              name = sample_name, 
              col = col_fun, 
              na_col = "lightblue",
              column_labels = row_labels,
              row_labels = column_labels,
              cluster_rows = F,
              cluster_columns = F,
              column_names_rot = 0,
              column_names_side = "top",
              column_title =sample_name,
              column_order = column_order,
              left_annotation = rowAnnotation(sample = anno_block(gp = gpar(fill = sample(x = 2:5, size = 1)),
                                                                  labels = sample_name, 
                                                                  labels_gp = gpar(col = "black", fontsize = 10))),
              row_km = 1)
  #draw(ht)
  
  return(ht)
  
}

prep_enrich_df_glob<-function(seu_obj,
                         overexpressed=TRUE,
                         name="whole_set") {
  
  
  Idents(seu_obj)<-"clusters"
  gene_mark<-FindAllMarkers(seu_obj, do.print = FALSE,
                            assay="SCT", slot="data",
                            logfc.threshold = 0.5, min.pct = 0.20,
                            only.pos = FALSE) %>% 
    dplyr::filter(p_val_adj < 0.05) 
  
  # filter over expressed or under expressed genes:
  
  list<-lapply(unique(seu_obj$clusters), function(x){
    markers<-gene_mark %>% 
      filter(cluster==x)%>% 
      filter(case_when(
        overexpressed ~ avg_log2FC > 0,
        !overexpressed ~ avg_log2FC < 0)) %>%
      pull(gene)
    return(markers)})
  names(list)<-unique(seu_obj$clusters)

  
  #- faire enrichr avec background sur la list
  # pour chaque enriched, filtre sur nombre de gènes >=1)
  
  bckg<-Features(seu_obj)
  
  enriched_list<-lapply(list, function(x) {
    input=x
    websiteLive <- getOption("enrichR.live")
    
    if (websiteLive) {
      listEnrichrSites()
      setEnrichrSite("Enrichr") # Human genes
    }
    #select databases
    dbs <- c("GO_Molecular_Function_2025", "GO_Cellular_Component_2025", 
             "GO_Biological_Process_2025","Reactome_Pathways_2024")
    bckg<-Features(merged_gene)
    
    #query
    
    if (websiteLive) {
      
      enriched_bg<-lapply(dbs, function(dbase){
        tryCatch(
          {
            enriched_dbase <- enrichr(input, dbase, background=bckg)
            return(enriched_dbase[[1]])
          }, error = function(e){
            enriched_dbase<-(data.frame(
              Term = character(),
              Rank=character(),
              P.value = double(),
              Adjusted.P.value = double(),
              Old.P.value = double(),
              Old.Adjusted.P.value = double(),
              Odds.Ratio= double(),
              Combined.Score= double(),
              Genes = character(),
              stringsAsFactors = FALSE))
            return(enriched_dbase)
          })
      })
      names(enriched_bg)<-dbs
    }
    #regroupement de toutes les databases en un df et filtre sur le nombre de gènes associés avec gene_number_treshold
    gene_number_treshold=1
    
    enriched_bg<-lapply(names(enriched_bg), function(x){
      return(enriched_bg[[x]]%>% mutate(dbs=x))
    })
    
    enriched_bg<-bind_rows(enriched_bg)%>% 
      mutate(gene_number = str_count(Genes, ";")+1 ) %>%
      filter(gene_number>=gene_number_treshold)
    
    return(enriched_bg)
    
  })
  names(enriched_list)<-names(list)
  
  #- compilation des enrichissements par clusters dans le même dataframe avec annotation par cluster
  
 
  enriched_list<-lapply(names(enriched_list), function(x) {
    df<-enriched_list[[x]]
    df<-df %>% mutate(cluster=x)
    return(df)
  })
  names(enriched_list)<-names(list)
  enrich_df<-bind_rows(enriched_list)
  if (overexpressed) {
    express_suffix="_overexpressed_genes"
  } else { express_suffix="_underexpressed_genes"}
  write_csv(enrich_df, file = paste0(rootdir,"/unsupervised/enrichment_per_cluster_",name,express_suffix,".csv"))
  return(enrich_df)
}

##########

#to prep enrich_df en splited object



prep_enrich_df_split<-function(seu_obj,
                     overexpressed=TRUE,
                     name="whole_set") {
  

sub_list<-SplitObject(seu_obj, split.by = "sample")
bckg<-Features(seu_obj)

enrich_df_list<-lapply(names(sub_list), function(sample_name) {
 # sample_name="S24" #(tests)
  
  seu_obj<-sub_list[[sample_name]]
  Idents(seu_obj)<-"clusters"
  gene_mark<-FindAllMarkers(seu_obj, do.print = FALSE,
                            assay="SCT", slot="data",
                            logfc.threshold = 0.5, min.pct = 0.10,
                            only.pos = FALSE) %>% 
    dplyr::filter(p_val_adj < 0.05) 
  
  #filter negative or positive markers in each cluster
  list<-lapply(unique(seu_obj$clusters), function(x){
    markers<-gene_mark %>% 
      filter(cluster==x)%>% 
      filter(case_when(
        overexpressed ~ avg_log2FC > 0,
        !overexpressed ~ avg_log2FC < 0)) %>%
      pull(gene)
    return(markers)})
  names(list)<-unique(seu_obj$clusters)
  
  #enrichment
  enriched_list<-lapply(list, function(x) {
    #x<-list[["1"]]#(tests)
    
    # si il y a suffisamment de genes DE dans un cluster:
    #an error happens in enrichr if there is only one gene or 0
    if (length(x)>1){
      input=x
      #websiteLive <- TRUE 
      websiteLive <-getOption("enrichR.live")
      if (websiteLive==FALSE) {
        message("bug enrichr")
        
      }
      
      if (websiteLive) {
        listEnrichrSites()
        setEnrichrSite("Enrichr") # Human genes
      }
      
      #select databases
      dbs <- c("GO_Molecular_Function_2025", "GO_Cellular_Component_2025", 
               "GO_Biological_Process_2025","Reactome_Pathways_2024")
      
      
      #query
      
      if (websiteLive) {
        
        enriched_bg<-lapply(dbs, function(dbase){
          tryCatch(
            {
          enriched_dbase <- enrichr(input, dbase, background=bckg)
          return(enriched_dbase[[1]])
          }, error = function(e){
            enriched_dbase<-(data.frame(
              Term = character(),
              Rank=character(),
              P.value = double(),
              Adjusted.P.value = double(),
              Old.P.value = double(),
              Old.Adjusted.P.value = double(),
              Odds.Ratio= double(),
              Combined.Score= double(),
              Genes = character(),
              stringsAsFactors = FALSE))
            return(enriched_dbase)
          })
        })
        names(enriched_bg)<-dbs
      }
      # merged result from different databases together and filter on gene number
      enriched_bg<-lapply(names(enriched_bg), function(x){
        return(enriched_bg[[x]]%>% mutate(dbs=x))})
      
      gene_number_treshold=1
      enriched_bg<-bind_rows(enriched_bg)%>% 
        mutate(gene_number = str_count(Genes, ";")+1 ) %>%
        filter(gene_number>=gene_number_treshold)
      
     
      
    } else { # si pas assez de gènes DE pour faire un enrichissment, on renvoit un enriched_bg vide. cela permet de ne pas stopper la fonction. 
      #les heatmap pourront ensuite être superposée grâce à l'ajout de champs vides dans la matrice d'heatmap dans la fonction qui fait les heatmaps
      enriched_bg<-(data.frame(
        Term = character(),
        Rank=character(),
        P.value = double(),
        Adjusted.P.value = double(),
        Old.P.value = double(),
        Old.Adjusted.P.value = double(),
        Odds.Ratio= double(),
        Combined.Score= double(),
        Genes = character(),
        dbs = character(),
        gene_number = double(),
        cluster=character(),
        stringsAsFactors = FALSE)) }
    return(enriched_bg)
    
  })
  # les noms des clusters:
  names(enriched_list)<-names(list)
  
  
  
  
  #- compilation dans le même dataframe avec annotation par cluster
  
  enriched_list<-lapply(names(enriched_list), function(x) {
    df<-enriched_list[[x]]
    df<-df %>% mutate(cluster=x)
    return(df)
  })
  names(enriched_list)<-names(list)
  
  enrich_df<-bind_rows(enriched_list)
  if (overexpressed) {
    express_suffix="_overexpressed_genes"
  } else { express_suffix="_underexpressed_genes"}
  write_csv(enrich_df, file = paste0(rootdir,"/unsupervised/enrichment_per_cluster_",sample_name,express_suffix,".csv"))
  return(enrich_df)
  
})
names(enrich_df_list)<-names(sub_list)
return(enrich_df_list)}


# focus on cluster and unsupervides analysis



filter_relevant_iso<-function(seu_obj, #!! mettre l'objet isoformes!!
                              expr_threshold = 1,
                              cell_threshold = 100,
                              more_than_one_iso=FALSE,
                              list, # a list genes of interest
                              output
){
  
  
  df<- data.frame(feature= Features(seu_obj))  
  # split transcript ids into gene and transcript id
  pseudobulk_data <- df %>% 
    mutate(
      gene_name = sub(".*--", "",feature))|> 
    mutate(transcript_id=sub("--.*", "",feature))
  #  Count the number of isoforms per gene
  isoform_count_per_gene <- pseudobulk_data %>%
    group_by(gene_name) %>%
    summarise(n_isoforms = n_distinct(transcript_id)) %>%
    filter(gene_name %in% list)
  # filter those with two or more isoforms
  list_filt1<- isoform_count_per_gene %>%
    filter(n_isoforms>=2)%>%
    pull(gene_name)
  # filter isoforms that are expressed above a specified threshold in a specified number of cells
  expr_matrix <- LayerData(seu_obj, assay = "SCT",layer  = "data")  
  feature_counts <- Matrix::rowSums(expr_matrix > expr_threshold) # calculate the number of cells where each feature is above the cell threshold
  df_expr <- data.frame( # and store the result in a dataframe
    feature = names(feature_counts),
    n_cells_expressed = feature_counts,
    row.names = NULL
  )
  # filter out the transcripts that express two isoforms and more and that are expressed in enough cells
  df_filt2<-df_expr %>% mutate(gene_name = sub(".*--", "",feature),
                               transcript_id=sub("--.*", "",feature)) %>%
    filter((gene_name %in% list_filt1) & (n_cells_expressed > cell_threshold )
    ) 
  print(df_filt2)
  # compute the proportion of isoforms that pass the filters for each gene (may be used later but nor yet)
  isoform_count_per_gene_pass <- df_filt2 %>%
    group_by(gene_name) %>%
    summarise(n_isoforms_pass = n_distinct(feature))
  pass_prop <-isoform_count_per_gene %>% full_join(isoform_count_per_gene_pass, by ="gene_name") %>%
    replace_na(replace=list(n_isoforms_pass=0)) %>%
    mutate(pass=n_isoforms_pass/n_isoforms)
  #print(pass_prop%>%filter(pass>0))
  write_csv(pass_prop, file = paste0(output,"/pass_prop.csv"))
  
  if (more_than_one_iso) { # filter the gene list in order to keep only the genes that have two features or more that pass the filters 
    # compute the proportion of isoforms that pass the filters for each gene (may be used later but nor yet)
    isoform_count_per_gene_pass <- df_filt2 %>%
      group_by(gene_name) %>%
      summarise(n_isoforms_pass = n_distinct(feature))
    gene_pass<-isoform_count_per_gene_pass %>%
      filter(n_isoforms_pass >1) %>%
      pull(gene_name)
    df_filt2<-df_filt2 %>%
      filter(gene_name %in% gene_pass)
  }
  write_csv(df_filt2, file=paste0(output, "/filtered_gene_iso.csv"))
  return(df_filt2)
  
}

dot_plot_switch<-function(seu_obj,
                          more_than_one_iso=FALSE,
                          filtered_genes_df,
                          output_dir) {# has to contain a gene_name column and a feature column with the transcripts
  #dotplot pour visualiser les switchs entre les conditions
  seu_obj$sample <- factor(seu_obj$sample , levels = c("ctrl", "S24", "R"))
  Idents(seu_obj)<-"sample"
  if (more_than_one_iso) { # filter the gene list in order to keep only the genes that have two features or more that pass the filters 
    # compute the proportion of isoforms that pass the filters for each gene (may be used later but nor yet)
    isoform_count_per_gene_pass <- filtered_genes_df %>%
      group_by(gene_name) %>%
      summarise(n_isoforms_pass = n_distinct(feature))
    gene_pass<-isoform_count_per_gene_pass %>%
      filter(n_isoforms_pass >1) %>%
      pull(gene_name)
  }
  for (gene in unique(gene_pass)) {
    plot_features_list <- grep(paste0("--",gene,"$"), filtered_genes_df$feature, value = TRUE)
    dplot<-DotPlot(object = seu_obj, #cols = c("green", "blue","red"), 
                   features =plot_features_list, #split.by = "sample"
    )+
      scale_x_discrete(guide = guide_axis(n.dodge=length(plot_features_list)))
    print(dplot)
    ggsave(plot = dplot,
           filename =paste0(output_dir,"/dotplot_",gene,".png") )
    
  }}

# a function to feat the feature plots for all the isoforms of a gene that are present in the specified marker_list

plot_feature_iso2<-function(seu_obj,
                            gene,
                            sample_name,
                            output_dir,
                            marker_list,
                            label=T) { # marker list est la liste de transcrits filtrés qu'on veut plot
  
  seu_obj<-merged_iso
  features <- Features(seu_obj)
  
  plot_features_list <- grep(paste0("--",gene,"$"), intersect(features,marker_list), value = TRUE)
  
  seu_obj$sample <- factor(seu_obj$sample , levels = c("ctrl", "S24", "R"))
  plot1<-FeaturePlot(seu_obj, features = plot_features_list,combine=TRUE,split.by="sample",label=label) +theme(strip.text = element_text(size = 1)) #+labs(title = paste0("Differential transcript expression per sample for ",gene))
  
  ggsave(filename = paste0(output_dir,"/featureplot",paste0(sample_name,"_",gene,".png")), 
         plot = plot1, width=8,height=5*length(plot_features_list),limitsize = FALSE)
  
  
  
  print(paste0("plotted in ",paste0(output_dir,"/featureplot",paste0(sample_name,"_",gene,".png"))))
  return(plot1)}


#sample_prop_per_cluster<-function
