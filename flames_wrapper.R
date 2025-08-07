# the purpose of these functions is to run flames and create seurat-ready count matrixe at genes and transcript levels
# the arguments passed to this script should be rootdir , sample_name and path_to_project

library(FLAMES)


#args = commandArgs(trailingOnly=TRUE)
#rootdir=args[1]
#sample_name=args[2]
#path_to_project=args[3]
setwd(rootdir)  # Set this to correct location


run_flames<-function(fastq,
                     sample_name,
                     annot=paste0(path_to_project,"./references/gencode.v48.annotation.gtf.gz"),
                     genome_fa=paste0(path_to_project,"./references/genome.fa") ){
  library(FLAMES)
  setwd(rootdir)
  dir.create(paste0(rootdir,"/flames/",sample_name), recursive = TRUE, showWarnings = FALSE)
  outdir<-paste0(rootdir,"/flames/",sample_name)
  
  
  config_file <- FLAMES::create_config(outdir, do_barcode_demultiplex = FALSE,type="sc_3end",multithread_isoform_identification=TRUE,oarfish_quantification=FALSE)
  config <- jsonlite::fromJSON(config_file)
  genome_bam <- rownames(minimap2_align(
    config = config, fa_file = genome_fa, fq_in = fastq, annot = annot,
    outdir = outdir))
  find_isoform(
    annotation = annot, genome_fa = genome_fa,
    genome_bam = genome_bam, outdir = outdir, config = config)
  minimap2_realign(
    config = config, fq_in = fastq,
    outdir = outdir)
  quantify_transcript(annotation = annot, outdir = outdir, config = config)
  sce <- create_sce_from_dir(outdir = outdir, annotation = annot)
  saveRDS(sce, file = paste0(outdir,"/sce.rds"))
  transcript_counts<-read.csv(paste0(outdir,"/transcript_count.csv.gz"))
  return(transcript_counts)
}


#Count matrices processing

import_ref_gtf<-function(reference_gtf="/home/labct/Documents/stage_eline_denis/gencode.v48.annotation.gtf.gz"){
  # Import reference gtf file with gene symbols
  gtf2 <- import(reference_gtf)
  gtf2_df <- as.data.frame(gtf2)
  
  # Select relevant columns from the second GTF
  selected_columns2 <- gtf2_df[, c("gene_name", "gene_id")]
  unique_gene_symbol <- unique(selected_columns2)
  print("recuperation de la référence")
  return(unique_gene_symbol)
}

make_transcript_level_counts <- function( sample_name,
                                          transcript_count_matrix,
                                          unique_gene_symbol,
                                          output_dir="./counts",
                                          make_gene_level_counts=TRUE) {
  
  
  
  df_c<-transcript_count_matrix
  
  df_c<-df_c |> left_join(unique_gene_symbol,by="gene_id")
  
  #create transcript level count matrix
  isoform_counts<-df_c|> dplyr::select(-c(gene_id))
  rownames(isoform_counts)<-paste0(isoform_counts$transcript_id, "--", isoform_counts$gene_name)
  isoform_counts<-isoform_counts |> dplyr::select(-c(transcript_id,gene_name))
  
  # Write the output to a CSV file
  output_path <- file.path(output_dir, paste0("gene_symbol_", sample_name, "transcript_level_counts.csv"))
  write.csv(isoform_counts, output_path)
  
  print(paste0("Processed sample: ", sample_name, " Output saved to: ", output_path))
  res<-list()
  res$iso<-isoform_counts
  
  if (make_gene_level_counts==TRUE) {
    #also create a gene level counts matrix
    test<-df_c|> dplyr::select(-c(transcript_id,gene_id)) |>dplyr::group_by(gene_name) %>%dplyr::summarise(across(everything(), sum, na.rm = TRUE))
    gene_counts<-as.data.frame(test)
    rownames(gene_counts)<-gene_counts$gene_name
    gene_counts$gene_name<-NULL
    # Write the output to a CSV file
    output_path <- file.path(output_dir, paste0("gene_symbol_", sample_name, "gene_level_counts.csv"))
    write.csv(gene_counts, output_path)
    print(paste0("Processed sample: ", sample_name, " Output saved to: ", output_path))
    res$gene<-gene_counts
  }
  
  return(res)
}


