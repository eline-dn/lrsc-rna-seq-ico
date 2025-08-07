#!/bin/bash
# Usage : ./run_sqanti3.sh <rootdir> <path to unzipped reference annotation gtf used in flames> <path to unzipped reference genome.fa used in flames> <tappas_ref>
#we have to merge the 3 isoform_annotated.gtf from flames into 1 using AGAT's #agat_sp_merge_annotations.pl

rootdir=$1
cd $rootdir
mkdir "$rootdir/sqanti"
conda create -c bioconda -n agat agat 
conda activate agat
mkdir "$rootdir/agat"
agat_sp_merge_annotations.pl --gff "$rootdir/flames/ctrl/isoform_annotated.gff3" --gff "$rootdir/flames/R/isoform_annotated.gff3" --gff "$rootdir/flames/S24/isoform_annotated.gff3" --out "$rootdir/agat/merged_not_filtered.gff"

agat_sp_statistics.pl --gff "$rootdir/agat/merged3.gff" -o "$rootdir/agat/stats" 

agat_sp_functional_statistics.pl --gff "$rootdir/agat/merged3.gff" -o "$rootdir/agat/stats_function"

#filter gtf file to remove annotation of unknown strands and keep only transcripts


#convert gff3file in gtf file for sqanti3 input

agat_convert_sp_gff2gtf.pl --gff "$rootdir/agat/merged3.gff" -o "$rootdir/agat/merged.gtf" 

# add gene_name flag to ref annotation (refGTF) for isoAnnotLite

#GTF="/home/labct/Documents/stage_eline_denis/gencode.v48.annotation.gtf"
#genome="/home/labct/Documents/stage_eline_denis/genome.fa"

#GTF=/media/inserm-root/P4/Single_cell_rna_seq_analyses_eline/ref/gencode.v48.annotation.gtf
#genome=/media/inserm-root/P4/Single_cell_rna_seq_analyses_eline/ref/genome.fa

GTF=$2
genome=$3


# warning: very ressource-intensive, run only once if possible
agat_sp_manage_attributes.pl -gff $GTF -att gene_id/gene_name --cp -o "$rootdir/agat/gencode.v48.annotation_with_gene_name.gtf"


conda deactivate

#install SQANTi3
conda update conda
cd "$rootdir/sqanti"
wget https://github.com/ConesaLab/SQANTI3/releases/download/v5.5/SQANTI3_v5.5.zip
mkdir sqanti3
unzip SQANTI3_v5.5.zip -d sqanti3

conda env create -f sqanti3/release_sqanti3/SQANTI3.conda_env.yml

#run SQANTIpython3 
conda activate sqanti3
echo "Launching SQANT3"
GTF="$rootdir/agat/gencode.v48.annotation_with_gene_name.gtf"



tappas_ref=$4
#tappas_ref=/media/inserm-root/P4/Single_cell_rna_seq_analyses_eline/ref/Homo_sapiens_Gencode_v39/human_tappas_gencode_annotation_file.gff3

#/home/labct/Documents/stage_eline_denis/Homo_sapiens_GRCh38_Ensembl_86.gff3

tappas_ref=/media/inserm-root/P4/Single_cell_rna_seq_analyses_eline/ref/Homo_sapiens_GRCh38_Ensembl_86.gff3

$rootdir/sqanti/sqanti3/release_sqanti3/sqanti3_qc.py --isoforms "$rootdir/agat/merged.gtf" --refGTF $GTF --refFasta $genome \
#--isoAnnotLite --gff3 $tappas_ref

echo "Workflow terminé avec succès."
