if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("IsoformSwitchAnalyzeR")
BiocManager::install("BiocGenerics")


library(IsoformSwitchAnalyzeR)





## préparer la pseudobulk matrix

seu_obj<-readRDS(paste0(rootdir,"/merge_integration/isoform_ctrl_S24_R_merged_unintegrated_seurat.rds"))


#lets aggeragte the expresstion data by cell type 
counts <- AggregateExpression(
  seu_obj, 
  assays = "SCT", 
  return.seurat = FALSE,
  group.by = "sample"
)
head(counts$SCT)

as.data.frame(counts) -> df
row.names(df) -> df$gene

#split transcript ids into gene and transcript id
pseudobulk_data <- df %>% separate(gene, into = c("transcript_id", "gene_id"), sep = "--",  extra = "merge") 

head(pseudobulk_data)



### Étape 0: Définir les chemins vers tes fichiers
# Adapter les chemins à ton système
count_matrix_file <- "/chemin/vers/mon_fichier_comptage.csv"         # output de FLAMES
gtf_file          <- "/chemin/vers/mon_annotation.gtf"              # GTF généré par FLAMES
fasta_file        <- "/chemin/vers/mes_sequences_transcripts.fa"    # Fasta des isoformes
output_dir        <- "/chemin/vers/output/"                         # Dossier de sortie

### Étape 1: Charger la matrice de comptage
# Elle doit être sous forme de tableau : isoform_id | sample1 | sample2 | ...
counts_raw <- read.csv(count_matrix_file, row.names = 1, check.names = FALSE)

# Supposons que tu n’as que les comptages (pas d’abundance en TPM)
# IsoformSwitchAnalyzeR peut calculer les TPM à partir des comptages
# mais attention à la qualité des longueurs si nécessaire

# Créer un design (adapter selon tes échantillons)
myDesign <- data.frame(
  sampleID = colnames(counts_raw),
  condition = c("cond1", "cond1", "cond2")  # <--- à adapter
)

### Étape 2: Créer un objet switchAnalyzeRlist
aSwitchList <- importRdata(
  isoformCountMatrix       = counts_raw,
  isoformRepExpression     = NULL,  # Pas de TPM fourni
  designMatrix             = myDesign,
  isoformExonAnnoation     = gtf_file,
  isoformNtFasta           = fasta_file,
  showProgress             = TRUE
)

summary(aSwitchList)

### Étape 3: Lancer l’analyse des isoform switches (Part 1)
aSwitchList <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = aSwitchList,
  outputSequences      = TRUE,
  pathToOutput         = output_dir,
  prepareForWebServers = TRUE   # génère des fichiers découpés pour webtools
)

extractSwitchSummary(aSwitchList)

### Étape 4: Effectuer les prédictions externes (FAIRE MANUELLEMENT)
# À ce stade, tu dois soumettre les fichiers générés dans `output_dir` :
# - output_dir/isoform_nt_sequences.fasta
# - output_dir/isoform_aa_sequences.fasta

# Voici les outils nécessaires (et comment les utiliser) :

#### 🧬 [OPTIONNEL] CPC2 (ou CPAT) – Potentiel codant
- Website: http://cpc2.cbi.pku.edu.cn
- Input: `isoform_nt_sequences.fasta`
- Output: fichier .txt avec colonne "coding potential"

#### 🧬 Pfam (domaines protéiques)
- Website: https://www.ebi.ac.uk/Tools/pfa/pfam/
  - Input: `isoform_aa_sequences.fasta`
- Output: fichier texte avec coordonnées des domaines

#### 🧬 SignalP (peptides signal)
- Website: https://services.healthtech.dtu.dk/service.php?SignalP
- Input: `isoform_aa_sequences.fasta`

#### 🧬 IUPred2A (régions désordonnées)
- Website: https://iupred2a.elte.hu
- Input: `isoform_aa_sequences.fasta`

⚠️ Regarde si tu peux lancer ces outils en ligne ou en local via leur version CLI si tu veux automatiser plus.

---
  
  ### Étape 5: Intégrer les résultats externes (Part 2)
  
  Une fois les fichiers résultats obtenus, indique les chemins :
  
  ```r
aSwitchList <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = aSwitchList,
  pathToCPC2resultFile      = "/chemin/vers/cpc2_output.txt",   # ou CPAT
  pathToPFAMresultFile      = "/chemin/vers/pfam_results.txt",
  pathToSignalPresultFile   = "/chemin/vers/signalp.txt",
  pathToIUPred2AresultFile  = "/chemin/vers/iupred2a.txt",
  removeNoncodinORFs        = TRUE,
  outputPlots               = TRUE,
  n                         = 10
)