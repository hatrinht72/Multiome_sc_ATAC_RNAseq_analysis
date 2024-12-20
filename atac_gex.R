#scATACseq
#1.Installing/Load
##1.1/installing needed package
###install.packages('hdf5r') #need to read h5 files
###install.packages("Signac") #seurat addon for analyzing chromatin
###install.packages('Seurat')
###if (!requireNamespace("TFBSTools", quietly = TRUE)) {
###install.packages("BiocManager")
###BiocManager::install("TFBSTools")
###}
###if (!requireNamespace("JASPAR2020", quietly = TRUE)) {
###BiocManager::install("JASPAR2020")
###}
###if (!requireNamespace("Biostrings", quietly = TRUE)) {
###BiocManager::install("Biostrings")
###}
###if (!require("BiocManager", quietly = TRUE))
###install.packages("BiocManager")
###BiocManager::install("EnsDb.Mmusculus.v79")#EnsDb.Hsapiens.v86 for human
###BiocManager::install("GenomeInfoDb") #translation between chromosome names
###BiocManager::install("biovizBase")

##1.2/loading needed libraries
library(Signac)
library(Seurat)
library(SeuratObject)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(tidyr)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)#EnsDb.Hsapiens.v86 for human 
library(GenomicRanges)


##1.3/loading needed files
###1.3.1/counts files
old <- Read10X_h5(filename = "old_filtered_feature_bc_matrix.h5")
old_rna <- old$'Gene Expression'
old_atac <- old$'Peaks'

young <- Read10X_h5(filename = "young_filtered_feature_bc_matrix.h5")
young_rna <- young$`Gene Expression`
young_atac <- young$Peaks

###1.3.2/meta files

old_meta<- read.csv(
  file = 'old_barcode_metrics.csv',
  header = TRUE,
  row.names = 1)
young_meta<- read.csv(
  file = 'young_per_barcode_metrics.csv',
  header = TRUE,
  row.names = 1)
#path to fragment file

old_fragment_path ="/home/ha/y1/projects/ATACseq_analysis/old_atac_fragments.tsv"
young_fragpath ="/home/ha/y1/projects/ATACseq_analysis/young_atac_fragments.tsv"
  


####1.3.2.1/select only column with atac
old_meta_atac <- old_meta %>% select(contains("atac"))
young_meta_atac <- dplyr::select(young_meta, contains("atac"))
###1.3.3/gunzip the file tbi if needed for each fragments files we will use : GSM5723631_Young_HSC_fragments.tsv.tbi.gz

#2.Data Processing with Seurat/Signac
##2.1/Create a Chromatin Assay: integrate and analyze chromatin accessibility data (ATAC-seq).  

so_young <-  CreateSeuratObject(
  counts = young_rna,
  assay = "RNA"
)

so_old <-  CreateSeuratObject(
  counts = old_rna,
  assay = "RNA"
)


atac_chrom_assay <- CreateChromatinAssay(
  counts = old_counts_atac,
  sep = c(":", "-"),
  genome = 'mm10',  # Ensure mm10 is registered or try another valid genome
  fragments = 'old_atac_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
young_atac_chrom_assay <- CreateChromatinAssay(
  counts = young_counts_atac,
  sep = c(":", "-"),
  genome = 'mm10',  # Ensure mm10 is registered or try another valid genome
  fragments = 'young_atac_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

##2.2/Create Seurat object 
old_data_atac <- CreateSeuratObject(
  counts = atac_chrom_assay,
  assay = "peaks",
  meta.data = old_meta_atac
)
young_data_atac <- CreateSeuratObject(
  counts = young_atac_chrom_assay,
  assay = "peaks",
  meta.data = young_meta_atac
)

##2.3/Create annotations file for each regions we found
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
###Convert annotation seqlevels to UCSC style (i.e., prefix with "chr")
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
###Set the genome to mm10
genome(annotations) <- "mm10"
###Set the annotation for each Seurat object  
Annotation(old_data_atac) <- annotations
Annotation(young_data_atac) <- annotations


so_young[["ATAC"]] <- CreateChromatinAssay(
  counts = young_atac,
  sep = c(":", "-"),
  fragments = young_fragpath,
  annotation = annotations
)


##2.4/Quality Control (to filtering data that we achieved)
###2.4.1/Assess the signal from nucleosome positioning based on the fragment length distribution from the ATAC-seq data
old_data_atac <- NucleosomeSignal(object = old_data_atac) #fragment ratio 147-294: <147
young_data_atac <- NucleosomeSignal(object = young_data_atac) #fragment ratio 147-294: <147
###2.4.2/assess the enrichment of transcription start sites (TSS) in ATAC-seq data, calculates TSS enrichment scores, which help determine the degree to which TSS regions are accessible in your chromatin accessibility data
old_data_atac <- TSSEnrichment(object = old_data_atac, fast = FALSE)
young_data_atac <- TSSEnrichment(object = young_data_atac, fast = FALSE)
###2.4.3/Access the proportion of reads that map to regions defined as "blacklisted.", consisting problematic genomic regions, often due to high background noise or sequencing artifacts, and should be considered with caution when interpreting results
####data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments
###2.4.4/Number of fragments (or reads) that map to specific accessible regions (peaks)
####old_data_atac$atac_peak_region_fragments 
###2.4.5/access the total number of reads that passed quality control filters during preprocessing
####pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100 
###2.4.6/Visualize
VlnPlot(
  object = old_data_atac,
  features = c('atac_peak_region_fragments', 'pct_reads_in_peaks', 
               'blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'),
  pt.size = 0.1,
  ncol = 5
)
VlnPlot(
  object = young_data_atac,
  features = c('atac_peak_region_fragments', 'pct_reads_in_peaks', 
               'blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'),
  pt.size = 0.1,
  ncol = 5
)
###2.4.7/FIltering the our data
####2.4.7.1/Identify the thresholds for each parameter
low_prf <- quantile(old_data_atac[["atac_peak_region_fragments"]]$atac_peak_region_fragments, probs = 0.02)
hig_prf <- quantile(old_data_atac[["atac_peak_region_fragments"]]$atac_peak_region_fragments, probs = 0.98)
#low_prp <- quantile(data[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02)
#high_blr <- quantile(data[["blacklist_ratio"]]$blacklist_ratio, probs = 0.98)
hig_ns <- quantile(old_data_atac[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)
low_ts <- quantile(old_data_atac[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)


low_prf <- quantile(young_data_atac[["atac_peak_region_fragments"]]$atac_peak_region_fragments, probs = 0.02)
hig_prf <- quantile(young_data_atac[["atac_peak_region_fragments"]]$atac_peak_region_fragments, probs = 0.98)
#low_prp <- quantile(data[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02)
#high_blr <- quantile(data[["blacklist_ratio"]]$blacklist_ratio, probs = 0.98)
hig_ns <- quantile(young_data_atac[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)
low_ts <- quantile(young_data_atac[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)

###2.4.7.2/Final filter
old_data_atac <- subset(
  x = old_data_atac,
  subset = atac_peak_region_fragments > low_prf &
    atac_peak_region_fragments < hig_prf &
    nucleosome_signal < hig_ns &
    TSS.enrichment > low_ts
)
young_data_atac <- subset(
  x = young_data_atac,
  subset = atac_peak_region_fragments > low_prf &
    atac_peak_region_fragments < hig_prf &
    nucleosome_signal < hig_ns &
    TSS.enrichment > low_ts
)

##2.5/Normalization, dimension reduction
###2.5.1/add dataset for each data file 
young_data_atac$dataset <- "young"
old_data_atac$dataset <- "old"
##2.5.2/Merge the two datasets
data <- merge(young_data_atac, old_data_atac, all = TRUE)
###Find top features which differenciated datasets  
data <- FindTopFeatures(data, min.cutoff = 'q0')
data <- RunTFIDF(data)
data <- RunSVD(data)
###Visualizing - like PCA of single cell
DepthCor(data)
### To identify which dimension will we take 
### only use the dimension thats we need here
data <- RunUMAP(object = data, reduction = 'lsi', dims = 2:30)
data <- FindNeighbors(object = data, reduction = 'lsi', dims = 2:30)
data <- FindClusters(object = data, verbose = FALSE, algorithm = 3, resolution = .4)
###UMAP 
DimPlot(object = data, label = TRUE) + NoLegend()
DimPlot(object = data, label = TRUE, group.by = "dataset") + NoLegend()

#3/Data analysis
##3.1/estimation expression of RNA based on ATAC-seq
gene.activities <- GeneActivity(data)
###data[['RNA']] <- CreateAssayObject(counts = gene.activities)
data <- NormalizeData(
  object = data,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(data$nCount_RNA)
)
DefaultAssay(data) <- 'RNA'
FeaturePlot(
  object = data,
  features = c('Kit', 'Pecam1', 'Itgam'),#gene marker of HSC  
  max.cutoff = 'q95'
)
##3.2/Peak analyzing
###3.2.1/Create the da_peaks file
DefaultAssay(data) <- 'peaks'
da_peaks <- FindMarkers(
  object = data,
  ident.1 = rownames(data[[]][data$dataset == "old",]),
  ident.2 = rownames(data[[]][data$dataset == "young",]),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments'
)
###features corresponds to each Pca 
da_peaks
###Find the closest gene with the chromatine region
da_peaks$closest_gene <-ClosestFeature(data, regions = rownames(da_peaks))$gene_name
###If not_matching_data
####Exclude regions with the problematic seqlevels from da_peaks
valid_regions <- rownames(da_peaks)[!grepl("GL456216.1|JH584304.1", rownames(da_peaks))]
####only with good regions
closest_genes <- ClosestFeature(data, regions = valid_regions)
###add new column for gene name
da_peaks$closest_gene <- NA
da_peaks[valid_regions, "closest_gene"] <- closest_genes$gene_name
da_peaks$distance <- ClosestFeature(data, regions = rownames(da_peaks))$distance
###add new column for distance
da_peaks$distance <- NA
da_peaks[valid_regions, "distance"] <- closest_genes$distance
###3.2.2/ Visualize with significant features 
CoveragePlot(
  object = data,
  region = rownames(da_peaks)[1],
  extend.upstream = 10000,
  extend.downstream = 5000,
  group.by = "dataset"
)
plot1 <- VlnPlot(
  object = data,
  features = rownames(da_peaks)[2],
  group.by = "dataset"
)
plot2 <- FeaturePlot(
  object = data,
  features = rownames(da_peaks)[2],
  max.cutoff = 'q95'
)

plot1 | plot2

##3.3/Motif analysis with JASPAR  
###3.3.1/convert peak region into bed file 
####Extract region data from rownames and bind it as columns
da_peaks_bed <- da_peaks %>%
  mutate(region = rownames(da_peaks)) %>%
  separate(region, into = c("chr", "start", "end"), sep = "-") %>%
  dplyr::select(chr, start, end, closest_gene, distance)
#####we have so many select method so we have to precise which select from what package that we want to use

###Ensure start and end are numeric for BED format
da_peaks_bed$start <- as.numeric(da_peaks_bed$start)
da_peaks_bed$end <- as.numeric(da_peaks_bed$end)

###Optionally save as a BED file
write.table(da_peaks_bed, "da_peaks.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

###3.3.2/Loading files 
####Load JASPAR motifs, sequences
opts <- list(species = 10090, all_versions = FALSE) # Adjust species (e.g., 9606 for human, 10090 for mouse)
pwm_list <- getMatrixSet(JASPAR2020, opts)
#pwm_list: contains transcription factor (TF) position weight matrices (PWMs)
###Extract motif IDs and associated TF names
motif_metadata <- lapply(pwm_list, function(pwm) {
  data.frame(
    motif_id = ID(pwm),
    tf_name = name(pwm),
    stringsAsFactors = FALSE
  )
})
# Combine metadata into a single data frame
motif_metadata_df <- do.call(rbind, motif_metadata)
####Load sequences from the .fa file
sequences <- readDNAStringSet("/home/ha/y1/tutorial/cellranger/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa")
####Create pwm matrix 
pwm_matrices <- lapply(pwm_list, function(pwm) {
  as.matrix(pwm)
})
###3.3.3/Binding Site Scan
# Loop through each PWM and match to each sequence
binding_results <- list()
for (i in seq_along(pwm_matrices)) {
  pwm <- pwm_matrices[[i]]
  pwm_name <- names(pwm_list)[i]  # Get the name for the PWM
  
  # Run motif matching for each sequence in the .fa file
  sites <- lapply(sequences, function(seq) matchPWM(pwm, seq, min.score = "80%"))
  
  # Store results with PWM name
  binding_results[[pwm_name]] <- sites
}




###3.3.4/Process and Integrate Binding Results
####Convert Binding results to data frame

top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])


enriched.motifs <- FindMotifs(
  object = data,
  features = top.da.peak
)
install.packages("motifmatchr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("motifmatchr")
n