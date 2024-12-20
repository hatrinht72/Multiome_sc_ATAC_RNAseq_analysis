#update packages
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
```
#old
young <- Read10X_h5("young_filtered_feature_bc_matrix.h5")
fragpath <- "young_atac_fragments.tsv.gz"
OUT_DIR = "/home/ha/y1/projects/Multiome_sc_ATAC_RNAseq_analysis/3_Secondary_Analysis/"
#annotations
saveRDS(annotations, "annotations.rds")
annotations = readRDS("annotations.rds")
# create a Seurat object containing the RNA adata

young_meta<- read.csv(
  file = 'young_per_barcode_metrics.csv',
  header = TRUE,
  row.names = 1)


so_young <- CreateSeuratObject(
  counts = young$`Gene Expression`,
  assay = "RNA",
  meta.data = young_meta
)
young_atac <- young$Peaks

# create ATAC assay and add it to the object
so_young[["ATAC"]] <- CreateChromatinAssay(
  counts = young$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations
)

#quality control
DefaultAssay(so_young) <- "ATAC"

so_young <- NucleosomeSignal(so_young)
so_young <- TSSEnrichment(so_young)

# add fraction of reads in peaks
so_young$pct_reads_in_peaks <- so_young$atac_peak_region_fragments / so_young$atac_fragments * 100
so_young$blacklist_fraction <- FractionCountsInRegion(
  object = so_young, 
  assay = 'ATAC',
  regions = blacklist_mm10
)

png(paste0(OUT_DIR, "3_1_ATAC_QC_metrics_young.png"), res = 150)

VlnPlot(
  object = so_young,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

saveRDS(so_young, "so_young.rds")
#relationship between variables 
png(paste0(OUT_DIR, "3_2_ATAC_relation_count_tss_young.png"), res = 150)
DensityScatter(so_young, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()


low_ncountATAC <- quantile(old[["nCount_ATAC"]]$nCount_ATAC, probs = 0.02)
hig_ncountATAC <- quantile(old[["nCount_ATAC"]]$nCount_ATAC, probs = 0.98)
hig_nCountRNA <- quantile(old[["nCount_RNA"]]$nCount_RNA, probs = 0.98)
low_nCountRNA  <- quantile(old[["nCount_RNA"]]$nCount_RNA, probs = 0.02)
hig_ns <- quantile(old[["nCount_RNA"]]$nCount_RNA, probs = 0.98)
low_ts <- quantile(old[["nCount_RNA"]]$nCount_RNA, probs = 0.02)

#afficher vinplot
library(ggplot2)

# Calcul des seuils pour 'nCount_ATAC' et 'nCount_RNA'
low_ncountATAC <- quantile(old[["nCount_ATAC"]]$nCount_ATAC, probs = 0.02)
hig_ncountATAC <- quantile(old[["nCount_ATAC"]]$nCount_ATAC, probs = 0.98)
low_nCountRNA  <- quantile(old[["nCount_RNA"]]$nCount_RNA, probs = 0.02)
hig_nCountRNA <- quantile(old[["nCount_RNA"]]$nCount_RNA, probs = 0.98)
low_ts <- quantile(old[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)
hig_ts <- quantile(old[["TSS.enrichment"]]$TSS.enrichment, probs = 0.98)

nCount_peaks < 100000 &
  pct_reads_in_peaks > 40 &
  blacklist_ratio < 0.01 &
  nucleosome_signal < 4 &
  TSS.enrichment > 4



old_test <- subset(
  x = old,
  subset = nCount_ATAC < hig_ncountATAC &
    nCount_RNA < hig_nCountRNA &
    nCount_ATAC > low_ncountATAC &
    nCount_RNA > low_nCountRNA &
    nucleosome_signal < hig_ns &
    TSS.enrichment > low_ts
)

#RNA
DefaultAssay(old_test) <- "RNA"
old_test <- SCTransform(old_test)
old_test <- RunPCA(old_test)

old_test <- NormalizeData(old_test)

# Mise à l'échelle des données
old_test <- ScaleData(old_test)

# Exécutez PCA
old_test <- RunPCA(old_test)

