### 1.1 Loading library 
Make sure all packages are installed
```
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
```
### 1.2 Loading files
#### 1.2.1 Loading count files
We will use directly filtered matrix created by cellranger arc.
For RNA data, its correspond to genes/cell matrix.
For ATAC data, its correspond to peak/cell matrix. 
```
INPUT_DIR <- "/path/to/directory/"
old <- Read10X_h5(filename = paste0(INPUT_DIR,file= "old_filtered_feature_bc_matrix.h5"))

young <- Read10X_h5(filename = paste0(INPUT_DIR,file="young_filtered_feature_bc_matrix.h5"))
```
#### 1.2.2 Files preparation
For ATAC seq data, we will need metadata file and the fragments files to be able using the other function of Signac. 
The fragment file must be stored in the same folder with .tbi file to be able map the position 

```
old_meta <- read.csv(
  paste0(INPUT_DIR,file = 'old_per_barcode_metrics.csv'),
  header = TRUE,
  row.names = 1)

young_meta <- read.csv(
  paste0(INPUT_DIR,file = 'young_per_barcode_metrics.csv'),
  header = TRUE,
  row.names = 1)
```

To identify the peak of open chromatin region, we need a annotation for each regions, so we have to prepare our annotations data

```
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
###Convert annotation seqlevels to UCSC style (i.e., prefix with "chr")
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
###Set the genome to mm10
genome(annotations) <- "mm10"
#Save the annotations dataset to be able to use it directly next time
save(annotations,paste0(INPUT_DIR, file="annotations_mm10.rds")) 
```

We can also store the path to our fragment files to facilitate the code writing after
```
old_fragpath = paste0(INPUT_DIR,file="old_atac_fragments.tsv.gz")
young_fragpath = paste0(INPUT_DIR,file="young_atac_fragments.tsv.gz")
```

#### 1.2.3 Create Seurat object
Create a Seurat object containing the RNA data first :
```
so_old<- CreateSeuratObject(
  counts = old$`Gene Expression`,
  assay = "RNA",
  meta.data = old_meta,
  project="Old_wk4"
)

so_young <- CreateSeuratObject(
  counts = young$`Gene Expression`,
  assay = "RNA",
  meta.data = young_meta,
  project="Young_wk4",
)
```
Then add ATAC assay to the object 
```
so_young[["ATAC"]] <- CreateChromatinAssay(
  counts = young$Peaks,
  sep = c(":", "-"),
  fragments = young_fragpath,
  annotation = annotations
)

so_old[["ATAC"]] <- CreateChromatinAssay(
  counts = old$Peaks,
  sep = c(":", "-"),
  fragments = old_fragpath,
  annotation = annotations
)
```
So now we can move to next part : Filtering 

### 1.3 Filtering ATAC and RNA
We can first clean our data then analyze independently ATAC and RNA.
We can start with ATAC 
#### 1.3.1 ATAC metrics
We have to set our default assay first
```
DefaultAssay(so_young) <- "ATAC"
DefaultAssay(so_old) <- "ATAC"
```
To filter scATACseq data, there is several parameters to look after :
- Nucleosome binding pattern : to make sure we have the right DNA fragment size, corresponding to the length of the DNA wrapped around a single nucleosome (147 pb) and to assess the signal from nucleosome positioning based on the fragment length distribution from the ATAC-seq data. The nucleosome_signal will store the value of the appromimate ratio of mononucleosomam and nucleosome-free fragments
- Transcriptional start site (TSS) enrichnment score : to assess the enrichment of transcription start sites (TSS) in ATAC-seq data, calculates TSS enrichment scores, which help determine the degree to which TSS regions are accessible in your chromatin accessibility data. Poor scATAC seq experiments would have low TSS enrichment score
- Total number of fragments in peaks : mesure depth sequencing; identify the number of fragments (or reads) that map to specific accessible regions (peaks)
- Fraction of reads in peaks : all fragments that fall within ATAC-seq peaks, access the total number of reads that passed quality control filters during preprocessing
- Ratio reads in genomic blacklist regions : reads which are often represent low-quality cells/technical artifacts. Access the proportion of reads that map to regions defined as "blacklisted.", consisting problematic genomic regions, often due to high background noise or sequencing artifacts, and should be considered with caution when interpreting results

```
so_young <- NucleosomeSignal(so_young) #fragment ratio 147-294: <147
so_old <- NucleosomeSignal(so_old) #fragment ratio 147-294: <147

so_young <- TSSEnrichment(object = so_young, fast = FALSE)
so_old <- TSSEnrichment(object = so_old, fast = FALSE)

so_young$blacklist_fraction <- FractionCountsInRegion(
  object = so_young, 
  assay = 'ATAC',
  regions = blacklist_mm10
)
so_old$blacklist_fraction <- FractionCountsInRegion(
  object = so_old, 
  assay = 'ATAC',
  regions = blacklist_mm10
)

##Number of fragments (or reads) that map to specific accessible regions (peaks)
old_data_atac$atac_peak_region_fragments 

so_young$pct_reads_in_peaks <- so_young$atac_peak_region_fragments / so_young$atac_fragments * 100
so_old$pct_reads_in_peaks <- so_old$atac_peak_region_fragments / so_old$atac_fragments * 100
```

we can then visualize all these parameter by violin plot.

```
png(paste0(OUT_DIR, "1_3_1_ATAC_QC_metrics_young.png"), res = 150, width = 2000, height = 1000 )
VlnPlot(
  object = so_young,
  features = c('atac_peak_region_fragments',  'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()
```
![2_1_ATAC_QC_metrics_young](https://github.com/user-attachments/assets/f76bab2a-a3ba-48e0-ae07-7730128a2c49)

```
png(paste0(OUT_DIR, "1_3_1_ATAC_QC_metrics_old.png"), res = 150, width = 2000, height = 1000 )
VlnPlot(
  object = so_old,
  features = c('atac_peak_region_fragments',  'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()
```
![1_3_1_ATAC_QC_metrics_old](https://github.com/user-attachments/assets/89af65be-3417-494e-a63b-c425addc1f2f)

we can also visualize the relationship between variables.
```
png(paste0(OUT_DIR, "1_3_1_ATAC_relation_count_tss_young.png"), res = 150,width = 1000, height = 1000 )
DensityScatter(so_young, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
```
![2_1_ATAC_relation_count_tss_young](https://github.com/user-attachments/assets/28ad481f-9973-4d1a-99a7-bb59721e3ab5)

```
png(paste0(OUT_DIR, "1_3_1_ATAC_relation_count_tss_old.png"), res = 150,width = 1000, height = 1000 )
DensityScatter(so_old, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
```
![1_3_1_ATAC_relation_count_tss_old](https://github.com/user-attachments/assets/0608d5e2-c540-418a-b480-763baf8f1568)

#### 1.3.2 RNA metrics
We can also visualize the RNA data.
```
png(paste0(OUT_DIR, "1_3_2_QC_RNA_young.png"), width = 1200, height = 1000, res = 150)
VlnPlot(so_young, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()
```
![1_3_2_RNA_QC_vinplot_young](https://github.com/user-attachments/assets/57bbd2c9-5d6a-4ec0-b931-e625bb1dfaa0)
```
png(paste0(OUT_DIR, "1_3_2_QC_RNA_old.png"), width = 1200, height = 1000, res = 150)
VlnPlot(so_old, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()
```
![1_3_2_RNA_QC_vinplot_old](https://github.com/user-attachments/assets/298b6dad-c53e-4e06-8d3e-53585d1a37c3)

```
png(paste0(OUT_DIR, "1_3_2_RNA_QC_scatter_young.png"), width = 1200, height = 1000, res = 150)
FeatureScatter(so_young, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
```
![1_3_2_RNA_QC_count_feature_young](https://github.com/user-attachments/assets/82bda713-9467-4404-851a-9a03392168e6)

```
png(paste0(OUT_DIR, "1_3_2_RNA_QC_scatter_old.png"), width = 1200, height = 1000, res = 150)
FeatureScatter(so_old, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
```
![1_3_2_RNA_QC_scatter_old](https://github.com/user-attachments/assets/1dcfa85d-c285-4709-a028-44566840bce2)

You can save now the object in case of you want to try different way to filter the cells
```
save(so_young, file="so_young")
save(so_old, file="so_old")

```

#### 1.3.3 Filtering 
##### 1.3.3.1 Set the thredshold
Then we can set some thredshold to filtering data by these parameters. I will keep only the most 98% abundant of these parameters.
```
#ATAC
low_prf <- quantile(so_young[["atac_peak_region_fragments"]]$atac_peak_region_fragments, probs = 0.02)
hig_prf <- quantile(so_young[["atac_peak_region_fragments"]]$atac_peak_region_fragments, probs = 0.98)
low_prp <- quantile(so_young[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02)
high_blr <- quantile(so_young[["blacklist_fraction"]]$blacklist_fraction, probs = 0.98)
hig_ns <- quantile(so_young[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)
low_ts <- quantile(so_young[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)


low_prf <- quantile(so_old[["atac_peak_region_fragments"]]$atac_peak_region_fragments, probs = 0.02)
hig_prf <- quantile(so_old[["atac_peak_region_fragments"]]$atac_peak_region_fragments, probs = 0.98)
low_prp <- quantile(so_old[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02)
high_blr <- quantile(so_old[["blacklist_fraction"]]$blacklist_fraction, probs = 0.98)
hig_ns <- quantile(so_old[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)
low_ts <- quantile(so_old[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)

#RNA
low_nCountRNA  <- quantile(so_young[["nCount_RNA"]]$nCount_RNA, probs = 0.02)
hig_nCountRNA <- quantile(so_young[["nCount_RNA"]]$nCount_RNA, probs = 0.98)

low_nCountRNA  <- quantile(so_old[["nCount_RNA"]]$nCount_RNA, probs = 0.02)
hig_nCountRNA <- quantile(so_old[["nCount_RNA"]]$nCount_RNA, probs = 0.98)
```
##### 1.3.3.2 Filter 
```
so_old_filtered <- subset(
  x = so_old,
  subset = atac_peak_region_fragments > low_prf &
    atac_peak_region_fragments < hig_prf &
    pct_reads_in_peaks > low_prp &
    blacklist_fraction < high_blr &
    nucleosome_signal < hig_ns &
    TSS.enrichment > low_ts &
    nCount_RNA < hig_nCountRNA &
    nCount_RNA > low_nCountRNA &
    nFeature_RNA < hig_nFeatureRNA &
    nFeature_RNA > low_nFeatureRNA
)

so_young_filtered <- subset(
  x = so_young,
  subset = atac_peak_region_fragments > low_prf &
    atac_peak_region_fragments < hig_prf &
    pct_reads_in_peaks > low_prp &
    blacklist_fraction < high_blr &
    nucleosome_signal < hig_ns &
    TSS.enrichment > low_ts &
    nCount_RNA < hig_nCountRNA &
    nCount_RNA > low_nCountRNA &
    nFeature_RNA < hig_nFeatureRNA &
    nFeature_RNA > low_nFeatureRNA
)
```
Save our filtered data
```
saveRDS(so_young_filtered, file="so_young_filtered.rds")
saveRDS(so_old_filtered, file="so_old_filtered.rds")
```
Now we can start with the ATAC analysis ! 
