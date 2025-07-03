# Tutorial of single-cell Multiome sequencing in R 
### Table of content 
 * [Section 1. Mono-modal data analysis of the scMultiome data](#section-1-mono-modal-data-analysis-of-the-scmultiome-data)
   * [0. Import the required packages](#0-import-the-required-packages)
   * [1. Load the data and create Seurat object](#1-load-the-data-and-create-seurat-object)
     * [1.1. RNA](#11-RNA)
     * [1.2. ATAC](#12-ATAC)
         * [1.2.1 Create a common peak set](#121-create-a-common-peak-set)
         * [1.2.2 Create fragment objects](#122-create-fragment-objects)
         * [1.2.3 Quantify peaks in each dataset](#123-quantify-peak-in-each-dataset)
         * [1.2.4 Create chromatin assay](#124-create-chromatin-assay)
  
### 0. Import the required packages 
We can use a variety of tools to analyze single cell sequencing data, depending on the language we choose and the type of data we wish to examine. Seurat and Signac on R were selected for this tutorial due of their improvements, extensive documentation, and beginner-friendly nature. 
```R
library(Signac)
library(Seurat)
library(GenomicRanges)
```
To make the command easier to use later, we may also set the input and output directory paths. 
```R
INPUT_DIR = "/path/to/your/data/"
OUTPUT_DIR ="/path/to/the/results/"
```
### 1 Load the data and create Seurat object
#### 1.1. RNA
Cellranger ARC immediately generates RNA data that correlates to the genes and cell matrix and ATAC data that correlates to the peaks and cell matrix. With low-quality cells and genes/peaks removed, they offer the raw and filtered matrix. 
We can do the filtering on our own by combining multiple tools, such as [DoubletFinder] (https://github.com/chris-mcginnis-ucsf/DoubletFinder) for gene expression level, [Amulet] (https://github.com/UcarLab/AMULET) for chromatin accessibility information, and [EmptyDropsMultiome] (https://github.com/MarioniLab/EmptyDropsMultiome) to identify low quality droplets (empty droplets) after sequencing.
The directly filtered matrix produced by Cellranger Arc will be used in this case.
In addition to importing the filtered matrix for each dataset—old and young—we will first load the meta data, which includes comprehensive information for each modality sequencing result.  
```R
# meta data
old_meta <- read.csv(
  paste0(INPUT_DIR,file = 'old_per_barcode_metrics.csv'),
  header = TRUE,
  row.names = 1)

young_meta <- read.csv(
  paste0(INPUT_DIR,file = 'young_per_barcode_metrics.csv'),
  header = TRUE,
  row.names = 1)

#Filtered matrix
old_raw <-Read10X_h5(paste0(INPUT_DIR,"old_filtered_feature_bc_matrix.h5"))
young_raw <- Read10X_h5(paste0(INPUT_DIR,"young_filtered_feature_bc_matrix.h5"))
```
Now we can use these matrix to create the RNA assay in our Seurat objects

```R
so_old<- CreateSeuratObject(
  counts = old_raw$`Gene Expression`,
  assay = "RNA",
  meta.data = old_meta,
  project="Old"
)
so_young<- CreateSeuratObject(
  counts = young_raw$`Gene Expression`,
  assay = "RNA",
  meta.data = old_meta,
  project="Young"
)
```
To keep the R environment as light as possible, it is preferable to get rid of the old_raw and young_raw since we won't be using them for our upcoming ATAC assay generation.  
```R
rm(young_raw)
rm(old_raw)  
gc()  
```

#### 1.2. ATAC
##### 1.2.1 Create a common peak set
As with the RNA assay, we may create our ATAC assay directly from our filtered matrix if we are working with a single dataset. 
However, ATAC data is more sparse than RNA data, which has a set of genes that are similar across all datasets. Among samples under the same experimental environment, as well as among various datasets, the precision of the Tn5 binding detected by sequencing may differ. In the end, if we wish to compare different conditions, the variance we recorded might come solely from technological noise rather than biological information. Therefore, in order to reach an agreement, we must first establish a shared set of peaks from which we can perform peak calling.  
Here, we will compare the peak between two conditions—young and old—with additional ATAC data produced by other studies.
We first import our peaks.
```R
peaks <- read.delim(paste0(INPUT_DIR,"/PMC10183972/GSE184851_all_blood_cells_masterpeaklist_anno.txt")[,2:5]

peaks.old <- read.table(
  file = paste0(INPUT_DIR,"old_atac_peaks.bed"),
  col.names = c("chr", "start", "end")
)

peaks.young <- read.table(
  file = paste0(INPUT_DIR,"young_atac_peaks.bed"),
  col.names = c("chr", "start", "end")
)
```
Then we convert the table we just imported into genomic ranges objects
```R
gr <- makeGRangesFromDataFrame(peaks)
gr.old <- makeGRangesFromDataFrame(peaks.old)
gr.young <- makeGRangesFromDataFrame(peaks.young)
```

How ever, one crucial point is the fragment file must be stored in the same folder with .tbi file to be able map the position.
```

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
