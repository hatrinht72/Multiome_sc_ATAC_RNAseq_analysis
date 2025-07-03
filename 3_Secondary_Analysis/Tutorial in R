# Tutorial of single-cell Multiome sequencing in R 
### Table of content 
 * [Section 1. Mono-modal data analysis](#section-1-mono-modal-data-analysis)
   * [0. Import the required packages](#0-import-the-required-packages)
   * [1. Load the data and create Seurat object](#1-load-the-data-and-create-seurat-object)
     * [1.1. RNA](#11-RNA)
     * [1.2. ATAC](#12-ATAC)
         * [1.2.1 Create a common peak set](#121-create-a-common-peak-set)
         * [1.2.2 Quantify peaks in each dataset](#122-quantify-peak-in-each-dataset)
         * [1.2.3 Create chromatin assay](#123-create-chromatin-assay)
   * [2. Data filtering](#2-data-filtering)
     * [2.1. ATAC metrics](#21-ATAC-metrics)
     * [2.2. RNA metrics](#22-RNA-metrics)
     * [2.3. Filtering](#23-Filtering)

## Section 1 : Mono-modal data analysis  
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
We can do the filtering on our own by combining multiple tools, such as [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) for gene expression level, [Amulet](https://github.com/UcarLab/AMULET) for chromatin accessibility information, and [EmptyDropsMultiome](https://github.com/MarioniLab/EmptyDropsMultiome) to identify low quality droplets (empty droplets) after sequencing.
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
Now we can combine them to get a unified set of peaks to quantify in each dataset later on, also process the filtering bad peaks based on length. 
```R
combined.peaks <- reduce(x = c(gr, gr.old, gr.young))
# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
```

##### 1.2.2 Quantify peaks object in each dataset
Now we have the common set of peak, we can generate new matrix count (peak/cells) for each dataset
First we gonna need to create the fragment objects to store fragment file information. One crucial point is the fragment file must be stored in the same folder with .tbi file to be able map the position.

```R
frags.old <- CreateFragmentObject(
  path = paste0(INPUT_DIR,"old_atac_fragments.tsv.gz"),
  cells = colnames(so_old)
)

frags.young <- CreateFragmentObject(
  path = paste0(INPUT_DIR,"young_atac_fragments.tsv.gz"),
  cells = colnames(so_young)
```
Then we can process the peak quantification. 
```R
old.counts <- FeatureMatrix(
  fragments = frags.old,
  features = combined.peaks,
  cells = colnames(so_old)
)

young.counts <- FeatureMatrix(
  fragments = frags.young,
  features = combined.peaks,
  cells = colnames(so_young)
)
```
Another important information is the annotation for each regions detected in our common peak. 
The annotations can contain the gene coded by the regions, or the non coding elements. 
```
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
###Convert annotation seqlevels to UCSC style (i.e., prefix with "chr")
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
###Set the genome to mm10
genome(annotations) <- "mm10"
#Save the annotations dataset to be able to use it directly next time
save(annotations,paste0(INPUT_DIR, file="annotations_mm10.rds")) 
```
##### 1.2.3 Create chromatin assay
Now we can add chromatin assay on our Seurat object that we cretaed before with RNA assay
```R
so_old[["ATAC"]] <- CreateChromatinAssay(
  counts = old.counts,
  sep = c(":", "-"),
  fragments = frags.old,
  annotation = annotations
)
so_young[["ATAC"]] <- CreateChromatinAssay(
  counts = young.counts,
  sep = c(":", "-"),
  fragments = frags.young,
  annotation = annotations
)
```
So now we can move to next part : Filtering 

### 2. Data filtering
We can first clean our datasets independently.

#### 2.1 ATAC metrics
Since now in our Seurat object, we have 2 assays (RNA and ATAC), so we have to set our default assay first
```R
DefaultAssay(so_young) <- "ATAC"
DefaultAssay(so_old) <- "ATAC"
```
To filter scATACseq data, there is several parameters to look after :
- Nucleosome binding pattern : to make sure we have the right DNA fragment size, corresponding to the length of the DNA wrapped around a single nucleosome (147 pb) and to assess the signal from nucleosome positioning based on the fragment length distribution from the ATAC-seq data. The nucleosome_signal will store the value of the appromimate ratio of mononucleosomam and nucleosome-free fragments
- Transcriptional start site (TSS) enrichnment score : to assess the enrichment of transcription start sites (TSS) in ATAC-seq data, calculates TSS enrichment scores, which help determine the degree to which TSS regions are accessible in your chromatin accessibility data. Poor scATAC seq experiments would have low TSS enrichment score
- Total number of fragments in peaks : mesure depth sequencing; identify the number of fragments (or reads) that map to specific accessible regions (peaks)
- Fraction of reads in peaks : all fragments that fall within ATAC-seq peaks, access the total number of reads that passed quality control filters during preprocessing
- Ratio reads in genomic blacklist regions : reads which are often represent low-quality cells/technical artifacts. Access the proportion of reads that map to regions defined as "blacklisted.", consisting problematic genomic regions, often due to high background noise or sequencing artifacts, and should be considered with caution when interpreting results
```R
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
so_young$pct_reads_in_peaks <- so_young$atac_peak_region_fragments / so_young$atac_fragments * 100
so_old$pct_reads_in_peaks <- so_old$atac_peak_region_fragments / so_old$atac_fragments * 100
```
#### 2.2 RNA metrics
For RNA data, we gonna add a new observation called mitochondrial percentage, to estimate the contamination of mitochondrial gene captured by the sequencing. The cells with high mitochondrial percent could be considered as low quality. We can also look into ribosome gene. 
```R
DefaultAssay(so_old) <-"RNA"
DefaultAssay(so_young) <-"RNA"
so_old[["percent.mt"]] <- PercentageFeatureSet(so_old, pattern = "^mt-")
so_young[["percent.mt"]] <- PercentageFeatureSet(so_young, pattern = "^mt-")
```
By this step, we can save our data to be able to try different filtering 
#### 2.3 Filtering
we can fist visualize the data before filtering 
```R
VlnPlot(
  object = so_old,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks','nCount_RNA','nFeature_RNA','percent.mt'),
  pt.size = 0.1,
  ncol = 5
)
VlnPlot(
  object = so_young,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks','nCount_RNA','nFeature_RNA','percent.mt'),
  pt.size = 0.1,
  ncol = 5
)
```
We can also check the number of cells before filtering 
```R
length(colnames(so_old))
[1] 20000
length(colnames(so_young))
[1] 20000

```
We can use some parameter proposed by Seurat, but we can also identify the threshold for each parameter on our own
Then we can set some thredshold to filtering data by these parameters. I will keep only the most 98% abundant of these parameters.
```R
#ATAC
low_prf <- quantile(so_young[["nCount_ATAC"]]$nCount_ATAC, probs = 0.02)
hig_prf <- quantile(so_young[["nCount_ATAC"]]$nCount_ATAC, probs = 0.98)
low_prp <- quantile(so_young[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02)
high_blr <- quantile(so_young[["blacklist_fraction"]]$blacklist_fraction, probs = 0.98)
hig_ns <- quantile(so_young[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)
low_ts <- quantile(so_young[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)

low_prf <- quantile(so_old[["nCount_ATAC"]]$nCount_ATAC, probs = 0.02)
hig_prf <- quantile(so_old[["nCount_ATAC"]]$nCount_ATAC, probs = 0.98)
low_prp <- quantile(so_old[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02)
high_blr <- quantile(so_old[["blacklist_fraction"]]$blacklist_fraction, probs = 0.98)
hig_ns <- quantile(so_old[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)
low_ts <- quantile(so_old[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)

#RNA
low_nCountRNA  <- quantile(so_young[["nCount_RNA"]]$nCount_RNA, probs = 0.02)
hig_nCountRNA <- quantile(so_young[["nCount_RNA"]]$nCount_RNA, probs = 0.98)
low_nFeatureRNA  <- quantile(so_young[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02)
hig_nFeatureRNA <- quantile(so_young[["nFeature_RNA"]]$nFeature_RNA, probs = 0.98)
low_mt  <- quantile(so_young[["percent.mt"]]$percent.mt, probs = 0.02)
hig_mt <- quantile(so_young[["percent.mt"]]$percent.mt, probs = 0.98)

low_nCountRNA  <- quantile(so_old[["nCount_RNA"]]$nCount_RNA, probs = 0.02)
hig_nCountRNA <- quantile(so_old[["nCount_RNA"]]$nCount_RNA, probs = 0.98)
low_nFeatureRNA  <- quantile(so_old[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02)
hig_nFeatureRNA <- quantile(so_old[["nFeature_RNA"]]$nFeature_RNA, probs = 0.98)
low_mt  <- quantile(so_old[["percent.mt"]]$percent.mt, probs = 0.02)
hig_mt <- quantile(so_old[["percent.mt"]]$percent.mt, probs = 0.98)
```
Now we can first vizualize our data after filtering 
```R
VlnPlot(
  object = subset(x=so_young,
                  subset = nCount_ATAC > low_prf &
                    nCount_ATAC < hig_prf &
                    pct_reads_in_peaks > low_prp &
                    blacklist_fraction < high_blr &
                    nucleosome_signal < hig_ns &
                    TSS.enrichment > low_ts &
                    nCount_RNA < hig_nCountRNA &
                    nCount_RNA > low_nCountRNA &
                    nFeature_RNA < hig_nFeatureRNA &
                    nFeature_RNA > low_nFeatureRNA &
                    percent.mt < hig_mt
  ),
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks','nCount_RNA','nFeature_RNA','percent.mt'),
  pt.size = 0.1,
  ncol = 5
)
```
Save our filtered data
```R
saveRDS(so_young_filtered, file="so_young_filtered.rds")
saveRDS(so_old_filtered, file="so_old_filtered.rds")
```
```
length(colnames(old))
[1] 17167
length(colnames(young))
[1] 17079
```
