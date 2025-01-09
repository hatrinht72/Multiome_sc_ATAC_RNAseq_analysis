## 2. Filtering ATAC and RNA
We can first clean our data then analyze independently ATAC and RNA.
We can start with ATAC 
### 2.1 ATAC metrics
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
VlnPlot(
  object = so_young,
  features = c('atac_peak_region_fragments',  'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)
```
![2_1_ATAC_QC_metrics_young](https://github.com/user-attachments/assets/f76bab2a-a3ba-48e0-ae07-7730128a2c49)

```
VlnPlot(
  object = so_old,
  features = c('atac_peak_region_fragments',  'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)
```
we can also visualize the relationship between variables.
```
png(paste0(OUT_DIR, "2_1_ATAC_relation_count_tss_young.png"), res = 150,width = 1000, height = 1000 )
DensityScatter(so_young, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
```
![2_1_ATAC_relation_count_tss_young](https://github.com/user-attachments/assets/28ad481f-9973-4d1a-99a7-bb59721e3ab5)

### 2.2 RNA metrics
We can also visualize the RNA data.
```
png(paste0(OUT_DIR, "2_2_QC_RNA_young.png"), width = 1200, height = 1000, res = 150)
VlnPlot(so_young, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()
```
![2_2_RNA_QC_young](https://github.com/user-attachments/assets/8ae6a0e1-4df1-4709-bf72-42f5a1c4b91a)

```
png(paste0(OUT_DIR, "2_2_QC_RNA_count_feature_young.png"), width = 1200, height = 1000, res = 150)
FeatureScatter(so_young, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
```
![2_2_RNA_QC_count_feature_young](https://github.com/user-attachments/assets/2f3de10f-35d2-49ba-b2e2-3f7ac146304b)

### 2.3 Filtering 
#### 2.3.1 Set the thredshold
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
#### 2.3.1 Filter 
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
