## 2. scATACseq analysis workflow
Since in our Seurat object now we have 2 assays so to do control quality for ATAC assay, we have to choose the default assay
```
DefaultAssay(so_young) <- "ATAC"
DefaultAssay(so_old) <- "ATAC"
```
### 2.1 QC metrics
To filter scATACseq data, there is several parameters to look after :
- Nucleosome bandng pattern : to make sure we have the right DNA fragment size, corresponding to the length of the DNA wrapped around a single nucleosome (147 pb) and to assess the signal from nucleosome positioning based on the fragment length distribution from the ATAC-seq data. The nucleosome_signal will store the value of the appromimate ratio of mononucleosomam and nucleosome-free fragments
- Transcriptional start site (TSS) enrichnment score : to assess the enrichment of transcription start sites (TSS) in ATAC-seq data, calculates TSS enrichment scores, which help determine the degree to which TSS regions are accessible in your chromatin accessibility data. Poor scATAC seq experiments would have low TSS enrichment score
- Total number of fragments in peaks : mesure depth sequencing; identify the number of fragments (or reads) that map to specific accessible regions (peaks)
- Fraction of reads in peaks : all fragments that fall within ATAC-seq peaks, access the total number of reads that passed quality control filters during preprocessing
- Ratio reads in genomic blacklist regions : reads which are often represent low-quality cells/technical artifacts. Access the proportion of reads that map to regions defined as "blacklisted.", consisting problematic genomic regions, often due to high background noise or sequencing artifacts, and should be considered with caution when interpreting results

```
so_young <- NucleosomeSignal(object = so_young) #fragment ratio 147-294: <147
so_old <- NucleosomeSignal(object = so_old) #fragment ratio 147-294: <147

so_young <- TSSEnrichment(object = so_young, fast = FALSE)
so_old <- TSSEnrichment(object = so_old, fast = FALSE)

so_young$blacklist_ratio <- so_young$blacklist_region_fragments / so_young$peak_region_fragments
so_old$blacklist_ratio <- so_old$blacklist_region_fragments / so_old$peak_region_fragments
old_data_atac$atac_peak_region_fragments 

so_young$pct_reads_in_peaks <- so_young$peak_region_fragments / so_young$passed_filters * 100 
so_oldg$pct_reads_in_peaks <- so_oldg$peak_region_fragments / so_oldg$passed_filters * 100 
```

we can then visualize all these parameter by violin plot.

```
VlnPlot(
  object = so_young,
  features = c('atac_peak_region_fragments', 'pct_reads_in_peaks', 
               'blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'),
  pt.size = 0.1,
  ncol = 5
)

VlnPlot(
  object = so_old,
  features = c('atac_peak_region_fragments', 'pct_reads_in_peaks', 
               'blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'),
  pt.size = 0.1,
  ncol = 5
)
```

