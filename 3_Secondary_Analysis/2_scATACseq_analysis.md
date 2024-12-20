## 2. scATACseq analysis workflow
Since in our Seurat object now we have 2 assays so to do control quality for ATAC assay, we have to choose the default assay
```
DefaultAssay(old) <- "ATAC"
```
### 2.1 QC metrics
To filter scATACseq data, there is several parameters to look after :
- Nucleosome bandng pattern : to make sure we have the right DNA fragment size, corresponding to the length of the DNA wrapped around a single nucleosome (72 pb). The nucleosome_signal will store the value of the appromimate ratio of mononucleosomam and nucleosome-free fragments
- Transcriptional start site (TSS) enrichnment score : poor ATAC seq experiments would have low TSS enrichment score
- Total number of fragments in peaks : mesure depth sequencing
- Fraction of reads in peaks : all fragments that fall within ATAC-seq peaks 
- Ratio reads in genomic blacklist regions : reads which are often represent low-quality cells/technical artifacts

  The three last data could be found in meta data
  
