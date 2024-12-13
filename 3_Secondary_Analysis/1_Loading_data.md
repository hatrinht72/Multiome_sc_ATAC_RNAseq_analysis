## 1. Loading library 
Make sure all packages are installed, 
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
## 2. Loading files
### 2.1 Loading count files
We will use directly filtered matrix created by cellranger arc.
For RNA data, its correspond to genes/cell matrix
For ATAC data, its correspond to peak/cell matrix 
```
old <- Read10X_h5(filename = "old_filtered_feature_bc_matrix.h5")
old_rna<- old$`Gene Expression`
old_atac <- old$Peaks
young <- Read10X_h5(filename = "young_filtered_feature_bc_matrix.h5")
young_rna <- young$`Gene Expression`
young_atac <- young$Peaks
```
### 2.2 Loading meta files
For ATAC seq data, we ll need cell metadata file and the fragments files to be able using the other functio of Signac. 
As we already talk in **2_Primary_Analysis**, this fragment file must be stored in the same folder with .tbi file to be able map the position 
```
old_meta<- read.csv(
  file = 'old_barcode_metrics.csv',
  header = TRUE,
  row.names = 1)
young_meta<- read.csv(
  file = 'young_per_barcode_metrics.csv',
  header = TRUE,
  row.names = 1)
```
We can also store the path to our fragment files to facilitate the code writing after
```
old_fragpath = "/path/to/old_atac_fragments.tsv"
young_fragpath = "/path/to/young_atac_fragments.tsv"
```
### 2.3 Create Assay
Create a Seurat object containing the RNA data first :
```
so_young <-  CreateSeuratObject(
  counts = young_rna,
  assay = "RNA"
)
so_old <-  CreateSeuratObject(
  counts = young_rna,
  assay = "RNA"
)
```
We also need annotation for the chromatin assay
```
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
###Convert annotation seqlevels to UCSC style (i.e., prefix with "chr")
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
###Set the genome to mm10
genome(annotations) <- "mm10"
```
Then add ATAC assay to the object 
```
so_young[["ATAC"]] <- CreateChromatinAssay(
  counts = young_atac,
  sep = c(":", "-"),
  fragments = young_fragpath,
  annotation = annotations
)

so_old[["ATAC"]] <- CreateChromatinAssay(
  counts = old_atac,
  sep = c(":", "-"),
  fragments = old_fragpath,
  annotation = annotations
)
```
So now we can move to next part : workflow for ATAC&RNA 

