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
old <- Read10X_h5(filename = "old_filtered_feature_bc_matrix.h5")

young <- Read10X_h5(filename = "young_filtered_feature_bc_matrix.h5")
```
#### 1.2.2 Files preparation
For ATAC seq data, we will need metadata file and the fragments files to be able using the other function of Signac. 
The fragment file must be stored in the same folder with .tbi file to be able map the position 

```
old_meta <- read.csv(
  file = 'old_per_barcode_metrics.csv',
  header = TRUE,
  row.names = 1)

young_meta <- read.csv(
  file = 'young_per_barcode_metrics.csv',
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
```

We can also store the path to our fragment files to facilitate the code writing after
```
old_fragpath = "/path/to/old_atac_fragments.tsv.gz"
young_fragpath = "/path/to/young_atac_fragments.tsv.gz"
```

#### 1.2.3 Create Seurat object
Create a Seurat object containing the RNA data first :
```
so_young <-  CreateSeuratObject(
  counts = young_rna,
  assay = "RNA",
  meta.data = young_meta
)

so_old <-  CreateSeuratObject(
  counts = young_rna,
  assay = "RNA",
  meta.data = old_meta
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
So now we can move to next part : workflow for ATAC&RNA 

