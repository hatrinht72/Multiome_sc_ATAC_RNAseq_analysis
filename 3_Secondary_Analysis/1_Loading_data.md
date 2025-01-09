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

