### 2. ATAC sequencing data analysis
I will follow the [workflow of Signac](https://stuartlab.org/signac/articles/pbmc_vignette)
#### 2.1 Normalization and linear dimensional reduction 
```
# we can change the cutoff to be faster (here is 100% peak will be use but if we set q75, only 25% will be use
so_young <- so_young_filtered
so_old <- so_old_filtered

DefaultAssay(so_young) <- "ATAC"
DefaultAssay(so_old) <- "ATAC"

so_young  <- RunTFIDF(so_young)
so_young <- FindTopFeatures(so_young, min.cutoff = 'q0') 
so_young <- RunSVD(so_young)

so_old  <- RunTFIDF(so_old)
so_old <- FindTopFeatures(so_old, min.cutoff = 'q0') 
so_old <- RunSVD(so_old)
```
We can then access the correlation between each LSI (you can consider as PCA) component
```
png(paste0(OUT_DIR, "2_1_ATAC_depthcor_old.png"), width = 1200, height = 1000, res = 150)
DepthCor(so_old)
dev.off()
```
![2_1_ATAC_depthcor_old](https://github.com/user-attachments/assets/86e3a00e-8603-4dfb-b6e6-0e934d1afb17)

```
png(paste0(OUT_DIR, "2_1_ATAC_depthcor_young.png"), width = 1200, height = 1000, res = 150)
DepthCor(so_young)
dev.off()
```
![2_1_ATAC_depthcor_young](https://github.com/user-attachments/assets/76f89c4f-3e23-41d4-9fda-b060ae895b02)
The purpose of this that we can choose the dimension that we want to perform our graph-based clustering and non-linear dimension reduction. 
Since the first component is strongly correlated with the total number of counts, we will not include it in our analysis.

#### 2.2  Non-linear dimension reduction and clustering 

We process similar steps in scRNA seq analysis 
When we process the runUMAP, we can set a seed that we will use to have the consistancy across 
When we process the FindClusters, we can choose the resolution that we need, the higher the resolution, the more clusters you will get.
```
so_young <- RunUMAP(object = so_young, reduction = 'lsi', seed.use = 1990, dims = 2:15, verbose = FALSE) %>% FindNeighbors(reduction = 'lsi', dims = 2:15)
so_young <- FindClusters(object = so_young, verbose = FALSE, algorithm = 3, resolution=c(0.2, 0.4,0.6,0.8))

png(paste0(OUT_DIR, "2_2_ATAC_umap_young.png"), width = 2500, height = 2000, res = 150)
p1=DimPlot(so_young,  reduction = "umap", group.by="ATAC_snn_res.0.2", label=TRUE)
p2=DimPlot(so_young,  reduction = "umap", group.by="ATAC_snn_res.0.4", label=TRUE) 
p3=DimPlot(so_young,  reduction = "umap", group.by="ATAC_snn_res.0.6", label=TRUE) 
p4=DimPlot(so_young,  reduction = "umap", group.by="ATAC_snn_res.0.8", label=TRUE) 
p1 + p2 + p3 + p4
dev.off()
```
![2_2_ATAC_umap_young](https://github.com/user-attachments/assets/e923e25c-302d-4083-a92d-05542c4feabe)
```
so_old <- RunUMAP(object = so_old, reduction = 'lsi', seed.use = 1990, dims = 2:15, verbose = FALSE, reduction.name = "ATAC_umap") %>% FindNeighbors(reduction = 'lsi', dims = 2:15)
so_old <- FindClusters(object = so_old, verbose = FALSE, algorithm = 3, resolution=c(0.2, 0.4,0.6,0.8))

png(paste0(OUT_DIR, "2_2_ATAC_umap_old.png"), width = 2500, height = 2000, res = 150)
p1=DimPlot(so_old,  reduction = "ATAC_umap", group.by="ATAC_snn_res.0.2", label=TRUE)
p2=DimPlot(so_old,  reduction = "ATAC_umap", group.by="ATAC_snn_res.0.4", label=TRUE) 
p3=DimPlot(so_old,  reduction = "ATAC_umap", group.by="ATAC_snn_res.0.6", label=TRUE) 
p4=DimPlot(so_old,  reduction = "ATAC_umap", group.by="ATAC_snn_res.0.8", label=TRUE) 
p1 + p2 + p3 + p4
dev.off()
```
![2_2_ATAC_umap_old](https://github.com/user-attachments/assets/19fcd3ce-7da5-49d3-a4c0-eec9dba93673)

#### 2.3 Integration 2 conditions 
```
so_young$condition <- young
so_old$condition <- old
```

We can now merge 2 data sets
```
so_merge <- merge(so_young, so_old, add.cell.ids = c("young", "old"))
so_merge <- 
```
We can now rerun the comined seurfile:///home/ha/y1/projects/Multiome_sc_ATAC_RNAseq_analysis/3_Secondary_Analysis/2_3_combined_ATAC_umap.png
at object to see if there are some difference between 2 conditions
```
so_merge <- FindTopFeatures(so_merge, min.cutoff = "q75")
so_merge <- RunTFIDF(so_merge)
so_merge <- RunSVD(so_merge)
so_merge <- RunUMAP(so_merge, reduction = "lsi", dims = 2:30, reduction.name = "combine_atac_umap")

png("2_3_combined_ATAC_umap.png"), width = 2500, height = 2000, res = 150)
DimPlot(so_merge, group.by = "condition")
dev.off()
```
![2_3_combined_ATAC_umap](https://github.com/user-attachments/assets/10084d05-dc8b-4fa3-aa3e-67494da8d82a)

So here we can see that there are a clear difference between 2 conditions old and young, that's why we have to process next the integration to correct batch effect 

