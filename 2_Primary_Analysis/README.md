# Primary analysis
## Raw data & Hierarchy 
Experiment protocol is the Chromium Next GEM Single Cell Multiome ATAC + Gene Expression, combining with HTO tag
After sequencing, normally, we'll have 2 folders, which is ATAC and GEX
In my **atac/**:

![image](https://github.com/user-attachments/assets/0e716b3e-a829-48a2-ae19-3d9088b2bc11)

In my **gex/**:

![image](https://github.com/user-attachments/assets/a0190beb-6ffa-4312-9bdd-f9b1fe97ea4f)

So I already have my file with fastq.gz file. Depend on which Cellranger pipeline you want to use, we can elaborate the correspond sample file
Since I decide to use directly the [Cellranger ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/tutorials/installation) pipeline, I will have 2 sample.csv 

Be sure to have the compatible reference **refdata-cellranger-arc-mm10-2020-A-2.0.0/** for mice 

Then now I can run my cellranger pipeline 
```
# Run cellranger-arc for the "old" sample
cellranger-arc count --id=old_output \
                     --reference=/beegfs/ttrinh/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                     --libraries=/beegfs/ttrinh/scatacseq/samplesheet_old.csv \
                   --localcores=16 \
                --localmem=990

# Run cellranger-arc for the "young" sample
cellranger-arc count --id=young_output \
                     --reference=/beegfs/ttrinh/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                     --libraries=/beegfs/ttrinh/scatacseq/samplesheet_young.csv \
                     --localcores=128 \
                     --localmem=990
```
This procedure will take about 24 hours of running and will take a enomous temporaly space so you can use some cluster to facilitate the task. 

after running, we will end up with a output folder, since I have 2 conditions, I have output_old and output_young 

![image](https://github.com/user-attachments/assets/f90b87cc-5e18-4ad7-bcda-fe65cfb87ac2)

We will only focus on the **outs/**

![image](https://github.com/user-attachments/assets/c69d960e-9f6b-4220-ab9c-dea80d0325cc)

So if we want to start from scratch for the downstream analysis, we can use **raw_feature_bc_matrix** but I will use directly filtered_feature_bc_matrix 

Another alternative to do is use **cellranger count** to preprocess each of our data atac and gex, then we continue separately analyse until we have the cluster, then we can perform the integration 

