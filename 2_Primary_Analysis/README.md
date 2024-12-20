# Primary analysis
## Raw data & Hierarchy 
Experiment protocol is the Chromium Next GEM Single Cell Multiome ATAC + Gene Expression, combining with HTO tag
After sequencing, normally, we'll have 2 folders, which is ATAC and GEX
In my **atac/**:
![image](https://github.com/user-attachments/assets/5dc6b3cf-905f-418d-a07e-13b6eca49c8a)
In my **gex/**:
![image](https://github.com/user-attachments/assets/a79992ac-6bef-4ed8-8850-6621e30b618c)
So I already have my file with fastq.gz file. Depend on which Cellranger pipeline you want to use, we can elaborate the correspond sample file
Since I decide to use directly the Cellranger ARC pipeline, I will have 2 sample_file correspond to 

they re 2 ways to do
cellranger arc
cellranger count
