This is a pipeline to process a specific dataset within the research group that employes the author. 

The main steps include:
- Merging vcfs and assessment (VQSR, custom R scripts and samtools vcfstats) 
- Additional variants filtering by VQSR and a set of hard filters 
- Variants annotation by VEP
- Export of annotated variants to plain text fils for downstream analysis in R.

This repository is intended for the author's pesonal use.  
However, as long as it is kept public, everyone is welcome to have a look. 