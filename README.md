# README

**By Sarah Salisbury and Pablo Carballo Pacoret**

&#x1F4D8; This repository contains the scripts used to analyze snRNAseq libraries generated from retina samples from the shark (_Scyliorhinus canicula_) for the manuscript titled "A single-nucleus RNA sequencing atlas of the postnatal retina of the shark (_Scyliorhinus canicula_) ".

## Citation

If these scripts are helpful to you please consider citing the accompanying manuscript:

(ADD ONCE AVAILABLE)

## Repository Scripts

This repository contains 4 scripts:
- &#128195; `1_STAR.md` <- outlines how to map libraries to genome using STAR and generate feature counts per cell matrices for Seurat analyses
  
- &#128195; `2_UP_TO_FIRST_UMAP.Rmd` - R Markdown script to do initial QC and cell clustering
  
- &#128195; `3_FIND_MARKERS.Rmd` - R Markdown script to identify marker genes for each cell type
  
- &#128195; `4_REMOVING_BAD_CLUSTERS.Rmd` - R Markdown script to remove dubious clusters after initial clustering and then recluster remaining cells
