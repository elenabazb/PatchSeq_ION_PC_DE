# PatchSeq differential expression analysis - developing Purkinje cells and Inferior Olivary neurons.

This repository contains the R code and functional enrichment results used to analyze 
age-associated transcriptional changes in mouse inferior olivary neurons (IONs) and Purkinje cells (PCs) profiled using Patch-seq.

The analysis compares two developmental groups, referred to as Young and Old in the metadata, separately within each cell type.

## Contents 
- `R script`: analysis script.
- `environment`: R version, operating system and package versions used for the analysis.
- `data`: data availability (raw files, processed count matrices and metadata).
- `results`: Gene Ontology enrichment tables exported from Enrichr for differentially expressed genes in IONs and PCs.


## Citation
This script was used to generate analysis and figures in the following publication: [under preparation].


## Analysis overview
The analysis workflow includes:

1. Import of raw gene count matrices and cell metadata.
2. Calculation of quality-control metrics.
3. Filtering of cells according to the quality-control thresholds defined in the analysis script.
4. Log-normalization of expression values using Seurat.
5. Visualization of ION and PC marker-gene expression.
6. Differential expression analysis between age groups using DESeq2.
7. Extraction of upregulated and downregulated gene sets.
8. Gene Ontology enrichment analysis using Enrichr.
9. Visualization of enriched pathways and selected differentially expressed genes.
10. Generation of dot plots, violin plots, heatmaps and volcano plots.



## Data availability 
The count matrices and cell metadata required to reproduce this analysis will be available from the associated data repository:

Data repository: E-MTAB-17525
