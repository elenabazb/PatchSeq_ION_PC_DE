# Input data

The data repository is currently under preparation.



The sequencing count matrices and cell metadata are maintained separately
from this analysis repository.

The files are available from:

- Data repository: [repository link]
- Release/version: [version or DOI]


## Required files

Place the following files in this directory before running the analysis:

- `ionraw_counts.RDS`
- `pcraw_counts.RDS`
- `metadata_ion_export.RDS`
- `metadata_pc_export.RDS`


## Expected structure

### ION count matrix

- Rows: genes
- Columns: cells
- Values: raw read counts

### PC count matrix

- Rows: genes
- Columns: cells
- Values: raw read counts


### Metadata

The metadata tables must contain at least:

- `Explant N`: sample ID.
- `Library`: cell ID.
- `Age_dis`: discrete age variable. "Young" corresponds to DIV12-13. "Old" corresponds to DIV28-30.
- `celltype`: Purkinke cell (PC) or Inferior olivary neuron (ION).