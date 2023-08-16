# ASD_STG_LCM_RNAseq_paper_code
This is the repository that contains code associated with the paper "**Neuron-specific transcriptomic signatures indicate neuroinflammation and altered neuronal activity in ASD temporal cortex**".

The raw RNAseq data is deposited in dbGAP with accession number phs003208.

The "script" folder contains code used to analyze the data.

The "Supplementary Datasets" folder contains all supplemental tables.

The "counts_data" folder contains gene-level counts data from RNAseq. The data is in a [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment) object. This object contains both counts data and counts per million (CPM) data in the *assays* slot, as well as metadata in the *colData* slot.
