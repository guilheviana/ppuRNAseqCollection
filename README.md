# A compilation of publicly available _P. putida_ expression datasets analyzed using DESeq2

This repository contains data derived from several _Pseudomonas putida_ expression studies available in NCBI's [SRA database](https://www.ncbi.nlm.nih.gov/sra) and scripts for processing them using an unified pipeline. 
This analysis is part of an upcoming publication comparing the usage of different _P. putida_ genome annotations across the literature.

## Main file structure
- `data/`: Contains information on the sequencing runs available in SRA for each Bioproject dataset analyzed in this study. This data is used for batch-processing files using scripts `run-rsem.R` and `run-trimmomatic.R` etc. 
- `metadata/`: Contains additional metadata retrieved from the Gene Expression Omnibus database for selected datasets, needed for the analysis where the contents in `data/`alone are not sufficient.
- `output/`: Contains tables for converting condition names to formatted condition names for replicating the manuscript figures.
- `rsem/`: Contains folders with transcript pseudocounts calculated by RSEM using either the GenBank or the RefSeq versions of the reference _P. putida_ KT2440 transcriptome. It is the main results folder for the differential expression analysis carried out in this project.

### Scripts and data

- `affected_newformat` and `affected_oldformat` are R vectors, stored as binary files, with the RefSeq and GenBank genomic loci codes of the 897 shifted genes described in the main text, respectively.
- `genbank_processed.RDS` and `refseq_processed.RDS` are the final expression datasets processed using GenBank and RefSeq annotations as references. More conveniently provided as text files in the supplemental material of the manuscript.
- `prepare-loci-info.R`, `process-rsem.R` are scripts for parsing RSEM output files and running them through DESeq2 for differential expression analysis.
- `figures.R`and `plot_utils.R` are scripts for generating the figures present in the final manuscript.
- `utils.R` contains handy functions used across different scripts in this project.
- `download-software.R`, `download-sra.R`, `run-trimmomatic.R`, and `run-rsem.R` are wrappers around bash calls for downloading and processing FASTQ files. Please see additional information below.
  
## Additional information and contact
Memory-intensive steps in the pipeline (notably FASTQ files download and processing, and alignments using RSEM) were carried out using a server running Ubuntu 22.04.5 LTS equipped with an Intel(R) Xeon(R) Bronze 3104 CPU and 512 GiB RAM memory.
If you have any questions or feedback, please contact the first author Guilherme (Gui) Viana de Siqueira at guiviana@proton.com.
