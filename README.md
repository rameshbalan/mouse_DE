# Differential Expression Analysis on Mouse Tissue Samples

## Basic Workflow

1. FastQC and MultiQC.
2. Index the Reference Transcriptome.
	a. It was downloaded from the this ncbi [link](ftp://ftp.ncbi.nih.gov/genomes/M_musculus/RNA/)
3. Quantify each sample
	a. I used salmon to quantify the samples. The script is available in src named "quantification.sh"
4. DE Analysis
	a. Create a metadata table. This table describes the experimental design.
	b. Run DESeq() in R