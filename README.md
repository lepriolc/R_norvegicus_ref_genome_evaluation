# Rattus norvegicus reference genome evaluation

## Content

This repository contains all the scripts to reproduce the results published in the article entitled "Rattus norvegicus reference genome evaluation for hippocampus RNA-seq data analysis".

## Code

### SRA to Fastq conversion

SRA files of the samples which compose the GSE136913, GSE140420 and GSE179101 datasets must be downloaded the NCBI Gene Expression Omnibus and placed in the subdirectory ./00-Raw_data/RNAseq/, one subdirectory per dataset. The 'SraRunTable.txt' and 'sra_result.csv' must also be downloaded.

The `sra2fastq_conversion.sh` script converts SRA files into Fastq files. SRA files are removed after conversion.

### Reference genomes

The Fasta and the GTF files of the reference genomes must be downloaded (see article for more details) and placed in the subdirectory ./00-Raw_data/RNAseq/reference_genomes/.

Genome indexes are built using the following command: STAR --runThreadN 1 --runMode genomeGenerate --genomeDir ./00-Raw_data/RNAseq/reference_genomes/ --genomeFastaFiles <FASTA_FILE> --sjdbGTFfile <GTF_FILE>

### Expression quantification

The `exp_quantif_pipeline_launcher.sh` script aligns the reads to a reference genome and counts the number of aligned reads in annotated genomic regions for all the samples of a dataset.

### Count matrix

The `create_exp_matrix.R` script reads the outputs of `exp_quantif_pipeline_launcher.sh` script and generates a count matrix stored in a .RData file.

### Consistency and stringency

The `consistency_stringency.R` script performs the comparison of all pairs of reference genomes using the consistency and stringency metrics.

### Differential expression analysis and Gene Ontology term enrichment analysis

The `DE_analysis.R` script performs differential expression analysis and Gene Ontology term enrichment analysis. Specific designs and contrasts are used depending on the analyzed dataset.

### Specific Gene Ontology term analysis

The `DE_analysis-specific_enriched_GO_terms.R` script performs the analysis of genome-
specific differentially expressed genes related to specifically enriched Gene Ontology terms depending on the
reference genome used for the differential expression analysis.

## Contact

For more details about the analysis, please refer to the article entitled 'Rattus norvegicus reference genome evaluation for hippocampus RNA-seq data analysis'

If you have any question, please contact chris.lepriol@gmail.com.

