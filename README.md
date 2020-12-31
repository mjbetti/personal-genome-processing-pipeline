# personal-genome-processing-pipeline
## Introduction
The inspiration for this processing workflow came from receiving a clinical-grade whole genome sequence (sequenced on the Illumina NovaSeq platform) through the company Dante Labs. While Dante Labs provides its customers with alignments (as a BAM file) and variant calls (in the form of SNP, indel, and CNV VCF files), all of these processed file are aligned to GRCh37 (hg19), not the most recent human genome build (GRCh38, or hg38). While one could perhaps more simply use UCSC LiftOver to lift the hg19 coordinates over to hg38, this workflow was written both as a learning exercise, as well as to provide a means of greater control over the WGS processing/analysis workflow (combining reads from multiple runs, etc.).

### Variant calling pipeline
As a basic overview, the included **personal_genome_processing_pipeline.sh** script takes in paired-end Illumina reads (in FASTQ format), runs intiial read QC metrics using FASTQC, and then proceeds to map the paired reads to the most current genome build (GRCh38) using BWA-MEM. After some downstream processing, variants (SNPs and indels) are called using GATK's HaplotypeCaller. The resulting single VCF can optionally be split into separate SNP and indel files.

### Variant Annotation
The workflow for setting up and running GenomeChronicler is based off of the turotials presented by the Personal Genome Project UK (https://github.com/PGP-UK/GenomeChronicler) and Singularity (https://sylabs.io/guides/3.1/user-guide/). This tool takes in aligned reads (as a BAM file) and generates a PDF report highlighting the clinical significance of a panel of SNPs found in the proband's WGS data.

### Future Additions
Detailed genetic ancestry reporting, as well as copy number variant (CNV) calling are areas of interest for adding greater functionality to this workflow.

## Getting started
