# personal-genome-processing-pipeline
## Introduction
The inspiration for this WGS processing/analysis workflow came from receiving a clinical-grade whole genome sequence (Illumina NovaSeq) through the company Dante Labs. While Dante Labs provides its customers with alignments (as a BAM file) and variant calls (in the form of SNP, indel, and CNV VCF files), all of these processed file are aligned to GRCh37 (hg19), not the most recent human genome build (GRCh38, or hg38). While one could perhaps more simply use UCSC LiftOver to lift the hg19 coordinates over to hg38, this workflow was written both as a learning exercise, as well as to provide a means of greater control over the WGS processing/analysis workflow (combining reads from multiple runs, etc.).

### Variant calling pipeline
As a basic overview, the included **personal_genome_processing_pipeline.sh** script takes in paired-end Illumina reads (in FASTQ format), runs intiial read QC metrics using FASTQC, and then proceeds to map the paired reads to the most current genome build (GRCh38) using BWA-MEM. After some downstream processing, variants (SNPs and indels) are called using GATK's HaplotypeCaller. The resulting single VCF can optionally be split into separate SNP and indel files.

### Variant Annotation
The workflow for setting up and running GenomeChronicler is based off of the turotials presented by the Personal Genome Project UK (https://github.com/PGP-UK/GenomeChronicler) and Singularity (https://sylabs.io/guides/3.1/user-guide/). This tool takes in aligned reads (as a BAM file) and generates a PDF report highlighting the clinical significance of a panel of SNPs found in the proband's WGS data.

### Future Additions
Detailed genetic ancestry reporting, as well as copy number variant (CNV) calling are areas of interest for adding greater functionality to this workflow.

## Getting started
### Installing dependencies
The easiest way to install all of the required tools is via a package manager such as Anaconda (https://docs.conda.io/en/latest/miniconda.html). In this particular case, Miniconda3 (py37_4.8.3, linux-64) was used to install all of the following dependencies:

* fastqc (0.11.9)
* bwa (0.7.17)
* picard (2.23.9)
* gatk (4.1.9.0)
* java (OpenJDK) (11.0.6)
* R (4.0.3)
* ggplot2 (3.1.1)
* r-gplots (from Anaconda) (3.0.1)
* r-gsalib (from Anaconda) (2.1)
* samtools (1.11)
* vcftools (0.1.16)

### Downloading reference files
Most of the required reference files can be downloaded from teh Broad Institute's Google Cloud Bucket (https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/).

In total, one will need the following files:
* GRCh38 (hg38) reference genome FASTA (resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta)
* Mills and 1000 Genomes gold standard indels (resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf)
* 1000 Genomes Phase 1 high-confidence SNPs (resources-broad-hg38-v0-1000G_phase1.snps.high_confidence.hg38.vcf)
* Omni reference panel (resources-broad-hg38-v0-1000G_omni2.5.hg38.vcf)
* HapMap reference panel (resources-broad-hg38-v0-hapmap_3.3.hg38.vcf)

There are variable assignments in the script for ClinVar and dbSNP reference files, as well, but these are not actually used in this workflow.

## Running the pipeline script
Before running the **personal_genome_processing_pipeline.sh** script, one should set preffered paths for all of the variables declared at the top:

* MAIN_DIR - The main root directory to which all sub-directories and output files will be written to
* FASTQ1 - The first paired-end FASTQ file (forward reads)
* FASTQ2 - The second paired-end FASTQ file (reverse reads)
* OUT_PREF - A string containing the prefix with which all output files will be named
* FASTQC_OUT - The output directory to which fastqc results will be written
* REF_GENOME - Directory containing the reference genome FASTA
* READ_GROUPS - A string containing all of the read groups within your sequencing data. If unknown, for Dante Labs data, at least, one can easily find read groups by using the following command on one's hg19-aligned BAM file:
```
samtools view -H aligned_hg19.bam | grep '@RG
```
* TMP_DIR - A directory that will be generated to temporarily store intermediate files over the course of the workflow
* INTER_DIR - A directory to which intermediate files will be saved (such as pre-filtered BAMs and VCFs)
* MILLS, SNPs_1000G, OMNI, HAPMAP - Paths to the required reference files described above
* DBSNP, CLINVAR, and CLINVAR_WITH_CHR - Paths to the non-required and unused DBSNP and CLINVAR reference files. These variable paths can simply be left as-is.
* THREADS - The number of CPU threads that one wishes to use (if unsure how many to use, 1 is likely the safest option)
* RAM - The amount of RAM (in GB) that one wishes to use (if unsure, 1 is likely the safest option)

Once all variable paths are properly specified, one should be able to run the script, either via submission to a job scheduler (SLURM or LSF) or by simply calling the script directly in a Linux terminal:
```
./personal_genome_processing_pipeline.sh
```

If one is unsure of the suitability of their FASTQ read depth/quality, it would be best to run only the fastqc command initially and evaluate those results before deciding if the quality is high enough to merit downstream read mapping:
```
mkdir $FASTQC_OUT
fastqc -o $FASTQC_OUT $FASTQ1 $FASTQ2
```

Once this script is completely finished running, the most significant files that can be used for downstream analysis will be **$MAIN_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.recalibrated.bam** (BAM containing your aligned reads), **$MAIN_DIR\/$OUT_PREF\.recal.snp.indel.vcf.gz** (VCF cotaining SNPs and indels), **$MAIN_DIR\/$OUT_PREF\.snp.vcf.gz** (VCF containing SNPs only), and **$MAIN_DIR\/$OUT_PREF\.indel.vcf.gz** (VCF containing indels only). 

## Initial variant annotation using GenomeChronicler
