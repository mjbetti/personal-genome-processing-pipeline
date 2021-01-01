#Installation of the following commandline tools is required: fastqc, bwa, picard, gatk, java (OpenJDK), ggplot2, r-gplots, r-gsalib, samtools, and vcftools. The easiest way to install all of these is via Anaconda.

#The GRCh38 Reference Sequence used here was downloaded from https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
#dbSNP file downloaded from dbSNP FTP (https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-All.vcf.gz)
#Mills, 1000G, omni, and HapMap files downloaded from Broad FTP (https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)
#Clinvar downloaded from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20200905.vcf.gz

#Declare all files and directories that will be used for inputs and outputs. Clinvar reference files can be used for downstream annotations with tools like ANNOVAR, but we will not use them here.
MAIN_DIR=/home/mbetti/personal_genome
FASTQ1=$MAIN_DIR\/dante_labs_downloads/60820188479382_SA_L001_R1_001.fastq.gz
FASTQ2=$MAIN_DIR\/dante_labs_downloads/60820188479382_SA_L001_R2_001.fastq.gz
OUT_PREF=mb_hg38_60820188479382
FASTQC_OUT=$MAIN_DIR\/fast_qc
REF_GENOME=$MAIN_DIR\/ref_files/hg38_ref/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta
READ_GROUPS="@RG\tID:1\tLB:LB0\tPL:PL0PU:PU0\tSM:60820188479382"
TMP_DIR=$MAIN_DIR\/tmp_dir
INTER_DIR=$MAIN_DIR\/intermediate_dir
DBSNP=$MAIN_DIR\/ref_files/dbsnp/00-All.vcf
MILLS=$MAIN_DIR\/ref_files/mills/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf
SNPS1000G=$MAIN_DIR\/ref_files/1000g/resources-broad-hg38-v0-1000G_phase1.snps.high_confidence.hg38.vcf
OMNI=$MAIN_DIR\/ref_files/omni/resources-broad-hg38-v0-1000G_omni2.5.hg38.vcf
HAPMAP=$MAIN_DIR\/ref_files/hapmap/resources-broad-hg38-v0-hapmap_3.3.hg38.vcf
CLINVAR=$MAIN_DIR\/ref_files/clinvar/clinvar_20200905.vcf
CLINVAR_WITH_CHR=$MAIN_DIR\/ref_files/clinvar/clinvar_20200905_with_chr.vcf

#Declare variables specifying thread and RAM usage (in GB)
THREADS=4
RAM=32

#First run fastqc on both raw FASTQ files
mkdir $FASTQC_OUT
fastqc -o $FASTQC_OUT $FASTQ1 $FASTQ2

#Index the downloaded refrence genome using BWA
bwa index -a bwtsw $REF_GENOME

#Align the FASTQ files using BWA-MEM
#For Dante Labs data, the read group(s) of your sample can be determined by running the following command with your hg19 aligned BAM file:
#   samtools view -H sample.bam | grep '@RG'
mkdir $INTER_DIR

bwa mem -t $THREADS -T 0 -R $READ_GROUPS $REF_GENOME $FASTQ1 $FASTQ2 | samtools view -Shb -o $INTER_DIR\/$OUT_PREF\.bam

#Sort the aligned BAM file using Picard tools
mkdir $TMP_DIR
ulimit -c unlimited
picard SortSam "-Xmx"$RAM\g  \
    CREATE_INDEX=true \
    INPUT=$INTER_DIR\/$OUT_PREF\.bam \
    OUTPUT=$INTER_DIR\/$OUT_PREF\.sorted.bam \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=STRICT \
    TMP_DIR=$TMP_DIR

rm $TMP_DIR\/*

#Merge reads from multiple runs
ulimit -c unlimited
picard MergeSamFiles "-Xmx"$RAM\g \
    ASSUME_SORTED=false \
    CREATE_INDEX=true \
    INPUT=$INTER_DIR\/$OUT_PREF\.sorted.bam \
    MERGE_SEQUENCE_DICTIONARIES=false \
    OUTPUT=$INTER_DIR\/$OUT_PREF\.sorted.merged.bam \
    SORT_ORDER=coordinate \
    USE_THREADING=true \
    VALIDATION_STRINGENCY=STRICT

samtools sort -m $(($RAM / 4))\g -@ $THREADS \
    -o $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.bam \
    $INTER_DIR\/$OUT_PREF\.sorted.merged.bam

#Identify duplicate reads originating from the same single fragment of DNA, i.e. PCR artifacts
export _JAVA_OPTIONS=-Djava.io.tmpdir=$TMP_DIR
ulimit -c unlimited
gatk MarkDuplicatesSpark \
    -I $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.bam \
    -O $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.bam \
    -M $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates_metrics.txt

#Generate an index and dictionary for the reference genome for use in gatk processing
gatk CreateSequenceDictionary \
    --REFERENCE $REF_GENOME

samtools faidx $REF_GENOME

#Add chromosome names to the clinvar reference VCF before indexing and output as a new VCF file
#awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $CLINVAR > $CLINVAR_WITH_CHR

#Before generating the table, the VCF files will need corresponding indices generated
gatk IndexFeatureFile --input $DBSNP
gatk IndexFeatureFile --input $MILLS
gatk IndexFeatureFile --input $SNPS1000G
#gatk IndexFeatureFile --input $CLINVAR_WITH_CHR
gatk IndexFeatureFile --input $OMNI
gatk IndexFeatureFile --input $HAPMAP

#Generate a recalibration table for Base Quality Score Recalibration
gatk BaseRecalibrator \
    -I $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.bam \
    -R $REF_GENOME \
    --known-sites $DBSNP \
    --known-sites $MILLS \
    --known-sites $SNPS1000G \
    -O $INTER_DIR\/recal_data.table

#Apply the calculated base quality score recalibration
gatk ApplyBQSR \
    -R $REF_GENOME \
    -I $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.bam \
    --bqsr-recal-file $INTER_DIR\/recal_data.table \
    -O $MAIN_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.recalibrated.bam

#Evaluate the base quality score recalibration tables
gatk AnalyzeCovariates -bqsr $INTER_DIR\/recal_data.table -plots $INTER_DIR\/AnalyzeCovariates_mb.pdf

#Somatic short variant discovery (SNVs + Indels)
#Use HaplotypeCaller to call germline SNPs and indels via local re-assemply of haplotypes. The output will be an intermediate GVCF file, containing raw, unfiltered SNP and indel calls
ulimit -c unlimited
export _JAVA_OPTIONS=-Djava.io.tmpdir=$TMP_DIR
gatk --java-options "-Xmx"$RAM\g HaplotypeCaller \
    -R $REF_GENOME \
    -I $MAIN_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.recalibrated.bam \
    -O $INTER_DIR\/$OUT_PREF\.g.vcf.gz \
    -ERC GVCF

#Perform genotyping on one or more samples pre-called with HaplotypeCaller
ulimit -c unlimited
export _JAVA_OPTIONS=-Djava.io.tmpdir=$TMP_DIR
gatk --java-options "-Xmx"$RAM\g GenotypeGVCFs -R $REF_GENOME -V $INTER_DIR\/$OUT_PREF\.g.vcf.gz -O $INTER_DIR\/$OUT_PREF\.vcf.gz

#Build a SNP recalibration model to score variant quality and then apply it to filter the SNPS in the generated VCF
#Referred to WGS recalibration code on (https://github.com/BD2KGenomics/gatk-whole-genome-pipeline/blob/master/HAPvariantCalling.sh)
gatk VariantRecalibrator \
    -R $REF_GENOME \
    -V $INTER_DIR\/$OUT_PREF\.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 $OMNI \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $SNPS1000G \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP\
    -an QD \
    -an DP \
    -an FS \
    -an ReadPosRankSum \
    --mode SNP \
    --output $INTER_DIR\/$OUT_PREF\.snp.recal \
    --tranches-file $INTER_DIR\/$OUT_PREF\.snp.tranches \
    --rscript-file $INTER_DIR\/$OUT_PREF\.snp.plots.R

gatk ApplyVQSR \
    -R $REF_GENOME \
    -V $INTER_DIR\/$OUT_PREF\.vcf.gz \
    --output $INTER_DIR\/$OUT_PREF\.recal.snp.vcf.gz \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file $INTER_DIR\/$OUT_PREF\.snp.tranches \
    --recal-file $INTER_DIR\/$OUT_PREF\.snp.recal \
    -mode SNP

#Perform a second round of recalibration, this time focusing on the indels
gatk VariantRecalibrator \
    -R $REF_GENOME \
    -V $INTER_DIR\/$OUT_PREF\.recal.snp.vcf.gz \
    --resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
    -an DP \
    -an FS \
    -an ReadPosRankSum \
    --mode INDEL \
    --output $INTER_DIR\/$OUT_PREF\.snp.indel.recal \
    --tranches-file $INTER_DIR\/$OUT_PREF\.snp.indel.tranches \
    --rscript-file $INTER_DIR\/$OUT_PREF\.snp.indel.plots.R

gatk ApplyVQSR \
    -R $REF_GENOME \
    -V $INTER_DIR\/$OUT_PREF\.vcf.gz \
    --output $MAIN_DIR\/$OUT_PREF\.recal.snp.indel.vcf.gz \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file $INTER_DIR\/$OUT_PREF\.snp.indel.tranches \
    --recal-file $INTER_DIR\/$OUT_PREF\.snp.indel.recal \
    -mode INDEL

#Split the final VCF output into two separate SNP and indel VCF files
vcftools \
--gzvcf $MAIN_DIR\/$OUT_PREF\.recal.snp.indel.vcf.gz \
--remove-indels \
--recode \
--recode-INFO-all \
--stdout | gzip -c >  $MAIN_DIR\/$OUT_PREF\.snp.vcf.gz

vcftools \
--gzvcf $MAIN_DIR\/$OUT_PREF\.recal.snp.indel.vcf.gz \
--keep-only-indels \
--recode \
--recode-INFO-all \
--stdout | gzip -c > $MAIN_DIR\/$OUT_PREF\.indel.vcf.gz
