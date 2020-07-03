#!/bin/bash

# Load modules
module load tools parallel/20200522 anaconda2/4.4.0 hisat2/2.2.0 subread/1.6.2 bwa/0.7.15
ID_PATH=/home/projects/ku_00015/data/lotus/ids/ids.txt

# QC(1) - fastqc
echo "Step 1 - QC"
module load perl/5.30.2 jdk/14 fastqc/0.11.8
cd /home/projects/ku_00015/data/lotus 
fastqc read/*.fastq -o qc/

# QC(2) - MultiQC
module unload fastqc/0.11.8 perl/5.30.2
module load anaconda3/4.4.0
multiqc qc/*_fastqc.zip -o qc/

# Mapping - by BWA MEM
echo "Step 2 - Mapping with BWA MEM"
module load tools parallel/20200522 bwa/0.7.15
mkdir -p /home/projects/ku_00015/data/lotus/mapping/bwa
bwa index /home/projects/ku_00015/data/lotus/reference/chr3.fasta
cat ${ID_PATH} | parallel "bwa mem -t 2 /home/projects/ku_00015/data/lotus/reference/chr3.fasta /home/projects/ku_00015/data/lotus/read/{}_R1.fastq /home/projects/ku_00015/data/lotus/read/{}_R2.fastq > /home/projects/ku_00015/data/lotus/mapping/bwa/{}.sam"

# Sort and output BWA MEM mapped results as BAM by samtools
mkdir -p /home/projects/ku_00015/data/lotus/mapping/bwa-sort-bam
module load parallel/20200522 samtools/1.9
cat ${ID_PATH} | parallel "samtools sort -o /home/projects/ku_00015/data/lotus/mapping/bwa-sort-bam/{}.bam /home/projects/ku_00015/data/lotus/mapping/bwa/{}.sam"

# Mapping --------- Hisat2
echo "Step 2 - Mapping with Hisat2" 
module load tools parallel/20200522 anaconda2/4.4.0 hisat2/2.2.0
mkdir -p /home/projects/ku_00015/data/lotus/reference/hisat2
mkdir -p /home/projects/ku_00015/data/lotus/mapping/hisat2
hisat2-build /home/projects/ku_00015/data/lotus/reference/chr3.fasta /home/projects/ku_00015/data/lotus/reference/hisat2/chr3
cat ${ID_PATH} | parallel "hisat2 -x /home/projects/ku_00015/data/lotus/reference/hisat2/chr3 -1 /home/projects/ku_00015/data/lotus/read/{}_R1.fastq -2 /home/projects/ku_00015/data/lotus/read/{}_R2.fastq -S /home/projects/ku_00015/data/lotus/mapping/hisat2/{}.sam"

# Sort and output hisat2 mapped results as BAM by samtools
mkdir -p /home/projects/ku_00015/data/lotus/mapping/hisat2-sort-bam
module load parallel/20200522 samtools/1.9
cat ${ID_PATH} | parallel "samtools sort -o /home/projects/ku_00015/data/lotus/mapping/hisat2-sort-bam/{}.bam /home/projects/ku_00015/data/lotus/mapping/hisat2/{}.sam"

# Mapping statistics by samtools flagstat
module load parallel/20200522 samtools/1.9
mkdir -p /home/projects/ku_00015/data/lotus/mapping/bwa-mappingstat
cat ${ID_PATH} | parallel "samtools flagstat /home/projects/ku_00015/data/lotus/mapping/bwa-sort-bam/{}.bam > /home/projects/ku_00015/data/lotus/mapping/bwa-mappingstat/{}.flagstat"

# Mapping stat (CollectInsertSizeMetrics) by picard tools
module load java/1.8.0 picard-tools/2.20.2 parallel/20200522 gcc/9.3.0 intel/perflibs/2019_update5 R/3.6.1
mkdir -p /home/projects/ku_00015/data/lotus/mapping/bwa-picard-CollectInsertSizeMetrics
cat ${ID_PATH} | parallel "java -jar /services/tools/picard-tools/2.20.2/picard.jar CollectInsertSizeMetrics \I=/home/projects/ku_00015/data/lotus/mapping/bwa-sort-bam/{}.bam \O=/home/projects/ku_00015/data/lotus/mapping/bwa-picard-CollectInsertSizeMetrics/{}.txt \H=/home/projects/ku_00015/data/lotus/mapping/bwa-picard-CollectInsertSizeMetrics/{}.pdf"

# Mapping stat (CollectWgsMetrics) by picard tools
module load java/1.8.0 picard-tools/2.20.2 parallel/20200522 gcc/9.3.0 intel/perflibs/2019_update5 R/3.6.1
mkdir -p /home/projects/ku_00015/data/lotus/mapping/bwa-picard-CollectWgsMetrics
cat ${ID_PATH} | parallel "java -jar /services/tools/picard-tools/2.20.2/picard.jar CollectWgsMetrics \I=/home/projects/ku_00015/data/lotus/mapping/bwa-sort-bam/{}.bam \O=/home/projects/ku_00015/data/lotus/mapping/bwa-picard-CollectWgsMetrics/{}.txt \R=/home/projects/ku_00015/data/lotus/reference/chr3.fasta"


# test picard tools
#java -jar /services/tools/picard-tools/2.20.2/picard.jar CollectInsertSizeMetrics \I=/home/projects/ku_00015/data/rna-seq/mapping/hisat2/NG-23829_02_AD_12_BI_LS_01_E12_lib390812_6747_2.bam \O=output.txt \H=hist.pdf \AS=false
#java -jar /services/tools/picard-tools/2.20.2/picard.jar CollectInsertSizeMetrics \I=test.bam \O=output.txt \H=hist.pdf

# Mapping PacBio
module load minimap2/2.17r941
mkdir -p /home/projects/ku_00015/data/lotus/mapping/minimap2-pacbio
minimap2 -ax map-pb /home/projects/ku_00015/data/lotus/reference/chr3.fasta /home/projects/ku_00015/data/lotus/read/Gifu_PacBio.fastq > /home/projects/ku_00015/data/lotus/mapping/minimap2-pacbio/Gifu_PacBio.sam

# Mapping PacBio - sort, flagstat, insert size histogram
mkdir -p /home/projects/ku_00015/data/lotus/mapping/minimap2-pacbio-sort-bam
module load samtools/1.9
samtools sort -o /home/projects/ku_00015/data/lotus/mapping/minimap2-pacbio-sort-bam/Gifu_PacBio.bam /home/projects/ku_00015/data/lotus/mapping/minimap2-pacbio/Gifu_PacBio.sam

module load parallel/20200522 samtools/1.9
mkdir -p /home/projects/ku_00015/data/lotus/mapping/minimap2-mappingstat
samtools flagstat /home/projects/ku_00015/data/lotus/mapping/minimap2-pacbio-sort-bam/Gifu_PacBio.bam > /home/projects/ku_00015/data/lotus/mapping/minimap2-mappingstat/Gifu_PacBio.flagstat

module load java/1.8.0 picard-tools/2.20.2 parallel/20200522 gcc/9.3.0 intel/perflibs/2019_update5 R/3.6.1
mkdir -p /home/projects/ku_00015/data/lotus/mapping/minimap2-picard-CollectInsertSizeMetrics
java -jar /services/tools/picard-tools/2.20.2/picard.jar CollectInsertSizeMetrics \I=/home/projects/ku_00015/data/lotus/mapping/minimap2-pacbio-sort-bam/Gifu_PacBio.bam \O=/home/projects/ku_00015/data/lotus/mapping/minimap2-picard-CollectInsertSizeMetrics/Gifu_PacBio.txt \H=/home/projects/ku_00015/data/lotus/mapping/minimap2-picard-CollectInsertSizeMetrics/Gifu_PacBio.pdf \VALIDATION_STRINGENCY=SILENT \MINIMUM_PCT=0

# Mapping PacBio - samtools coverage
module load samtools/1.9
samtools coverage /home/projects/ku_00015/data/lotus/mapping/minimap2-pacbio-sort-bam/Gifu_PacBio.bam

# Merge multiple BAM files
module load samtools/1.9
samtools merge /home/projects/ku_00015/data/lotus/mapping/bwa-sort-bam/MG020_merged.bam /home/projects/ku_00015/data/lotus/mapping/bwa-sort-bam/MG020*.bam 

# Varient calling - bcf tools
module load bcftools/1.9
mkdir -p /home/projects/ku_00015/data/lotus/vcf
bcftools mpileup -Ov -f /home/projects/ku_00015/data/lotus/reference/chr3.fasta /home/projects/ku_00015/data/lotus/mapping/bwa-sort-bam/*[^bp].bam | \
bcftools call -Ov -mv > /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv.vcf

# Varient calling - splitting
egrep '^#|GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]0/0' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv.vcf > /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00.vcf

egrep '^#|GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]0/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv.vcf > /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-01.vcf

egrep '^#|GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]1/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv.vcf > /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-11.vcf

grep -v "#" /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00.vcf | wc -l
grep -v "#" /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-01.vcf | wc -l
grep -v "#" /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-11.vcf | wc -l

# SNPs Filtering - vcftools
module load perl/5.30.2 vcftools/0.1.16
cat /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00.vcf | vcf-annotate --filter d=10/Q=70 >> /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcftools.vcf
grep -v "#" /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter.vcf | wc -l

# SNPs Filtering - vcflib (vcffilter)
module load vcflib/1.0.0-rc2
rm /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf
vcffilter -f "DP > 50 & MQB > 0.8" /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00.vcf >> /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf
grep -v "#" /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l

# Burttii - polymorphism statistics
egrep 'GT:PL[[:blank:]]0/0' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_REF
egrep 'GT:PL[[:blank:]]1/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_ALT
egrep 'GT:PL[[:blank:]]0/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_ALT

# Gifu - polymorphism statistics
egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]0/0' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_REF

egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]1/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_ALT

egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]0/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_ALT

# MG012 - polymorphism statistics
egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]0/0' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_REF

egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]1/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_ALT

egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]0/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_ALT

# MG020 - polymorphism statistics
egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]0/0' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_REF

egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]1/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_ALT

egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]0/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_ALT

# MG042 - polymorphism statistics
egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]0/0' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_REF

egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]1/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_ALT

egrep 'GT:PL[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]][[:digit:]]/[[:digit:]]:[[:digit:]]{1,},[[:digit:]]{1,},[[:digit:]]{1,}[[:blank:]]0/1' /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf | wc -l # Calculate HOM_ALT

# RNA-seq build reference index
module load bowtie2/2.3.4.1
bowtie2-build /home/projects/ku_00015/data/lotus/reference/chr3.fasta /home/projects/ku_00015/data/lotus/reference/chr3

# RNA-seq mapping-tophat
module load parallel/20200522 tophat/2.1.1
mkdir /home/projects/ku_00015/data/lotus/mapping/tophat2
cat /home/projects/ku_00015/data/lotus/ids/ids_rnaseq | parallel "tophat -o /home/projects/ku_00015/data/lotus/mapping/tophat2/{} /home/projects/ku_00015/data/lotus/reference/chr3 /home/projects/ku_00015/data/lotus/read/{}_s1.fastq /home/projects/ku_00015/data/lotus/read/{}_s2.fastq"

# RNA-seq cufflinks
module load cufflinks/2.2.1
mkdir -p /home/projects/ku_00015/data/lotus/gene_model

# Varient calling from mRNA
module load bcftools/1.9
mkdir -p /home/projects/ku_00015/data/lotus/vcf
bcftools mpileup -Ov -f /home/projects/ku_00015/data/lotus/reference/chr3.fasta /home/projects/ku_00015/data/lotus/mapping/tophat2/Gifu_rnaSeq/accepted_hits.bam /home/projects/ku_00015/data/lotus/mapping/tophat2/mg20_rnaSeq/accepted_hits.bam | \
bcftools call -Ov -mv > /home/projects/ku_00015/data/lotus/vcf/vcf-f-mRNA-callmv.vcf

grep -v "#" /home/projects/ku_00015/data/lotus/vcf/vcf-f-mRNA-callmv.vcf | wc -l

# Common SNPs gDNA & mRNA
module load bedtools/2.28.0
bedtools intersect -a /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf -b /home/projects/ku_00015/data/lotus/vcf/vcf-f-mRNA-callmv.vcf| wc -l # common SNPs
bedtools intersect -v -a /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf -b /home/projects/ku_00015/data/lotus/vcf/vcf-f-mRNA-callmv.vcf | wc -l # in gDNA while not in mRNA
bedtools intersect -v -b /home/projects/ku_00015/data/lotus/vcf/vcf-f-bwa-callmv-MG20-00-filter-vcffilter.vcf -a /home/projects/ku_00015/data/lotus/vcf/vcf-f-mRNA-callmv.vcf | wc -l # in mRNA while not in gDNA

# Visualization - igv
module load jdk/14 igv/2.8.0
java -Xmx8000m  -jar igvtools.jar -g
