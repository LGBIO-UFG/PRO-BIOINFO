#!/bin/bash
# leonardo Corvalan
# lcjcorvalan@gmail.com
# 10/09/2025
# 15/09/2025 - last update

# this is script is just adptation of arima mapping pipeline to lg.bio.br server
# https://github.com/ArimaGenomics/mapping_pipeline/tree/master

# execute
# bash arima_map.sh <ref.fa> <PE_trimmed_1.fa> <PE_trimmed_2.fa> <prefix>


reference=$1
PE_1=$2 # reads trimmed
PE_2=$3
PREFIX=$4 # a name for your files
PE_alg_1=$PREFIX.1.trimmed.bam
PE_alg_2=$PREFIX.2.trimmed.bam
PE_filt_1=$PREFIX.1.filt.bam
PE_filt_2=$PREFIX.2.filt.bam
FAI=$reference.fai
PAIR=$PREFIX.paired.bam
Hicbam=$PREFIX.Hic.bam

# in put your trimmed reds
# the number is seted for 16

export PATH=$PATH:/media/lgbio-nas1/lcorvalan/Hic-data/Mfa

echo "### Step 0: Index reference"

bwa index -a bwtsw -p $PREFIX $reference

# alinhando
echo "### Step 1.A: FASTQ to BAM (1st)"

bwa mem -t 16 $PREFIX $PE_1| samtools view -@ 16 -Sb - > $PE_alg_1
echo "### Step 1.B: FASTQ to BAM (2nd)"

echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t 16 $PREFIX $PE_2| samtools view -@ 16 -Sb - > $PE_alg_2

echo "### Step 2.A: Filter 5' end (1st)"
samtools view -h $PE_alg_1 | perl /media/lgbio-nas1/lcorvalan/Hic-data/Mfa/filter_five_end.pl| samtools view -Sb - > $PE_filt_1

echo "### Step 2.B: Filter 5' end (2nd)"
samtools view -h $PE_alg_2 | perl /media/lgbio-nas1/lcorvalan/Hic-data/Mfa/filter_five_end.pl| samtools view -Sb - > $PE_filt_2

echo "### Step 3A: Pair reads & mapping quality filter"
perl /media/lgbio-nas1/lcorvalan/Hic-data/Mfa/two_read_bam_combiner.pl $PE_filt_1 $PE_filt_2 samtools 10 | samtools view -bS -t $FAI - | samtools sort -@ 18 -o tmp.bam -

PicardCommandLine AddOrReplaceReadGroups INPUT=tmp.bam OUTPUT=$PAIR ID=$PREFIX LB=$PREFIX SM=$PREFIX.hic PL=ILLUMINA PU=none

echo "### Step 4: Mark duplicates"

PicardCommandLine MarkDuplicates INPUT=$PAIR  OUTPUT=$Hicbam METRICS_FILE=$PREFIX.hic_METRICS_FILE.txt TMP_DIR=temp ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index $Hicbam

echo "### Step 5: Hic reads distance stats"

perl /media/lgbio-nas1/lcorvalan/Hic-data/Mfa/get_stats.pl $Hicbam > $Hicbam.stats
