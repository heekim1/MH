#!/bin/sh
reference=$1
fastq=$2
jobname=$3
scriptdir=$4
#bwa index $reference
bwa aln $reference $fastq > ${jobname}_intermediate.sai
bwa samse $reference ${jobname}_intermediate.sai $fastq > ${jobname}_intermediate.sam
samtools view -b ${jobname}_intermediate.sam -o ${jobname}_intermediate.bam
samtools sort ${jobname}_intermediate.bam -o ${jobname}.bam
samtools index ${jobname}.bam
rm ${jobname}_intermediate.sai ${jobname}_intermediate.sam ${jobname}_intermediate.bam
python $scriptdir/adjustMQ.py --bam ${jobname}.bam --min 65 --max 75 --seed 1000
samtools index ${jobname}_adjusted_mapq.bam
python $scriptdir/editAlignments.py --bam ${jobname}_adjusted_mapq.bam --pos chr5:2448145,chr13:54060892
samtools index ${jobname}_adjusted_mapq_adjusted_align.bam