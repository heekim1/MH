#!/bin/sh
reference=$1
fastq=$2
jobname=$3
scriptdir=$4
#tmap index -v -f $reference
tmap map2 -v -f $reference -r $fastq -o 2 -s ${jobname}_intermediate.bam
samtools sort ${jobname}_intermediate.bam -o ${jobname}.bam
samtools index ${jobname}.bam
rm ${jobname}_intermediate.bam
python $scriptdir/adjustMQ.py --bam ${jobname}.bam --min 60 --max 60 --seed 1000
samtools index ${jobname}_adjusted_mapq.bam