#!/bin/bash
while read acc
do
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 10 -phred33 /workdir/data/lettuce/DArT_raw_seq/validation_plates/${acc}.FASTQ.gz ./trimmed/${acc}.trim.FASTQ.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
bwa mem -t 20 Lsativa_467_v8.fa ./trimmed/${acc}.trim.FASTQ.gz > ./samfile/${acc}.sam
samtools view -@ 20 -bS -F 4 ./samfile/${acc}.sam > ./sam2bam/${acc}.bam
samtools sort ./sam2bam/${acc}.bam -o ./sam2bam/${acc}.sorted.bam
samtools view -b -q 50 ./sam2bam/${acc}.sorted.bam > ./filteredbam/${acc}.filtered.bam
samtools index ${acc}.filtered.bam 
done <accnew.list
