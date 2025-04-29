#!/bin/bash
while read acc
do
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 10 -phred33 /workdir/data/lettuce/DArT_raw_seq/validation_plates/batch2_29/${acc}.FASTQ.gz ./trimmed/${acc}.trim.FASTQ.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
done <accnew.list
