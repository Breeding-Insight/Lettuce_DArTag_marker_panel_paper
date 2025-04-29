#!/bin/bash
while read acc
do
./gatk-4.3.0.0/gatk --java-options "-Xmx20g -XX:ParallelGCThreads=2" HaplotypeCaller --native-pair-hmm-threads 2 -R Lsativa_467_v8.fa -I ./filteredbam/${acc}.filtered.RG.bam -ERC GVCF -O ./gatk_gvcf/${acc}_gvcf.vcf.gz
done <accnew.list
