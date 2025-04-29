#!/bin/bash
java -Xmx100g  -jar /workdir/ms3743/lettuce/realignment/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar CombineGVCFs -R /workdir/ms3743/lettuce/realignment/Lsativa_467_v8.fa --variant gvcfsnew.list -O /workdir/ms3743/lettuce/realignment/varcall/combined_new_gvcf.vcf.gz
