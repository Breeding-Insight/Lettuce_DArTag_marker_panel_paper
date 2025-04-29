#!/bin/bash
java -Xmx100g  -jar /workdir/ms3743/lettuce/realignment/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar CombineGVCFs -R /workdir/ms3743/lettuce/realignment/Lsativa_467_v8.fa --variant gvcfsnew.list -O /workdir/ms3743/lettuce/realignment/varcall/combined_new_gvcf.vcf.gz
java -Xmx100g  -jar /workdir/ms3743/lettuce/realignment/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar GenotypeGVCFs -R /workdir/ms3743/lettuce/realignment/Lsativa_467_v8.fa --variant combined_gvcf.vcf.gz -O combined_genotype_gvcf.vcf.gz

#Extracting SNPs
java -Xmx100g  -jar /workdir/ms3743/lettuce/realignment/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar SelectVariants -R /workdir/ms3743/lettuce/realignment/Lsativa_467_v8.fa -V combined_genotype_gvcf.vcf.gz --select-type-to-include SNP -O combined_snp_raw.vcf.gz
java -Xmx100g  -jar /workdir/ms3743/lettuce/realignment/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar VariantFiltration -R /workdir/ms3743/lettuce/realignment/Lsativa_467_v8.fa -V combined_snp_raw.vcf.gz -O combined_snp_raw_tagged.vcf.gz --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum  < -12.5 || ReadPosRankSum < -8.0" --filter-name "Default_recommended"
vcftools --gzvcf combined_snp_raw_tagged.vcf.gz --remove-filtered-all --recode -c > final_combined_snps.vcf


#Extracting INDELs
java -Xmx100g  -jar /workdir/ms3743/lettuce/realignment/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar SelectVariants -R /workdir/ms3743/lettuce/realignment/Lsativa_467_v8.fa -V combined_genotype_gvcf.vcf.gz --select-type-to-include INDEL -O combined_indel_raw.vcf.gz
java -Xmx100g  -jar /workdir/ms3743/lettuce/realignment/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar VariantFiltration -R /workdir/ms3743/lettuce/realignment/Lsativa_467_v8.fa -V combined_indel_raw.vcf.gz -O combined_indel_raw_tagged.vcf.gz --filter-expression "QD < 2.0 || FS > 200.0 || MQ < 40.0 || InbreedingCoeff < -0.8 || ReadPosRankSum < -20.0" --filter-name "Default_recommended"
vcftools --gzvcf combined_indel_raw_tagged.vcf.gz --remove-filtered-all --recode -c > final_combined_indel.vcf
