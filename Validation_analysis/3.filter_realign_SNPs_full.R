## 1. change fastq ID to sample ID
## 2. filter for SNP quality & read depth
## 3. filter for MAF (overall PCA)
## 4. from 2, seperate by population, filter for 
## (1) line-based missing rate & heterozygosity (2) SNP-based missing rate & heterozygosity (3) MAF
## 5. Check how many targeted loci got captured by realignment

## 1. change fastq ID to sample ID
cd /workdir/ml2498/Lettuce/Genotyping_2022/Realign_validation2/
snps=final_combined_new_snps
indel=final_combined_new_indel

chgtaxa=/workdir/ml2498/Lettuce/script/validation_plates/chgTaxa_vcf.py
python3 ${chgtaxa} ${snps}.vcf ${snps}_sampleID.vcf ../Realign_validation/Fastq_sampleID_map.txt
python3 ${chgtaxa} ${indel}.vcf ${indel}_sampleID.vcf ../Realign_validation/Fastq_sampleID_map.txt

sed -i '/^#/! s/Lsat_1_v8_lg_/Chr/g' ${snps}_sampleID.vcf
sed -i '/^#/! s/Lsat_1_v8_lg_/Chr/g' ${indel}_sampleID.vcf

grep "^[^#]" ${snps}_sampleID.vcf | wc -l # 147,043
grep "^[^#]" ${indel}_sampleID.vcf | wc -l # 3,893


## nano ${snps}_sampleID.vcf ## manually delete the last tab in the header line
bcftools annotate --set-id +'%CHROM\_%POS' ${snps}_sampleID.vcf > ${snps}_sampleID_snpID.vcf
bcftools annotate --set-id +'%CHROM\_%POS' ${indel}_sampleID.vcf > ${indel}_sampleID_snpID.vcf

## 2. filter for SNP quality & read depth
vcftools --vcf ${snps}_sampleID_snpID.vcf  --site-depth ## out.ldepth

## investigate cutoff for read depth (DP)
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")
site_depth<-fread("out.ldepth",data.table=F)
hist(site_depth$SUM_DEPTH,breaks=50)
quantile(site_depth$SUM_DEPTH,0.90) #5042; 5042/350=14.4
quantile(site_depth$SUM_DEPTH,0.10) #6; 10/350=0.029
length(site_depth$SUM_DEPTH[which(site_depth$SUM_DEPTH<500 & site_depth$SUM_DEPTH>10)]) # 72,831 SNPs
##

bcftools view -e '%QUAL<30' ${snps}_sampleID_snpID.vcf | bcftools view -e '%MAX(FORMAT/AD[:0])<1 || %MAX(FORMAT/AD[:1])<2' > ${snps}_qual_AD.vcf
grep "^[^#]" ${snps}_qual_AD.vcf | wc -l #130,279

#bcftools view -e '%QUAL<30' ${snps}_sampleID_snpID.vcf > ${snps}_qual.vcf
#grep "^[^#]" ${snps}_qual.vcf | wc -l # 143,299

#bcftools view -e '%MAX(FORMAT/AD[:0])<1 || %MAX(FORMAT/AD[:1])<2' ${snps}_sampleID_snpID.vcf > ${snps}_AD.vcf
#grep "^[^#]" ${snps}_AD.vcf | wc -l # 126,706

#bcftools view -e '%MIN(FORMAT/DP)<5 || %MAX(FORMAT/DP)>1000' ${snps}_sampleID_snpID.vcf > ${snps}_DP.vcf
#grep "^[^#]" ${snps}_DP.vcf | wc -l # 87,853

## filter for mean DP; ~350 samples; bi-allelic
vcftools --vcf ${snps}_qual_AD.vcf --max-alleles 2 --min-alleles 2 --out ${snps}_qual_AD_BIallele --recode --recode-INFO-all
## retained 128,572

## *** potentially add this
vcftools --vcf ${snps}_qual_AD_BIallele.recode.vcf --minDP 5 --not-chr Chr0 --out ${snps}_qual_AD_DP_BIallele_chr0 --recode --recode-INFO-all ## cutoff>=5 was based on mimic F1 error vs good SNP read depth distribution
#125,852

## 3. filter for MAF (overall PCA)****

vcftools --vcf ${snps}_qual_AD_DP_BIallele_chr0.recode.vcf --maf 0.05 --out ${snps}_qual_AD_DP_BIallele_chr0_MAF --recode --recode-INFO-all
## retained 72,428

vcftools --vcf ${snps}_qual_AD_DP_BIallele_chr0.recode.vcf --maf 0.05 --max-missing 0.5 --out ${snps}_qual_AD_DP_BIallele_chr0_MAF_MISS --recode --recode-INFO-all
## retained 5,852

## 4. from 2, seperate by population, filter for 
## (1) line-based missing rate & heterozygosity (2) SNP-based missing rate & heterozygosity (3) MAF
# vcftools --vcf final_combined_snps_qual_AD_DP_BIallele.recode.vcf --not-chr Chr0 --out final_combined_snps_qual_AD_DP_BIallele_chr0 --recode --recode-INFO-all
bcftools sort -o ${snps}_qual_AD_DP_BIallele_sort.recode.vcf ${snps}_qual_AD_DP_BIallele_chr0.recode.vcf
vcftools --vcf ${snps}_qual_AD_DP_BIallele_sort.recode.vcf --keep ../Realign_validation/div_lines.txt --out ${snps}_div --recode --recode-INFO-all
vcftools --vcf ${snps}_qual_AD_DP_BIallele_sort.recode.vcf --keep ../Realign_validation/RIL_lines.txt --out ${snps}_RIL --recode --recode-INFO-all
vcftools --vcf ${snps}_qual_AD_DP_BIallele_sort.recode.vcf --keep ../Realign_validation/OtherLines.txt --out ${snps}_other --recode --recode-INFO-all

/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -vcf  ${snps}_div.recode.vcf -genotypeSummary taxa -export ${snps}_div_TaxaSum

##
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")
TaxaSum<-fread("final_combined_new_snps_div_TaxaSum.txt",data.table=F)

pdf("./plot/TaxaSum_prefilter_div.pdf")
hist(TaxaSum$`Proportion Heterozygous`,breaks=50)
hist(TaxaSum$`Proportion Missing`,breaks=50)
dev.off()
rm_div_line<-TaxaSum$`Taxa Name`[which(TaxaSum$`Proportion Missing`>0.9)]
write.table(rm_div_line,"remove_div_lines_highMissing.txt",row.names=F,col.names=F,sep="\t",quote=F)
###
vcftools --vcf ${snps}_div.recode.vcf --remove remove_div_lines_highMissing.txt --out ${snps}_div_IndvMiss --recode --recode-INFO-all ## 12 lines removed, 132 retained

/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -vcf  ${snps}_div_IndvMiss.recode.vcf -genotypeSummary site -export ${snps}_div_SiteSum

##
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")
SiteSum<-fread("final_combined_new_snps_div_SiteSum.txt",data.table=F)

pdf("./plot/SiteSum_prefilter_div.pdf")
hist(SiteSum$`Proportion Heterozygous`,breaks=30)
hist(SiteSum$`Proportion Missing`,breaks=20)
hist(SiteSum$`Minor Allele Frequency`,breaks=20)
dev.off()

## check statistics for 3K target sites
pos_3K<-fread("3K_targetLoci_pos_info.txt",data.table=F)
pos_3K$SNP<-paste(pos_3K$Chrom,pos_3K$ChromPos,sep="_")
cm_snp<-intersect(pos_3K$SNP,SiteSum$`Site Name`) # 2,554 target sites captured
SiteSum_3Kcm<-SiteSum[which(SiteSum$`Site Name` %in% cm_snp),]
pdf("./plot/SiteSum_prefilter_div_3Ktarget.pdf")
hist(SiteSum_3Kcm$`Proportion Heterozygous`,breaks=30)
hist(SiteSum_3Kcm$`Proportion Missing`,breaks=20)
hist(SiteSum_3Kcm$`Minor Allele Frequency`,breaks=20)
dev.off()


MAF<-0.05
missing<-0.5
het<-0.15
keep_site<-SiteSum$`Site Name`[which(SiteSum$`Proportion Missing`<=missing & SiteSum$`Proportion Heterozygous`<=het &SiteSum$`Minor Allele Frequency`>=MAF)]
## 5346 SNPs retained
keep_site_3K<-keep_site[which(keep_site %in% pos_3K$SNP)]
## 2046 SNPs in 3K target retained
write.table(keep_site,"keep_SNPs_div.txt",row.names=F,col.names=F,sep="\t",quote=F)

##
vcftools --vcf ${snps}_div_IndvMiss.recode.vcf --snps keep_SNPs_div.txt --out ${snps}_div_final --recode --recode-INFO-all

/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -vcf  ${snps}_div_final.recode.vcf -genotypeSummary taxa -export ${snps}_div_TaxaSum2
/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -vcf  ${snps}_div_final.recode.vcf -genotypeSummary site -export ${snps}_div_SiteSum2
##
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")
TaxaSum<-fread("final_combined_new_snps_div_TaxaSum2.txt",data.table=F)

pdf("./plot/TaxaSum_filtered_div.pdf")
hist(TaxaSum$`Proportion Heterozygous`,breaks=50)
hist(TaxaSum$`Proportion Missing`,breaks=50)
dev.off()

SiteSum<-fread("final_combined_new_snps_div_SiteSum2.txt",data.table=F)
pdf("./plot/SiteSum_filter_div.pdf")
hist(SiteSum$`Proportion Heterozygous`,breaks=30)
hist(SiteSum$`Proportion Missing`,breaks=20)
hist(SiteSum$`Minor Allele Frequency`,breaks=20)
dev.off()

###
## get the pos and allele information of filtered SNPs
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${snps}_div_final.recode.vcf > ${snps}_div_final.info



#########################
#### 1M window plot ####
#########################
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")
info<-fread("final_combined_new_snps_div_final.info",data.table=F)
colnames(info)<-c("Chrom","ChromPos","MarkerID","RefAllele","AltAllele")

setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/from115K")
REGION<-read.table("SNP_counts_115K_1MWin_INT30_MAF_missing.txt",header=T)
REGION$Chr<-gsub("LG_","Chr",REGION$Chr,fixed=T)

REGION$counts_realign<-0
i=1
for (i in 1:nrow(REGION)){
  chr<-REGION$Chr[i]
  ST<-REGION$ST[i]
  END<-REGION$END[i]
  
  temp<-info[which(info$Chrom==chr & info$ChromPos>=ST & info$ChromPos<END),]
  if(nrow(temp)>0){
    ct<-nrow(temp)
    REGION$counts_realign[i]<-ct
  }
}

library(ggplot2)
p1<-ggplot(REGION,aes(END,counts_realign))+
  geom_bar(stat="identity")+
  facet_grid(rows = vars(Chr))+
  labs(x="Physical position",y="SNP counts (1Mb window; re-alignment)")+
  theme_light()+
  theme(strip.background = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        strip.text.y = element_text(size=10,colour = 'black'),
        strip.text.x = element_text(size=5,colour = 'black'),
        strip.placement = "outside",
        #axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=7),
        text = element_text(size=10),
        #legend.title=element_blank(),
        legend.text.align = 0,
        legend.text=element_text(size=rel(1)))

setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/plot")
pdf("SNP_counts_realign_1M_window.pdf",height=7,width=7)
print(p1)
dev.off()
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")
write.table(REGION,"SNP_counts_5.3K_realign.txt",quote=F,row.names=F,sep="\t")

### subset the same SNPs for OtherLines and check quality

vcftools --vcf ${snps}_other.recode.vcf --snps keep_SNPs_div.txt --out ${snps}_other_final --recode --recode-INFO-all
vcftools --vcf ${snps}_RIL.recode.vcf --snps keep_SNPs_div.txt --out ${snps}_RIL_final --recode --recode-INFO-all

/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -vcf  ${snps}_other_final.recode.vcf -genotypeSummary taxa -export ${snps}_other_TaxaSum
/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -vcf  ${snps}_other_final.recode.vcf -genotypeSummary site -export ${snps}_other_SiteSum


library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")
TaxaSum<-fread("final_combined_new_snps_other_TaxaSum.txt",data.table=F)

pdf("./plot/TaxaSum_filtered_other.pdf")
hist(TaxaSum$`Proportion Heterozygous`,breaks=10) ## mimic Fn's among these lines, expect to see higher heterozygosity
hist(TaxaSum$`Proportion Missing`,breaks=10)
dev.off()


SiteSum<-fread("final_combined_new_snps_other_SiteSum.txt",data.table=F)

pdf("./plot/SiteSum_filtered_other.pdf")
hist(SiteSum$`Proportion Heterozygous`,breaks=30)
hist(SiteSum$`Proportion Missing`,breaks=20)
hist(SiteSum$`Minor Allele Frequency`,breaks=20)
dev.off()


