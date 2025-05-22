cd /workdir/ml2498/Lettuce/Genotyping_2022/Realign_validation2/
snps=final_combined_new_snps
  
grep "^[^#]" ${snps}_RIL.recode.vcf | wc -l # 125,852

vcftools --vcf ${snps}_RIL.recode.vcf --max-alleles 2 --min-alleles 2 --max-missing 0.5 --out ${snps}_RIL_qual_AD_BIallele_MISS --recode --recode-INFO-all
# 27,182

vcftools --vcf ${snps}_RIL_qual_AD_BIallele_MISS.recode.vcf --out ${snps}_RIL_qual_AD_BIallele_MISS --012

#########################

library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")

RIL<-fread("final_combined_new_snps_RIL_qual_AD_BIallele_MISS.012",data.table=F) ## dosage of alt alleles
RIL_ind<-read.table("final_combined_new_snps_RIL_qual_AD_BIallele_MISS.012.indv",header=F)
RIL_pos<-fread("final_combined_new_snps_RIL_qual_AD_BIallele_MISS.012.pos",data.table=F)
colnames(RIL_pos)<-c("Chr","Pos")
RIL_pos$SNP<-paste(RIL_pos$Chr,RIL_pos$Pos,sep="_")
RIL<-RIL[,-1]
rownames(RIL)<-RIL_ind[,1]
colnames(RIL)<-RIL_pos$SNP

## missing data -1 to NA
for(i in 1:ncol(RIL)){
  RIL[which(RIL[,i]==-1),i]<-NA
}

#### select parents and code progenies based on parents (A B H)
## Lac_0038 is Eruption (P1); Lac_0044 is Reine des Glaces (P2)
RIL.1<-data.frame(t(RIL))


Parents<-RIL.1[,which(colnames(RIL.1) %in% c("Lac_0038","Lac_0044"))]
rils<-RIL.1[,-which(colnames(RIL.1) %in% c("Lac_0038","Lac_0044"))]

temp<-strsplit(rownames(RIL.1),"_")
temp<-data.frame(matrix(unlist(temp),ncol=2,byrow=T))
colnames(temp)<-c("sequence","sequence_pos")
temp$sequence<-gsub("Chr","",temp$sequence)
RIL.2<-cbind.data.frame(rownames(RIL.1),Parents,temp,rils)
colnames(RIL.2)[1:3]<-c("snp_name","P1","P2")
write.table(RIL.2,"Mappoly_formatted_BxE_RILs_012.csv",row.names=F,sep=",",quote=F)
length(RIL.2$P1[which(RIL.2$P1==1)])/nrow(RIL.2) # 0.03642116
length(RIL.2$P2[which(RIL.2$P2==1)])/nrow(RIL.2) # 0.03785593

######### convert to A,B,H format
## remove P1=P2=0 or 2
## remove P1=NA and /or P2=NA

rm_index<-c()
for(i in 1:nrow(Parents)){
  #if((Parents$Lac_0038[i]==0 & Parents$Lac_0044[i]==0)|(Parents$Lac_0038[i]==2 & Parents$Lac_0044[i]==2)|is.na(Parents$Lac_0038[i])|is.na(Parents$Lac_0044[i])){
  if((Parents$Lac_0038[i]==Parents$Lac_0044[i])|is.na(Parents$Lac_0038[i])|is.na(Parents$Lac_0044[i])){
    rm_index<-c(rm_index,i)
  }
}

Parents.f<-Parents[-rm_index,]  
rils.f<-rils[-rm_index,] 
temp.f<-temp[-rm_index,]
## 2,446 SNPs retained
length(Parents.f$Lac_0038[which(Parents.f$Lac_0038==1)])/nrow(Parents.f) # heterozygosity 0.1548767
length(Parents.f$Lac_0044[which(Parents.f$Lac_0044==1)])/nrow(Parents.f) # heterozygosity 0.1618435


Parents.abh<-Parents.f
rils.abh<-rils.f
for(j in 1:nrow(Parents.f)){
  if((Parents.f[j,1]==0 & Parents.f[j,2]==2) | (Parents.f[j,1]==2 & Parents.f[j,2]==0)){
    Parents.abh[j,1]<-"A"
    Parents.abh[j,2]<-"B"
    rils.abh[j,which(rils.abh[j,]==Parents.f[j,1])]<-"A"
    rils.abh[j,which(rils.abh[j,]==Parents.f[j,2])]<-"B"
    rils.abh[j,which(rils.abh[j,]==1)]<-"H"
  }else if(Parents.f[j,1]==1 & Parents.f[j,2]==2){
    Parents.abh[j,1]<-"H"
    Parents.abh[j,2]<-"B"
    rils.abh[j,which(rils.abh[j,]==Parents.f[j,1])]<-"H"
    rils.abh[j,which(rils.abh[j,]==Parents.f[j,2])]<-"B"
    #rils.abh[j,which(rils.abh[j,] %in% c(0,2))]<-"A"
    rils.abh[j,which(rils.abh[j,] %in% c(0))]<-"A"
  }else if(Parents.f[j,1]==1 & Parents.f[j,2]==0){
    Parents.abh[j,1]<-"H"
    Parents.abh[j,2]<-"B"
    rils.abh[j,which(rils.abh[j,]==Parents.f[j,1])]<-"H"
    rils.abh[j,which(rils.abh[j,]==Parents.f[j,2])]<-"B"
    #rils.abh[j,which(rils.abh[j,] %in% c(0,2))]<-"A"
    rils.abh[j,which(rils.abh[j,] %in% c(2))]<-"A"
  }else if(Parents.f[j,1]==0 & Parents.f[j,2]==1){
    Parents.abh[j,1]<-"A"
    Parents.abh[j,2]<-"H"
    rils.abh[j,which(rils.abh[j,]==Parents.f[j,1])]<-"A"
    rils.abh[j,which(rils.abh[j,]==Parents.f[j,2])]<-"H"
    rils.abh[j,which(rils.abh[j,] %in% c(2))]<-"B"
  }else if(Parents.f[j,1]==2 & Parents.f[j,2]==1){
    Parents.abh[j,1]<-"A"
    Parents.abh[j,2]<-"H"
    rils.abh[j,which(rils.abh[j,]==Parents.f[j,1])]<-"A"
    rils.abh[j,which(rils.abh[j,]==Parents.f[j,2])]<-"H"
    rils.abh[j,which(rils.abh[j,] %in% c(0))]<-"B"
  }
}
#rils.f[1,]
#rils.f[7,] 
#rils.abh[7,]
RIL.abh<-cbind.data.frame(rownames(Parents.abh),Parents.abh,temp.f,rils.abh)
colnames(RIL.abh)[1:3]<-c("snp_name","P1","P2")
write.table(RIL.abh,"Mappoly_formatted_BxE_RILs_ABH.csv",row.names=F,sep=",",quote=F)



