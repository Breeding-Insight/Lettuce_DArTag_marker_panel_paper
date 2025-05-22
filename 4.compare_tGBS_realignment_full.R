### Overall PCA
cd /workdir/ml2498/Lettuce/Genotyping_2022/Realign_validation2/
snps=final_combined_new_snps
#/programs/tassel-5-standalone/run_pipeline.pl -importGuess ${snps}_qual_AD_DP_BIallele_chr0_MAF.recode.vcf -KinshipPlugin -method Centered_IBS -endPlugin -export Lettuce_realign_kinship_validation_full.txt -exportType SqrMatrix
/programs/tassel-5-standalone/run_pipeline.pl -importGuess ${snps}_qual_AD_DP_BIallele_chr0_MAF_MISS.recode.vcf -KinshipPlugin -method Centered_IBS -endPlugin -export Lettuce_realign_kinship_validation_full_2025.txt -exportType SqrMatrix


### PCA for 14K SNPs
library(data.table)
library("RColorBrewer")
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")
#setwd("/workdir/ml2498/Lettuce/Genotyping_2022/Realign_validation2/")
#kin<-fread("Lettuce_realign_kinship_validation_full.txt",data.table=F)
kin<-fread("Lettuce_realign_kinship_validation_full_2025.txt",data.table=F)
info<-fread("../Validation/passport/passport_PCA.txt",data.table=F)

#length(intersect(kin$`72`,info$AMID))
#commonLine<-intersect(kin$`72`,info$AMID)
#info[-which(info$AMID %in% commonLine),]
# GBS-329 Lactuca indica (sexually uncompatible)

rownames(kin)<-kin[,1]
kin<-kin[,-1]
colnames(kin)<-rownames(kin)
#kin<-kin[which(rownames(kin) %in% commonLine),which(colnames(kin) %in% commonLine)]
which(rownames(kin)!=colnames(kin))

eigen_res<-eigen(kin)
eigen_factor<-eigen_res$vectors[,1:10]
eigen_factor<-cbind.data.frame(rownames(kin),eigen_factor)
colnames(eigen_factor)<-c("Taxa",paste("PC",1:10,sep=""))
eigen_factor<-merge(info,eigen_factor,by.x="Sample_ID",by.y="Taxa",all.y=T)
#eigen_factor$PC2.new<-eigen_factor$PC2*(-1)

eigen_val<-eigen_res$values
var_explained<-round(eigen_val*100/sum(eigen_val),2)
col_platte<-brewer.pal(n = 10, name = "Paired") # maximum 12 colors
#"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A"
col_platte<-col_platte[c(2,4,6,8,10)]
#"#1F78B4" "#33A02C" "#E31A1C" "#FF7F00" "#6A3D9A"
dup.sub<-eigen_factor[which(eigen_factor$dup_sample_pair %in% c(1:18)),]

eigen_factor$population_MS<-factor(eigen_factor$population_MS,levels=c("Diverse panel","BxE RILs","Artificially created F1",
                                                                       "Breeding lines","Mutation" ))
write.table(eigen_factor,"PCA_passport_5852SNP_376indv.txt",row.names=F,quote=F,sep="\t")

library(ggplot2)
library(ggrepel)
p<-ggplot(eigen_factor,aes(PC1,PC2,color=population_MS,label=label))+
  geom_point(aes(shape=If_duplicated))+ 
  #geom_line(data=dup.sub,mapping=aes(group=dup_sample_pair))+
  scale_shape_manual(values=c(1,2))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#33A02C",  "#FF7F00", "#6A3D9A"))+
  geom_text_repel(color="black")+
  theme_light()+
  labs(x=paste("PC1 (",var_explained[1],"%)",sep=""),y=paste("PC2 (",var_explained[2],"%)",sep=""),
       color="Population",shape="Duplicated sample")
p

#pdf("./plot/PCA_376x72K_realign.pdf",width=8,height=5)
pdf("./plot/PCA_376x72K_realign_2025.pdf",width=12,height=7,family="sans")
print(p)
dev.off()

p2<-ggplot(eigen_factor,aes(PC1,PC2,color=population_MS,label=label))+
  geom_point(aes(shape=If_duplicated))+ 
  geom_line(data=dup.sub,mapping=aes(group=dup_sample_pair))+
  scale_shape_manual(values=c(1,2))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#33A02C",  "#FF7F00", "#6A3D9A"))+
  geom_text_repel(color="black")+
  theme_light()+
  labs(x=paste("PC1 (",var_explained[1],"%)",sep=""),y=paste("PC2 (",var_explained[2],"%)",sep=""),
       color="Population",shape="Duplicated sample")
p2
pdf("./plot/PCA_376x72K_realign_2025_wLines.pdf",width=12,height=7,family="sans")
print(p2)
dev.off()


### convert vcf to 0, 1, 2 format
cd /workdir/ml2498/Lettuce/Genotyping_2022/Realign_validation2/
vcftools --vcf ${snps}_div_final.recode.vcf --out final_combined_snps_div_final --012
vcftools --vcf ${snps}_other_final.recode.vcf --out final_combined_snps_other_final --012

####
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")
div<-fread("final_combined_snps_div_final.012",data.table=F) ## dosage of alt alleles
div_ind<-read.table("final_combined_snps_div_final.012.indv",header=F)
div_pos<-fread("final_combined_snps_div_final.012.pos",data.table=F)
colnames(div_pos)<-c("Chr","Pos")
div_pos$SNP<-paste(div_pos$Chr,div_pos$Pos,sep="_")
div<-div[,-1]
rownames(div)<-div_ind[,1]
colnames(div)<-div_pos$SNP

## missing data -1 to NA
for(i in 1:ncol(div)){
  div[which(div[,i]==-1),i]<-NA
}
## change from alt dosage to ref dosage
div.adj<-2-div


#### format tGBS data
reseq_snp012<-fread("/Users/ml2498/Desktop/Labserver/Lettuce/tGBS_SNPs/raw_genotype/Lettuce.all.snps.reanalysis.genotype.txt",data.table=F)
reseq_snp012$Chr<-gsub("Lsat_1_v8_lg_","Chr",reseq_snp012$Chr)
reseq_snp012$snpID<-paste(reseq_snp012$Chr,reseq_snp012$Pos,sep="_")

reseq_snp012.adj<-as.matrix(reseq_snp012[,5:518])
reseq_snp012.adj[which(reseq_snp012.adj==0)]<-NA
reseq_snp012.adj[which(reseq_snp012.adj==2)]<-0
reseq_snp012.adj[which(reseq_snp012.adj==1)]<-2
reseq_snp012.adj[which(reseq_snp012.adj==3)]<-1
row.names(reseq_snp012.adj)<-reseq_snp012[,519]
reseq_snp012.adj<-t(reseq_snp012.adj)
rownames(reseq_snp012.adj)<-gsub("A","",rownames(reseq_snp012.adj))
rownames(reseq_snp012.adj)<-gsub("MIX","",rownames(reseq_snp012.adj))
## make accession name consistant
info<-fread("../Validation/passport/passport_PCA.txt",data.table=F)
info$GBS<-gsub("-","_",info$GBS,fixed=T)
info.sub<-info[which(info$Population=="tGBS.1"),]
for(i in 1:nrow(reseq_snp012.adj)){
  gbsID<-rownames(reseq_snp012.adj)[i]
  sampleID<-info.sub$Sample_ID[which(info.sub$GBS==gbsID)]
  if(length(sampleID)>0){
    rownames(reseq_snp012.adj)[i]<-sampleID
  }
  
}
#rownames(reseq_snp012.adj)

##### start the comparison
cm_snp<-intersect(colnames(reseq_snp012.adj),colnames(div.adj))
cm_taxa<-intersect(rownames(reseq_snp012.adj),rownames(div.adj))

geno.gbs<-reseq_snp012.adj[which(rownames(reseq_snp012.adj) %in% cm_taxa),which(colnames(reseq_snp012.adj) %in% cm_snp)]
geno.dart<-div.adj[which(rownames(div.adj) %in% cm_taxa),which(colnames(div.adj) %in% cm_snp)]

which(colnames(geno.dart)!=colnames(geno.gbs))
which(rownames(geno.dart)!=rownames(geno.gbs))

geno.dart<-t(geno.dart)
geno.gbs<-t(geno.gbs)

IBS_pairs<-c()
NM_markers<-c()
ONE_dosage_diff<-c()
TWO_dosage_diff<-c()
for(i in 1:length(cm_taxa)){
  geno.1<-geno.dart[,i]
  
  geno.2<-geno.gbs[,i] 
  geno.11<-geno.1[which(!is.na(geno.1)&!is.na(geno.2))]
  geno.22<-geno.2[which(!is.na(geno.1)&!is.na(geno.2))]
  ibs<-length(which(geno.11==geno.22))/length(geno.11)
  IBS_pairs<-c(IBS_pairs,ibs)
  
  nm_markers<-length(geno.11)
  NM_markers<-c(NM_markers,nm_markers)
  
  dosage_diff<-abs(geno.11-geno.22)
  one_dosage_diff<-length(dosage_diff[which(dosage_diff==1)])
  two_dosage_diff<-length(dosage_diff[which(dosage_diff==2)])
  ONE_dosage_diff<-c(ONE_dosage_diff,one_dosage_diff)
  TWO_dosage_diff<-c(TWO_dosage_diff,two_dosage_diff)
}

colnames(IBS_pairs)[1]<-"SampleID"
write.table(IBS_pairs,"IndIBS_wtGBS_lines_commonSNP.txt",quote=FALSE,col.names=T,row.names=T,sep="\t")

####
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")
IBS_pairs<-fread("IndIBS_wtGBS_lines_commonSNP.txt",data.table=F)

pdf("./plot/Hist_IndIBS_wtGBS_lines_commonSNP.pdf",width=6,height=5,family="sans")
hist(IBS_pairs$IBS_pairs,col="cyan4",xlab="Line-based identical-by-state genotype ratio",cex.lab=1.3,main=NULL,breaks=30)
dev.off()

## Lac_0017 is not aligned well 

### pairwise IBS
colnames(geno.dart)<-paste(colnames(geno.dart),"dart",sep="_")
colnames(geno.gbs)<-paste(colnames(geno.gbs),"gbs",sep="_")
test_geno<-cbind(geno.dart,geno.gbs)

nm_sample<-ncol(test_geno)
IBS.reseq<-matrix(0,nrow=nm_sample,ncol=nm_sample)
colnames(IBS.reseq)<-colnames(test_geno)
rownames(IBS.reseq)<-colnames(test_geno)
for(i in 1:nm_sample){
  geno.1<-test_geno[,i]
  for(j in 1:nm_sample){
    geno.2<-test_geno[,j] 
    geno.11<-geno.1[which(!is.na(geno.1)&!is.na(geno.2))]
    geno.22<-geno.2[which(!is.na(geno.1)&!is.na(geno.2))]
    ibs<-length(which(geno.11==geno.22))/length(geno.11)
    IBS.reseq[i,j]<-ibs
  }
}
getwd()
#write.table(IBS.reseq,"pwIndIBS_47tGBS_lines_2.5KSNP.txt",quote=FALSE,col.names=T,row.names=T,sep="\t")
#write.table(t(test_geno),"47_tGBS_DArT_refDosage_2.5K.txt",quote=FALSE,col.names=T,row.names=T,sep="\t")
write.table(IBS.reseq,"pwIndIBS_wtGBS_lines_commonSNP.txt",quote=FALSE,col.names=T,row.names=T,sep="\t")
write.table(t(test_geno),"tGBS_DArT_refDosage_commonSNP_geno.txt",quote=FALSE,col.names=T,row.names=T,sep="\t")

#pdf("../plot/heatmap_pwIndIBS_47tGBS_lines_2.5KSNP.pdf",width=16,height=16)
pdf("./plot/heatmap_pwIndIBS_wtGBS_lines_commonSNP.pdf",width=16,height=16)
par(mar = c(16, 4, 4, 16))
heatmap(IBS.reseq)
dev.off()

########################
# neibourgh joining tree
########################
library(ape)
dist<-as.matrix(1-IBS.reseq)
tr<-nj(dist)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2")
#write.tree(tr, file = "pwIndIBS_47tGBS_lines_2.5KSNP_tree.txt", append = FALSE, digits = 10, tree.names = T)
write.tree(tr, file = "pwIndIBS_wtGBS_lines_commonSNP_tree.txt", append = FALSE, digits = 10, tree.names = T)

########### compare repeated samples ######
info<-fread("../Validation/passport/passport_PCA.txt",data.table=F)
info$GBS<-gsub("-","_",info$GBS,fixed=T)
info$tempID<-paste(info$Sample_ID,info$GBS,sep=":")
tGBS.1<-info[which(info$Population=="tGBS.1"),]
#tGBS.1$GBS<-paste(tGBS.1$GBS,"A",sep="")
#tGBS.1<-tGBS.1[-which(tGBS.1$GBS=="GBS-330A"),] ## this sample was filtered out

tGBS.2<-info[which(info$Population=="tGBS.2"),]
#tGBS.2$GBS<-paste(tGBS.2$GBS,"A",sep="")

tGBS.1$GBS[-which(tGBS.1$GBS %in% rownames(reseq_snp012.adj))]

tGBS.1.sub<-tGBS.1[which(tGBS.1$GBS %in% tGBS.2$GBS),]
tGBS.1.geno<-div.adj[which(rownames(div.adj) %in% tGBS.1.sub$Sample_ID),]
tGBS.2.geno<-div.adj[which(rownames(div.adj) %in% tGBS.2$Sample_ID),]




GBSID.1<-c()
for(i in 1:nrow(tGBS.1.geno)){
  gbsID<-tGBS.1.sub$GBS[which(tGBS.1.sub$Sample_ID==rownames(tGBS.1.geno)[i])]
  tempID<-tGBS.1.sub$tempID[which(tGBS.1.sub$Sample_ID==rownames(tGBS.1.geno)[i])]
  rownames(tGBS.1.geno)[i]<-tempID
  GBSID.1<-c(GBSID.1,gbsID)
}
GBSID.2<-c()
for(i in 1:nrow(tGBS.2.geno)){
  gbsID<-tGBS.2$GBS[which(tGBS.2$Sample_ID==rownames(tGBS.2.geno)[i])]
  tempID<-tGBS.2$tempID[which(tGBS.2$Sample_ID==rownames(tGBS.2.geno)[i])]
  rownames(tGBS.2.geno)[i]<-tempID
  GBSID.2<-c(GBSID.2,gbsID)
}

which(colnames(tGBS.1.geno)!=colnames(tGBS.2.geno))
rep_taxa<-intersect(GBSID.1,GBSID.2)
rep_taxa_temp.1<-tGBS.1.sub$temp[which(tGBS.1.sub$GBS %in% rep_taxa)]
rep_taxa_temp.2<-tGBS.2$temp[which(tGBS.2$GBS %in% rep_taxa)] 
tGBS.1.geno<-tGBS.1.geno[which(rownames(tGBS.1.geno) %in% rep_taxa_temp.1),]
tGBS.2.geno<-tGBS.2.geno[which(rownames(tGBS.2.geno) %in% rep_taxa_temp.2),]

test_geno<-rbind.data.frame(tGBS.1.geno,tGBS.2.geno)
test_geno<-t(test_geno)
dim(test_geno) #5346 SNPs

nm_sample<-length(rep_taxa)
IBS_pairs<-c()
for(i in 1:nm_sample){
  geno.1<-test_geno[,i]
  
  geno.2<-test_geno[,(i+nm_sample)] 
  geno.11<-geno.1[which(!is.na(geno.1)&!is.na(geno.2))]
  geno.22<-geno.2[which(!is.na(geno.1)&!is.na(geno.2))]
  ibs<-length(which(geno.11==geno.22))/length(geno.11)
  IBS_pairs<-c(IBS_pairs,ibs)
  
}
#pdf("../plot/Hist_IndIBS_18repeated_tGBS_lines_9KSNP.pdf",width=16,height=16)
hist(IBS_pairs)
dev.off()
IBS_pairs<-cbind.data.frame(rep_taxa,IBS_pairs)
#write.table(IBS_pairs,"IndIBS_18repeated_tGBS_lines_9KSNP.txt",quote=FALSE,col.names=T,row.names=T,sep="\t")
write.table(IBS_pairs,"IndIBS_repeated_tGBS_lines_commonSNP.txt",quote=FALSE,col.names=T,row.names=T,sep="\t")
## GBS-539 seems to be contaminated? With heterozygosity > 0.3.

IBS_pairs<-fread("IndIBS_repeated_tGBS_lines_commonSNP.txt",data.table=F)

pdf("./plot/Hist_IndIBS_repeated_tGBS_lines_commonSNP.pdf",width=5,height=5)
hist(IBS_pairs$IBS_pairs,main="Alignment: Rep1 vs Rep2",xlab="IBS of taxa pairs",breaks=10)
dev.off()



nm_sample<-ncol(test_geno)
IBS.reseq<-matrix(0,nrow=nm_sample,ncol=nm_sample)
colnames(IBS.reseq)<-colnames(test_geno)
rownames(IBS.reseq)<-colnames(test_geno)
for(i in 1:nm_sample){
  geno.1<-test_geno[,i]
  for(j in 1:nm_sample){
    geno.2<-test_geno[,j] 
    geno.11<-geno.1[which(!is.na(geno.1)&!is.na(geno.2))]
    geno.22<-geno.2[which(!is.na(geno.1)&!is.na(geno.2))]
    ibs<-length(which(geno.11==geno.22))/length(geno.11)
    IBS.reseq[i,j]<-ibs
  }
}
getwd()
#write.table(IBS.reseq,"pwIndIBS_18repeated_tGBS_lines_9KSNP.txt",quote=FALSE,col.names=T,row.names=T,sep="\t")
#write.table(t(test_geno),"18repeated_tGBS_lines_DArT_refDosage_9K.txt",quote=FALSE,col.names=T,row.names=T,sep="\t")
write.table(IBS.reseq,"pwIndIBS_repeated_tGBS_lines_commonSNP.txt",quote=FALSE,col.names=T,row.names=T,sep="\t")
#write.table(t(test_geno),"18repeated_tGBS_lines_DArT_refDosage_2.3K.txt",quote=FALSE,col.names=T,row.names=T,sep="\t")
########################
# neibourgh joining tree
########################
library(ape)
dist<-as.matrix(1-IBS.reseq)
tr<-nj(dist)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2")
#write.tree(tr, file = "pwIndIBS_18repeated_tGBS_lines_9KSNP_tree.txt", append = FALSE, digits = 10, tree.names = T)
write.tree(tr, file = "pwIndIBS_repeated_tGBS_lines_commonSNP_tree.txt", append = FALSE, digits = 10, tree.names = T)



#### mimic F1's
library(stringr)
### diversity lines
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/")
#setwd("/workdir/ml2498/Lettuce/Genotyping_2022/Realign_validation2/")
div<-fread("final_combined_snps_div_final.012",data.table=F) ## dosage of alt alleles
div_ind<-read.table("final_combined_snps_div_final.012.indv",header=F)
div_pos<-fread("final_combined_snps_div_final.012.pos",data.table=F)
colnames(div_pos)<-c("Chr","Pos")
div_pos$SNP<-paste(div_pos$Chr,div_pos$Pos,sep="_")
div<-div[,-1]
rownames(div)<-div_ind[,1]
colnames(div)<-div_pos$SNP

## missing data -1 to NA
for(i in 1:ncol(div)){
  div[which(div[,i]==-1),i]<-NA
}
## change from alt dosage to ref dosage
workingset<-2-div

###
others<-fread("final_combined_snps_other_final.012",data.table=F)
others_ind<-read.table("final_combined_snps_other_final.012.indv",header=F)
others_pos<-fread("final_combined_snps_other_final.012.pos",data.table=F)
colnames(others_pos)<-c("Chr","Pos")
others_pos$SNP<-paste(others_pos$Chr,others_pos$Pos,sep="_")
others<-others[,-1]
rownames(others)<-others_ind[,1]
colnames(others)<-others_pos$SNP

## missing data -1 to NA
for(i in 1:ncol(others)){
  others[which(others[,i]==-1),i]<-NA
}
## change from alt dosage to ref dosage
workingset2<-2-others


####
pop_str<-fread("../Validation/passport/passport_PCA.txt",data.table=F)

pop_F1<-pop_str[which(pop_str$Population=="mimicF1"),]
parents<-strsplit(pop_F1$GBS,"\\+")
parents<-as.data.frame(matrix(unlist(parents),ncol=2,byrow=T))
colnames(parents)<-c("P1","P2")
pop_F1<-cbind.data.frame(pop_F1,parents)
pop_F1$P1<-str_trim(pop_F1$P1, side = c("both"))
pop_F1$P2<-str_trim(pop_F1$P2, side = c("both"))
#str(pop_str)
#str(pop_F1)

### process accession names to remove GBS ID

GBSID<-strsplit(pop_str$Accession,split="\\[")
gbsID<-c()
for(i in 1:length(GBSID)){
  temp<-GBSID[[i]][1]
  temp<-str_trim(temp, side = c("both"))
  gbsID<-c(gbsID,temp)
}
pop_str$gbsID<-gbsID

uniq_parents<-unique(c(pop_F1$P1,pop_F1$P2)) #40 parents
parent_sampleID.1<-pop_str[which(pop_str$gbsID %in% uniq_parents),] ## non-tGBS panel
parent_sampleID.2<-pop_str[which(pop_str$GBS %in% uniq_parents),] ## tGBS panel

IDtable1<-parent_sampleID.1[,c(3,7)]
colnames(IDtable1)<-c("SampleID","gbsID")
IDtable2<-parent_sampleID.2[,c(3,2)]
colnames(IDtable2)<-c("SampleID","gbsID")
IDtable<-rbind.data.frame(IDtable1,IDtable2)
i=3
i=1
for(i in 1:nrow(IDtable)){
  working_index<-which(rownames(workingset)==IDtable$SampleID[i])
  if(length(working_index)>0){
    rownames(workingset)[working_index]<-IDtable$gbsID[i]
    #print(paste(IDtable$SampleID[i],":", IDtable$gbsID[i]))
  }else{
    working_index2<-which(rownames(workingset2)==IDtable$SampleID[i])
    if(length(working_index2)>0){
      rownames(workingset2)[working_index2]<-IDtable$gbsID[i]
    }else{
      print(paste("Did not find accession for ",IDtable$gbsID[i],sep=""))
    }
  }
}
# [1] "Did not find accession for IVT 280"
# [1] "Did not find accession for PI 253597"
# [1] "Did not find accession for PI 261653"
# [1] "Did not find accession for PI 490999"
# [1] "Did not find accession for PI 491204"
# [1] "Did not find accession for SAL 012"
# [1] "Did not find accession for GBS-329"
# [1] "Did not find accession for GBS-330"

pop_F1.test<-pop_F1[-which(pop_F1$P2 %in% c("IVT 280","PI 253597","PI 261653","PI 490999","PI 491204","GBS-329","SAL 012","GBS-330")),]

cm_snp<-intersect(colnames(workingset),colnames(workingset2))
WS1<-workingset[which(rownames(workingset) %in% IDtable$gbsID),which(colnames(workingset) %in% cm_snp)]
WS2<-workingset2[which(rownames(workingset2) %in% IDtable$gbsID),which(colnames(workingset2) %in% cm_snp)]
WS3<-workingset2[which(rownames(workingset2) %in% pop_F1.test$Sample_ID),which(colnames(workingset2) %in% cm_snp)]
which(colnames(WS1)!=colnames(WS2))
WS<-rbind.data.frame(WS1,WS2,WS3)

## check heterozygosity
HET<-c()
i=1
for(i in 1:nrow(WS)){
  temp<-WS[i,]
  nm_snp<-length(temp[!is.na(temp[1,])])
  het<-length(temp[1,which(temp[1,]==1)])
  HET<-c(HET,het/nm_snp)
}
HET<-cbind.data.frame(rownames(WS),HET)
hist(HET$HET,breaks=10) # ***Lac_0327 Lac_0354 heterozygosity > 0.25
#### start exam genotype:
## 1. extract P1 $ P2 genotype
## 2. extract F1 genotype
## 3. for each loci that P1!=P2, check if F1 is 1 (het)
i=1
#NM_snp<-c()
#SEG_het<-c()
#SEG_homo<-c()
ALL_samples<-c()
i=9
for(i in 1:nrow(pop_F1.test)){
  sampleID<-pop_F1.test$Sample_ID[i]
  P1<-pop_F1.test$P1[i]
  P2<-pop_F1.test$P2[i]
  F1<-pop_F1.test$Sample_ID[i]
  
  genoP1<-WS[which(rownames(WS)==P1),]
  genoP2<-WS[which(rownames(WS)==P2),]
  genoF1<-WS[which(rownames(WS)==F1),]
  geno<-rbind.data.frame(genoP1,genoP2,genoF1)
  
  rm_snp<-c()
  #j=19
  for (j in 1:ncol(geno)){
    temp_geno<-geno[,j]
    #if(anyNA(temp_geno) | temp_geno[1]==temp_geno[2]){ # if any missing data or no polymorphysm
    if(anyNA(temp_geno)){ # if any missing data
      rm_snp<-c(rm_snp,colnames(geno)[j])
    }
  }
  
  ## need to have P1, P2 and F1 all exist
  if(nrow(genoP1)>0&nrow(genoP2)>0&nrow(genoF1)>0){
    
    geno.test<-geno[,-which(colnames(geno) %in% rm_snp)]
    nm_snp<-ncol(geno.test)
    
    # P1==P2
    no_seg<-length(which(geno.test[1,]==geno.test[2,]))
    # 0,2->1; 1,0->1; 1,2->1
    seg02_1<-length(which((geno.test[1,]==0&geno.test[2,]==2&geno.test[3,]==1)|(geno.test[1,]==2&geno.test[2,]==0&geno.test[3,]==1)))
    seg01_1<-length(which((geno.test[1,]==0&geno.test[2,]==1&geno.test[3,]==1)|(geno.test[1,]==1&geno.test[2,]==0&geno.test[3,]==1)))
    seg12_1<-length(which((geno.test[1,]==1&geno.test[2,]==2&geno.test[3,]==1)|(geno.test[1,]==2&geno.test[2,]==1&geno.test[3,]==1)))
    # 1,0->0; 1,2->2
    seg01_0<-length(which((geno.test[1,]==0&geno.test[2,]==1&geno.test[3,]==0)|(geno.test[1,]==1&geno.test[2,]==0&geno.test[3,]==0)))
    seg12_2<-length(which((geno.test[1,]==1&geno.test[2,]==2&geno.test[3,]==2)|(geno.test[1,]==2&geno.test[2,]==1&geno.test[3,]==2)))
    ## likely errors
    rest_cases<-nm_snp-(no_seg+seg02_1+seg01_1+seg12_1+seg01_0+seg12_2)

    one_sample<-c(sampleID,nm_snp,no_seg,seg02_1,seg01_1,seg12_1,seg01_0,seg12_2,rest_cases)
    ALL_samples<-rbind.data.frame(ALL_samples,one_sample)
  }

  #NM_snp<-c(NM_snp,nm_snp)
  #SEG_het<-c(SEG_het,seg_het)
  #SEG_homo<-c(SEG_homo,seg_homo)
}
colnames(ALL_samples)<-c("Sample_ID","nm_snp","P1=P2","seg02_1","seg01_1","seg12_1","seg01_0","seg12_2","rest_cases")
pop_F1.test<-merge(pop_F1.test,ALL_samples,by="Sample_ID",all.x=T)

write.table(pop_F1.test,"mimicF1_genotype_comparison.txt",row.names=F,sep="\t",quote=F)

####### plot the mimic F1 results
library(reshape2)
library(ggplot2)
library(data.table)
pop_F1.test<-fread("mimicF1_genotype_comparison.txt",data.table=F)
pop_F1.test<-pop_F1.test[!is.na(pop_F1.test$nm_snp),]

pop_F1.plot0<-pop_F1.test[,c(3,9:16)]

pop_F1.plot<-pop_F1.plot0
for (i in 3:ncol(pop_F1.plot)){
  pop_F1.plot[,i]<-pop_F1.plot[,i]/pop_F1.plot[,2]
}
pop_F1.plot<-pop_F1.plot[,-2]
pop_F1.plot1<-melt(data=pop_F1.plot,id.var=c("GBS"),variable.name="Genotype_Group")
pop_F1.plot1$plot_label<-as.character(pop_F1.plot1$Genotype_Group)
pop_F1.plot1$plot_label[which(pop_F1.plot1$plot_label=="P1=P2")]<-"Parent1 = Parent2"
pop_F1.plot1$plot_label[which(pop_F1.plot1$plot_label=="seg02_1")]<-"Parents 0 or 2; Progenies 1"
pop_F1.plot1$plot_label[which(pop_F1.plot1$plot_label=="seg01_1")]<-"Parents 0 or 1; Progenies 1"
pop_F1.plot1$plot_label[which(pop_F1.plot1$plot_label=="seg12_1")]<-"Parents 1 or 2; Progenies 1"
pop_F1.plot1$plot_label[which(pop_F1.plot1$plot_label=="seg01_0")]<-"Parents 0 or 1; Progenies 0"
pop_F1.plot1$plot_label[which(pop_F1.plot1$plot_label=="seg12_2")]<-"Parents 1 or 2; Progenies 2"
pop_F1.plot1$plot_label[which(pop_F1.plot1$plot_label=="rest_cases")]<-"Other cases"

pop_F1.plot1$plot_label<-factor(pop_F1.plot1$plot_label,
                                levels=c("Parent1 = Parent2","Parents 0 or 2; Progenies 1","Parents 0 or 1; Progenies 1",
                                         "Parents 1 or 2; Progenies 1","Parents 0 or 1; Progenies 0",
                                         "Parents 1 or 2; Progenies 2","Other cases"))

#pop_F1.plot1$Genotype_Group<-factor(pop_F1.plot1$Genotype_Group,levels=c("P1=P2","seg02_1","seg01_1","seg12_1","seg01_0","seg12_2","rest_cases"))

pop_F1.plot1$GBS<-as.factor(pop_F1.plot1$GBS)

p<-ggplot(data = pop_F1.plot1, aes(x=GBS, y=value, fill = plot_label)) + 
  geom_bar(stat='identity')+
  scale_fill_brewer(direction = -1)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        legend.title = element_text( size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.5, 'cm'))+
  labs(x ="Mimic F1 progenies", y = "Proportion of genotype groups",fill="Genotype group")
p

setwd("/Users/ml2498/Desktop/Labserver/Lettuce/Genotyping_2022/Realign_validation2/plot")
pdf("Genotype_Proportion_mimicF1.pdf",width=15,height=9,family="sans")
print(p)
dev.off()

