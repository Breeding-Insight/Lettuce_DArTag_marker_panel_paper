library(data.table)
#### DArT 4K selection
Dart_pre_sel<-fread("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/Prelim_from_DArT/Lactuca_PreliminarySelection.txt",data.table=F)

#### 8K list ####
lettuce_8K<-fread("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/format_file/Lettuce_unique_alignment_tGBS_MAS_8.2K_07062022.txt",data.table=F)
lettuce_RLL<-c("Chr3_129893362","Chr4_049498375","Chr5_336147078","Chr9_063962630","Chr9_152765187")

#### 1M window counts #####
REGION<-fread("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/from115K/SNP_counts_115K_1MWin_INT30_MAF_missing.txt",data.table=F)

#### proportional sub-selection ######
lettuce_MAS<-c(lettuce_RLL,lettuce_8K$MarkerName[1:41])
lettuce_MAS.1<-lettuce_MAS[which(lettuce_MAS %in% Dart_pre_sel$MarkerName)]

Dart_pre_sel.tGBS<-Dart_pre_sel[-which(Dart_pre_sel$MarkerName %in% lettuce_MAS.1),]
Dart_pre_sel.tGBS$chr<-gsub("Chr0","",Dart_pre_sel.tGBS$Chrom)
Dart_pre_sel.tGBS$chr<-gsub("Chr","",Dart_pre_sel.tGBS$Chrom)

REGION$counts_4K<-0
#i=1
for (i in 1:nrow(REGION)){
  chr<-REGION$Chr[i]
  chr<-gsub("LG_","",chr)
  ST<-REGION$ST[i]
  END<-REGION$END[i]
  
  temp.i<-Dart_pre_sel.tGBS[which(Dart_pre_sel.tGBS$chr==chr & Dart_pre_sel.tGBS$ChromPosPhysical>=ST & Dart_pre_sel.tGBS$ChromPosPhysical<END),]
  if(nrow(temp.i)>0){
    ct.i<-nrow(temp.i)
    REGION$counts_4K[i]<-ct.i
  }
}

library(ggplot2)
p1<-ggplot(REGION,aes(END,counts_4K))+
  geom_bar(stat="identity")+
  facet_grid(rows = vars(Chr))+
  labs(x="Physical position",y="SNP counts (1Mb window; 4K)")+
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

p1

hist(REGION$counts_4K,breaks=20)

######
SEEDS<-c(2000:(1999+nrow(REGION)))
KEEP_SNP<-c()
#for(i in 1:2){
for(i in 1:nrow(REGION)){
    
  nm_4K<-REGION$counts_4K[i]
  
  seed<-SEEDS[i]
  start<-REGION$ST[i]
  end<-REGION$END[i]
  chr<-REGION$Chr[i]
  chr<-gsub("LG_","",chr)
  
  sub_snp<-Dart_pre_sel.tGBS[which(Dart_pre_sel.tGBS$chr==chr & Dart_pre_sel.tGBS$ChromPosPhysical>=start & Dart_pre_sel.tGBS$ChromPosPhysical<end),]
  snp_nm<-nrow(sub_snp)
  
  if(nm_4K!=snp_nm){break}
  
  if(nm_4K<=3&nm_4K>0){
    keep_snp<-sub_snp
    keep_snp$window_index<-i
    keep_snp$snp_nm<-snp_nm
    KEEP_SNP<-rbind.data.frame(KEEP_SNP,keep_snp)
  }else if(nm_4K>3){
    set.seed(seed)
    pick<-sample(1:snp_nm,3,replace=F)
    keep_snp<-sub_snp[pick,]
    keep_snp$window_index<-i
    keep_snp$snp_nm<-snp_nm
    KEEP_SNP<-rbind.data.frame(KEEP_SNP,keep_snp)
  }
}

table(REGION$counts_4K)

#### remove additional markers from windows based on distance

chroms<-unique(KEEP_SNP$Chrom)
KEEP_SNP.2<-c()
for(ch in chroms){
  temp<-KEEP_SNP[which(KEEP_SNP$Chrom==ch),]
  temp$dist<-0
  temp<-temp[order(temp$ChromPosPhysical),]
  for (i in 2:nrow(temp)){
    temp$dist[i]<-temp$ChromPosPhysical[i]-temp$ChromPosPhysical[(i-1)]
  }
  KEEP_SNP.2<-rbind.data.frame(KEEP_SNP.2,temp)
}

#### remove the SNP with smallest distance with other SNPs in each 3-SNP window
KEEP_SNP.3_window<-KEEP_SNP.2[which(KEEP_SNP.2$snp_nm==3),]
min_dist_snps<-c()
Windows<-unique(KEEP_SNP.3_window$window_index)

for(w in Windows){
  temp<-KEEP_SNP.3_window[which(KEEP_SNP.3_window$window_index==w),]
  temp<-temp[order(temp$ChromPosPhysical),]
  temp.sub<-temp[2:3,]
  
  min_dist_snp<-temp.sub[which.min(temp.sub$dist),]
  min_dist_snps<-rbind.data.frame(min_dist_snps,min_dist_snp)
  
}

####
KEEP_SNP.2_window<-KEEP_SNP.2[which(KEEP_SNP.2$snp_nm==2),]
min_dist_snps_2<-c()
Windows<-unique(KEEP_SNP.2_window$window_index)

for(w in Windows){
  temp<-KEEP_SNP.2_window[which(KEEP_SNP.2_window$window_index==w),]
  temp<-temp[order(temp$ChromPosPhysical),]
  temp.sub<-temp[2,]
  
  min_dist_snps_2<-rbind.data.frame(min_dist_snps_2,temp.sub)
  
}

nm_to_remove<-nrow(KEEP_SNP)+length(lettuce_MAS.1)-3000-nrow(min_dist_snps) # nrow(min_dist_snps) is number of SNPs removed from 3-SNP windows
## still need to remove 60 SNPs from 2-SNP windows
min_dist_snps_2<-min_dist_snps_2[order(min_dist_snps_2$dist),]
min_dist_snps_2_rm<-min_dist_snps_2[1:nm_to_remove,]


#### sub-set the 4K list
Dart_pre_sel.tGBS.keep<-Dart_pre_sel.tGBS[-which(Dart_pre_sel.tGBS$MarkerName %in% min_dist_snps$MarkerName),]
Dart_pre_sel.tGBS.keep<-Dart_pre_sel.tGBS.keep[-which(Dart_pre_sel.tGBS.keep$MarkerName %in% min_dist_snps_2_rm$MarkerName),]
Dart_pre_sel.tGBS.keep<-Dart_pre_sel.tGBS.keep[which(Dart_pre_sel.tGBS.keep$MarkerName %in% KEEP_SNP$MarkerName),]

Dart_pre_sel.mas<-Dart_pre_sel[which(Dart_pre_sel$MarkerName %in% lettuce_MAS.1),]

DArT3K<-rbind.data.frame(Dart_pre_sel.mas,Dart_pre_sel.tGBS.keep[,1:10])
length(unique(DArT3K$MarkerName))

### final counts in 1M window ##

DArT3K$chr<-gsub("Chr0","",DArT3K$Chrom)
DArT3K$chr<-gsub("Chr","",DArT3K$Chrom)

REGION$counts_3K<-0
#i=1
for (i in 1:nrow(REGION)){
  chr<-REGION$Chr[i]
  chr<-gsub("LG_","",chr)
  ST<-REGION$ST[i]
  END<-REGION$END[i]
  
  temp.i<-DArT3K[which(DArT3K$chr==chr & DArT3K$ChromPosPhysical>=ST & DArT3K$ChromPosPhysical<END),]
  if(nrow(temp.i)>0){
    ct.i<-nrow(temp.i)
    REGION$counts_3K[i]<-ct.i
  }
}
sum(REGION$counts_3K)

library(ggplot2)
p1<-ggplot(REGION,aes(END,counts_3K))+
  geom_bar(stat="identity")+
  facet_grid(rows = vars(Chr))+
  labs(x="Physical position",y="SNP counts (1Mb window; 3K)")+
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

p1

setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/DArT3K")
write.table(REGION,"SNP_counts_3K_1MWin.txt",quote=F,row.names=F,sep="\t")
write.table(DArT3K,"Lettuce_DArT3K_08172022.txt",quote=F,row.names=F,sep="\t")

pdf("SNP_counts_DArT3K_1M_window.pdf",height=7,width=7)
print(p1)
dev.off()
#### for a improved figure ####
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/DArT3K")
REGION<-fread("SNP_counts_3K_1MWin.txt",data.table=F)
REGION$Chr<-gsub("LG_","Chr ",REGION$Chr)

library(ggplot2)
p1<-ggplot(REGION,aes(END,counts_3K))+
  geom_bar(stat="identity",color="cyan4")+
  facet_grid(rows = vars(Chr))+
  labs(x="Physical position (bp)",y="SNP count (1Mb window)")+
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
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        #legend.title=element_blank(),
        legend.text.align = 0,
        legend.text=element_text(size=rel(1)))

p1
pdf("SNP_counts_DArT3K_1M_window_PAG.pdf",height=5.5,width=8,family="sans")
print(p1)
dev.off()
################################################

#### final PCA plot ####
cd /Users/ml2498/Desktop/Labserver/Lettuce/10K_selection
vcftools --vcf Lettuce.MCR50.biallelic.77lines.recode.vcf --positions Lettuce_3K_positions.txt --out Lettcue_3K_genotype --recode



library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/from115K/")
let_115sub<-fread("Lettuce_115KSNP_selected_lines.txt",data.table=F)
let_115sub[1:3,1:13]
let_115sub$ID<-paste(let_115sub$Chr,let_115sub$Pos,sep="-")
Let_3Kpos<-fread("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/DArT3K/Lettuce_3K_positions.txt",data.table=F)
Let_3Kpos$ID<-paste(Let_3Kpos$V1,Let_3Kpos$V2,sep="-")
let_3K<-let_115sub[which(let_115sub$ID %in% Let_3Kpos$ID),]
## make up a hmp file
let_3K$alleles<-paste(let_3K$REF,let_3K$ALT,sep="/")
let_3K$strand<-"+"
let_3K$assembly<-NA
let_3K$center<-NA
let_3K$protLSID<-NA	
let_3K$assayLSID<-NA	
let_3K$panelLSID<-NA	
let_3K$QCcode<-NA
let_3K<-let_3K[,c(81,82,1:2,83:89,5:80)]
colnames(let_3K)[1:6]<-c("rs#","alleles","chrom","pos","strand","assembly#")
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/DArT3K/")
write.table(let_3K,"Lettcue_3K_genotype.hmp.txt",row.names=F,sep="\t",quote=F)

### make up phenotype ##

Taxa<-colnames(let_3K)[12:87]
value<-sample(1:3,length(Taxa),replace=T)
pheno<-cbind.data.frame(Taxa,value)
write.table(pheno,"Lettcue_76lines_makeup_pheno.txt",row.names=F,sep="\t",quote=F)
## generate kinship using GAPIT
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

myG<-read.table("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/DArT3K/Lettcue_3K_genotype.hmp.txt",header=F, sep="\t", comment.char="")
myY<-pheno
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/DArT3K/GAPIT_kinship")
myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=6,
  Major.allele.zero=F,
  Model.selection=TRUE
)

#### plotting PCA ###
PCA<-fread("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/DArT3K/GAPIT_kinship/GAPIT.Genotype.PCA.csv",data.table=F)
PCA$taxa<-gsub("A","",PCA$taxa)
PCA$taxa<-gsub("_","-",PCA$taxa)
PCA$taxa<-gsub("MIX","",PCA$taxa)
PCA$taxa[which(PCA$taxa=="GBS-373-1")]<-"GBS-373"
PCA$taxa[which(PCA$taxa=="GBS-385-B")]<-"GBS-385B"
PCA$taxa[which(PCA$taxa=="GBS-385-W")]<-"GBS-385W"

info<-fread("/Users/ml2498/Desktop/Labserver/Lettuce/tGBS_SNPs/accession_ID_group.txt",data.table=F)
length(intersect(PCA$taxa,info$GBS_ID))
commonLine<-intersect(PCA$taxa,info$GBS_ID)
eigen_factor<-merge(info,PCA,by.x="GBS_ID",by.y="taxa",all.y=T)

library(ggplot2)
p<-ggplot(eigen_factor,aes(PC1,PC2))+
  geom_point(aes(color=Type))+ 
  theme_light()+
  labs(x=paste("PC1 (23.82%)",sep=""),y=paste("PC2 (6.75%)",sep=""),color="Type")
p

setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/DArT3K/")
pdf("PCA_76lines_3K.pdf",width=8,height=5)
print(p)
dev.off()


###### physical positions of 3K ############
library(RColorBrewer)
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/DArT3K")
dart_3K<-fread("Lettuce_DArT3K_08172022.txt",data.table=F)

map<-dart_3K[order(dart_3K$Chrom,dart_3K$ChromPosPhysical),]

#x_lim<-max(map$ChromPosPhysical)
x_lim<-4e8
## color panel, no y axis, 
color_panel<-brewer.pal(n = 8, name = "Set2")
display.brewer.pal(n = 8, name = "Set2")
color_panel<-c(color_panel,color_panel)


{pdf("Lettuce_3K_Physical_position.pdf",width=10,height=6)
  plot(NULL,xlim=c(0,x_lim),ylim=c(0,10),yaxt="n",bty="n",
       xlab="Physical position (bp)",ylab="Chromosome")
  axis(2, at=c(1:9),labels=c(1:9), las=1)
  CHR<-unique(map$chr)
  #chr=CHR[1]
  for(chr in CHR){
    #chr_nm<-as.numeric(substr(chr,4,4))
    temp<-map[which(map$chr==chr),]
    xrange<-range(temp$ChromPosPhysical)
    rect(xleft=0,xright=xrange[2],ybottom=(chr-0.35),
         ytop=(chr+0.35),xlim=c(0,max),col=color_panel[chr])
    
  }
  
  for (i in 1:nrow(map)){
    chr<-map$chr[i]
    segments(x0=map$ChromPosPhysical[i],x1=map$ChromPosPhysical[i],y0=(chr-0.35),y1=(chr+0.35),lwd=0.05)
  }
  dev.off()
}

