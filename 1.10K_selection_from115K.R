#Lettuce.all.snps.reanalysis.genotype.txt: 0 - missing; 1 & 2 - homozygous; 3 - heterozygous

library(data.table)

let_115<-fread("/Users/ml2498/Desktop/Labserver/Lettuce/tGBS_SNPs/Lettuce.all.snps.reanalysis.genotype.txt",data.table=F)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection")
taxa<-fread("Core_selection_forBreeding_77_tGBS_ID.txt",data.table=F)
let_115sub<-let_115[,c(1:4,which(colnames(let_115) %in% taxa[,1]))]
#taxa[-which(taxa[,1] %in% colnames(let_115sub)),]

write.table(let_115sub,"./from115K/Lettuce_115KSNP_selected_lines.txt",row.names=F,sep="\t",quote=F)
##########################
setwd("/workdir/ml2498/Lettuce/10K_selection/from115K")
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/from115K")
let_115sub<-fread("Lettuce_115KSNP_selected_lines.txt",data.table = F)
nm_snp<-nrow(let_115sub)
missing_ct<-colSums(let_115sub==0)
homo_ct1<-colSums(let_115sub==1)
homo_ct2<-colSums(let_115sub==2)
het_ct<-colSums(let_115sub==3)

stat_115<-rbind.data.frame(missing_ct,homo_ct1,homo_ct2,het_ct)
colnames(stat_115)<-colnames(let_115sub)
stat_115<-stat_115[,-(1:4)]
rownames(stat_115)<-c("missing_ct","homo_ct1","homo_ct2","het_ct")
stat_115<-as.data.frame(t(stat_115))
stat_115$missing_rate<-stat_115$missing_ct/nm_snp
stat_115$het_rate<-stat_115$het_ct/nm_snp

pdf("Taxa_stats_115K.pdf")
hist(stat_115$missing_rate)
hist(stat_115$het_rate)
dev.off()
###################################
rownames(let_115sub)<-paste(let_115sub[,1],let_115sub[,2],sep="_")
snp_info<-let_115sub[,1:4]
let_115sub<-let_115sub[,-(1:4)]
let_115sub<-as.data.frame(t(let_115sub))
nm_taxa<-nrow(let_115sub)

missing_ct<-colSums(let_115sub==0)
homo_ct1<-colSums(let_115sub==1)
homo_ct2<-colSums(let_115sub==2)
het_ct<-colSums(let_115sub==3)

stat_115.snp<-cbind.data.frame(missing_ct,homo_ct1,homo_ct2,het_ct)
rownames(stat_115.snp)<-colnames(let_115sub)

stat_115.snp$missing_rate<-stat_115.snp$missing_ct/nm_taxa
stat_115.snp$het_rate<-stat_115.snp$het_ct/nm_taxa

which(stat_115.snp$homo_ct2 > stat_115.snp$homo_ct1)
stat_115.snp$MAF<-0
i=1
i=3
for(i in 1:nrow(stat_115.snp)){
  temp<-stat_115.snp[i,2:3]
  minor<-which.min(temp)
  
  minor_allele_c<-(2*temp[minor]+stat_115.snp[i,4])/(76*2)
  stat_115.snp$MAF[i]<-minor_allele_c[1,1]
  
}
stat_115.snp<-cbind.data.frame(snp_info,stat_115.snp)
write.table(stat_115.snp,"SNP_stats_115K.txt",row.names=T,sep="\t",quote=F)
write.table(stat_115,"Taxa_stats_115K.txt",row.names=F,sep="\t",quote=F)

pdf("SNP_stats_115K.pdf")
hist(stat_115.snp$missing_rate)
hist(stat_115.snp$het_rate)
hist(stat_115.snp$MAF)

dev.off()
##############################################

####################################
# keep SNPs in gene or CDS regions 
####################################
library(data.table) 

setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/from115K")
setwd("/workdir/ml2498/Lettuce/10K_selection/from115K")

stat_115.snp<-fread("SNP_stats_115K.txt",data.table=F)
colnames(stat_115.snp)[1]<-"ID"

mask_anno<-"/workdir/data/lettuce/Ref/Phytozome/PhytozomeV12/early_release/Lsativa_467_v5/annotation/Lsativa_467_v5.gene_exons.gff3"
gff<-read.delim(mask_anno,header=F,sep="\t")
REGION<-fread("../SNP_counts_50K_1M_window.txt",data.table=F)

gff$V1<-gsub("Lsat_1_v8_lg","LSAT_1_V8_LG",gff$V1)
gff.gene<-gff[which(gff$V3=="gene"),]
gff.CDS<-gff[which(gff$V3=="CDS"),]
## V1=chr, V4=ST, V5=END

#### plot gene models ####
library(ggplot2)
gff.gene$y0<-0
gff.gene$y1<-1
gff.gene1<-gff.gene
gff.gene1$V1<-gsub("LSAT_1_V8_","",gff.gene1$V1,fixed=T)
gff.gene1<-gff.gene1[which(gff.gene1$V1 %in% paste("LG_",1:9,sep="")),]

p<-ggplot()+
  geom_rect(data=gff.gene1, mapping=aes(xmin=V4, xmax=V5, ymin=y0, ymax=y1), fill="blue", color=NA) +
  facet_grid(rows = vars(V1))+
  labs(x="Physical position",y="",title="Genic regions (V8)")+
  theme_light()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_text(size=10,colour = 'black'),
        strip.text.x = element_text(size=5,colour = 'black'),
        strip.placement = "outside",
        #axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=7),
        text = element_text(size=10),
        #legend.title=element_blank(),
        legend.text.align = 0,
        legend.text=element_text(size=rel(1)))
  
pdf("V8_gene_model_regions.pdf",height=8,width=12)
print(p)
dev.off()

#############################

stat_115.snp$Chr<-gsub("Lsat_1_v8_lg_","LSAT_1_V8_LG_",stat_115.snp$Chr,fixed=T)
keep.gene<-c()
keep.CDS<-c()
for (i in 1:9){
  chr<-paste("LSAT_1_V8_LG_",i,sep="")
  print(chr)
  let.chr<-stat_115.snp[which(stat_115.snp$Chr==chr),1:3]
  gff.gene.chr<-gff.gene[which(gff.gene$V1==chr),]
  gff.CDS.chr<-gff.CDS[which(gff.CDS$V1==chr),]
  
  #for(j in 1:1000){
  for(j in 1:nrow(let.chr)){
    pos<-let.chr[j,3]
    temp.gene<-which(gff.gene.chr$V4<=pos & gff.gene.chr$V5>=pos)
    temp.CDS<-which(gff.CDS.chr$V4<=pos & gff.CDS.chr$V5>=pos)
    
    if(length(temp.gene)>0){
      #print(paste(j,": gene",sep=""))
      keep.gene<-rbind.data.frame(keep.gene,let.chr[j,])
      
    }
    if(length(temp.CDS)>0){
      #print(paste(j,": CDS",sep=""))
      keep.CDS<-rbind.data.frame(keep.CDS,let.chr[j,])
    }
  }
}

dim(keep.gene) # 4332,3
dim(keep.CDS) # 2298,3
getwd()
write.table(keep.gene,"Retained_SNPs_GenicRegion.txt",row.names=F,sep="\t",quote=F)

nm_SNP<-table(keep.gene$Chr)
pdf("NM_SNPs_byChr_inGenicRegion.pdf")
barplot(nm_SNP)
dev.off()
#######################
# no SNPs within 50bp up/down-stream a SNP
#######################
library(data.table)
## Assuming the 115K is the full SNP set, to calculate distance between SNPs
stat_115.snp<-fread("SNP_stats_115K.txt",data.table=F)
colnames(stat_115.snp)[1]<-"ID"

stat_115.snp$info<-paste(stat_115.snp$Chr,stat_115.snp$Pos,sep="-")
stat_115.snp[,2]<-gsub("Lsat_1_v8_lg_","LSAT_1_V8_LG_",stat_115.snp[,2],fixed=T)
head(stat_115.snp)

####### get info for extra variants ######
index<-which(stat_115.snp$ID=="Lsat_1_v8_lg_1_121673947")
stat_115.snp[(index-3):index,]

index<-which(stat_115.snp$ID=="Lsat_1_v8_lg_2_1081821")
stat_115.snp[(index-3):index,]

index<-which(stat_115.snp$ID=="Lsat_1_v8_lg_4_105718090")
stat_115.snp[(index-3):index,]

index<-which(stat_115.snp$ID=="Lsat_1_v8_lg_5_80705220")
stat_115.snp[(index-3):index,]

index<-which(stat_115.snp$ID=="Lsat_1_v8_lg_5_291392586")
stat_115.snp[(index-3):(index+3),]

index<-which(stat_115.snp$ID=="Lsat_1_v8_lg_7_161895696")
stat_115.snp[(index-3):index,]

index<-which(stat_115.snp$ID=="Lsat_1_v8_lg_9_151974232")
stat_115.snp[(index-3):index,]


### get distance for adjacent SNPs
stat_115.snp<-stat_115.snp[order(stat_115.snp$Chr,stat_115.snp$Pos),]

stat_115.snp$Dist<-0
stat_115.snp.1<-c()

for (i in 1:9){
  chr<-paste("LSAT_1_V8_LG_",i,sep="")
  print(chr)
  temp<-stat_115.snp[which(stat_115.snp$Chr==chr),]
  for(j in 2:nrow(temp)){
    temp$Dist[j]<-temp$Pos[j]-temp$Pos[j-1]
  }
  stat_115.snp.1<-rbind.data.frame(stat_115.snp.1,temp)
}
###############################
### selection based on distance
###############################

interval<-c(30,50,70)
stat_115.snp.1$interval30bp<-"N"
stat_115.snp.1$interval50bp<-"N"
stat_115.snp.1$interval70bp<-"N"

stat_115.snp.2<-c()

for (c in 1:9){
  chr<-paste("LSAT_1_V8_LG_",c,sep="")
  print(chr)
  temp<-stat_115.snp.1[which(stat_115.snp.1$Chr==chr),]
  
  for(i in interval){
    COL_nm<-which(colnames(temp)==paste("interval",i,"bp",sep=""))
    
    for(j in 1:(nrow(temp)-1)){
      ### first row
      if(j==1){
        if(temp$Pos[j]>=i & temp$Dist[j+1]>=i){
          temp[j,COL_nm]<-"Y"
        }
      }else{
        if(temp$Dist[j]>=i & temp$Dist[j+1]>=i){
          temp[j,COL_nm]<-"Y"
        }
      }
    }
  }
  stat_115.snp.2<-rbind.data.frame(stat_115.snp.2,temp)
}

length(which(stat_115.snp.2$interval30bp=="Y")) #23633
length(which(stat_115.snp.2$interval50bp=="Y")) # 15895
length(which(stat_115.snp.2$interval70bp=="Y")) # 12414
write.table(stat_115.snp.2,"Retained_SNP_bySNPinterval_115K.txt",row.names=F,sep="\t",quote=F)

################################################
# selection based on MAF and missing rate
################################################
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/from115K")
stat_115.snp.2 <- fread("Retained_SNP_bySNPinterval_115K.txt", data.table=F)
#taxaSummaryTable <- fread("TaxaSum_for50K.txt", sep="\t", data.table=F)

MAF<-0.05
#77*0.05=3.85
Missing<-0.75

keep2<-stat_115.snp.2[which(stat_115.snp.2$missing_rate<=Missing & stat_115.snp.2$MAF>=MAF),] #36,614
length(keep2$info[which(keep2$interval30bp=="Y")]) #8599
keep2.int30<-keep2[which(keep2$interval30bp=="Y"),]
write.table(keep2.int30,"Retained_8.5KSNP_30bpFlanking_MAF_missing.txt",row.names=F,sep="\t",quote=F)

####################
# check MAS markers
#####################
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/MAS")
MAS<-fread("MAS_in_115K.txt",data.table=F)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/from115K")
keep2.int30<-fread("Retained_8.5KSNP_30bpFlanking_MAF_missing.txt",data.table=F)

MAS$If_in_8.5K<-""
MAS$If_in_8.5K[which(MAS$SNP_ID %in% keep2.int30$ID)]<-"Y"
MAS$SNP_115K<-""
MAS$SNP_115K[which(MAS$SNP_ID %in% stat_115.snp.2$ID)]<-"Y"

MAS$interval_30<-""
MAS$interval_30[which(MAS$SNP_ID %in% stat_115.snp.2$ID[which(stat_115.snp.2$interval30bp=="Y")])]<-"Y"

MAS$Chr<-""
MAS$Pos<-0
MAS$REF<-"N"
MAS$ALT<-"N"
for(i in 1:nrow(MAS)){
  ID<-MAS$SNP_ID[i]
  if(ID!="Lsat_1_v8_lg_4_109726844"){
    MAS$REF[i]<-stat_115.snp.2$REF[which(stat_115.snp.2$ID==ID)]
    MAS$ALT[i]<-stat_115.snp.2$ALT[which(stat_115.snp.2$ID==ID)]
    MAS$Chr[i]<-stat_115.snp.2$Chr[which(stat_115.snp.2$ID==ID)]
    MAS$Pos[i]<-stat_115.snp.2$Pos[which(stat_115.snp.2$ID==ID)]
    
  }
}
MAS$Chr<-gsub("LSAT_1_V8_LG_","Lsat_1_v8_lg_",MAS$Chr)
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/MAS")
write.table(MAS,"MAS_in115K_for_key_file.txt",row.names=F,sep="\t",quote=F)

all_to_check<-c()
for (i in 1:nrow(MAS)){
  if(MAS$interval_30[i]==""&MAS$SNP_ID[i]!="Lsat_1_v8_lg_4_109726844"){
    ID<-MAS$SNP_ID[i]
    row_index<-which(stat_115.snp.2$ID==ID)
    check<-stat_115.snp.2[(row_index-5):(row_index+5),]
    all_to_check<-rbind.data.frame(all_to_check,check)
  }
}

setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/MAS")
write.table(all_to_check,"MAS_marker_dist_withOtherSNPs.txt",row.names=F,sep="\t",quote=F)


#########################
#### 1M window plot ####
#########################
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection")
REGION<-read.table("SNP_counts_50K_1M_window.txt",header=T)
#REGION$Chr<-gsub("LSAT_1_V8_","",REGION$chr)

# keep2: MAF 0.05; missing 0.6 (28,472 SNPs)
INT30<-stat_115.snp.2[which(stat_115.snp.2$interval30bp=="Y"),]

REGION$counts_INT30<-0
REGION$counts_MAF_missing<-0
REGION$counts_INT30_MAF_missing<-0
#i=1
for (i in 1:nrow(REGION)){
  chr<-REGION$chr[i]
  ST<-REGION$ST[i]
  END<-REGION$END[i]
  
  temp.i<-INT30[which(INT30$Chr==chr & INT30$Pos>=ST & INT30$Pos<END),]
  if(nrow(temp.i)>0){
    ct.i<-nrow(temp.i)
    REGION$counts_INT30[i]<-ct.i
  }
  
  temp.mm<-keep2[which(keep2$Chr==chr & keep2$Pos>=ST & keep2$Pos<END),]
  if(nrow(temp.mm)>0){
    ct.mm<-nrow(temp.mm)
    REGION$counts_MAF_missing[i]<-ct.mm
  }
  
  temp.i.mm<-keep2.int30[which(keep2.int30$Chr==chr & keep2.int30$Pos>=ST & keep2.int30$Pos<END),]
  if(nrow(temp.i.mm)>0){
    ct.i.mm<-nrow(temp.i.mm)
    REGION$counts_INT30_MAF_missing[i]<-ct.i.mm
  }
}

library(ggplot2)
p1<-ggplot(REGION,aes(END,counts_INT30))+
  geom_bar(stat="identity")+
  facet_grid(rows = vars(Chr))+
  labs(x="Physical position",y="SNP counts (1Mb window; >=30bp between SNPs)")+
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

p2<-ggplot(REGION,aes(END,counts_MAF_missing))+
  geom_bar(stat="identity")+
  facet_grid(rows = vars(Chr))+
  labs(x="Physical position",y="SNP counts (1Mb window; MAF>5% & Missing rate<0.75)")+
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

p3<-ggplot(REGION,aes(END,counts_INT30_MAF_missing))+
  geom_bar(stat="identity")+
  facet_grid(rows = vars(Chr))+
  labs(x="Physical position",y="SNP counts (1Mb window; 30bp interval; MAF>5% & Missing rate<0.75)")+
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


pdf("SNP_counts_INT30bp_MAF_missing_1M_window.pdf",height=7,width=7)
print(p1)
print(p2)
print(p3)
dev.off()
setwd("/Users/ml2498/Desktop/Labserver/Lettuce/10K_selection/from115K")
write.table(REGION,"SNP_counts_115K_1MWin_INT30_MAF_missing.txt",quote=F,row.names=F,sep="\t")

getwd()

