# source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
#        local = TRUE)
# install_phyloseq(branch = "devel")
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(agricolae)
sharedfile = "1stning16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared"
taxfile = "1stning16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy"
mothur_data <- import_mothur(mothur_shared_file = sharedfile, mothur_constaxonomy_file = taxfile)
mothur_meta_df<-read.csv("mothur_meta.csv",sep=",",header=TRUE,row.names = 1)
head(mothur_meta_df)
mothur_meta_df <- sample_data(mothur_meta_df)
moth_merge<-merge_phyloseq(mothur_data,mothur_meta_df)
colnames(tax_table(moth_merge))
colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
sample_names(moth_merge)
moth_merge_up<-prune_samples(sample_names(moth_merge),moth_merge)
sample_names(moth_merge_up)

# write.csv(otu_table(moth_merge),file = "May1816snifa_otutable.csv")
# write.csv(tax_table(moth_merge),file = "May1816snifa_taxtable.csv")

### bar plot of read depth
Read_depth<-sample_sums(moth_merge_up)
summary(Read_depth)
sample_names(mothur_meta_df)
mothur_meta<-prune_samples(sample_names(mothur_meta_df),mothur_meta_df)
mothur_meta_up<-prune_samples(sample_names(mothur_meta_df)[-c(45,46,47,36)],mothur_meta_df)#remove VCNO_1-3 and VNTN60_4
dim(mothur_meta_up)#number of rows and columns
sample_ID<- sample_names(mothur_meta_up)
mothur_meta<-data.frame(mothur_meta_up,Read_depth=sample_sums(moth_merge_up))
tiff('Sequencing depth barplot.tiff', units="in", width=16, height=8, res=300)
ggplot(data=mothur_meta_up, aes(x=sample_ID,y=Read_depth))+geom_bar(stat = "identity") + xlab("Sample_ID")+ylab("Sequencing depth") + scale_y_continuous(labels = scales::comma)+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(family="Arial"), plot.title = element_text(family='Arial', hjust = .5,vjust = 3,size = 24,face="bold"), axis.text.x = element_text(family="Arial",size=8,angle=90),axis.text.y = element_text(family="Arial",size=14),axis.title.x = element_text(family="Arial",size=20,face="bold"),axis.title.y = element_text(family="Arial",size=20,face="bold"))
dev.off()
sort(Read_depth,decreasing = TRUE)

###rarefaction plot
library(reshape2)
library(vegan)
library(phyloseq)
write.csv(t(otu_table(moth_merge_up)),file="moth_merge_up_OTU.csv")
write.csv(otu_table(moth_merge_up),file='moth_merge_up_OTU_t.csv')
moth_merge_up_otu_table<-read.csv("moth_merge_up_OTU.csv",row.names = 1)
dim(moth_merge_up_otu_table) # correct OTUs number
rare<-rarefy(moth_merge_up_otu_table,sample = seq(0,max(sample_sums(moth_merge_up))))
rownames(rare)<-rownames(moth_merge_up_otu_table)
long_rare<-t(rare)
long_rare<-data.frame(subsample=seq(0,max(sample_sums(moth_merge_up)),by=100),long_rare)
head(long_rare)
rare_melt<-melt(long_rare,id.vars = "subsample")# reshape2 package
head(rare_melt)
## 1) Rarefaction curve using maximum read depth

ggplot(rare_melt,aes(x=subsample,y=value,colour=variable))+geom_point()+geom_line(linetype = "dashed")+theme(legend.position="top")+ylab("Bactria OTU numbers")+theme(title=(element_text(size=15,family = "Times",face="bold")),text=element_text(family = "Times",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Times"), axis.title=element_text(size = 15,face="bold",family = "Times"),legend.title = element_text(size=13,face="bold",family = "Times"),legend.text = (element_text(size=10,family = "Times")),legend.position="right")+ggtitle("Rarefaction curve with maximum read depth")

## 2) Rarefaction curve using minimum read depth

rare_min<-rarefy(moth_merge_up_otu_table,sample = seq(0,min(sample_sums(moth_merge_up)),by=100))
rownames(rare_min)<-rownames(moth_merge_up_otu_table)
long_rare_min<-t(rare_min)
long_rare_min<-data.frame(subsample=seq(0,min(sample_sums(moth_merge_up)),by=100),long_rare_min)
head(long_rare_min)
rare_melt_min<-melt(long_rare_min,id.vars = "subsample")
head(rare_melt_min)

#step2 singleton remove()
rms_mothur_phyloseq<-filter_taxa(moth_merge_up,function (x) sum(x > 0) > (0.1*length(x)), TRUE)
rms_mothur_phyloseq#21298 taxa

#step3 rarefied data (1.rarefied--relative abundance calculate--category(phylum, class,...)--merge(replications)--prune(relative abundance>0.02)
rarefied_mothur_phyloseq<-rarefy_even_depth(rms_mothur_phyloseq,sample.size = min(sample_sums(rms_mothur_phyloseq)),rngseed = 1013,replace = TRUE,trimOTUs = TRUE,verbose = TRUE)#normalize the sequences based on the least number of sequences
rarefied_mothur_phyloseq #20008 taxa
write.csv(tax_table(rarefied_mothur_phyloseq),"tax.csv")
write.csv(otu_table(rarefied_mothur_phyloseq),file = "afterrarefysnifa_otutable.csv")

### relative abundance
r_rarefied_mothur_phyloseq<-transform_sample_counts(rarefied_mothur_phyloseq,function(x)x/sum(x))
r_rarefied_mothur_phyloseq #20008 taxa

##phylum(relative abundance)
phylum_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[2],NArm=TRUE)
phylum_r_rarefied_mothur_phyloseq #34 phylum
#write.csv(phylum_r_rarefied_mothur_phyloseq, file = "44_sample")
merged_r_phylum_rarefied_mothur_phyloseq<-merge_samples(phylum_r_rarefied_mothur_phyloseq,"Treat")#merge
merged_r_phylum_rarefied_mothur_phyloseq
r_merged_phylum_r_rarefied_mothur_phyloseq<-transform_sample_counts(merged_r_phylum_rarefied_mothur_phyloseq,function(x) x/sum(x))
r_merged_phylum_r_rarefied_mothur_phyloseq#24
prune_r_merged_phylum_r_rarefied_mothur_phyloseq<-prune_taxa(taxa_sums(r_merged_phylum_r_rarefied_mothur_phyloseq)>=0.02,r_merged_phylum_r_rarefied_mothur_phyloseq)
prune_r_merged_phylum_r_rarefied_mothur_phyloseq#18 

phylum_combine<-cbind((tax_table(prune_r_merged_phylum_r_rarefied_mothur_phyloseq)[,2]),otu_table(t(prune_r_merged_phylum_r_rarefied_mothur_phyloseq)))
write.csv(phylum_combine,file='phylum_combine.csv')
phylum_combine_sub<-read.csv('phylum_combine_sub.csv',sep=",",header=TRUE,row.names = 1)#manully select relative abundance >1%
phylum_melt<-melt(phylum_combine_sub,id.vars = "Phylum")
phylum_melt# totally there have 12 phylum have been selected
head(phylum_melt)
#phylum_combine_sub_transpose<-read.csv('phylum_combinesub_transpose.csv',sep=",",header=TRUE)#manully select relative abundance >1%

#lm_phylum=lm(phylum_combine_sub_transpose$Proteobacteria ~ phylum_combine_sub_transpose$group,data = phylum_combine_sub_transpose)#ra is relative abundance
#(HSD.test(lm_phylum,"phylum_combine_sub_transpose$group")) # no sig
#anova(lm_phylum)
###figure1
cols_phylum <- c("#1B9E77", "#7570B3", "#f2a7cd", "#66A61E", "#E6AB02", "#A6761D","#666666", "#56B4E9","#ABB065", "#999999","#F0E442", "#0072B2" ,"#00AFBB", "#E69F00","#D55E00", "#D95F02", "#CC79A7","#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#009E73","#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#E7B800", "#FC4E07")
tiff('phylum_barplot', units="in", width=16, height=8, res=300)
ggplot(phylum_melt, aes(x =variable, y = value, fill=Phylum))+
  geom_bar(stat='identity',colour="grey",size=0.1)+
  labs(title="",x="Treatment",y="Relative Abundance")+
  scale_fill_manual(values=cols_phylum)+
  theme(title=(element_text(size=15,family = "Times",face="bold")),
        plot.margin = unit(c(0.5,0.1,0.5,0.1),"cm"),
        axis.text.x = element_text(angle=90,hjust = 1),
        text=element_text(family = "Times",face="plain"),
        panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        panel.border = element_rect(colour = NA, fill=NA),
        axis.text=element_text(size=20,family="Times"),
        axis.title=element_text(size = 20,face="bold",family = "Times"),
        legend.title = element_text(size=20,face="bold",family = "Times"),
        legend.text = (element_text(size=20,family = "Times")))
dev.off()

#figure2
Diversity<-data.frame(estimate_richness(rarefied_mothur_phyloseq,split=TRUE,measures = NULL),sample_data(rarefied_mothur_phyloseq))
Diversity
identical(rownames(Diversity),sample_names(rarefied_mothur_phyloseq)) # if true, two objects are equal.
write.csv(Diversity,file = "alpha_diversity.csv")
library(tidyverse)
library(lubridate)
library(ggplot2)
library(agricolae)
library(dplyr)#for mutate function
head(Diversity)

d_diversity <- read.csv(file = "alpha_diversity_rev.csv", sep = ",", row.names = 1)
head(d_diversity)
d = as.data.frame(d_diversity)
df_shannon=d %>% group_by(Treat) %>% mutate(sem_shannon=sd(Shannon)/sqrt(n())) #stdev
mean_shannon = df_shannon %>% group_by(Treat) %>% mutate(avg_shannon=sum(Shannon)/(n()))#combine mean and stdev in "mean_shannon"
write.csv(mean_shannon, "shannon_meananddev.csv")

##shannon_tillage
#mutiple barplot
tiff('tillage_shannon diversity.tiff', units="in", width=7, height=5, res=300)
shannon_tillage <- ggplot(df_shannon, aes(x = d_diversity$Cover_nitro, y = mean_shannon$avg_shannon, fill = Tillage))+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,8)+
  scale_fill_manual(values = c("#E6AB02","aquamarine4"),name="Treatment")+
  labs(title="",x="Treatment",y="Shannon")+
  geom_errorbar(aes(ymin =mean_shannon$avg_shannon - mean_shannon$sem_shannon, ymax = mean_shannon$avg_shannon + mean_shannon$sem_shannon),width = 0.2, colour = "black", position = position_dodge(.9))+
  facet_wrap(~Tillage)+
  theme_bw()
dev.off()

tiff('cover_shannon diversity.tiff', units="in", width=7, height=5, res=300)
shannon_cover <- ggplot(df_shannon, aes(x = d_diversity$Nitro_till, y = mean_shannon$avg_shannon, fill = Cover))+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,8)+
  scale_fill_manual(values = c("#E6AB02","aquamarine4", "#665e8f"), name ="Treatment")+
  labs(title="",x="Treatment",y="Shannon")+
  geom_errorbar(aes(ymin =mean_shannon$avg_shannon - mean_shannon$sem_shannon, ymax = mean_shannon$avg_shannon + mean_shannon$sem_shannon),width = 0.2, colour = "black", position = position_dodge(.9))+
  facet_wrap(~Cover)+
  theme_bw()
dev.off()

tiff('nitro_shannon diversity.tiff', units="in", width=7, height=5, res=300)
shannon_nitro <- ggplot(df_shannon, aes(x = d_diversity$Cover_till, y = mean_shannon$avg_shannon, fill = Nitrogen))+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,8)+
  scale_fill_manual(values = c("#E6AB02","#665e8f"), name ="Treatment")+
  labs(title="",x="Treatment",y="Shannon")+
  geom_errorbar(aes(ymin =mean_shannon$avg_shannon - mean_shannon$sem_shannon, ymax = mean_shannon$avg_shannon + mean_shannon$sem_shannon),width = 0.2, colour = "black", position = position_dodge(.9))+
  facet_wrap(~Nitrogen)+
  theme_bw()
dev.off()
#arrange multiple figures in one figure
library("cowplot")
tiff('arrange_shannon diversity.tiff', units="in", width=13, height=5, res=300)
plot_grid(shannon_nitro,shannon_tillage,shannon_cover,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
dev.off()

#multiple boxplot

#figure2
#chao1
tiff('boxTillage_Chao1 diversity.tiff', units="in", width=7, height=5, res=300)
box_Chao1_tillage <- ggplot(d_diversity,aes(x=d_diversity$Cover_nitro,y=Chao1,fill=Tillage))+
  geom_boxplot(size=0.5)+
  geom_jitter(shape=16, position=position_jitter(0),aes(fill=Tillage))+
  scale_fill_manual(values = c("#E6AB02","aquamarine4"),name="Treatment")+
  labs(title="",x="Treatment",y="Chao1 Index")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),
        axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="bold",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+
  facet_wrap(~Tillage)+
  theme_bw()
dev.off()

tiff('boxcover_Chao1 diversity.tiff', units="in", width=7, height=5, res=300)
box_Chao1_cover <- ggplot(d_diversity,aes(x=d_diversity$Nitro_till,y=Chao1,fill=Cover))+
  geom_boxplot(size=0.5)+
  geom_jitter(shape=16, position=position_jitter(0),aes(fill=Cover))+
  scale_fill_manual(values = c("#E6AB02","aquamarine4","#665e8f"),name="Treatment")+
  labs(title="",x="Treatment",y="Chao1 Index")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),
        axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="bold",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+
  facet_wrap(~Cover)+
  theme_bw()
dev.off()

tiff('boxnitrogen_Chao1 diversity.tiff', units="in", width=7, height=5, res=300)
box_Chao1_nitro <- ggplot(d_diversity,aes(x=d_diversity$Cover_till,y=Chao1,fill=Nitrogen))+
  geom_boxplot(size=0.5)+
  geom_jitter(shape=16, position=position_jitter(0),aes(fill=Nitrogen))+
  scale_fill_manual(values = c("#E6AB02","#665e8f"),name="Treatment")+
  labs(title="",x="Treatment",y="Chao1 Index")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),
        axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="bold",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+
  facet_wrap(~Nitrogen)+
  theme_bw()
dev.off()

treatlm_chao1=lm(d_diversity$Chao1~d_diversity$Treat,data =d_diversity)
(HSD.test(treatlm_chao1,"d_diversity$Treat")) # no sig

tillagelm_chao1=lm(d_diversity$Chao1~d_diversity$Cover_nitro,data =d_diversity)
(HSD.test(tillagelm_chao1,"d_diversity$Cover_nitro")) # no sig

coverlm_chao1=lm(d_diversity$Chao1~d_diversity$Nitro_till,data =d_diversity)
(HSD.test(coverlm_chao1,"d_diversity$Nitro_till")) # no sig

nitrolm_chao1=lm(d_diversity$Chao1~d_diversity$Cover_till,data =d_diversity)
(HSD.test(nitrolm_chao1,"d_diversity$Cover_till")) # no sig

#numbers of OTUs in each sample 
write.csv(otu_table(rarefied_mothur_phyloseq), file = "rarefied_otutable.csv" )
write.csv(sample_data(rarefied_mothur_phyloseq), file = "rarefied_sampledata.csv" )
otu_sampledata <- read.csv("rarefied_otutable_sampledata.csv", header = TRUE, sep = ",",row.names = 1)

tiff('box_otus.diversity.tiff', units="in", width=10, height=5, res=300)
otu_num <- ggplot(otu_sampledata,aes(x=otu_sampledata$Treat,y=otu_sampledata$OTUs,fill=Treat))+
  geom_boxplot(size=0.5)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  scale_fill_manual(values = c("#00798c","#d1495b","#edae49","#66a182","#2e4057","#8d96a3","#a2d729","#E6AB02","aquamarine4","#665e8f", "lightsalmon1","lightskyblue"),name="Treatment")+
  labs(title="",x="Treatment",y="Phylotype richness (OTUs)")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=5,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y = element_text(size=5,family="Arial"),
        axis.text.x = element_text(family="Arial",size=5,angle=0,hjust = 1),
        axis.title=element_text(size = 5,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=5,face="plain",family = "Arial"),
        legend.text = (element_text(size=5,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+theme_bw() 
dev.off()

otu_till <- ggplot(otu_sampledata,aes(x=Tillage,y=otu_sampledata$OTUs,color=Tillage))+
  geom_boxplot(size=0.5,width=0.5,alpha=0.5,notch = FALSE,na.rm=TRUE,outlier.shape=NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  scale_color_manual(values = c("#d1495b","#edae49","#66a182","#8d96a3","#a2d729","#E6AB02","aquamarine4","#665e8f", "lightsalmon1","lightskyblue"),name="Treatment")+
  labs(title="",x="Treatment",y="Phylotype richness (OTUs)")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=5,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y = element_text(size=5,family="Arial"),
        axis.text.x = element_text(family="Arial",size=5,angle=0,hjust = 1),
        axis.title=element_text(size = 5,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=5,face="plain",family = "Arial"),
        legend.text = (element_text(size=5,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+theme_bw() 
otu_nitrogen <- ggplot(otu_sampledata,aes(x=Nitrogen,y=otu_sampledata$OTUs,color=Nitrogen))+
  geom_boxplot(size=0.5,width=0.5,alpha=0.5,notch = FALSE,na.rm=TRUE,outlier.shape=NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  scale_color_manual(values = c("#d1495b","#2e4057","#edae49","#66a182","#8d96a3","#a2d729","#E6AB02","aquamarine4","#665e8f", "lightsalmon1","lightskyblue"),name="Treatment")+
  labs(title="",x="Treatment",y="Phylotype richness (OTUs)")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=5,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y = element_text(size=5,family="Arial"),
        axis.text.x = element_text(family="Arial",size=5,angle=0,hjust = 1),
        axis.title=element_text(size = 5,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=5,face="plain",family = "Arial"),
        legend.text = (element_text(size=5,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+theme_bw()
otu_cover <- ggplot(otu_sampledata,aes(x=Cover,y=otu_sampledata$OTUs,color=Cover))+
  geom_boxplot(size=0.5,width=0.5,alpha=0.5,notch = FALSE,na.rm=TRUE,outlier.shape=NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  scale_color_manual(values = c("#d1495b","#edae49","#2e4057","#d1495b","#2e4057","#a2d729","#E6AB02","aquamarine4","#665e8f", "lightsalmon1","lightskyblue"),name="Treatment")+
  labs(title="",x="Treatment",y="Phylotype richness (OTUs)")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=5,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y = element_text(size=5,family="Arial"),
        axis.text.x = element_text(family="Arial",size=5,angle=0,hjust = 1),
        axis.title=element_text(size = 5,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=5,face="plain",family = "Arial"),
        legend.text = (element_text(size=5,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+theme_bw()

otu.mod <- lmer(OTUs~Nitrogen*Tillage*Cover+(1|block),data =otu_sampledata )
anova(otu.mod)
Anova(otu.mod)
treatlm_otu=lm(otu_sampledata$OTUs~otu_sampledata$Treat,data =otu_sampledata)
(HSD.test(treatlm_otu,"otu_sampledata$Treat")) # no sig

#legend won't show up if theme_bw() not include
#arrange multiple figures in one figure
library("cowplot")
tiff('arrange_Chao1_otu diversity.tiff', units="in", width=17, height=7, res=300)
plot_grid(box_Chao1_nitro ,box_Chao1_tillage,box_Chao1_cover,otu_num,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()

#figure 3
tiff('boxTillage_shannon diversity.tiff', units="in", width=7, height=5, res=300)
box_shannon_tillage <- ggplot(d_diversity,aes(x=Tillage,y=Shannon,color=Tillage))+
  geom_boxplot(size=0.5,width=0.5,alpha=0.5,notch = FALSE,na.rm=TRUE,outlier.shape=NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  scale_color_manual(values = c("#E6AB02","aquamarine4"),name="Treatment")+
  labs(title="",x="Treatment",y="Shannon Index (H')")+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),
        axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="bold",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+
#  facet_wrap(~Tillage)+
  theme_bw()
dev.off()

tiff('boxcover_shannon diversity.tiff', units="in", width=7, height=5, res=300)
box_shannon_cover <- ggplot(d_diversity,aes(x=Cover,y=Shannon,color=Cover))+
  geom_boxplot(size=0.5,width=0.5,alpha=0.5,notch = FALSE,na.rm=TRUE,outlier.shape=NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  scale_color_manual(values = c("#E6AB02","aquamarine4","#665e8f"),name="Treatment")+
  labs(title="",x="Treatment",y="Shannon Index (H')")+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),
        axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="bold",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+
#  facet_wrap(~Cover)+
  theme_bw()
dev.off()

tiff('boxnitrogen_shannon diversity.tiff', units="in", width=7, height=5, res=300)
box_shannon_nitro <- ggplot(d_diversity,aes(x=Nitrogen,y=Shannon,color=Nitrogen))+
  geom_boxplot(size=0.5,width=0.5,alpha=0.5,notch = FALSE,na.rm=TRUE,outlier.shape=NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  scale_color_manual(values = c("#E6AB02","#665e8f"),name="Treatment")+
  labs(title="",x="Treatment",y="Shannon Index (H')")+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),
        axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="bold",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+
#  facet_wrap(~Nitrogen)+
  theme_bw()
dev.off()

head(d_diversity)
shan.mod <- lmer(Shannon~Nitrogen*Tillage*Cover+(1|block),data = d_diversity)
Anova(shan.mod)
anova(shan.mod)

library("ggpubr")
tiff('figures/alpha_diversity.tiff', units="in", width=11, height=4, res=300)
ggarrange(box_shannon_nitro,box_shannon_tillage,box_shannon_cover,
          labels = c("A","B","C"),ncol = 3, nrow = 1,common.legend = F,legend = "bottom")
dev.off()

treatlm_Shannon=lm(d_diversity$Shannon~d_diversity$Treat,data =d_diversity)
(HSD.test(treatlm_Shannon,"d_diversity$Treat")) # no sig

tillagelm_Shannon=lm(d_diversity$Shannon~d_diversity$Cover_nitro,data =d_diversity)
(HSD.test(tillagelm_Shannon,"d_diversity$Cover_nitro")) # no sig

coverlm_Shannon=lm(d_diversity$Shannon~d_diversity$Nitro_till,data =d_diversity)
(HSD.test(coverlm_Shannon,"d_diversity$Nitro_till")) # no sig

nitrolm_Shannon=lm(d_diversity$Shannon~d_diversity$Cover_till,data =d_diversity)
(HSD.test(nitrolm_Shannon,"d_diversity$Cover_till")) # no sig


#Pielou's measure, ref:https://www.rdocumentation.org/packages/pez/versions/1.0-0/topics/evenness
library(vegan)
H <- diversity(t(otu_table(rarefied_mothur_phyloseq)))
Pielou <- H/log(specnumber(t(otu_table(rarefied_mothur_phyloseq))))
head(Pielou)
write.csv(Pielou,file = "Pielou.csv")
d_diversity_Pielou <- cbind(d_diversity,Pielou)
head(d_diversity_Pielou)
treatlm_Pielou=lm(d_diversity_Pielou$Pielou~d_diversity_Pielou$Treat,data = d_diversity_Pielou)
(HSD.test(treatlm_Pielou,"d_diversity_Pielou$Treat")) # no sig

tiff('boxTillage_Pielou diversity.tiff', units="in", width=7, height=5, res=300)
box_Pielou_tillage <- ggplot(d_diversity_Pielou,aes(x=d_diversity_Pielou$Cover_nitro,y=Pielou,fill=Tillage))+
  geom_boxplot(size=0.5)+
  geom_jitter(shape=16, position=position_jitter(0),aes(fill=Tillage))+
  scale_fill_manual(values = c("#00798c","#d1495b"),name="Treatment")+
  labs(title="",x="Treatment",y="Pielou's Index (J')")+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),
        axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="bold",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+
  facet_wrap(~Tillage)+
  theme_bw()
dev.off()

tiff('boxcover_Pielou diversity.tiff', units="in", width=7, height=5, res=300)
box_Pielou_cover <- ggplot(d_diversity_Pielou,aes(x=d_diversity_Pielou$Nitro_till,y=Pielou,fill=Cover))+
  geom_boxplot(size=0.5)+
  geom_jitter(shape=16, position=position_jitter(0),aes(fill=Cover))+
  scale_fill_manual(values = c("#00798c","#d1495b","#edae49"),name="Treatment")+
  labs(title="",x="Treatment",y="Pielou's Index (J')")+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),
        axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="bold",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+
  facet_wrap(~Cover)+
  theme_bw()
dev.off()

tiff('boxnitrogen_Pielou diversity.tiff', units="in", width=7, height=5, res=300)
box_Pielou_nitro <- ggplot(d_diversity_Pielou,aes(x=d_diversity_Pielou$Cover_till,y=Pielou,fill=Nitrogen))+
  geom_boxplot(size=0.5)+
  geom_jitter(shape=16, position=position_jitter(0),aes(fill=Nitrogen))+
  scale_fill_manual(values = c("#d1495b","#edae49"),name="Treatment")+
  labs(title="",x="Treatment",y="Pielou's Index (J')")+
  theme(legend.position="top",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),
        axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="bold",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")),
        plot.margin = unit(c(1,1,1,1), "cm"))+
  facet_wrap(~Nitrogen)+
  theme_bw()
dev.off()

treatlm_Pielou=lm(d_diversity_Pielou$Pielou~d_diversity_Pielou$Treat,data =d_diversity_Pielou)
(HSD.test(treatlm_Pielou,"d_diversity_Pielou$Treat")) # no sig

tillagelm_Pielou=lm(d_diversity_Pielou$Pielou~d_diversity_Pielou$Cover_nitro,data =d_diversity_Pielou)
(HSD.test(tillagelm_Pielou,"d_diversity_Pielou$Cover_nitro")) # no sig


coverlm_Pielou=lm(d_diversity_Pielou$Pielou~d_diversity_Pielou$Nitro_till,data =d_diversity_Pielou)
(HSD.test(coverlm_Pielou,"d_diversity_Pielou$Nitro_till")) # no sig

nitrolm_Pielou=lm(d_diversity_Pielou$Pielou~d_diversity_Pielou$Cover_till,data =d_diversity_Pielou)
(HSD.test(nitrolm_Pielou,"d_diversity_Pielou$Cover_till")) # no sig

tillagelm_Pielou=lm(d_diversity_Pielou$Pielou~d_diversity_Pielou$Tillage,data =d_diversity_Pielou)
(HSD.test(tillagelm_Pielou,"d_diversity_Pielou$Tillage")) # no sig

coverlm_Pielou=lm(d_diversity_Pielou$Pielou~d_diversity_Pielou$Cover,data =d_diversity_Pielou)
(HSD.test(coverlm_Pielou,"d_diversity_Pielou$Cover")) # no sig

nitrolm_Pielou=lm(d_diversity_Pielou$Pielou~d_diversity_Pielou$Nitrogen,data =d_diversity_Pielou)
(HSD.test(nitrolm_Pielou,"d_diversity_Pielou$Nitrogen")) # no sig

#arrange multiple figures in one figure
library("cowplot")
tiff('arrange_shannon_Pielou.tiff', units="in", width=21, height=7, res=300)
#plot_align_shannon <- align_plots(box_shannon_nitro, box_shannon_tillage, box_shannon_cover, align = 'v', axis = 'l')
#plot_align_Pielou <- align_plots(box_Pielou_nitro, box_Pielou_tillage, box_Pielou_cover, align = 'v', axis = 'l')
#plot_grid(plot_align_shannon, plot_align_Pielou,
#         labels = c("A", "B", "C", "D", "E", "F"),
#         ncol = 2, nrow = 3, align = "hv",rel_widths = c(1, 1))
plot_grid(box_shannon_nitro, box_shannon_tillage, box_shannon_cover, box_Pielou_nitro,  box_Pielou_tillage, box_Pielou_cover,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2, align = "hv",rel_widths = c(1, 1),axis = 'l')
dev.off()
write.csv(d_diversity_Pielou, file = "d_diversity_pielou.csv")

##figure 4
#CCA
library(vegan)
library(ggplot2)
r_rarefied_mothur_phyloseq_bray<-vegdist(t(otu_table(r_rarefied_mothur_phyloseq)),method="bray",binary = FALSE)
r_rarefied_mothur_phyloseq_CCA<-ordinate(r_rarefied_mothur_phyloseq,method = "CCA", r_rarefied_mothur_phyloseq_bray, ~Nitrogen+Tillage+Cover)
# Calculate bray curtis distance matrix
bray_not_na <- phyloseq::distance(rarefied_mothur_phyloseq, method = "bray")
cap_ord <- ordinate(physeq = rarefied_mothur_phyloseq, method = "CAP",distance = bray_not_na, formula = ~VBR+Temp+Virus_abund+Bacteria_abund)
# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, yend = CAP2, x = 0,  y = 0,  shape = NULL,color = NULL, label = labels)
label_map <- aes(x = 1.3 * CAP1, y = 1.3 * CAP2, shape = NULL, color = NULL,label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))
#block effect
tiff('CCA blockeffect.tiff', units="in", width=7, height=5, res=300)
dev.off()
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_CCA,type="samples",color="block")+geom_point(size=4)+ggtitle("Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c","gold","aquamarine"))+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=13,family = "Arial")))
#nitrogen&cover--good
r_rarefied_mothur_phyloseq_CCA_nitrogen<-ordinate(r_rarefied_mothur_phyloseq,method = "CCA", r_rarefied_mothur_phyloseq_bray, ~Nitrogen)
library(ggforce)
library(italic)
tiff('CCA nitrogenandcover.tiff', units="in", width=7, height=5, res=300)
CCA_nitro <- plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_CCA,type="samples",color="Nitrogen", shape = "Cover")+
  geom_point(size=4)+ggtitle("")+scale_color_manual(values = c("#2c5491","#d1495b","#edae49"))+
  labs(x="CCA1",y="CCA2")+
  geom_mark_ellipse(aes(color = Nitrogen, group=Nitrogen),expand=0)+
  geom_vline(xintercept = c(0), size=0.2,linetype=2)+
  geom_hline(yintercept = c(0), size=0.2,linetype=2)+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=10,family="Arial"),
        axis.title=element_text(size = 10,face="plain",family = "Arial"),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  annotate(geom="text", x=1.6, y=2.5, size=8,label="p",color="black", fontface="bold.italic")+
  annotate(geom="text", x=2.5, y=2.5, size=8,label="< 0.0001",color="black", fontface="bold")


dev.off()

#cover&nitrogen--not good
dev.off()
tiff('CCA coverandnitrogen.tiff', units="in", width=7, height=5, res=300)
CCA_cover <- plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_CCA,type="samples",color="Cover", shape = "Nitrogen")+
  geom_point(size=4)+ggtitle("")+
  scale_color_manual(values = c("#2c5491","#d1495b","#edae49"))+
  labs(x="CCA1",y="CCA2")+
  geom_mark_ellipse(aes(color = Cover, group=Cover),expand=0)+
  geom_vline(xintercept = c(0), size=0.2,linetype=2)+
  geom_hline(yintercept = c(0), size=0.2,linetype=2)+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=10,family="Arial"),
        axis.title=element_text(size = 10,face="plain",family = "Arial"),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  annotate(geom="text", x=1.8, y=2.5, size=8,label="p",color="black", fontface="bold.italic")+
  annotate(geom="text", x=2.5, y=2.5, size=8,label="< 0.05",color="black", fontface="bold")
dev.off()

#tillage&nitrogen--good
dev.off()
tiff('CCA tillageandnitrogen.tiff', units="in", width=7, height=5, res=300)
CCA_tillage <- plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_CCA,type="samples",color="Tillage", shape = "Nitrogen")+
  geom_point(size=4)+ggtitle("")+
  scale_color_manual(values = c("#2c5491","#d1495b","#edae49"))+
  labs(x="CCA1",y="CCA2")+
  geom_mark_ellipse(aes(color = Tillage, group=Tillage),expand=0)+
  geom_vline(xintercept = c(0), size=0.2,linetype=2)+
  geom_hline(yintercept = c(0), size=0.2,linetype=2)+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=10,family="Arial"),
        axis.title=element_text(size = 10,face="plain",family = "Arial"),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
    annotate(geom="text", x=1.6, y=2.5, size=8,label="p",color="black", fontface="bold.italic")+
    annotate(geom="text", x=2.5, y=2.5, size=8,label="< 0.0001",color="black", fontface="bold")
dev.off()

tiff('arrownitro_cca.tiff', units="in", width=7, height=5, res=300)
arrownitro_cca <- CCA_nitro + geom_segment(mapping = arrow_map, size = .5,data = arrowdf,color = "gray",arrow = arrowhead )+
  geom_text(mapping = label_map, size = 4,  data = arrowdf,  show.legend = FALSE)
dev.off()

tiff('arrowtillage_cca.tiff', units="in", width=7, height=5, res=300)
arrowtillage_cca <- CCA_tillage + geom_segment(mapping = arrow_map, size = .5,data = arrowdf,color = "gray",arrow = arrowhead )+
  geom_text(mapping = label_map, size = 4,  data = arrowdf,  show.legend = FALSE)
dev.off()

tiff('arrowcover_cca.tiff', units="in", width=7, height=5, res=300)
arrowcover_cca <- CCA_cover+ geom_segment(mapping = arrow_map, size = .5,data = arrowdf,color = "gray",arrow = arrowhead )+
  geom_text(mapping = label_map, size = 4,  data = arrowdf,  show.legend = FALSE)
dev.off()


#arrange multiple figures in one figure
library("cowplot")
tiff('arrwoarrange_cca.tiff', units="in", width=16, height=5, res=300)
plot_grid(arrownitro_cca,arrowtillage_cca,arrowcover_cca,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1, align = "hv",rel_widths = c(1, 1),axis = 'l')
dev.off()

tiff('newcolor_cca.tiff', units="in", width=16, height=7, res=300)
plot_grid(CCA_nitro,CCA_tillage,CCA_cover,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1, align = "hv",rel_widths = c(1, 1),axis = 'l')
dev.off()

r_rarefied_mothur_phyloseq_PCoA_1<-ordinate(r_rarefied_mothur_phyloseq,method = "PCoA", r_rarefied_mothur_phyloseq_bray, ~Nitrogen)

pcoa_nitro <- plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="biplot",color="Nitrogen", shape = "Cover")+
  geom_point(size=4)+ggtitle("")+scale_color_manual(values = c("#2c5491","#d1495b","#edae49"))+
  #labs(x="CCA1",y="CCA2")+
  geom_mark_ellipse(aes(color = Nitrogen, group=Nitrogen),expand=0)+
  geom_vline(xintercept = c(0), size=0.2,linetype=2)+
  geom_hline(yintercept = c(0), size=0.2,linetype=2)+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=10,family="Arial"),
        axis.title=element_text(size = 10,face="plain",family = "Arial"),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))
  # annotate(geom="text", x=1.6, y=2.5, size=8,label="p",color="black", fontface="bold.italic")+
  # annotate(geom="text", x=2.5, y=2.5, size=8,label="< 0.0001",color="black", fontface="bold")


dev.off()

#cover&nitrogen--not good
dev.off()
tiff('CCA coverandnitrogen.tiff', units="in", width=7, height=5, res=300)
pcoa_cover <- plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA_1,type="samples",color="Cover", shape = "Nitrogen")+
  geom_point(size=4)+ggtitle("")+
  scale_color_manual(values = c("#2c5491","#d1495b","#edae49"))+
  labs(x="CCA1",y="CCA2")+
  geom_mark_ellipse(aes(color = Cover, group=Cover),expand=0)+
  geom_vline(xintercept = c(0), size=0.2,linetype=2)+
  geom_hline(yintercept = c(0), size=0.2,linetype=2)+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=10,family="Arial"),
        axis.title=element_text(size = 10,face="plain",family = "Arial"),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  annotate(geom="text", x=1.8, y=2.5, size=8,label="p",color="black", fontface="bold.italic")+
  annotate(geom="text", x=2.5, y=2.5, size=8,label="< 0.05",color="black", fontface="bold")
dev.off()

#tillage&nitrogen--good
dev.off()
tiff('CCA tillageandnitrogen.tiff', units="in", width=7, height=5, res=300)
CCA_tillage <- plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA_1,type="samples",color="Tillage", shape = "Nitrogen")+
  geom_point(size=4)+ggtitle("")+
  scale_color_manual(values = c("#2c5491","#d1495b","#edae49"))+
  labs(x="CCA1",y="CCA2")+
  geom_mark_ellipse(aes(color = Tillage, group=Tillage),expand=0)+
  geom_vline(xintercept = c(0), size=0.2,linetype=2)+
  geom_hline(yintercept = c(0), size=0.2,linetype=2)+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=10,family="Arial"),
        axis.title=element_text(size = 10,face="plain",family = "Arial"),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  annotate(geom="text", x=1.6, y=2.5, size=8,label="p",color="black", fontface="bold.italic")+
  annotate(geom="text", x=2.5, y=2.5, size=8,label="< 0.0001",color="black", fontface="bold")
dev.off()


#table2
phylum_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[2],NArm=TRUE)
phylum_r_rarefied_mothur_phyloseq

class_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[3],NArm=TRUE)
class_r_rarefied_mothur_phyloseq
write.csv(tax_table(class_r_rarefied_mothur_phyloseq), file = "classtaxa.csv ")
order_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[4],NArm=TRUE)
order_r_rarefied_mothur_phyloseq
write.csv(tax_table(order_r_rarefied_mothur_phyloseq), file = "ordertaxa.csv ")

###Simper analysis--phylum level # the average/the sum of average=contribution of each species to the total dissmilarity between groups
library(vegan)
sampledata <- sample_data(phylum_r_rarefied_mothur_phyloseq)
sr_phylum_r_rarefied_mothur_phyloseq <- t(sqrt(otu_table(phylum_r_rarefied_mothur_phyloseq)))#square root normalization
simper_output_tillage <- simper(sr_phylum_r_rarefied_mothur_phyloseq,sampledata$Tillage,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
tillage_phylum <- summary(simper_output_tillage,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_tillage
write.csv(as.data.frame(tillage_phylum$Conventional_No_tillage),file = "1stsimper_output_tillage.csv")


simper_output_nitrogen <- simper(sr_phylum_r_rarefied_mothur_phyloseq,sampledata$Nitrogen,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
nitrogen_simper <- summary(simper_output_nitrogen,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_nitrogen
write.csv(as.data.frame(nitrogen_simper$N0_N60),file = "1stsimper_output_nitrogen.csv")


simper_output_cover <- simper(sr_phylum_r_rarefied_mothur_phyloseq,sampledata$Cover,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
cover_simper <- summary(simper_output_cover,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_cover
write.csv(as.data.frame(cover_simper$No_cover_Vetch),file = "1stsimper_output_cover_No_cover_Vetch.csv")
write.csv(as.data.frame(cover_simper$No_cover_Wheat),file = "1stsimper_output_cover_No_cover_Wheat.csv")
write.csv(as.data.frame(cover_simper$Vetch_Wheat),file = "1stsimper_output_cover_Vetch_Wheat.csv")

#average:Average contribution to overall dissimilarity.overall:The overall between-group dissimilarity.sd:Standard deviation of contribution.ratio:Average to sd ratio.ava, avb/ava:Average abundances per group.ord:An index vector to order vectors by their contribution or order cusum back to the original data order.cusum:Ordered cumulative contribution.p:Permutation p-value. Probability of getting a larger or equal average contribution in random permutation of the group factor.

##Simper analysis--class level
sampledata_class <- sample_data(class_r_rarefied_mothur_phyloseq)
sr_class_r_rarefied_mothur_phyloseq <- t(sqrt(otu_table(class_r_rarefied_mothur_phyloseq)))#square root normalization
simper_output_tillage_class <- simper(sr_class_r_rarefied_mothur_phyloseq,sampledata_class$Tillage,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
tillage_simper_class <- summary(simper_output_tillage_class,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_tillage_class
write.csv(as.data.frame(tillage_simper_class$Conventional_No_tillage),file = "1stsimper_output_tillage_class.csv")
summary(simper_output_tillage_class)
simper_output_nitrogen_class <- simper(sr_class_r_rarefied_mothur_phyloseq,sampledata_class$Nitrogen,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
nitrogen_simper_class <-summary(simper_output_nitrogen_class,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_nitrogen_class
write.csv(as.data.frame(nitrogen_simper_class$N0_N60),file = "1stsimper_output_nitrogen_class.csv")

simper_output_cover_class <- simper(sr_class_r_rarefied_mothur_phyloseq,sampledata_class$Cover,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
cover_simper_class <- summary(simper_output_cover_class,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_cover_class
write.csv(as.data.frame(cover_simper_class$No_cover_Vetch),file = "1stsimper_output_cover_class_No_cover_Vetch.csv")
write.csv(as.data.frame(cover_simper_class$No_cover_Wheat),file = "1stsimper_output_cover_class_No_cover_Wheat.csv")
write.csv(as.data.frame(cover_simper_class$Vetch_Wheat),file = "1stsimper_output_cover_class_Vetch_Wheat.csv")

##Simper analysis--order level
sampledata_order <- sample_data(order_r_rarefied_mothur_phyloseq)
sr_order_r_rarefied_mothur_phyloseq <- t(sqrt(otu_table(order_r_rarefied_mothur_phyloseq)))#square root normalization
simper_output_tillage_order <- simper(sr_order_r_rarefied_mothur_phyloseq,sampledata_order$Tillage,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
tillage_simper_order <-summary(simper_output_tillage_order,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_tillage_order
write.csv(as.data.frame(tillage_simper_order$Conventional_No_tillage),file = "1stsimper_output_tillage_order.csv")

simper_output_nitrogen_order <- simper(sr_order_r_rarefied_mothur_phyloseq,sampledata_order$Nitrogen,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
nitrogen_simper_order <- summary(simper_output_nitrogen_order,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_nitrogen_order
write.csv(as.data.frame(nitrogen_simper_order$N0_N60),file = "1stsimper_output_nitrogen_order.csv")

simper_output_cover_order <- simper(sr_order_r_rarefied_mothur_phyloseq,sampledata_order$Cover,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
cover_simper_order <- summary(simper_output_cover_order,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_cover_order
write.csv(as.data.frame(cover_simper_order$No_cover_Vetch),file = "1stsimper_output_cover_order_No_cover_Vetch.csv")
write.csv(as.data.frame(cover_simper_order$No_cover_Wheat),file = "1stsimper_output_cover_order_No_cover_Wheat.csv")
write.csv(as.data.frame(cover_simper_order$Vetch_Wheat),file = "1stsimper_output_cover_order_Vetch_Wheat.csv")


#phylum relative abudance stats
phylum_r_stats<-data.frame(t(otu_table(phylum_r_rarefied_mothur_phyloseq)),sample_data(phylum_r_rarefied_mothur_phyloseq))
phylum_r_stats
identical(rownames(phylum_r_stats),sample_names(phylum_r_rarefied_mothur_phyloseq)) # if true, two objects are equal.
write.csv(phylum_r_stats,file = "phylum_r_stats.csv")
write.csv(tax_table(phylum_r_rarefied_mothur_phyloseq),file = "phylum_taxa_all.csv")


#class relative abudance stats
class_r_stats<-data.frame(t(otu_table(class_r_rarefied_mothur_phyloseq)),sample_data(class_r_rarefied_mothur_phyloseq))
class_r_stats
identical(rownames(class_r_stats),sample_names(class_r_rarefied_mothur_phyloseq)) # if true, two objects are equal.
write.csv(class_r_stats,file = "class_r_stats.csv")
write.csv(tax_table(class_r_rarefied_mothur_phyloseq),file = "class_taxa_all.csv")

phylum_r_stats_new <- read.csv(file = "phylum_r_stats_new.csv", header = TRUE, sep = ",", row.names = 1 )
class_r_stats_new <- read.csv(file = "class_r_stats_new.csv", header = TRUE, sep = ",", row.names = 1 )
head(phylum_r_stats_new)
head(class_r_stats_new)



##plylum_barplot
library(dplyr)
df_Proteobacteria=phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Proteobacteria=sd(Proteobacteria)/sqrt(n())) #stdev
mean_Proteobacteria = df_Proteobacteria %>% group_by(Treat) %>% mutate(avg_Proteobacteria=sum(Proteobacteria)/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Proteobacteria, file = "df_Proteobacteria.csv")

tiff('Proteobacteria.tiff', units="in", width=7, height=5, res=300)
Proteobacteria <- ggplot(df_Proteobacteria, aes(x = phylum_r_stats_new$Treat, y = mean_Proteobacteria$avg_Proteobacteria, fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Proteobacteria",x="",y="Relative abundance")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Proteobacteria$avg_Proteobacteria - mean_Proteobacteria$sem_Proteobacteria, ymax = mean_Proteobacteria$avg_Proteobacteria + mean_Proteobacteria$sem_Proteobacteria),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()

df_Acidobacteria=phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Acidobacteria=sd(Acidobacteria)/sqrt(n())) #stdev
mean_Acidobacteria = df_Acidobacteria %>% group_by(Treat) %>% mutate(avg_Acidobacteria=sum(Acidobacteria)/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Acidobacteria, file = "df_Acidobacteria.csv")

tiff('Acidobacteria.tiff', units="in", width=7, height=5, res=300)
Acidobacteria <- ggplot(df_Acidobacteria, aes(x = phylum_r_stats_new$Treat, y = mean_Acidobacteria$avg_Acidobacteria, fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Acidobacteria",x="",y="")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Acidobacteria$avg_Acidobacteria - mean_Acidobacteria$sem_Acidobacteria, ymax = mean_Acidobacteria$avg_Acidobacteria + mean_Acidobacteria$sem_Acidobacteria),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()

df_Actinobacteria=phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Actinobacteria=sd(Actinobacteria)/sqrt(n())) #stdev
mean_Actinobacteria = df_Actinobacteria %>% group_by(Treat) %>% mutate(avg_Actinobacteria=sum(Actinobacteria)/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Actinobacteria, file = "df_Actinobacteria.csv")

tiff('Actinobacteria.tiff', units="in", width=7, height=5, res=300)
Actinobacteria <- ggplot(df_Actinobacteria, aes(x = phylum_r_stats_new$Treat, y = mean_Actinobacteria$avg_Actinobacteria, fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Actinobacteria",x="",y="")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Actinobacteria$avg_Actinobacteria - mean_Actinobacteria$sem_Actinobacteria, ymax = mean_Actinobacteria$avg_Actinobacteria + mean_Actinobacteria$sem_Actinobacteria),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()


df_Planctomycetes=phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Planctomycetes=sd(Planctomycetes)/sqrt(n())) #stdev
mean_Planctomycetes = df_Planctomycetes %>% group_by(Treat) %>% mutate(avg_Planctomycetes=sum(Planctomycetes)/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Planctomycetes, file = "df_Planctomycetes.csv")

tiff('Planctomycetes.tiff', units="in", width=7, height=5, res=300)
Planctomycetes <- ggplot(df_Planctomycetes, aes(x = phylum_r_stats_new$Treat, y = mean_Planctomycetes$avg_Planctomycetes, fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Planctomycetes",x="",y="Relative abundance")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Planctomycetes$avg_Planctomycetes - mean_Planctomycetes$sem_Planctomycetes, ymax = mean_Planctomycetes$avg_Planctomycetes + mean_Planctomycetes$sem_Planctomycetes),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()

df_Bacteroidetes=phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Bacteroidetes=sd(Bacteroidetes)/sqrt(n())) #stdev
mean_Bacteroidetes = df_Bacteroidetes %>% group_by(Treat) %>% mutate(avg_Bacteroidetes=sum(Bacteroidetes)/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Bacteroidetes, file = "df_Bacteroidetes.csv")

tiff('Bacteroidetes.tiff', units="in", width=7, height=5, res=300)
Bacteroidetes <- ggplot(df_Bacteroidetes, aes(x = phylum_r_stats_new$Treat, y = mean_Bacteroidetes$avg_Bacteroidetes, fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Bacteroidetes",x="",y="")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Bacteroidetes$avg_Bacteroidetes - mean_Bacteroidetes$sem_Bacteroidetes, ymax = mean_Bacteroidetes$avg_Bacteroidetes + mean_Bacteroidetes$sem_Bacteroidetes),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()

df_Chloroflexi=phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Chloroflexi=sd(Chloroflexi)/sqrt(n())) #stdev
mean_Chloroflexi = df_Chloroflexi %>% group_by(Treat) %>% mutate(avg_Chloroflexi=sum(Chloroflexi)/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Chloroflexi, file = "df_Chloroflexi.csv")

tiff('Chloroflexi.tiff', units="in", width=7, height=5, res=300)
Chloroflexi <- ggplot(df_Chloroflexi, aes(x = phylum_r_stats_new$Treat, y = mean_Chloroflexi$avg_Chloroflexi, fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Chloroflexi",x="",y="")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Chloroflexi$avg_Chloroflexi - mean_Chloroflexi$sem_Chloroflexi, ymax = mean_Chloroflexi$avg_Chloroflexi + mean_Chloroflexi$sem_Chloroflexi),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()

df_Firmicutes=phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Firmicutes=sd(Firmicutes)/sqrt(n())) #stdev
mean_Firmicutes = df_Firmicutes %>% group_by(Treat) %>% mutate(avg_Firmicutes=sum(Firmicutes)/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Firmicutes, file = "df_Firmicutes.csv")

tiff('Firmicutes.tiff', units="in", width=7, height=5, res=300)
Firmicutes <- ggplot(df_Firmicutes, aes(x = phylum_r_stats_new$Treat, y = mean_Firmicutes$avg_Firmicutes, fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Firmicutes",x="",y="")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Firmicutes$avg_Firmicutes - mean_Firmicutes$sem_Firmicutes, ymax = mean_Firmicutes$avg_Firmicutes + mean_Firmicutes$sem_Firmicutes),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()

df_Gemmatimonadetes=phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Gemmatimonadetes=sd(Gemmatimonadetes)/sqrt(n())) #stdev
mean_Gemmatimonadetes = df_Gemmatimonadetes %>% group_by(Treat) %>% mutate(avg_Gemmatimonadetes=sum(Gemmatimonadetes)/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Gemmatimonadetes, file = "df_Gemmatimonadetes.csv")

tiff('Gemmatimonadetes.tiff', units="in", width=7, height=5, res=300)
Gemmatimonadetes <- ggplot(df_Gemmatimonadetes, aes(x = phylum_r_stats_new$Treat, y = mean_Gemmatimonadetes$avg_Gemmatimonadetes, fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Gemmatimonadetes",x="Treatment",y="Relative abundance")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=90,hjust = 1),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Gemmatimonadetes$avg_Gemmatimonadetes - mean_Gemmatimonadetes$sem_Gemmatimonadetes, ymax = mean_Gemmatimonadetes$avg_Gemmatimonadetes + mean_Gemmatimonadetes$sem_Gemmatimonadetes),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()

df_Nitrospirae=phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Nitrospirae=sd(Nitrospirae)/sqrt(n())) #stdev
mean_Nitrospirae = df_Nitrospirae %>% group_by(Treat) %>% mutate(avg_Nitrospirae=sum(Nitrospirae)/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Nitrospirae, file = "df_Nitrospirae.csv")

tiff('Nitrospirae.tiff', units="in", width=7, height=5, res=300)
Nitrospirae <- ggplot(df_Nitrospirae, aes(x = phylum_r_stats_new$Treat, y = mean_Nitrospirae$avg_Nitrospirae, fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Nitrospirae",x="Treatment",y="")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=90,hjust = 1),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Nitrospirae$avg_Nitrospirae - mean_Nitrospirae$sem_Nitrospirae, ymax = mean_Nitrospirae$avg_Nitrospirae + mean_Nitrospirae$sem_Nitrospirae),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()

df_Patescibacteria  =phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Patescibacteria  =sd(Patescibacteria  )/sqrt(n())) #stdev
mean_Patescibacteria   = df_Patescibacteria   %>% group_by(Treat) %>% mutate(avg_Patescibacteria  =sum(Patescibacteria  )/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Patescibacteria  , file = "df_Patescibacteria  .csv")

tiff('Patescibacteria.tiff', units="in", width=7, height=5, res=300)
Patescibacteria   <- ggplot(df_Patescibacteria  , aes(x = phylum_r_stats_new$Treat, y = mean_Patescibacteria  $avg_Patescibacteria  , fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Patescibacteria",x="Treatment",y="")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=90,hjust = 1),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Patescibacteria  $avg_Patescibacteria   - mean_Patescibacteria  $sem_Patescibacteria  , ymax = mean_Patescibacteria  $avg_Patescibacteria   + mean_Patescibacteria  $sem_Patescibacteria  ),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()

df_Verrucomicrobia =phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Verrucomicrobia =sd(Verrucomicrobia )/sqrt(n())) #stdev
mean_Verrucomicrobia  = df_Verrucomicrobia  %>% group_by(Treat) %>% mutate(avg_Verrucomicrobia =sum(Verrucomicrobia )/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Verrucomicrobia , file = "df_Verrucomicrobia .csv")

tiff('Verrucomicrobia.tiff', units="in", width=7, height=5, res=300)
Verrucomicrobia  <- ggplot(df_Verrucomicrobia , aes(x = phylum_r_stats_new$Treat, y = mean_Verrucomicrobia $avg_Verrucomicrobia , fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Verrucomicrobia",x="",y="")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Verrucomicrobia $avg_Verrucomicrobia  - mean_Verrucomicrobia $sem_Verrucomicrobia , ymax = mean_Verrucomicrobia $avg_Verrucomicrobia  + mean_Verrucomicrobia $sem_Verrucomicrobia ),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()

df_Rokubacteria=phylum_r_stats_new %>% group_by(Treat) %>% mutate(sem_Rokubacteria=sd(Rokubacteria)/sqrt(n())) #stdev
mean_Rokubacteria = df_Rokubacteria %>% group_by(Treat) %>% mutate(avg_Rokubacteria=sum(Rokubacteria)/(n()))#combine mean and stdev in "mean_shannon"
write.csv(df_Rokubacteria, file = "df_Rokubacteria.csv")

tiff('Rokubacteria.tiff', units="in", width=7, height=5, res=300)
Rokubacteria <- ggplot(df_Rokubacteria, aes(x = phylum_r_stats_new$Treat, y = mean_Rokubacteria$avg_Rokubacteria, fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Rokubacteria",x="Treatment",y="")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=90,hjust = 1),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Rokubacteria$avg_Rokubacteria - mean_Rokubacteria$sem_Rokubacteria, ymax = mean_Rokubacteria$avg_Rokubacteria + mean_Rokubacteria$sem_Rokubacteria),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()


#arrange multiple figures in one figure
library("cowplot")
tiff('arrange_phylum.tiff', units="in", width=20, height=16, res=300)
plot_grid(Proteobacteria,Acidobacteria,Actinobacteria,Verrucomicrobia,Planctomycetes,Chloroflexi,Firmicutes,Bacteroidetes,Gemmatimonadetes,Patescibacteria,Rokubacteria,Nitrospirae,
          labels = c("A", "B", "C","D", "E", "F", "G", "H", "I","J", "K", "L"),
          ncol = 4, nrow = 3, rel_heights = c(1,1,1.2), rel_widths = c(1,1,1))
dev.off()

#fertilization--publication#
fer_phyla_dat <- phylum_r_stats_new[,c(1,2,10,12,22,35:42)]
fer_phyla_dat_long <- pivot_longer(fer_phyla_dat,cols = c("Proteobacteria","Acidobacteria","Bacteroidetes","Patescibacteria","Rokubacteria"))

df_fer_phyla_dat=fer_phyla_dat_long %>% group_by(name,Nitrogen) %>% mutate(sem=sd(value)/sqrt(n())) #stdev
mean_df_fer_phyla_dat = df_fer_phyla_dat %>% group_by(name,Nitrogen) %>% mutate(avg=sum(value)/(n()))#combine mean and stdev in "mean_shannon"

df_fer_phyla <- ggplot(mean_df_fer_phyla_dat, aes(x = name, y = avg, fill = Nitrogen))+
  geom_bar(position = 'dodge',stat="identity", colour='black',width = 0.3)+
  scale_fill_manual(values = c('#999999',"#ab6d63"),name="Treatment")+
  labs(title="",x="Treatment",y="")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=0,hjust = 0.5),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_df_fer_phyla_dat$avg - mean_df_fer_phyla_dat$sem, 
                    ymax = mean_df_fer_phyla_dat$avg + mean_df_fer_phyla_dat$sem),
                width = 0.2, colour = "black", position = position_dodge(.3))

fer_phyla_dat <- phylum_r_stats_new[,c(1,2,10,12,22,35:42)]
fer_phyla_dat_long <- pivot_longer(fer_phyla_dat,cols = c("Proteobacteria","Acidobacteria","Bacteroidetes","Patescibacteria","Rokubacteria"))

df_fer_phyla_dat=fer_phyla_dat_long %>% group_by(name,Nitrogen) %>% mutate(sem=sd(value)/sqrt(n())) #stdev
mean_df_fer_phyla_dat = df_fer_phyla_dat %>% group_by(name,Nitrogen) %>% mutate(avg=sum(value)/(n()))#combine mean and stdev in "mean_shannon"

df_fer_phyla <- ggplot(mean_df_fer_phyla_dat, aes(x = name, y = avg, fill = Nitrogen))+
  geom_bar(position = 'dodge',stat="identity", colour='black',width = 0.3)+
  scale_fill_manual(values = c('#999999',"#ab6d63"),name="Treatment")+
  labs(title="",x="Treatment",y="")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.text.x = element_text(family="Arial",size=12,angle=0,hjust = 0.5),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_df_fer_phyla_dat$avg - mean_df_fer_phyla_dat$sem, 
                    ymax = mean_df_fer_phyla_dat$avg + mean_df_fer_phyla_dat$sem),
                width = 0.2, colour = "black", position = position_dodge(.3))

#stats
library(agricolae)
treatlm_Proteobacteria=lm(phylum_r_stats_new$Proteobacteria~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE))
#phylum_r_stats_new$Proteobacteria groups
#VNTN60                          0.3256137      a
#WCN60                           0.3081989     ab
#VCN60                           0.3058665     ab
#NCCN60                          0.3051943     ab
#WCN0                            0.2906103     ab
#VNTN0                           0.2889639     ab
#WNTN0                           0.2885523     ab
#NCNTN60                         0.2856163     ab
#WNTN60                          0.2838190     ab
#NCCN0                           0.2774668      b
#NCNTN0                          0.2745994      b
treatlm_Proteobacteria=lm(phylum_r_stats_new$Acidobacteria~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE))
#phylum_r_stats_new$Acidobacteria groups
#VNTN0                          0.2299830      a
#WNTN0                          0.2163045     ab
#WCN0                           0.2074964     ab
#NCNTN60                        0.2040391     ab
#WNTN60                         0.2019674     ab
#NCNTN0                         0.1935161     ab
#NCCN0                          0.1928301     ab
#VNTN60                         0.1925511     ab
#NCCN60                         0.1846669      b
#WCN60                          0.1839946      b
#VCN60                          0.1776698      b
treatlm_Proteobacteria=lm(phylum_r_stats_new$Actinobacteria~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE)) # no sig

treatlm_Proteobacteria=lm(phylum_r_stats_new$Verrucomicrobia~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE))
#phylum_r_stats_new$Verrucomicrobia groups
#VNTN0                           0.08041104      a
#WNTN0                           0.07308473     ab
#NCNTN60                         0.07278290     ab
#VNTN60                          0.07008012     ab
#NCNTN0                          0.06899627     ab
#WNTN60                          0.06756942     ab
#WCN0                            0.06526452     ab
#WCN60                           0.05824004     ab
#NCCN60                          0.05651136     ab
#NCCN0                           0.05572934      b
#VCN60                           0.05555098      b
treatlm_Proteobacteria=lm(phylum_r_stats_new$Planctomycetes~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE)) # no sig

treatlm_Proteobacteria=lm(phylum_r_stats_new$Chloroflex~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE)) # no sig

treatlm_Proteobacteria=lm(phylum_r_stats_new$Firmicutes~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE)) 
#phylum_r_stats_new$Firmicutes groups
#NCCN0                      0.05578422      a
#WCN60                      0.05015915     ab
#VCN60                      0.04923993     ab
#NCCN60                     0.04864998     ab
#WNTN60                     0.03744101     ab
#NCNTN60                    0.03733125     ab
#NCNTN0                     0.03379157     ab
#WCN0                       0.02837230     ab
#WNTN0                      0.02524421     ab
#VNTN0                      0.01721820     ab
#VNTN60                     0.01584166      b
treatlm_Proteobacteria=lm(phylum_r_stats_new$Bacteroidetes~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE)) # no sig

treatlm_Proteobacteria=lm(phylum_r_stats_new$Gemmatimonadetes~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE)) # no sig

treatlm_Proteobacteria=lm(phylum_r_stats_new$Patescibacteria~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE)) 
#phylum_r_stats_new$Patescibacteria groups
#VNTN60                          0.02925036      a
#VCN60                           0.02184173     ab
#NCNTN60                         0.01901548     bc
#WCN60                           0.01869992     bc
#WNTN60                          0.01797278     bc
#WNTN0                           0.01455658     bc
#VNTN0                           0.01436450     bc
#NCCN60                          0.01376084     bc
#WCN0                            0.01249863      c
#NCCN0                           0.01096202      c
#NCNTN0                          0.01005653      c
treatlm_Proteobacteria=lm(phylum_r_stats_new$Rokubacteria~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE)) 
#phylum_r_stats_new$Rokubacteria groups
#VNTN0                       0.014021512      a
#WCN0                        0.012690704     ab
#NCCN0                       0.012676984     ab
#NCNTN0                      0.012224235     ab
#WNTN0                       0.011702887     ab
#WNTN60                      0.011030623     ab
#NCCN60                      0.010413237     ab
#WCN60                       0.009452859     ab
#NCNTN60                     0.009178465     ab
#VCN60                       0.007477226     ab
#VNTN60                      0.004372004      b

treatlm_Proteobacteria=lm(phylum_r_stats_new$Nitrospirae~phylum_r_stats_new$Treat,data =phylum_r_stats_new)
(HSD.test(treatlm_Proteobacteria,"phylum_r_stats_new$Treat",unbalanced=TRUE)) # no sig



#---------------------------------------------may 11-------------------















###pca ggbiplot
install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
##PCA of phylum
phylum_r_rarefied_mothur_phyloseq43 <- subset_samples(phylum_r_rarefied_mothur_phyloseq,sample_names(phylum_r_rarefied_mothur_phyloseq)!="VCN0_4")
phylum <- t(otu_table(phylum_r_rarefied_mothur_phyloseq43))
tax <- as.matrix(tax_table(phylum_r_rarefied_mothur_phyloseq43))
write.csv(tax,"dataset/tax.csv")
write.csv(phylum,file = "dataset/phylum.csv")
#write.csv(sample_data(phylum_r_rarefied_mothur_phyloseq),file = "sample_data.csv")
phylum.rev.1 <- read.csv(file = "dataset/phylum_rev.csv", sep = ",", row.names = 1)#phylum_rev is the file been revised accroding to the simper results, we selected first several phylum which cumulation explaination are over 70%

#phylum.rev <- phylum.rev.1[,-9]
pca.phylum <- prcomp(as.data.frame(phylum.rev.1),center = TRUE, scale. = TRUE)
summary(pca.phylum)
str(pca.phylum)
phylum.group_tillage <- c(rep("Tillage",8),rep("No_tillage", 8),rep("Tillage",4),rep("No_tillage", 7),rep("Tillage", 8),rep("No_tillage", 8))
ggbiplot_tillage <- 
  ggbiplot(pca.phylum,ellipse=TRUE,labels.size = 1.5, 
                             circle=FALSE,obs.scale = 1, var.scale = 1,
                             groups=phylum.group_tillage, 
                             labels = row.names(phylum),
                             var.axes=TRUE,
                             varname.size = 2.5,
                             varname.adjust =2,
                             varname.abbrev = TRUE,
                             arrow.color="blue")+
  scale_colour_manual(name="Treatment", 
                      values= c("forest green", "coral", "dark blue"))#+
  ggtitle("")+
  theme_minimal(base_size = 11)+
  theme(legend.position = "bottom")+xlim(-6, 6) + ylim(-5, 5)

phylum.group_nitrogen <- c(rep("N0",4),rep("N60", 4),rep("N0",4),rep("N60", 4),rep("N60", 4),rep("N0",4),rep("N60", 3),rep("N0",4),rep("N60", 4),rep("N0",4),rep("N60", 4))
ggbiplot_nitro <- ggbiplot(pca.phylum,ellipse=TRUE,
                           labels.size = 1.5, 
                           circle=FALSE,
                           obs.scale = 1, 
                           var.scale = 1,
                           groups=phylum.group_nitrogen, 
                           labels = row.names(phylum),
                           var.axes=TRUE,
                           varname.size = 2.5,
                           varname.adjust =2,
                           varname.abbrev = TRUE,arrow.color="blue")+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue"))#+
  ggtitle("")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-6, 6) + ylim(-5, 5)
 
phylum.group_cover <- c(rep("No_cover",16),rep("Vetch", 11),rep("Wheat",16))
ggbiplot_cover <- ggbiplot(pca.phylum,ellipse=TRUE,labels.size = 1.5, 
                           circle=FALSE,obs.scale = 1, 
                           var.scale = 1,groups=phylum.group_cover, 
                           labels = row.names(phylum),
                           var.axes=TRUE,
                           varname.size = 2.5,
                           varname.adjust =2,
                           varname.abbrev = TRUE,arrow.color="blue")+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue"))#+
  ggtitle("")+
  theme_minimal()+
  xlim(-6, 6) + ylim(-5, 5)

adonis(t(otu_table(phylum_r_rarefied_mothur_phyloseq))~Nitrogen*Cover+Tillage*Nitrogen+Cover*Tillage,
       data=data.frame(sample_data(phylum_r_rarefied_mothur_phyloseq)),
       strata=data.frame(sample_data(phylum_r_rarefied_mothur_phyloseq))$block
       ,permutations=9999,by="margin")


library("ggpubr")
tiff('figures/pca_all.tiff', units="in", width=11, height=5, res=300)
ggarrange(ggbiplot_tillage,ggbiplot_nitro,ggbiplot_cover,
          labels = c("A","B","C"),ncol = 3, nrow = 1,
          common.legend = F,legend = "bottom",hjust = -0.5,vjust =4 )
dev.off()







#CCA
library(vegan)
library(ggplot2)


#block effect
tiff('CCA blockeffect.tiff', units="in", width=7, height=5, res=300)
dev.off()
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_CCA,type="samples",color="block")+geom_point(size=4)+ggtitle("Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c","gold","aquamarine"))+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=13,family = "Arial")))
#nitrogen&cover--good
tiff('CCA nitrogenandcover.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_CCA,type="samples",color="Nitrogen", shape = "Cover")+
  geom_point(size=4)+ggtitle("")+scale_color_manual(values = c("#00798c","#d1495b","#edae49"))+
  labs(x="CCA1",y="CCA2")+
  theme(panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=13,family="Arial"),
        axis.title=element_text(size = 15,face="plain",family = "Arial"),
        legend.title = element_text(size=13,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))
dev.off()

#cover&nitrogen--not good
dev.off()
tiff('CCA coverandnitrogen.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_CCA,type="samples",color="Cover", shape = "Nitrogen")+
  geom_point(size=4)+ggtitle("")+
  scale_color_manual(values = c("#00798c","#d1495b","#edae49"))+
  labs(x="CCA1",y="CCA2 [9%]")+
  theme(panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=13,family="Arial"),
        axis.title=element_text(size = 15,face="plain",family = "Arial"),
        legend.title = element_text(size=13,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))
dev.off()

#tillage&nitrogen--good
dev.off()
tiff('CCA tillageandnitrogen.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_CCA,type="samples",color="Tillage", label="Treat", shape = "Nitrogen")+
  geom_point(size=4)+ggtitle("")+
  scale_color_manual(values = c("#d1495b","#edae49","yellow"))+
  labs(x="CCA1",y="CCA2")+
  theme(panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=10,family="Arial"),
        axis.title=element_text(size = 10,face="plain",family = "Arial"),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))
dev.off()

# Calculate bray curtis distance matrix
bray_not_na <- phyloseq::distance(rarefied_mothur_phyloseq, method = "bray")
cap_ord <- ordinate(physeq = rarefied_mothur_phyloseq, method = "CAP",distance = bray_not_na, formula = ~pH+Moisture+Temperature+NH4+NO3+PO4+MBC+MBN+TOC+TON)
# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, yend = CAP2, x = 0,  y = 0,  shape = NULL,color = NULL, label = labels)
label_map <- aes(x = 1.3 * CAP1, y = 1.3 * CAP2, shape = NULL, color = NULL,label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))

geom_segment(mapping = arrow_map, size = .5,data = arrowdf,color = "gray",arrow = arrowhead )+
  geom_text(mapping = label_map, size = 2,  data = arrowdf,  show.legend = FALSE)
dev.off()

#table1
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~Nitrogen*Cover*Tillage,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=9999,by="margin")
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~Nitrogen*Cover+Tillage*Nitrogen+Cover*Tillage,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),strata=data.frame(sample_data(r_rarefied_mothur_phyloseq))$block,permutations=9999,by="margin")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Nitrogen          1    0.4350 0.43497  4.7403 0.09294 0.0001 ***
#  Cover             2    0.2626 0.13131  1.4311 0.05612 0.0219 *  
#  Tillage           1    0.3469 0.34685  3.7800 0.07411 0.0001 ***
#  Nitrogen:Cover    2    0.2409 0.12047  1.3129 0.05148 0.0535 .  
#Nitrogen:Tillage  1    0.0940 0.09403  1.0248 0.02009 0.3532    
#Cover:Tillage     2    0.1809 0.09047  0.9859 0.03866 0.4406    
#Residuals        34    3.1198 0.09176         0.66660           
#Total            43    4.6801                 1.00000   
#Table 2

adonis(t(otu_table(r_rarefied_mothur_phyloseq))~Temp+VBR+Bacteria_abund+Virus_abund,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=5000,by="margin")







#------------------May10-------------------

##class(relative abundance)
class_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[3],NArm=TRUE)
class_r_rarefied_mothur_phyloseq#116

class_prune_r_rarefied_mothur_phyloseq<-prune_taxa(taxa_sums(class_r_rarefied_mothur_phyloseq)>=0.02,class_r_rarefied_mothur_phyloseq)
class_prune_r_rarefied_mothur_phyloseq#62 taxa

merged_r_class_rarefied_mothur_phyloseq<-merge_samples(class_prune_r_rarefied_mothur_phyloseq,"Treat")#merge
merged_r_class_rarefied_mothur_phyloseq
merged_r_class_rarefied_mothur_phyloseq<-merge_samples(class_r_rarefied_mothur_phyloseq,"Treat")#merge
merged_r_class_rarefied_mothur_phyloseq

final_prune_merged_r_class_rarefied_mothur_phyloseq<-transform_sample_counts(merged_r_class_rarefied_mothur_phyloseq,function(x) x/sum(x))
final_prune_merged_r_class_rarefied_mothur_phyloseq

sample_sums(final_prune_merged_r_class_rarefied_mothur_phyloseq)
#NCCN0  NCCN60  NCNTN0 NCNTN60   VCN60   VNTN0  VNTN60    WCN0   WCN60   WNTN0  WNTN60 
#1       1       1       1       1       1       1       1       1       1       1
class_combine<-cbind(tax_table(final_prune_merged_r_class_rarefied_mothur_phyloseq)[,3],otu_table(t(final_prune_merged_r_class_rarefied_mothur_phyloseq)))
write.csv(class_combine,file='class_combine.csv')

#then sort the "class_combine" excelfile from large to small value, subset file the species large than 0.1,new file name"phylum_combine_subset" 
class_combine_sub<-read.csv('class_combine_subset.csv',row.names = 1)
class_melt<-melt(class_combine_sub,id.vars = "Class")
class_melt# totally there have 11 phylum have been selected

library(colorspace)
library(ggplot2)
pal <- choose_palette()
cols_class = pal(12)
tiff('class_barplot', units="in", width=16, height=8, res=300)
ggplot(class_melt, 
       aes(x =variable, 
           y = value, fill=Class))+
  geom_bar(stat='identity',colour="grey",size=0.1)+
  labs(title="May18 stacked barplot_class",x="Treatment",y="Relative Abundance")+
  scale_fill_manual(values=cols_class)+
  theme(title=(element_text(size=20,family = "Times",face="bold")),
        plot.margin = unit(c(0.5,0.1,0.5,0.1),"cm"),
        axis.text.x = element_text(angle=90,hjust = 1),
        text=element_text(family = "Times",face="plain"),
        panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        panel.border = element_rect(colour = NA, fill=NA),
        axis.text=element_text(size=20,family="Times"),
        axis.title=element_text(size = 20,face="bold",family = "Times"),
        legend.title = element_text(size=20,face="bold",family = "Times"),
        legend.text = (element_text(size=20,family = "Times")))#pdf save:5x6
dev.off()

###  %>% -diversity

Diversity<-data.frame(Chaos<-estimate_richness(rarefied_mothur_phyloseq,split=TRUE,measures = NULL),sample_data(rarefied_mothur_phyloseq))
Diversity
identical(rownames(Diversity),sample_names(rarefied_mothur_phyloseq)) # if true, two objects are equal.
write.csv(Diversity,file = "alpha_diversity.csv")
##Shannon diversity
library(colorspace)
library(ggplot2)
library(agricolae)
install.packages("agricolae")

pal <- choose_palette()
cols_diversity = pal(3)
pal <- choose_palette()
cols_diversity1=pal(2)
pal <- choose_palette()
cols_treat=pal(14)
tiff('treat_shannon diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Treat,y=Shannon,fill=Treat))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Treat))+scale_fill_manual(values = cols_treat,name="Treatment")+labs(title="May18Shannon Diversity",x="Treatment",y="Shannon Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

treatlm_shannon=lm(Diversity$Shannon~Diversity$Treat,data = Diversity)
(HSD.test(treatlm_shannon,"Diversity$Treat")) # no sig

tiff('Cover_shannon diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Cover,y=Shannon,fill=Cover))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Cover))+scale_fill_manual(values = cols_diversity,name="Treatment")+labs(title="May18Shannon Diversity",x="Treatment",y="Shannon Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

coverlm_shannon = lm(Diversity$Shannon~Diversity$Cover,data=Diversity)
(HSD.test(coverlm_shannon,"Diversity$Cover"))#no sig

tiff('Tillage_shannon diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Tillage,y=Shannon,fill=Tillage))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Tillage))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="May18Shannon Diversity",x="Treatment",y="Shannon Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

tillage_shannon_ttest <- t.test(Diversity$Shannon ~ Diversity$Tillage, data = Diversity)
tillage_shannon_ttest# p-value = 0.5951

tiff('Nitrogen_shannon diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Nitrogen,y=Shannon,fill=Nitrogen))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Nitrogen))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="May18Shannon Diversity",x="Treatment",y="Shannon Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
head(Diversity)

nitrogen_shannon_ttest <- t.test(Diversity$Shannon ~ Diversity$Nitrogen, data = Diversity)
nitrogen_shannon_ttest#p-value = 0.6604

##Simpson diversity

tiff('treat_Simpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Treat,y=Simpson,fill=Treat))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Treat))+scale_fill_manual(values = cols_treat,name="Treatment")+labs(title="May18Simpson Diversity",x="Treatment",y="Simpson Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

treatlm_simpson=lm(Diversity$Simpson~Diversity$Treat,data = Diversity)
(HSD.test(treatlm_simpson,"Diversity$Treat"))#no sig

tiff('Cover_Simpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Cover,y=Simpson,fill=Cover))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Cover))+scale_fill_manual(values = cols_diversity,name="Treatment")+labs(title="May18Simpson Diversity",x="Treatment",y="Simpson Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
coverlm_simpson=lm(Diversity$Simpson~Diversity$Cover,data = Diversity)
(HSD.test(coverlm_simpson,"Diversity$Cover"))

tiff('tillage_Simpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Tillage,y=Simpson,fill=Tillage))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Tillage))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="May18Simpson Diversity",x="Treatment",y="Simpson Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
tiff('Nitrogen_Simpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Nitrogen,y=Simpson,fill=Nitrogen))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Nitrogen))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="May18Simpson Diversity",x="Treatment",y="Simpson Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

##InvSimpson diversity

tiff('treat_InvSimpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Treat,y=InvSimpson,fill=Treat))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Treat))+scale_fill_manual(values = cols_treat,name="Treatment")+labs(title="May18InvSimpson Diversity",x="Treatment",y="InvSimpson Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

treatlm_InvSimpson=lm(Diversity$InvSimpson~Diversity$Treat,data = Diversity)
(HSD.test(treatlm_InvSimpson,"Diversity$Treat"))

tiff('Cover_InvSimpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Cover,y=InvSimpson,fill=Cover))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Cover))+scale_fill_manual(values = cols_diversity,name="Treatment")+labs(title="May18InvSimpson Diversity",x="Treatment",y="InvSimpson Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

coverlm_InvSimpson = lm(Diversity$InvSimpson~Diversity$Cover,data=Diversity)
(HSD.test(coverlm_InvSimpson,"Diversity$Cover"))#no sig

tiff('tillage_InvSimpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Tillage,y=InvSimpson,fill=Tillage))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Tillage))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="May18InvSimpson Diversity",x="Treatment",y="InvSimpson Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

tillage_InvSimpson_ttest <- t.test(Diversity$InvSimpson ~ Diversity$Tillage, data = Diversity)
tillage_InvSimpson_ttest#p-value = 0.9841

tiff('Nitrogen_InvSimpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Nitrogen,y=InvSimpson,fill=Nitrogen))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Nitrogen))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="May18InvSimpson Diversity",x="Treatment",y="InvSimpson Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

nitrogen_InvSimpson_ttest <- t.test(Diversity$InvSimpson ~ Diversity$Nitrogen, data = Diversity)
nitrogen_InvSimpson_ttest #p-value = 0.5472

##Chao1 diversity
tiff('treat_Chao1 diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Treat,y=Chao1 ,fill=Treat))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Treat))+scale_fill_manual(values = cols_treat,name="Treatment")+labs(title="May18Chao1 Diversity",x="Treatment",y="Chao1 Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

treatlm_Chao1=lm(Diversity$Chao1~Diversity$Treat,data = Diversity)
(HSD.test(treatlm_Chao1,"Diversity$Treat"))

tiff('Cover_Chao1 diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Cover,y=Chao1 ,fill=Cover))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Cover))+scale_fill_manual(values = cols_diversity,name="Treatment")+labs(title="May18Chao1 Diversity",x="Treatment",y="Chao1 Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

coverlm_Chao1 = lm(Diversity$Chao1~Diversity$Cover,data=Diversity)
(HSD.test(coverlm_Chao1,"Diversity$Cover"))#no sig


tiff('tillage_Chao1 diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Tillage,y=Chao1 ,fill=Tillage))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Tillage))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="May18Chao1 Diversity",x="Treatment",y="Chao1 Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

tillage_Chao1_ttest <- t.test(Diversity$Chao1 ~ Diversity$Tillage, data = Diversity)
tillage_Chao1_ttest#p-value = 0.3719

tiff('Nitrogen_Chao1diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Nitrogen,y=Chao1 ,fill=Nitrogen))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Nitrogen))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="May18Chao1 Diversity",x="Treatment",y="Chao1 Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

nitrogen_Chao1_ttest <- t.test(Diversity$Chao1 ~ Diversity$Nitrogen, data = Diversity)
nitrogen_Chao1_ttest#p-value = 0.1862

###PCoA_r_rarefied_mothur_phyloseq at otu level
## PCoA
library(vegan)
library(ggplot2)
r_rarefied_mothur_phyloseq_bray<-vegdist(t(otu_table(r_rarefied_mothur_phyloseq)),method="bray",binary = FALSE)
r_rarefied_mothur_phyloseq_PCoA<-ordinate(r_rarefied_mothur_phyloseq,method = "PCoA", r_rarefied_mothur_phyloseq_bray)
#block effect
tiff('PCoA blockeffect.tiff', units="in", width=7, height=5, res=300)
dev.off()
#plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="samples",color="block")+geom_point(size=4)+ggtitle("Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c","gold","aquamarine"))+labs(x="PCoA1 [25.5%]",y="PCoA2 [13.3%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="samples",color="block")+geom_point(size=4)+ggtitle("Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c","gold","aquamarine"))+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=13,family = "Arial")))
#p4=plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="split",color="phylum",shape = "Nitrogen")+geom_point(size=4)+ggtitle("Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c","gold","aquamarine"))+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=13,family = "Arial")))

#nitrogen&cover--good
tiff('PCoA nitrogenandcover.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="samples",color="Nitrogen", shape = "Cover")+geom_point(size=4)+ggtitle("May18Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c"))+labs(x="PCoA1 [14.8%]",y="PCoA2 [9%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))
dev.off()

#cover&nitrogen--not good
dev.off()
tiff('PCoA coverandnitrogen.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="samples",color="Cover", shape = "Nitrogen")+geom_point(size=4)+ggtitle("May18Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c","aquamarine"))+labs(x="PCoA1 [14.8%]",y="PCoA2 [9%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))
dev.off()

#nitrogen&tillage--good
dev.off()
tiff('PCoA nitrogenandtillage.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="samples",color="Nitrogen", shape = "Tillage")+geom_point(size=4)+ggtitle("May18Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c"))+labs(x="PCoA1 [14.8%]",y="PCoA2 [9%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))
dev.off()

#tillage&nitrogen--good
dev.off()
tiff('PCoA tillageandnitrogen.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="samples",color="Tillage", shape = "Nitrogen")+geom_point(size=4)+ggtitle("May18Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c"))+labs(x="PCoA1 [14.8%]",y="PCoA2 [9%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))
dev.off()

#tillage&cover--good
dev.off()
tiff('PCoA tillageandcover.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="samples",color="Tillage", shape = "Cover")+geom_point(size=4)+ggtitle("May18Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c"))+labs(x="PCoA1 [14.8%]",y="PCoA2 [9%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))
dev.off()

#cover&tillage--not good
dev.off()
tiff('PCoA coverandtillage.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="samples",color="Cover", shape = "Tillage")+geom_point(size=4)+ggtitle("May18Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c","aquamarine"))+labs(x="PCoA1 [14.8%]",y="PCoA2 [9%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))
dev.off()

#betadisper
phylum.group_tillage <- c(rep("Tillage",8),rep("No_tillage", 8),rep("Tillage",4),rep("No_tillage", 7),rep("Tillage", 8),rep("No_tillage", 8))
phylum.group_nitrogen <- c(rep("N0",4),rep("N60", 4),rep("N0",4),rep("N60", 4),rep("N60", 4),rep("N0",4),rep("N60", 3),rep("N0",4),rep("N60", 4),rep("N0",4),rep("N60", 4))
phylum.group_cover <- c(rep("No_cover",16),rep("Vetch", 11),rep("Wheat",16))

mod_tillage <- betadisper(r_rarefied_mothur_phyloseq_bray,phylum.group_tillage,bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
plot(mod_tillage,hull = FALSE, ellipse = TRUE)
anova(mod_tillage)
permutest(mod_tillage)

mod_nitrogen <- betadisper(r_rarefied_mothur_phyloseq_bray,phylum.group_nitrogen,bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
plot(mod_nitrogen,hull = FALSE, ellipse = TRUE)
anova(mod_nitrogen)
permutest(mod_nitrogen)
(mod.HSD_nitrogen <- TukeyHSD(mod_nitrogen))
plot(mod.HSD_nitrogen)


mod_cover <- betadisper(r_rarefied_mothur_phyloseq_bray,phylum.group_cover,bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
plot(mod_cover,hull = FALSE, ellipse = TRUE)
anova(mod_cover)
permutest(mod_cover)

###NMDS
library(dplyr)
library(scales)
prune_r_rarefied_mothur_phyloseq 
set.seed(12345)
prune_r_rarefied_mothur_phyloseq_ordinate <- ordinate(prune_r_rarefied_mothur_phyloseq,"NMDS", "bray")
dev.off()
tiff('NMDS_cover.tiff', units="in", width=7, height=5, res=300)
plot_ordination(prune_r_rarefied_mothur_phyloseq, prune_r_rarefied_mothur_phyloseq_ordinate, type="samples", color="Cover", title="May18Cover_nmds")
dev.off()
tiff('NMDS_tillage.tiff', units="in", width=7, height=5, res=300)
plot_ordination(prune_r_rarefied_mothur_phyloseq, prune_r_rarefied_mothur_phyloseq_ordinate, type="samples", color="Tillage", title="May18Tillage_nmds")
dev.off()
tiff('NMDS_nitrogen.tiff', units="in", width=7, height=5, res=300)
plot_ordination(prune_r_rarefied_mothur_phyloseq, prune_r_rarefied_mothur_phyloseq_ordinate, type="samples", color="Nitrogen", title="May18Nitrogen_nmds")
dev.off()

###pca ggbiplot
install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
##PCA of phylum
phylum <- t(otu_table(phylum_r_rarefied_mothur_phyloseq))
write.csv(phylum,file = "phylum.csv")
write.csv(sample_data(phylum_r_rarefied_mothur_phyloseq),file = "sample_data.csv")
phylum.rev <- read.csv(file = "phylum_rev.csv", sep = ",", row.names = 1)#phylum_rev is the file been revised accroding to the simper results, we selected first several phylum which cumulation explaination are over 70%
phylum.rev
pca.phylum <- prcomp(as.data.frame(phylum.rev),center = TRUE, scale. = TRUE)
summary(pca.phylum)
str(pca.phylum)
phylum.group_tillage <- c(rep("Tillage",8),rep("No_tillage", 8),rep("Tillage",4),rep("No_tillage", 7),rep("Tillage", 8),rep("No_tillage", 8))
tiff('pca.phylum_tillage.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.phylum,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=phylum.group_tillage, labels = row.names(phylum),var.axes=TRUE,varname.size = 2.5,varname.adjust =1.5
         ,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue"))+
  ggtitle("May18_PCA of phylum")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-4, 4)
dev.off()

phylum.group_nitrogen <- c(rep("N0",4),rep("N60", 4),rep("N0",4),rep("N60", 4),rep("N60", 4),rep("N0",4),rep("N60", 3),rep("N0",4),rep("N60", 4),rep("N0",4),rep("N60", 4))
tiff('pca.phylum_nitrogen.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.phylum,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=phylum.group_nitrogen, labels = row.names(phylum),var.axes=TRUE,varname.size = 2.5,varname.adjust =1.5
         ,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue"))+
  ggtitle("May18_PCA of phylum")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-4, 4)
dev.off()

phylum.group_cover <- c(rep("No_cover",16),rep("Vetch", 11),rep("Wheat",16))
tiff('pca.phylum_cover.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.phylum,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=phylum.group_cover, labels = row.names(phylum),var.axes=TRUE,varname.size = 2.2,varname.adjust =1.5
         ,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue"))+
  ggtitle("May18_PCA of phylum")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-4, 4)
dev.off()

##pca of class
class <- t(otu_table(class_r_rarefied_mothur_phyloseq))
write.csv(class,file = "class.csv")
class.rev <- read.csv(file = "class_rev.csv", sep = ",", row.names = 1)#phylum_rev is the file been revised accroding to the simper results, we selected first several phylum which cumulation explaination are over 70%
pca.class <- prcomp(as.data.frame(class.rev),center = TRUE, scale. = TRUE)
summary(pca.class)
str(pca.class)
class.group_tillage <- c(rep("Tillage",8),rep("No_tillage", 8),rep("Tillage",4),rep("No_tillage", 7),rep("Tillage", 8),rep("No_tillage", 8))
tiff('pca.class_tillage.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.class,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=class.group_tillage, labels = row.names(class),var.axes=TRUE,varname.size = 2.2,varname.adjust =1.2,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue"))+
  ggtitle("PCA of class")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-4, 4)
dev.off()

class.group_nitrogen <- c(rep("N0",4),rep("N60", 4),rep("N0",4),rep("N60", 4),rep("N60", 4),rep("N0",4),rep("N60", 3),rep("N0",4),rep("N60", 4),rep("N0",4),rep("N60", 4))
tiff('pca.class_nitrogen.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.class,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=class.group_nitrogen, labels = row.names(class),var.axes=TRUE,varname.size = 2.2,varname.adjust =1.2,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue"))+
  ggtitle("PCA of class")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-4, 4)
dev.off()

class.group_cover <- c(rep("No_cover",16),rep("Vetch", 11),rep("Wheat",16))
tiff('pca.class_cover.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.class,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=class.group_cover, labels = row.names(class),var.axes=TRUE,varname.size = 2.2,varname.adjust =1.2,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue"))+
  ggtitle("PCA of class")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-4, 4)
dev.off()

###CAP plot(linear)?
library(ggplot2)
rarefied_mothur_phyloseq_cap<-ordinate(rarefied_mothur_phyloseq,'CAP','bray',~Treat+Nitrogen+Tillage+Cover)
dev.off()
tiff('CAP rarefied_mothur_phyloseq.tiff', units="in", width=7, height=5, res=300)
cap_plot<-plot_ordination(rarefied_mothur_phyloseq,rarefied_mothur_phyloseq_cap,type="samples",color="Tillage",shape = "Cover")+ geom_point(size=4)+scale_shape_manual(values = c(17,16,15))+scale_color_manual(values = c("#838B8B","#e83535","#10cc8e","#0A0A0A","#020de2","#ab0ef9","#75891b","#f9830c"))+ ggtitle("CAP plot of soil microbiome \n ~Soil_Treatment+Tillage+Cover") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
plot_ordination(rarefied_mothur_phyloseq,rarefied_mothur_phyloseq_cap,type="samples",color="Cover",shape = "Tillage")+ geom_point(size=4)+scale_shape_manual(values = c(17,16,15))+scale_color_manual(values = c("#838B8B","#e83535","#10cc8e","#0A0A0A","#020de2","#ab0ef9","#75891b","#f9830c"))+ ggtitle("CAP plot of soil microbiome \n ~Soil_Treatment+Tillage+Cover") +theme(plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="bold",family = "Arial"),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
ls(rarefied_mothur_phyloseq_cap)
rarefied_mothur_phyloseq_cap$adjust
# Calculate bray curtis distance matrix
bray_not_na <- phyloseq::distance(physeq = rarefied_mothur_phyloseq, method = "bray")
#cap_ord <- ordinate(physeq = rarefied_mothur_phyloseq, method = "CAP",distance = bray_not_na, formula = ~ Cover+Nitrogen+Tillage+Temp+pH+WHC+POXC+NH4.N+NO3.N+WEOC+Virus_abund+Bacteria_abund+VBR+Read_depth)
cap_ord <- ordinate(physeq = rarefied_mothur_phyloseq, method = "CAP",distance = bray_not_na, formula = ~ Cover+Nitrogen+Tillage)
#cap_ord <- ordinate(physeq = rarefied_mothur_phyloseq, method = "CAP",distance = bray_not_na, formula = ~ Cover+Nitrogen+Tillage)

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, yend = CAP2, x = 0,  y = 0,  shape = NULL,color = NULL, label = labels)
label_map <- aes(x = 1.3 * CAP1, y = 1.3 * CAP2, shape = NULL, color = NULL,label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))
# Make a new graphic
dev.off()
tiff('CAP_NIFA.tiff', units="in", width=7, height=5, res=300) #+ ggtitle("Agriculture soil ~ Treatment \n CAP plot ")#save as tiff
cap_plot + geom_segment( mapping = arrow_map, size = .5,data = arrowdf,color = "gray",arrow = arrowhead )+geom_text(mapping = label_map, size = 4,  data = arrowdf,  show.legend = FALSE)
dev.off()
## you can also do cover and tillage cap
#~Cover
rarefied_mothur_phyloseq_cap_cover<-ordinate(rarefied_mothur_phyloseq,'CAP','bray',~Cover)
tiff('CAP_NIFA_cover.tiff', units="in", width=7, height=5, res=300) #+ ggtitle("Agriculture soil ~ Treatment \n CAP plot ")#save as tiff
plot_ordination(rarefied_mothur_phyloseq,rarefied_mothur_phyloseq_cap_cover,type="samples",color="Cover")+ggtitle("CAP_Cover")+geom_point(size=4)+scale_color_manual(values = c("#838B8B","#e83535","#10cc8e","#0A0A0A","#020de2","#ab0ef9","#75891b","#f9830c")) +theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()
#~Tillage
rarefied_mothur_phyloseq_cap_tillage<-ordinate(rarefied_mothur_phyloseq,'CAP','bray',~Tillage)
tiff('CAP_NIFA_tillage.tiff', units="in", width=7, height=5, res=300) #+ ggtitle("Agriculture soil ~ Treatment \n CAP plot ")#save as tiff
plot_ordination(rarefied_mothur_phyloseq,rarefied_mothur_phyloseq_cap_tillage,type="samples",color="Tillage")+ggtitle("CAP_tillage")+geom_point(size=4)+scale_color_manual(values = c("#838B8B","#e83535","#10cc8e","#0A0A0A","#020de2","#ab0ef9","#75891b","#f9830c")) +theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()# or dev.set(dev.prev())
#~Nitrogen
rarefied_mothur_phyloseq_cap_nitrogen<-ordinate(rarefied_mothur_phyloseq,'CAP','bray',~Nitrogen)
tiff('CAP_NIFA_nitrogen.tiff', units="in", width=7, height=5, res=300) #+ ggtitle("Agriculture soil ~ Treatment \n CAP plot ")#save as tiff
plot_ordination(rarefied_mothur_phyloseq,rarefied_mothur_phyloseq_cap_nitrogen,type="samples",color="Nitrogen")+ggtitle("CAP_nitrogen")+geom_point(size=4)+scale_color_manual(values = c("#838B8B","#e83535","#10cc8e","#0A0A0A","#020de2","#ab0ef9","#75891b","#f9830c")) +theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="bold",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title.x=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=10,b=0,l=0,r=0)),axis.title.y=element_text(size = 15,face="bold",family = "Arial",margin=margin(t=0,b=0,l=0,r=10)),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()# or dev.set(dev.prev())


###PCA
library(stats)
library(ggplot2)
pca <- prcomp(t(otu_table(prune_r_rarefied_mothur_phyloseq)),scale = TRUE) #prcomp() returns three things.1.x 2.sedv 3. rotation
pca.data <- data.frame(sample = rownames(pca$x), X=pca$x[,1],Y=pca$x[,2])
pca.data
plot(pca$x[,1],pca$x[,2])#we use first two columns in x to draw a 2d plot, x contains PCs for drawing a graph,the first PC accounts for the most variation in the original data(the otu across all the samples)
pca.var <- pca$sdev^2 #square of sdev which stands for "standard deviation" to calculate how much variation in the original data each principal component accounts for
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per,main = "Scree plot",xlab = "principle component",ylab = "percentage variation")
tiff('PCA.tiff', units="in", width=7, height=5, res=300)
ggplot(data = pca.data,aes(x=X,y=Y,label = sample))+geom_text()+xlab(paste("PC1 -", pca.var.per[1],"%",sep = ""))+ylab(paste("PC2 -", pca.var.per[2],"%",sep = ""))+theme_bw()+ggtitle("PCA graph")
dev.off()
#e.g. 1 unit long vector, consisting of 0.97 parts of otu1 and 0.242 otu2 is called singular vector, or the Eigenvector for PC1, to make pc1 mix 0,97 parts of otu1 with 0.242 parts of otu2... and the proportions of each otu are called "loading scores", it is simlar for PC2
loading_scores <- pca$rotation[,1]
otu_scores <- abs(loading_scores)
otu_scores_ranked <- sort(otu_scores,decreasing = TRUE)
top10_otus <- names(otu_scores_ranked[1:10])
top10_otus
pca$rotation[top10_otus,1]#show the scores(and +/- sign)


#####VBR ration barplot
head(Diversity)
library(car)
library(MASS)
library(fBasics)
library(lme4)
library(tidyverse)
normalTest(log10(Diversity$VBR),"sw")
vbr.lm <- lm(log10(VBR)~ Cover*Nitrogen*Tillage,data = Diversity)
vbr.mod <- lmer(log10(VBR)~ Cover*Nitrogen*Tillage +(1|block), data = Diversity)
Anova(vbr.mod)#
anova(vbr.mod)
summary(vbr.mod)
coef(vbr.mod)

normalTest(sqrt(Diversity$Virus_abund),"sw")
vbr.lm <- lm(sqrt(Virus_abund)~ Cover*Nitrogen*Tillage,data = Diversity)
vbr.mod <- lmer(sqrt(Virus_abund)~ Cover*Nitrogen*Tillage +(1|block), data = Diversity)
Anova(vbr.mod)#
anova(vbr.mod)
summary(vbr.mod)
coef(vbr.mod)

normalTest(log10(Diversity$Bacteria_abund),"sw")
vbr.lm <- lm(log10(Virus_abund)~ Cover*Nitrogen*Tillage,data = Diversity)
vbr.mod <- lmer(log10(Bacteria_abund)~ Cover*Nitrogen*Tillage +(1|block), data = Diversity)
Anova(vbr.mod)#
anova(vbr.mod)
summary(vbr.mod)
coef(vbr.mod)


ggplot(Diversity,aes(x=Treat,y=VBR))+
  geom_boxplot(size=0.5,alpha=0.8)+
  geom_jitter(shape=16, position=position_jitter(0),aes(fill=Treat))+
#  scale_fill_manual(values = c("#4D897C", "#C6B7EC", "#ba892f","#c76d70"),name="Treatment")+
  labs(title="",x="Treatment",y="VBR")+
  theme_bw() 
dev.off()



tiff('Proteobacteria.tiff', units="in", width=7, height=5, res=300)
ggplot(df_Proteobacteria, aes(x = phylum_r_stats_new$Treat, y = mean_Proteobacteria$avg_Proteobacteria, fill = Treat))+
  geom_bar(position = 'dodge',stat="identity", colour='black')+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.5)+
  scale_fill_manual(values = cols_phylum,name="Treatment")+
  labs(title="Proteobacteria",x="",y="Relative abundance")+
  theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=10,family="Arial"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size = 10,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))+
  geom_errorbar(aes(ymin=mean_Proteobacteria$avg_Proteobacteria - mean_Proteobacteria$sem_Proteobacteria, ymax = mean_Proteobacteria$avg_Proteobacteria + mean_Proteobacteria$sem_Proteobacteria),width = 0.2, colour = "black", position = position_dodge(.9))
dev.off()


##########pcoa biplot#########
head(phylum.rev.1)
# phylum_r_rarefied_mothur_phyloseq43
# phylum.rev.otu_table <- otu_table(t(phylum.rev.1), taxa_are_rows = T)
# phylum.rev.phyloseq <- merge_phyloseq(phylum.rev.otu_table, sample_data(phylum_r_rarefied_mothur_phyloseq43))
colnames(phylum.rev.1)
colnames(phylum.rev.1) <- c("Prt." , "Acd." ,  "Act." , "Frm.",   "Ntr.",  "Vrr.",
                            "Chl." , "Gmm."  , "Bct.", "Pln.", "Rok." ,   "Pts.")
phylum.rev_bray<-vegdist((phylum.rev.1),method="bray",binary = FALSE)
phylum.rev.pcoa <- cmdscale(phylum.rev_bray,
                            k=(nrow((phylum.rev.1)))-1, eig = T)

pcoa_score <- scores(phylum.rev.pcoa, choices = c(1, 2))
colnames(pcoa_score) <- c("PC1 (46.3%)", "PC2 (31.2%)")
env.data <- as.data.frame(sample_data(phylum_r_rarefied_mothur_phyloseq43))
phylum.wa <- wascores(phylum.rev.pcoa$points[, 1:2]*5, phylum.rev.1)[-c(5,8,10),]

tiff('figures/PCoA_bioplot.tiff', units="in", width=10, height=3.5, res=300)
par(mfrow= c(1,3))
ordiplot(pcoa_score,
         type = "n",
         main = "",xlim = c(-0.1,0.1),ylim = c(-0.1,0.1))
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
points(pcoa_score, 
       col = env.data$Tillage,pch=19)
# Add weighted average projection of species
text(phylum.wa, rownames(phylum.wa), cex = 0.85, col = "darkcyan")
#phylum.wa.arrow <- wascores(phylum.rev.pcoa$points[, 1:2]*4.5, phylum.rev.1)
legend("bottomleft", legend = c("Tillage","No-tillage"), pch = 19,
       col = c("black","red"))
#arrows(0,0,phylum.wa.arrow[,1],phylum.wa.arrow[,2], col = "darkcyan",code=2, length = 0.05)

ordiplot(pcoa_score,
         type = "n",
         main = "",xlim = c(-0.1,0.1),ylim = c(-0.1,0.1))
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
points(pcoa_score, 
       col = env.data$Nitrogen,pch=19)
# Add weighted average projection of species
text(phylum.wa, rownames(phylum.wa), cex = 0.85, col = "darkcyan")
#phylum.wa.arrow <- wascores(phylum.rev.pcoa$points[, 1:2]*4.5, phylum.rev.1)
legend("bottomleft", legend = c("N0","N60"), pch = 19,
       col = c("black","red"))
#arrows(0,0,phylum.wa.arrow[,1],phylum.wa.arrow[,2], col = "darkcyan",code=2, length = 0.05)

ordiplot(pcoa_score,
         type = "n",
         main = "",xlim = c(-0.1,0.1),ylim = c(-0.1,0.1))
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
points(pcoa_score, 
       col = env.data$Cover,pch=19)
# Add weighted average projection of species
text(phylum.wa, rownames(phylum.wa), cex = 0.85, col = "darkcyan")
#phylum.wa.arrow <- wascores(phylum.rev.pcoa$points[, 1:2]*4.5, phylum.rev.1)
legend("bottomleft", legend = c("No-cover","Wheat","Vetch"), pch = 19,
       col = c("black","green","red"))
#arrows(0,0,phylum.wa.arrow[,1],phylum.wa.arrow[,2], col = "darkcyan",code=2, length = 0.05)
dev.off()



