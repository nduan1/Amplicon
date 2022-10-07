#-----------------May.07.2020-----------------------
#another way to merge mothur(more complicated)
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)
install_phyloseq(branch = "devel")
install_phyloseq(branch = "github")
#change the tsv file to csv file, and then revise the sample names and read in r
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)
library(reshape2)
abund_table<- read.csv("samplegreenland_OTUs_counts.csv", row.names = 1)
taxonomy <- read.csv("greenland_OTUs_taxonomy.csv", row.names =1)
meta_table <- read.csv("glmothur_meta.csv", row.names =1)
OTU = otu_table(as.data.frame(abund_table), taxa_are_rows = TRUE)
TAX = tax_table(t(as.data.frame(taxonomy)))
SAM = sample_data(meta_table)
otutax <- phyloseq(OTU,t(TAX))
phyloseq <- merge_phyloseq(otutax,SAM)
phyloseq
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 15465 taxa and 48 samples ]
#sample_data() Sample Data:       [ 48 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 15465 taxa by 6 taxonomic ranks ]
#phyloseq_up<-prune_samples(sample_names(phyloseq)[-c(6,42,30)],phyloseq)#remove the samples have seqs less than 2000, here is C_20_30V_1,T2_20_30V-4,T2_20_30V_3
#sample_sums(phyloseq_up)
sort(sample_sums(phyloseq),decreasing = TRUE)
### bar plot of read depth
Read_depth<-sample_sums(phyloseq)
summary(Read_depth)
sample_names(SAM)
mothur_meta<-prune_samples(sample_names(SAM),SAM)
dim(mothur_meta)#number of rows and columns
mothur_meta_rev<-as.data.frame(Read_depth,Read_depth=sample_sums(phyloseq))
sample_ID<- rownames(mothur_meta_rev)

tiff('Sequencing depth barplot.tiff', units="in", width=16, height=8, res=300)
ggplot(data=as.data.frame(mothur_meta), aes(x=sample_ID,y=Read_depth))+geom_bar(stat = "identity") + xlab("Sample_ID")+ylab("Sequencing depth") + scale_y_continuous(labels = scales::comma)+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(family="Arial"), plot.title = element_text(family='Arial', hjust = .5,vjust = 3,size = 24,face="bold"), axis.text.x = element_text(family="Arial",size=8,angle=90),axis.text.y = element_text(family="Arial",size=14),axis.title.x = element_text(family="Arial",size=20,face="bold"),axis.title.y = element_text(family="Arial",size=20,face="bold"))
dev.off()
sort(Read_depth,decreasing = TRUE)

phyloseq#15465 taxa
sort(sample_sums(phyloseq),decreasing = FALSE)

#step2 singleton remove(very important to remove the singleton, remove the otu number is less than 1, you can also set the number is 5 or 10 according to your results)
rms_mothur_phyloseq<-filter_taxa(phyloseq,function (x) {sum(x > 0) > 1}, prune=TRUE)
#rms_mothur_phyloseq<-prune_taxa(taxa_sums(phyloseq_up)>16, phyloseq_up) do not remove singleton, just remove the total number of sequences which less than 16
rms_mothur_phyloseq
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3573 taxa and 48 samples ]
#sample_data() Sample Data:       [ 48 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 3573 taxa by 6 taxonomic ranks ]
#write.csv(otu_table(rms_mothur_phyloseq), file = "beforerare.csv") 
phylum_mothur_phyloseq <-tax_glom(rms_mothur_phyloseq,taxrank = rank_names(rms_mothur_phyloseq)[2],NArm=TRUE)
phylum_mothur_phyloseq#25
class_mothur_phyloseq<-tax_glom(rms_mothur_phyloseq,taxrank = rank_names(rms_mothur_phyloseq)[3],NArm=TRUE)
class_mothur_phyloseq#54
order_mothur_phyloseq<-tax_glom(rms_mothur_phyloseq,taxrank = rank_names(rms_mothur_phyloseq)[4],NArm=TRUE)
order_mothur_phyloseq#95
family_mothur_phyloseq<-tax_glom(rms_mothur_phyloseq,taxrank = rank_names(rms_mothur_phyloseq)[5],NArm=TRUE)
family_mothur_phyloseq#104
genus_mothur_phyloseq<-tax_glom(rms_mothur_phyloseq,taxrank = rank_names(rms_mothur_phyloseq)[6],NArm=TRUE)
genus_mothur_phyloseq#137


#step3 rarefied data (1.rarefied--relative abundance calculate--category(phylum, class,...)--merge(replications)--prune(relative abundance>0.02)
rarefied_mothur_phyloseq<-rarefy_even_depth(rms_mothur_phyloseq,sample.size = min(sample_sums(rms_mothur_phyloseq)),rngseed = 1013,replace = TRUE,trimOTUs = TRUE,verbose = TRUE)#normalize the sequences based on the least number of sequences
rarefied_mothur_phyloseq#3559 taxa
#write.csv(otu_table(rarefied_mothur_phyloseq), file = "afterrare.csv")
phylum_rarefied_mothur_phyloseq <-tax_glom(rarefied_mothur_phyloseq,taxrank = rank_names(rarefied_mothur_phyloseq)[2],NArm=TRUE)
phylum_rarefied_mothur_phyloseq#25
class_rarefied_mothur_phyloseq<-tax_glom(rarefied_mothur_phyloseq,taxrank = rank_names(rarefied_mothur_phyloseq)[3],NArm=TRUE)
class_rarefied_mothur_phyloseq#54
order_rarefied_mothur_phyloseq<-tax_glom(rarefied_mothur_phyloseq,taxrank = rank_names(rarefied_mothur_phyloseq)[4],NArm=TRUE)
order_rarefied_mothur_phyloseq#95
family_rarefied_mothur_phyloseq<-tax_glom(rarefied_mothur_phyloseq,taxrank = rank_names(rarefied_mothur_phyloseq)[5],NArm=TRUE)
family_rarefied_mothur_phyloseq#111
genus_rarefied_mothur_phyloseq<-tax_glom(rarefied_mothur_phyloseq,taxrank = rank_names(rarefied_mothur_phyloseq)[6],NArm=TRUE)
genus_rarefied_mothur_phyloseq#137
###otu counts in each sample
#before rarefaction
#genfac = factor(tax_table(phyloseq_up)[, "Genus"])
#gentab = apply(otu_table(phyloseq_up), MARGIN = 2, function(x) {
#  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
#})
#head(gentab)[, 1:10]
#observationThreshold = 1
#counts_otu_before <- apply(gentab > observationThreshold, 2, sum)
#counts_otu_before_tab <- as.data.frame(counts_otu_before)
#after
#genfac_after = factor(tax_table(rarefied_mothur_phyloseq)[, "Genus"])
#gentab_after = apply(otu_table(rarefied_mothur_phyloseq), MARGIN = 2, function(x) {
#  tapply(x, INDEX = genfac_after, FUN = sum, na.rm = TRUE, simplify = TRUE)
#})
#head(gentab_after)[, 1:10]
#observationThreshold = 1
#counts_otu_after <- apply(gentab_after > observationThreshold, 2, sum)
#counts_otu_after_tab <- as.data.frame(counts_otu_after)
#row.names(counts_otu_after_tab)
#dim(counts_otu_after_tab)#[1] 170   1
#dim(counts_otu_before_tab)#[1] 170   1
#otucounts <- cbind(counts_otu_before_tab,counts_otu_after_tab)
#otucounts$colnames <- c("sample","before","after")
#write.csv(otucounts,file = "otucounts_rare.csv")
#otu_counts <- read.csv("otucounts_rare_new.csv",sep=",",header=TRUE, row.names = 1)#add colnames
#ggplot(otu_counts, aes(fill=otu_counts$counts_otu_before, y=value, x=otu_counts$sample)) + 
#  geom_bar(position="dodge", stat="identity")
#ggplot(otu_counts, aes(x = Fer_cover2, y = mean_shannon$avg_shannon, fill = Tillage))+
#  geom_col(alpha = 0.5, position = 'dodge')+
#  ylim(0,8)+
#  scale_fill_manual(values = c("yellow","aquamarine4"),name="Treatment")+
#  labs(title="Diversity index",x="Treatment",y="Shannon")

### relative abundance
r_rarefied_mothur_phyloseq<-transform_sample_counts(rarefied_mothur_phyloseq,function(x)x/sum(x))
r_rarefied_mothur_phyloseq#1093 taxa
#prune_r_rarefied_mothur_phyloseq<-prune_taxa(taxa_sums(r_rarefied_mothur_phyloseq)>=0.02,r_rarefied_mothur_phyloseq)
#prune_r_rarefied_mothur_phyloseq#319
colnames(tax_table(r_rarefied_mothur_phyloseq))<-c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
rank_names(r_rarefied_mothur_phyloseq)
##phylum(relative abundance)
phylum_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[2],NArm=TRUE)
phylum_r_rarefied_mothur_phyloseq #7 phylum
#phylum_prune_r_rarefied_mothur_phyloseq<-prune_taxa(taxa_sums(phylum_r_rarefied_mothur_phyloseq)>=0.02,phylum_r_rarefied_mothur_phyloseq)
#write.csv(otu_table(phylum_prune_r_rarefied_mothur_phyloseq),file = "before.csv")

merged_r_phylum_rarefied_mothur_phyloseq<-merge_samples(phylum_r_rarefied_mothur_phyloseq,"Treat")#merge
merged_r_phylum_rarefied_mothur_phyloseq
depth_merged_r_phylum_rarefied_mothur_phyloseq<-merge_samples(phylum_r_rarefied_mothur_phyloseq,"Depth")
depth_merged_r_phylum_rarefied_mothur_phyloseq
sample_merged_r_phylum_rarefied_mothur_phyloseq<-merge_samples(phylum_r_rarefied_mothur_phyloseq,"Sample")
sample_merged_r_phylum_rarefied_mothur_phyloseq

tax_table(merged_r_phylum_rarefied_mothur_phyloseq)[,2]
r_merged_r_phylum_rarefied_mothur_phyloseq<-transform_sample_counts(merged_r_phylum_rarefied_mothur_phyloseq,function(x) x/sum(x))
r_merged_r_phylum_rarefied_mothur_phyloseq#by treat
prune_r_merged_r_phylum_rarefied_mothur_phyloseq<-prune_taxa(taxa_sums(r_merged_r_phylum_rarefied_mothur_phyloseq)>=0.02,r_merged_r_phylum_rarefied_mothur_phyloseq)
prune_r_merged_r_phylum_rarefied_mothur_phyloseq

tax_table(depth_merged_r_phylum_rarefied_mothur_phyloseq)[,2]
r_depth_merged_r_phylum_rarefied_mothur_phyloseq<-transform_sample_counts(depth_merged_r_phylum_rarefied_mothur_phyloseq,function(x) x/sum(x))
r_depth_merged_r_phylum_rarefied_mothur_phyloseq#by depth
prune_r_depth_merged_r_phylum_rarefied_mothur_phyloseq<-prune_taxa(taxa_sums(r_depth_merged_r_phylum_rarefied_mothur_phyloseq)>=0.02,r_depth_merged_r_phylum_rarefied_mothur_phyloseq)
prune_r_depth_merged_r_phylum_rarefied_mothur_phyloseq

r_sample_merged_r_phylum_rarefied_mothur_phyloseq<-transform_sample_counts(sample_merged_r_phylum_rarefied_mothur_phyloseq,function(x) x/sum(x))
r_sample_merged_r_phylum_rarefied_mothur_phyloseq#by sample
prune_r_sample_merged_r_phylum_rarefied_mothur_phyloseq<-prune_taxa(taxa_sums(r_sample_merged_r_phylum_rarefied_mothur_phyloseq)>=0.02,r_sample_merged_r_phylum_rarefied_mothur_phyloseq)
prune_r_sample_merged_r_phylum_rarefied_mothur_phyloseq

#treat
phylum_combine<-cbind(tax_table(r_merged_r_phylum_rarefied_mothur_phyloseq)[,2],otu_table(t(_r_merged_r_phylum_rarefied_mothur_phyloseq)))
write.csv(phylum_combine,file='phylum_combine.csv')
phylum_combine<-read.csv('phylum_combine.csv',row.names = 1)
#depth
depth_phylum_combine<-cbind(tax_table(r_depth_merged_r_phylum_rarefied_mothur_phyloseq)[,2],otu_table(t(r_depth_merged_r_phylum_rarefied_mothur_phyloseq)))
write.csv(depth_phylum_combine,file='depth_phylum_combine.csv')
depth_phylum_combine<-read.csv('depth_phylum_combine.csv',row.names = 1)
#sample
sample_phylum_combine<-cbind(tax_table(r_sample_merged_r_phylum_rarefied_mothur_phyloseq)[,2],otu_table(t(r_sample_merged_r_phylum_rarefied_mothur_phyloseq)))
write.csv(sample_phylum_combine,file='sample_phylum_combine.csv')
sample_phylum_combine<-read.csv('sample_phylum_combine.csv',row.names = 1)

##then sort the "phylum_combine" excelfile from large to small value, subset file the species large than 0.1,new file name"phylum_combine_subset" 
phylum_combine_sub<-read.csv('phylum_combine.csv',row.names = 1)
phylum_melt<-melt(phylum_combine_sub,id.vars = "Phylum")
head(phylum_melt)

depth_phylum_combine_sub<-read.csv('depth_phylum_combine.csv',row.names = 1)
depth_phylum_melt<-melt(depth_phylum_combine_sub,id.vars = "Phylum")
head(depth_phylum_melt)

sample_phylum_combine_sub<-read.csv('sample_phylum_combine.csv',row.names = 1)
sample_phylum_melt<-melt(sample_phylum_combine_sub,id.vars = "Phylum")
head(sample_phylum_melt)

###rarefaction plot
library(reshape2)
library(vegan)
library(ggplot2)

moth_merge_otu_table<-t(otu_table(phyloseq))
dim(moth_merge_otu_table) # correct OTUs number
rare<-rarefy(moth_merge_otu_table,sample = seq(0,max(sample_sums(phyloseq)),by=300))
rownames(rare)<-rownames(moth_merge_otu_table)
long_rare<-t(rare)
long_rare<-data.frame(subsample=seq(0,max(sample_sums(phyloseq)),by=300),long_rare)
head(long_rare)
rare_melt<-melt(long_rare,id.vars = "subsample")# reshape2 package
head(rare_melt)
## 1) Rarefaction curve using maximum read depth

ggplot(rare_melt,aes(x=subsample,y=value,colour=variable))+geom_point()+geom_line(linetype = "dashed")+theme(legend.position="top")+ylab("OTU numbers")+theme(title=(element_text(size=15,family = "Times",face="bold")),text=element_text(family = "Times",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Times"), axis.title=element_text(size = 15,face="bold",family = "Times"),legend.title = element_text(size=13,face="bold",family = "Times"),legend.text = (element_text(size=10,family = "Times")),legend.position="right")+ggtitle("Rarefaction curve with maximum read depth")

## 2) Rarefaction curve using minimum read depth

rare_min<-rarefy(moth_merge_otu_table,sample = seq(0,min(sample_sums(phyloseq)),by=300))
rownames(rare_min)<-rownames(moth_merge_otu_table)
long_rare_min<-t(rare_min)
long_rare_min<-data.frame(subsample=seq(0,min(sample_sums(phyloseq)),by=300),long_rare_min)
head(long_rare_min)
rare_melt_min<-melt(long_rare_min,id.vars = "subsample")
head(rare_melt_min)

tiff('Rarefaction_curve', units="in", width=16, height=8, res=300)
ggplot(rare_melt_min,aes(x=subsample,y=value,colour=variable))+geom_point()+geom_line(linetype = "dashed")+theme(legend.position="top")+ylab("OTU numbers")+theme(title=(element_text(size=15,family = "Times",face="bold")),text=element_text(family = "Times",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Times"))
dev.off()

#stacked barplot                                                                                                                                                                                                                                                                      
library(colorspace)
library(ggplot2)
cols_phylum <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D","#666666","#00AFBB", "#E7B800", "#FC4E07","#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7","#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#ABB065")

tiff('phylum_barplot_treat', units="in", width=16, height=8, res=300)
ggplot(phylum_melt, aes(x =variable, y = value, fill=Phylum))+geom_bar(stat='identity',colour="grey",size=0.1)+labs(title="16S stacked barplot",x="Treatment",y="Relative Abundance")+scale_fill_manual(values=cols_phylum)+theme(title=(element_text(size=25,family = "Times",face="bold")),plot.margin = unit(c(0.5,0.1,0.5,0.1),"cm"),axis.text.x = element_text(angle=90,hjust = 1),text=element_text(family = "Times",face="plain",size = 25),panel.background = element_rect(fill = "grey90"),panel.grid.major = element_line(color = "white"),panel.grid.minor = element_line(color = "white"),panel.border = element_rect(colour = NA, fill=NA),axis.text=element_text(size=25,family="Times"),axis.title=element_text(size = 25,face="bold",family = "Times"),legend.title = element_text(size=25,face="bold",family = "Times"),legend.text = (element_text(size=20,family = "Times")))#pdf save:5x6
dev.off()

tiff('phylum_barplot_depth', units="in", width=16, height=8, res=300)
ggplot(depth_phylum_melt, aes(x =variable, y = value, fill=Phylum))+geom_bar(stat='identity',colour="grey",size=0.1)+labs(title="16S stacked barplot",x="Treatment",y="Relative Abundance")+scale_fill_manual(values=cols_phylum)+theme(title=(element_text(size=25,family = "Times",face="bold")),plot.margin = unit(c(0.5,0.1,0.5,0.1),"cm"),axis.text.x = element_text(angle=90,hjust = 1),text=element_text(family = "Times",face="plain",size = 25),panel.background = element_rect(fill = "grey90"),panel.grid.major = element_line(color = "white"),panel.grid.minor = element_line(color = "white"),panel.border = element_rect(colour = NA, fill=NA),axis.text=element_text(size=25,family="Times"),axis.title=element_text(size = 25,face="bold",family = "Times"),legend.title = element_text(size=25,face="bold",family = "Times"),legend.text = (element_text(size=20,family = "Times")))#pdf save:5x6
dev.off()

tiff('phylum_barplot_sample', units="in", width=16, height=8, res=300)
ggplot(sample_phylum_melt, aes(x =variable, y = value, fill=Phylum))+geom_bar(stat='identity',colour="grey",size=0.1)+labs(title="16S stacked barplot",x="Treatment",y="Relative Abundance")+scale_fill_manual(values=cols_phylum)+theme(title=(element_text(size=25,family = "Times",face="bold")),plot.margin = unit(c(0.5,0.1,0.5,0.1),"cm"),axis.text.x = element_text(angle=90,hjust = 1),text=element_text(family = "Times",face="plain",size = 25),panel.background = element_rect(fill = "grey90"),panel.grid.major = element_line(color = "white"),panel.grid.minor = element_line(color = "white"),panel.border = element_rect(colour = NA, fill=NA),axis.text=element_text(size=25,family="Times"),axis.title=element_text(size = 25,face="bold",family = "Times"),legend.title = element_text(size=25,face="bold",family = "Times"),legend.text = (element_text(size=20,family = "Times")))#pdf save:5x6
dev.off()

##class(relative abundance)
class_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[3],NArm=TRUE)
class_r_rarefied_mothur_phyloseq 

merged_r_class_rarefied_mothur_phyloseq<-merge_samples(class_r_rarefied_mothur_phyloseq,"Treat")#merge
merged_r_class_rarefied_mothur_phyloseq
depth_merged_r_class_rarefied_mothur_phyloseq<-merge_samples(class_r_rarefied_mothur_phyloseq,"Depth")
depth_merged_r_class_rarefied_mothur_phyloseq
sample_merged_r_class_rarefied_mothur_phyloseq<-merge_samples(class_r_rarefied_mothur_phyloseq,"Sample")
sample_merged_r_class_rarefied_mothur_phyloseq

tax_table(merged_r_class_rarefied_mothur_phyloseq)[,3]
r_merged_r_class_rarefied_mothur_phyloseq<-transform_sample_counts(merged_r_class_rarefied_mothur_phyloseq,function(x) x/sum(x))
r_merged_r_class_rarefied_mothur_phyloseq#by treat
prune_r_merged_r_class_rarefied_mothur_phyloseq<-prune_taxa(taxa_sums(r_merged_r_class_rarefied_mothur_phyloseq)>=0.02,r_merged_r_class_rarefied_mothur_phyloseq)
prune_r_merged_r_class_rarefied_mothur_phyloseq

tax_table(depth_merged_r_class_rarefied_mothur_phyloseq)[,3]
r_depth_merged_r_class_rarefied_mothur_phyloseq<-transform_sample_counts(depth_merged_r_class_rarefied_mothur_phyloseq,function(x) x/sum(x))
r_depth_merged_r_class_rarefied_mothur_phyloseq#by depth
prune_r_depth_merged_r_class_rarefied_mothur_phyloseq<-prune_taxa(taxa_sums(r_depth_merged_r_class_rarefied_mothur_phyloseq)>=0.02,r_depth_merged_r_class_rarefied_mothur_phyloseq)
prune_r_depth_merged_r_class_rarefied_mothur_phyloseq

r_sample_merged_r_class_rarefied_mothur_phyloseq<-transform_sample_counts(sample_merged_r_class_rarefied_mothur_phyloseq,function(x) x/sum(x))
r_sample_merged_r_class_rarefied_mothur_phyloseq#by sample
prune_r_sample_merged_r_class_rarefied_mothur_phyloseq<-prune_taxa(taxa_sums(r_sample_merged_r_class_rarefied_mothur_phyloseq)>=0.02,r_sample_merged_r_class_rarefied_mothur_phyloseq)
prune_r_sample_merged_r_class_rarefied_mothur_phyloseq

#treat
class_combine<-cbind(tax_table(r_merged_r_class_rarefied_mothur_phyloseq)[,3],otu_table(t(r_merged_r_class_rarefied_mothur_phyloseq)))
write.csv(class_combine,file='class_combine.csv')
class_combine<-read.csv('class_combine.csv',row.names = 1)
#depth
depth_class_combine<-cbind(tax_table(r_depth_merged_r_class_rarefied_mothur_phyloseq)[,3],otu_table(t(r_depth_merged_r_class_rarefied_mothur_phyloseq)))
write.csv(depth_class_combine,file='depth_class_combine.csv')
depth_class_combine<-read.csv('depth_class_combine.csv',row.names = 1)
#sample
sample_class_combine<-cbind(tax_table(r_sample_merged_r_class_rarefied_mothur_phyloseq)[,3],otu_table(t(r_sample_merged_r_class_rarefied_mothur_phyloseq)))
write.csv(sample_class_combine,file='sample_class_combine.csv')
sample_class_combine<-read.csv('sample_class_combine.csv',row.names = 1)

##then sort the "class_combine" excelfile from large to small value, subset file the species large than 0.1,new file name"class_combine_subset" 
class_combine_sub<-read.csv('class_combine.csv',row.names = 1)
class_melt<-melt(class_combine_sub,id.vars = "Class")
head(class_melt)

depth_class_combine_sub<-read.csv('depth_class_combine.csv',row.names = 1)
depth_class_melt<-melt(depth_class_combine_sub,id.vars = "Class")
head(depth_class_melt)# totally there have 7 class have been selected

sample_class_combine_sub<-read.csv('sample_class_combine.csv',row.names = 1)
sample_class_melt<-melt(sample_class_combine_sub,id.vars = "Class")
head(sample_class_melt)

library(colorspace)
library(ggplot2)
pal <- choose_palette()
library(colorspace)
library(ggplot2)
cols_phylum <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D","#666666","#00AFBB", "#E7B800", "#FC4E07","#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7","#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#ABB065")

tiff('phylum_barplot_treat', units="in", width=16, height=8, res=300)
ggplot(phylum_melt, aes(x =variable, y = value, fill=Phylum))+geom_bar(stat='identity',colour="grey",size=0.1)+labs(title="16S stacked barplot",x="Treatment",y="Relative Abundance")+scale_fill_manual(values=cols_phylum)+theme(title=(element_text(size=25,family = "Times",face="bold")),plot.margin = unit(c(0.5,0.1,0.5,0.1),"cm"),axis.text.x = element_text(angle=90,hjust = 1),text=element_text(family = "Times",face="plain",size = 25),panel.background = element_rect(fill = "grey90"),panel.grid.major = element_line(color = "white"),panel.grid.minor = element_line(color = "white"),panel.border = element_rect(colour = NA, fill=NA),axis.text=element_text(size=25,family="Times"),axis.title=element_text(size = 25,face="bold",family = "Times"),legend.title = element_text(size=25,face="bold",family = "Times"),legend.text = (element_text(size=20,family = "Times")))#pdf save:5x6
dev.off()

tiff('phylum_barplot_depth', units="in", width=16, height=8, res=300)
ggplot(depth_phylum_melt, aes(x =variable, y = value, fill=Phylum))+geom_bar(stat='identity',colour="grey",size=0.1)+labs(title="16S stacked barplot",x="Treatment",y="Relative Abundance")+scale_fill_manual(values=cols_phylum)+theme(title=(element_text(size=25,family = "Times",face="bold")),plot.margin = unit(c(0.5,0.1,0.5,0.1),"cm"),axis.text.x = element_text(angle=90,hjust = 1),text=element_text(family = "Times",face="plain",size = 25),panel.background = element_rect(fill = "grey90"),panel.grid.major = element_line(color = "white"),panel.grid.minor = element_line(color = "white"),panel.border = element_rect(colour = NA, fill=NA),axis.text=element_text(size=25,family="Times"),axis.title=element_text(size = 25,face="bold",family = "Times"),legend.title = element_text(size=25,face="bold",family = "Times"),legend.text = (element_text(size=20,family = "Times")))#pdf save:5x6
dev.off()

tiff('phylum_barplot_sample', units="in", width=16, height=8, res=300)
ggplot(sample_phylum_melt, aes(x =variable, y = value, fill=Phylum))+geom_bar(stat='identity',colour="grey",size=0.1)+labs(title="16S stacked barplot",x="Treatment",y="Relative Abundance")+scale_fill_manual(values=cols_phylum)+theme(title=(element_text(size=25,family = "Times",face="bold")),plot.margin = unit(c(0.5,0.1,0.5,0.1),"cm"),axis.text.x = element_text(angle=90,hjust = 1),text=element_text(family = "Times",face="plain",size = 25),panel.background = element_rect(fill = "grey90"),panel.grid.major = element_line(color = "white"),panel.grid.minor = element_line(color = "white"),panel.border = element_rect(colour = NA, fill=NA),axis.text=element_text(size=25,family="Times"),axis.title=element_text(size = 25,face="bold",family = "Times"),legend.title = element_text(size=25,face="bold",family = "Times"),legend.text = (element_text(size=20,family = "Times")))#pdf save:5x6
dev.off()

### alpha-diversity
library(plotly)
Diversity<-data.frame(Chaos<-estimate_richness(rarefied_mothur_phyloseq,split=TRUE,measures = NULL),sample_data(rarefied_mothur_phyloseq))
Diversity
identical(rownames(Diversity),sample_names(rarefied_mothur_phyloseq)) # if true, two objects are equal.
write.csv(Diversity,file = "alpha_diversity.csv")

##Shannon diversity
library(colorspace)
library(ggplot2)
library(agricolae)
pal <- choose_palette()
cols_diversity = pal(3)
pal <- choose_palette()
cols_diversity1=pal(4)
pal <- choose_palette()
cols_diversity2=pal(12)
tiff('Sample_shannon diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Sample,y=Shannon,fill=Sample))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Sample))+scale_fill_manual(values = cols_diversity2,name="Treatment")+labs(title="ITS Shannon",x="Treatment",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
samplelm_Shannon = lm(Diversity$Shannon~Diversity$Sample,data=Diversity)
(HSD.test(samplelm_Shannon,"Diversity$Sample"))
#$groups
#Diversity$Shannon groups
#W_0_2V              4.558566      a
#T2_0_2V             4.297728     ab
#C_0_2V              4.293205     ab
#W_2_10V             4.088126     ab
#T2W_0_2V            4.031083     ab
#T2_20_30V           3.894755     ab
#C_2_10V             3.848406     ab
#T2_2_10V            3.807944     ab
#T2W_2_10V           3.728413     ab
#T2W_20_30V          3.538769     ab
#C_20_30V            3.406319     ab
#W_20_30V            3.338620      b
tiff('Sample_InvSimpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Sample,y=InvSimpson,fill=Sample))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Sample))+scale_fill_manual(values = cols_diversity2,name="Treatment")+labs(title="ITS InvSimpson",x="Treatment",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
samplelm_InvSimpson = lm(Diversity$InvSimpson~Diversity$Sample,data=Diversity)
(HSD.test(samplelm_InvSimpson,"Diversity$Sample")) #all a, no sig


tiff('Sample_Chao1 diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Sample,y=Chao1 ,fill=Sample))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Sample))+scale_fill_manual(values = cols_diversity2,name="Treatment")+labs(title="ITS Chao1",x="Treatment",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
samplelm_Chao1 = lm(Diversity$Chao1~Diversity$Sample,data=Diversity)
(HSD.test(samplelm_Chao1,"Diversity$Sample"))
#Diversity$Chao1 groups
#W_0_2V           238.02971      a
#T2W_0_2V         180.31607     ab
#T2_0_2V          175.76062     ab
#C_0_2V           167.38830     ab
#W_2_10V          162.39073     ab
#T2_20_30V        153.10357     ab
#T2W_2_10V        137.75446      b
#T2_2_10V         118.78125      b
#C_2_10V          113.40206      b
#C_20_30V          93.53030      b
#T2W_20_30V        93.44048      b
#W_20_30V          87.17946      b

tiff('treat_shannon diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Treat,y=Shannon,fill=Treat))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Treat))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="ITS Shannon",x="Treatment",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
treatlm_shannon = lm(Diversity$Shannon~Diversity$Treat,data=Diversity)
anova(treatlm_shannon)
#Analysis of Variance Table

#Response: Diversity$Shannon
#Df Sum Sq Mean Sq F value Pr(>F)
#Diversity$Treat  3  0.463 0.15433  0.5469 0.6531
#Residuals       41 11.571 0.28221    

#Tukey and other multiple comparison tests can be performed with a handful of functions.  The functions TukeyHSD, HSD.test, and LSD.test are probably not appropriate for cases where there are unbalanced data or unequal variances among levels of the factor, though TukeyHSD does make an adjustment for mildly unbalanced data.  It is my understanding that the multcomp and lsmeans packages are more appropriate for unbalanced data.  Another alternative is the DTK package that performs mean separation tests on data with unequal sample sizes and no assumption of equal variances.
install.packages('agricolae', repos="http://cran.rstudio.com/")#if there have something goes wrong remeber try to install the packrat package for dependency
library(agricolae)
library(vegan)
#install.packages("packrat")
#library(packrat)
(HSD.test(treatlm_shannon,"Diversity$Treat"))#here we use Turkey test

#$statistics
#MSerror Df     Mean       CV
#0.28221 41 3.914042 13.57253

#$parameters
#test          name.t ntr StudentizedRange alpha
#Tukey Diversity$Treat   4         3.786726  0.05

#$means
#Diversity$Shannon       std  r      Min      Max      Q25      Q50      Q75
#Control          3.889582 0.3984774 11 3.236011 4.506328 3.601657 3.840392 4.186500
#T2               4.021220 0.4004932 10 3.474255 4.556617 3.650333 4.086982 4.246662
#T2W              3.766088 0.6578211 12 2.499788 4.764983 3.374425 3.776249 4.286304
#W                3.995104 0.5861429 12 2.769083 4.705526 3.671944 4.092737 4.401383

#$comparison
#NULL

#$groups
#Diversity$Shannon groups
#T2               4.021220      a
#W                3.995104      a
#Control          3.889582      a
#T2W              3.766088      a

#attr(,"class")
#[1] "group"


tiff('depth_shannon diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Depth,y=Shannon,fill=Depth))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Depth))+scale_fill_manual(values = cols_diversity,name="Depth")+labs(title="ITS Shannon",x="Depth",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
#tillage_shannon_ttest <- t.test(Diversity$Shannon ~ Diversity$Tillage, data = Diversity)
#tillage_shannon_ttest#p-value = 0.01317
depthlm_shannon = lm(Diversity$Shannon~Diversity$Depth,data=Diversity)
(HSD.test(depthlm_shannon,"Diversity$Depth"))#here we use Turkey test, significant!
#$groups
#Diversity$Shannon groups
#V0_2            4.295146      a
#V2_10           3.868222      b
#V20_30          3.501386      b

#attr(,"class")
#[1] "group"

##Simpson diversity
tiff('treat_Simpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Treat,y=Simpson,fill=Treat))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Treat))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="ITS Simpson",x="Treatment",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
treatlm_Simpson = lm(Diversity$Simpson~Diversity$Treat,data=Diversity)
(HSD.test(treatlm_Simpson,"Diversity$Treat"))

tiff('depth_Simpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Depth,y=Simpson,fill=Depth))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Depth))+scale_fill_manual(values = cols_diversity,name="Depth")+labs(title="ITS Simpson",x="Depth",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
depthlm_Simpson = lm(Diversity$Simpson~Diversity$Depth,data=Diversity)
(HSD.test(depthlm_Simpson,"Diversity$Depth"))


tiff('Sample_Simpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Sample,y=Simpson,fill=Sample))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Sample))+scale_fill_manual(values = cols_diversity2,name="Treatment")+labs(title="ITS Simpson",x="Treatment",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
samplelm_Simpson = lm(Diversity$Simpson~Diversity$Sample,data=Diversity)
(HSD.test(samplelm_Simpson,"Diversity$Sample")) #all a, no sig



##InvSimpson diversity
tiff('treat_InvSimpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Treat,y=InvSimpson,fill=Treat))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Treat))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="ITS InvSimpson",x="Treatment",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
treatlm_InvSimpson = lm(Diversity$InvSimpson~Diversity$Treat,data=Diversity)
(HSD.test(treatlm_InvSimpson,"Diversity$Treat"))#no significant differences.

tiff('depth_InvSimpson diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Depth,y=InvSimpson,fill=Depth))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Depth))+scale_fill_manual(values = cols_diversity,name="Depth")+labs(title="ITS InvSimpson",x="Depth",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
depthlm_InvSimpson = lm(Diversity$InvSimpson~Diversity$Depth,data=Diversity)
(HSD.test(depthlm_InvSimpson,"Diversity$Depth"))#significant differences
#$groups
#Diversity$InvSimpson groups
#V0_2               45.78293      a
#V2_10              27.64449      b
#V20_30             19.18709      b

#attr(,"class")
#[1] "group"

##Chao1 diversity
tiff('treat_Chao1 diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Treat,y=Chao1 ,fill=Treat))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Treat))+scale_fill_manual(values = cols_diversity1,name="Treatment")+labs(title="ITS Chao1",x="Treatment",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

treatlm_Chao1 = lm(Diversity$Chao1~Diversity$Treat,data=Diversity)
(HSD.test(treatlm_Chao1,"Diversity$Treat")) #no sig

tiff('depth_Chao1 diversity.tiff', units="in", width=7, height=5, res=300)
ggplot(Diversity,aes(x=Depth,y=Chao1,fill=Depth))+geom_boxplot(size=0.5)+geom_jitter(shape=16, position=position_jitter(0),aes(fill=Depth))+scale_fill_manual(values = cols_diversity,name="Depth")+labs(title="ITS Chao1",x="Depth",y="Index")+theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=17,family = "Arial",face="bold",vjust = 3,hjust=0.5)),text=element_text(family = "Arial",face="bold"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text.y=element_text(size=12,family="Arial"),axis.text.x = element_text(family="Arial",size=12,angle=60,hjust = 1),axis.title=element_text(size = 15,face="bold",family = "Arial",vjust = 1),legend.title = element_text(size=13,face="bold",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")),plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
depthlm_Chao1 = lm(Diversity$Chao1~Diversity$Depth,data=Diversity)
(HSD.test(depthlm_Chao1,"Diversity$Depth"))#significant differences
#$groups
#Diversity$Chao1 groups
#V0_2          190.3737      a
#V2_10         133.0821      b
#V20_30        100.7137      b

###PCoA_r_rarefied_mothur_phyloseq at otu level
## PCoA
library(vegan)
library(ggplot2)
r_rarefied_mothur_phyloseq_bray<-vegdist(t(otu_table(r_rarefied_mothur_phyloseq)),method="bray",binary = FALSE)
r_rarefied_mothur_phyloseq_PCoA<-ordinate(r_rarefied_mothur_phyloseq,method = "PCoA", r_rarefied_mothur_phyloseq_bray)
#treat and depth
tiff('PCoA depthandtreat.tiff', units="in", width=7, height=5, res=300)
#plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="samples",color="block",label="Sample")+geom_point(size=4)+ggtitle("Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c","gold","aquamarine"))+labs(x="PCoA1 [25.5%]",y="PCoA2 [13.3%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="samples",color="Depth",shape = "Treat")+geom_point(size=4)+ggtitle("Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c","gold","aquamarine"))+labs(x="PCoA1 [19.7%]",y="PCoA2 [8.2%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=13,family = "Arial")))
dev.off()


tiff('PCoA treatanddepth.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq,r_rarefied_mothur_phyloseq_PCoA,type="samples",color="Treat", shape = "Depth")+geom_point(size=4)+ggtitle("Principal Coordinates Analysis")+scale_color_manual(values = c("#0A0A0A","#f9830c","gold","aquamarine"))+labs(x="PCoA1 [19.7%]",y="PCoA2 [8.2%]")+theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))
dev.off()

###pca ggbiplot
install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
##pca of phylum
write.csv(tax_table(r_rarefied_mothur_phyloseq),file = "taxa_table.csv")
phylum <- t(otu_table(phylum_r_rarefied_mothur_phyloseq))
write.csv(phylum,file = "phylum.csv")
write.csv(sample_data(phylum_r_rarefied_mothur_phyloseq),file = "sample_data.csv")
phylum.rev <- read.csv(file = "phylum_rev.csv", sep = ",", row.names = 1)#phylum_rev is the file been revised accroding to the simper results, we selected first several phylum which cumulation explaination are over 70%
pca.phylum <- prcomp(as.data.frame(phylum.rev),center = TRUE, scale. = TRUE)
summary(pca.phylum)
str(pca.phylum)
phylum.group_treat <- c(rep("W",3),rep("C", 2),rep("T2W",6),rep("T2",3),rep("W",3),rep("T2",3),rep("W",3),rep("C", 3),rep("T2",2),rep("T2W",3),rep("W",3),rep("C", 3),rep("T2",2),rep("T2W",3),rep("C", 3))
tiff('pca.phylum_treat.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.phylum,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=phylum.group_treat, labels = row.names(phylum.rev),var.axes=TRUE,varname.size = 2.5,varname.adjust =1.5
         ,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue","gold"))+
  ggtitle("PCA of phylum")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-5, 5)
dev.off()
#ggbiplot(pca.phylum,ellipse=TRUE,labels.size = 2, circle=FALSE,obs.scale = 1, var.scale = 1,groups=phylum.group_tillage, labels = row.names(phylum),var.axes=TRUE,varname.size = 2.5,varname.adjust =1.8
#         ,varname.abbrev = FALSE,arrow.color="blue")+
#  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
#  scale_colour_manual(name="Treatment", values= c("forest green", "goldenrod3", "dark blue"))+
#  ggtitle("PCA of phylum")+
#  theme_minimal()+
#  theme(legend.position = "bottom")+xlim(-6, 6) + ylim(-4, 4)
#dev.off()
#cap_plot + geom_segment( mapping = arrow_map, size = .5,data = arrowdf,color = "gray",arrow = arrowhead )+geom_text(mapping = label_map, size = 4,  data = arrowdf,  show.legend = FALSE)

phylum.group_depth <- c("V0_2","V2_10","V20_30","V0_2","V2_10","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30")
tiff('pca.phylum_depth.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.phylum,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=phylum.group_depth, labels = row.names(phylum.rev),var.axes=TRUE,varname.size = 2.5,varname.adjust =1.5
         ,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue"))+
  ggtitle("PCA of phylum")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-5, 5)
dev.off()

##pca of class
class <- t(otu_table(class_r_rarefied_mothur_phyloseq))
write.csv(class,file = "class.csv")
class.rev <- read.csv(file = "class_rev.csv", sep = ",", row.names = 1)#phylum_rev is the file been revised accroding to the simper results, we selected first several phylum which cumulation explaination are over 70%
pca.class <- prcomp(as.data.frame(class.rev),center = TRUE, scale. = TRUE)
summary(pca.class)
str(pca.class)
class.group_treat <- c(rep("W",3),rep("C", 2),rep("T2W",6),rep("T2",3),rep("W",3),rep("T2",3),rep("W",3),rep("C", 3),rep("T2",2),rep("T2W",3),rep("W",3),rep("C", 3),rep("T2",2),rep("T2W",3),rep("C", 3))
tiff('pca.class_treat.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.class,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=class.group_treat, labels = row.names(class),var.axes=TRUE,varname.size = 2.2,varname.adjust =1.2,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue","gold"))+
  ggtitle("PCA of class")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-4, 4)
dev.off()

class.group_depth <- c("V0_2","V2_10","V20_30","V0_2","V2_10","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30")
tiff('pca.class_depth.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.class,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=class.group_depth, labels = row.names(class),var.axes=TRUE,varname.size = 2.2,varname.adjust =1.2,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Depth", values= c("forest green", "coral", "dark blue"))+
  ggtitle("PCA of class")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-4, 4)
dev.off()

###pca of genus
genus_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[6],NArm=TRUE)
genus_r_rarefied_mothur_phyloseq

genus <- t(otu_table(genus_r_rarefied_mothur_phyloseq))
write.csv(genus,file = "genus.csv")
genus.rev <- read.csv(file = "genus_rev.csv", sep = ",", row.names = 1)#phylum_rev is the file been revised accroding to the simper results, we selected first several phylum which cumulation explaination are over 70%
pca.genus <- prcomp(as.data.frame(genus.rev),center = TRUE, scale. = TRUE)
summary(pca.genus)
str(pca.genus)
genus.group_treat <- c(rep("W",3),rep("C", 2),rep("T2W",6),rep("T2",3),rep("W",3),rep("T2",3),rep("W",3),rep("C", 3),rep("T2",2),rep("T2W",3),rep("W",3),rep("C", 3),rep("T2",2),rep("T2W",3),rep("C", 3))
tiff('pca.genus_treat.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.genus,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=genus.group_treat, labels = row.names(genus),var.axes=TRUE,varname.size = 2.2,varname.adjust =1.2,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue","gold"))+
  ggtitle("PCA of genus")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-4, 4)
dev.off()

genus.group_depth <- c("V0_2","V2_10","V20_30","V0_2","V2_10","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30","V0_2","V2_10","V0_2","V2_10","V20_30","V0_2","V2_10","V20_30")
tiff('pca.genus_depth.tiff', units="in", width=7, height=5, res=300)
ggbiplot(pca.genus,ellipse=TRUE,labels.size = 1.5, circle=FALSE,obs.scale = 1, var.scale = 1,groups=genus.group_depth, labels = row.names(genus),var.axes=TRUE,varname.size = 2.2,varname.adjust =1.2,varname.abbrev = FALSE,arrow.color="blue")+
  theme(panel.background = element_rect(fill = "white",color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=(element_text(size=15,family = "Arial",face="plain",hjust = 0.5,vjust = 3)),text=element_text(family = "Arial",face="plain"),panel.border = element_rect(colour = "black", fill=NA, size=0.5),axis.text=element_text(size=13,family="Arial"),axis.title=element_text(size = 15,face="plain",family = "Arial"),legend.title = element_text(size=13,face="plain",family = "Arial"),legend.text = (element_text(size=10,family = "Arial")))+
  scale_colour_manual(name="Treatment", values= c("forest green", "coral", "dark blue"))+
  ggtitle("PCA of genus")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-5, 5) + ylim(-4, 4)
dev.off()


###NMDS
library(dplyr)
library(scales)
r_rarefied_mothur_phyloseq 
set.seed(12345)
r_rarefied_mothur_phyloseq_ordinate <- ordinate(r_rarefied_mothur_phyloseq,"NMDS", "bray")
dev.off()
tiff('NMDS_treat.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq, r_rarefied_mothur_phyloseq_ordinate, type="samples", color="Treat", title="ITS_treat_nmds")
dev.off()
tiff('NMDS_depth.tiff', units="in", width=7, height=5, res=300)
plot_ordination(r_rarefied_mothur_phyloseq, r_rarefied_mothur_phyloseq_ordinate, type="samples", color="Depth", title="ITS_depth_nmds")
dev.off()

#betadisper
library(vegan)
group_depth <- c(rep("V0_2",1),rep("V2_10", 1),rep("V20_30",1),rep("V0_2",1),rep("V2_10", 1),
                 rep("V0_2",1),rep("V2_10", 1),rep("V20_30",1),rep("V0_2",1),rep("V2_10", 1),
                 rep("V20_30",1),rep("V0_2",1),rep("V2_10", 1),rep("V20_30",1),rep("V0_2",1),
                 rep("V2_10", 1),rep("V20_30",1),rep("V0_2",1),rep("V2_10", 1),rep("V20_30",1),
                 rep("V0_2",1),rep("V2_10", 1),rep("V20_30",1),rep("V0_2",1),rep("V2_10", 1),
                 rep("V20_30",1),rep("V0_2",1),rep("V2_10", 1),rep("V0_2",1),rep("V2_10", 1),
                 rep("V20_30",1),rep("V0_2",1),rep("V2_10", 1),rep("V20_30",1),rep("V0_2",1),
                 rep("V2_10", 1),rep("V20_30",1),rep("V0_2",1),rep("V2_10", 1),rep("V0_2",1),
                 rep("V2_10", 1),rep("V20_30",1),rep("V0_2",1),rep("V2_10", 1),rep("V20_30",1))
group_treat <- c(rep("W",3),rep("C", 2),rep("T2W",6),rep("T2", 3),rep("W",3),rep("T2", 3),rep("W",3),rep("C", 3),rep("T2", 2),rep("T2W",3),rep("W",3),rep("C", 3),rep("T2", 2),rep("T2W", 3),rep("C", 3))

mod_depth <- betadisper(r_rarefied_mothur_phyloseq_bray,group_depth,bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
tiff('depth_disper.tiff', units="in", width=12, height=8, res=300)
plot(mod_depth,hull = FALSE, ellipse = TRUE,cex=2)
dev.off()
anova(mod_depth)
permutest(mod_depth)
#Response: Distances

#Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)  
#Groups     2 0.06191 0.030953 2.989    999  0.049 *
#  Residuals 42 0.43494 0.010356                      

mod_treat <- betadisper(r_rarefied_mothur_phyloseq_bray,group_treat,bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
tiff('treat_disper.tiff', units="in", width=12, height=8, res=300)
plot(mod_treat,hull = FALSE, ellipse = TRUE,cex=2)
dev.off()
anova(mod_treat)
permutest(mod_treat)#no sig