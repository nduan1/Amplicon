#install.packages("indicspecies")
library(indicspecies)
library(phyloseq)
library(tidyverse)
#setwd("/Users/Ning/Desktop/R/New project/greenland_16s/")
#setwd("/Users/Ning/Desktop/R/New project/greenland_16s/sparcc_analysis/")
#at otu_level
sample_sums(rarefied_mothur_phyloseq)
filter_r_rarefied_mothur_phyloseq<-prune_taxa(taxa_sums(r_rarefied_mothur_phyloseq)>=0.02,r_rarefied_mothur_phyloseq)
r_rarefied_mothur_phyloseq_prune <- filter_taxa(filter_r_rarefied_mothur_phyloseq,function(x) sum(x > 0) > (4/48*length(x)), TRUE)#Remove taxa not seen more than 0 times in at least 4/48 of the samples. 
rarefied_mothur_phyloseq_prune<-transform_sample_counts(r_rarefied_mothur_phyloseq_prune,function(x) x*3262)
head(otu_table(rarefied_mothur_phyloseq_prune))
its.phyloseq.v0_2 <- subset_samples(rarefied_mothur_phyloseq_prune, Depth=="V0_2")
its.phyloseq.v0_2 
its.phyloseq.v2_10 <- subset_samples(rarefied_mothur_phyloseq_prune, Depth=="V2_10")
its.phyloseq.v2_10
its.phyloseq.v20_30 <- subset_samples(rarefied_mothur_phyloseq_prune, Depth=="V20_30")
its.phyloseq.v20_30
#0_2
its.otu.0_2.indi.data <- as.data.frame(t(otu_table(its.phyloseq.v0_2)))
head(its.otu.0_2.indi.data[1:5])
row.names(its.otu.0_2.indi.data)
nrow(its.otu.0_2.indi.data)
ncol(its.otu.0_2.indi.data)
its.otu.0_2.indi<- as.data.frame(mutate(its.otu.0_2.indi.data,
                                        treat =c("W","C","T2W","T2W","T2","W","T2","W",
                                                 "C","T2","T2W","W","C","T2","T2W","C")))
colnames(its.otu.0_2.indi.data)
treat <- its.otu.0_2.indi$treat
its.otu.0_2.indicator <-  multipatt(its.otu.0_2.indi.data, its.otu.0_2.indi$treat, func = "IndVal.g",duleg = TRUE,control = how(nperm=1000))
summary(its.otu.0_2.indicator,indvalcomp=TRUE)
names(its.otu.0_2.indicator)
# [1] "call"    "func"    "cluster" "comb"    "str"    
# [6] "A"       "B"       "sign"
# its.otu.0_2.indicator$sign
# its.otu.0_2.indicator$cluster
# its.otu.0_2.indicator$comb
# its.otu.0_2.indicator$str
its.otu.0_2.indicator.dat <- as.matrix(summary(its.otu.0_2.indicator,indvalcomp=TRUE))

#2_10
its.otu.2_10.indi.data <- as.data.frame(t(otu_table(its.phyloseq.v2_10)))
head(its.otu.2_10.indi.data[1:5])
row.names(its.otu.2_10.indi.data)
nrow(its.otu.2_10.indi.data)
ncol(its.otu.2_10.indi.data)
its.otu.2_10.indi<- as.data.frame(mutate(its.otu.2_10.indi.data,
                                         treat =c("W","C","T2W","T2W","T2","W","T2","W",
                                                  "C","T2","T2W","W","C","T2","T2W","C")))
colnames(its.otu.2_10.indi.data)
treat <- its.otu.2_10.indi$treat
its.otu.2_10.indicator <-  multipatt(its.otu.2_10.indi.data, its.otu.2_10.indi$treat, func = "IndVal.g",duleg = TRUE, control = how(nperm=1000))
summary(its.otu.2_10.indicator,indvalcomp=TRUE)

#20_30
its.otu.20_30.indi.data <- as.data.frame(t(otu_table(its.phyloseq.v20_30)))[-c(4,6),]
head(its.otu.20_30.indi.data[1:5])
row.names(its.otu.20_30.indi.data)
nrow(its.otu.20_30.indi.data)
ncol(its.otu.20_30.indi.data)
its.otu.20_30.indi<- as.data.frame(mutate(its.otu.20_30.indi.data,
                                          treat =c("W","T2W","T2W","W","W",
                                                   "C","T2W","W","C","T2W","C")))
colnames(its.otu.20_30.indi.data)
treat <- its.otu.20_30.indi$treat
its.otu.20_30.indicator <-  multipatt(its.otu.20_30.indi.data, its.otu.20_30.indi$treat, func = "IndVal.g",duleg = TRUE, control = how(nperm=1000))
summary(its.otu.20_30.indicator,indvalcomp=TRUE)
