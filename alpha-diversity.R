###beta-diversity-hill numbers####
#ref:https://github.com/daijiang/hillR
library(vegan)
library(hillR)
library(phyloseq)
library(tidyverse)
new.phyloseq_up.no
#rows are sites and 
bd_data <- as.data.frame(t(otu_table(new.phyloseq_up.no)))
rownames(bd_data)
nccn0_bd_data <- bd_data[1:4,]
nccn0_beta <- hill_taxa_parti(nccn0_bd_data,q=0)
nccn60_bd_data <- bd_data[5:8,]
nccn60_beta <- hill_taxa_parti(nccn60_bd_data,q=0)
ncntn0_bd_data <- bd_data[9:12,]
ncntn0_beta <- hill_taxa_parti(ncntn0_bd_data,q=0)
ncntn60_bd_data <- bd_data[13:16,]
ncntn60_beta <- hill_taxa_parti(ncntn60_bd_data,q=0)
vcn60_bd_data <- bd_data[17:20,]
vcn60_beta <- hill_taxa_parti(vcn60_bd_data,q=0)
vntn0_bd_data <- bd_data[21:24,]
vntn0_beta <- hill_taxa_parti(vntn0_bd_data,q=0)
vntn60_bd_data <- bd_data[25:27,]
vntn60_beta <- hill_taxa_parti(vntn60_bd_data,q=0)
wcn0_bd_data <- bd_data[28:31,]
wcn0_beta <- hill_taxa_parti(wcn0_bd_data,q=0)
wcn60_bd_data <- bd_data[32:35,]
wcn60_beta <- hill_taxa_parti(wcn60_bd_data,q=0)
wntn0_bd_data <- bd_data[36:39,]
wntn0_beta <- hill_taxa_parti(wntn0_bd_data,q=0)
wntn60_bd_data <- bd_data[40:43,]
wntn60_beta <- hill_taxa_parti(wntn60_bd_data,q=0)
beta <- rbind(nccn0_beta$TD_beta,
            nccn60_beta$TD_beta,
            ncntn0_beta$TD_beta,
            ncntn60_beta$TD_beta,
            vcn60_beta$TD_beta,
            vntn0_beta$TD_beta,
            vntn60_beta$TD_beta,
            wcn0_beta$TD_beta,
            wcn60_beta$TD_beta,
            wntn0_beta$TD_beta,
            wntn60_beta$TD_beta)
row.names(beta) <- c("NCCN0","NCCN60","NCNTN0","NCNTN60","VCN60","VNTN0","VNTN60","WCN0",
                     "WCN60","WNTN0","WNTN60")

colnames(beta) <- "beta_diversity"
head(beta)
beta_1 <-mutate(as.data.frame(beta),Treat=rownames(as.data.frame(beta)))
head(beta_1 )
dong_data <- read.csv("Diversity.updated.dong.csv",sep = ",")
head(dong_data)
update_beta <- left_join(dong_data,beta_1,by="Treat")
head(update_beta)
write.csv(update_beta,"/Users/Ning/Desktop/R/New project/1stnifa_16snewreference/dataset/with_beta_diversity.csv")
