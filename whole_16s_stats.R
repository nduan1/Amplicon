adonis(t(otu_table(r_rarefied_mothur_phyloseq))~block+Season+Cover*Nitrogen*Tillage,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=999,by="margin")#at least 1000
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~block,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=999,by="margin")#at least 1000
#adonis(t(otu_table(prune_r_rmsmay_rarefied_mothur_phyloseq))~block+Season+Cover*Nitrogen*Tillage,data=data.frame(sample_data(prune_r_rmsmay_rarefied_mothur_phyloseq)),permutations=999,by="margin")#at least 1000
adonis(t(otu_table(prune_r_rmsmay_rarefied_mothur_phyloseq))~Season+Cover*Nitrogen*Tillage,data=data.frame(sample_data(prune_r_rmsmay_rarefied_mothur_phyloseq)),strata=block,permutations=999,by="margin")#at least 1000
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~Season,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=999,by="margin")#at least 1000

#tillage
library(permute)
set.seed(1013)
perm<-how(nperm = 999,plots = Plots(strata=paste(data.frame(sample_data(r_rarefied_mothur_phyloseq))$Nitrogen,data.frame(sample_data(r_rarefied_mothur_phyloseq))$Cover,sep="_" )),within=Within(type = "free"))
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~Tillage,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=perm,by="margin")#at least 1000

#cover
library(permute)
set.seed(1013)
perm<-how(nperm = 999,plots = Plots(strata=paste(data.frame(sample_data(r_rarefied_mothur_phyloseq))$Tillage,data.frame(sample_data(r_rarefied_mothur_phyloseq))$Cover,sep="_" )),within=Within(type = "free"))
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~Cover,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=perm,by="margin")#at least 1000

#Nitrogen
library(permute)
set.seed(1013)
perm<-how(nperm = 999,plots = Plots(strata=paste(data.frame(sample_data(r_rarefied_mothur_phyloseq))$Tillage,data.frame(sample_data(r_rarefied_mothur_phyloseq))$Cover,sep="_" )),within=Within(type = "free"))
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~Nitrogen,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=perm,by="margin")#at least 1000
