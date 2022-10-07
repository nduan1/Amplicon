r_rarefied_mothur_phyloseq
phylum_r_rarefied_mothur_phyloseq#35 taxa
#class(relative abundance)
class_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[3],NArm=TRUE)
class_r_rarefied_mothur_phyloseq#116 
##order(relative abundance)
order_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[4],NArm=TRUE)
order_r_rarefied_mothur_phyloseq#298
##family(relative abundance)
family_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[5],NArm=TRUE)
family_r_rarefied_mothur_phyloseq#501
##genus(relative abundance)
genus_r_rarefied_mothur_phyloseq<-tax_glom(r_rarefied_mothur_phyloseq,taxrank = rank_names(r_rarefied_mothur_phyloseq)[6],NArm=TRUE)
genus_r_rarefied_mothur_phyloseq

###PERMANOVA
## at OTU level - cover crop impact
library(permute)
library(vegan)
set.seed(1013)
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~block+Nitrogen*Cover*Tillage,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=9999,by="margin")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#block                   3    0.3371 0.11238  1.2330 0.07204 0.0634 .  
#Nitrogen                1    0.4244 0.42443  4.6568 0.09069 0.0001 ***
 # Cover                   2    0.2630 0.13150  1.4428 0.05620 0.0180 *  
#  Tillage                 1    0.3408 0.34076  3.7388 0.07281 0.0001 ***
#  Nitrogen:Cover          2    0.2299 0.11497  1.2615 0.04913 0.0750 .  
#Nitrogen:Tillage        1    0.0945 0.09452  1.0371 0.02020 0.3522    
#Cover:Tillage           2    0.1810 0.09051  0.9930 0.03868 0.4450    
#Nitrogen:Cover:Tillage  2    0.1662 0.08309  0.9116 0.03551 0.6717    
#Residuals              29    2.6431 0.09114         0.56475           
#Total                  43    4.6801                 1.00000           

#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~Nitrogen*Cover*Tillage,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=9999,by="margin")

perm<-how(nperm = 9999,plots = Plots(strata=paste(data.frame(sample_data(r_rarefied_mothur_phyloseq))$Nitrogen,data.frame(sample_data(r_rarefied_mothur_phyloseq))$Tillage,sep="_" )),within=Within(type = "free"))
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~Cover,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=perm,by="margin")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Cover      2    0.2768 0.13842  1.2889 0.05915 0.0023 **
 # Residuals 41    4.4033 0.10740         0.94085          
#Total     43    4.6801                 1.00000          

#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 

## at OTU level - tillage impact
set.seed(1013)
perm<-how(nperm = 9999,plots = Plots(strata=paste(data.frame(sample_data(r_rarefied_mothur_phyloseq))$Nitrogen,data.frame(sample_data(r_rarefied_mothur_phyloseq))$Cover,sep="_" )),within=Within(type = "free"))
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~Tillage,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=perm,by="margin")#at least 1000
#Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
#Tillage    1    0.3365 0.33647  3.1798 0.07197 0.0004998 ***
#Residuals 41    4.3384 0.10582         0.92803              
#Total     42    4.6749                 1.00000 

## at OTU level - nitrogen impact
set.seed(1013)
perm<-how(nperm = 9999,plots = Plots(strata=paste(data.frame(sample_data(r_rarefied_mothur_phyloseq))$Tillage,data.frame(sample_data(r_rarefied_mothur_phyloseq))$Cover,sep="_" )),within=Within(type = "free"))
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~Nitrogen,data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),permutations=perm,by="margin")
#Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
#Nitrogen   1    0.4382 0.43816  4.2402 0.09373 0.0004998 ***
#Residuals 41    4.2367 0.10334         0.90627              
#Total     42    4.6749                 1.00000 

## at OTU level-  block and cover nitrogen and tillage
set.seed(1013)
adonis(t(otu_table(r_rarefied_mothur_phyloseq))~block+Cover*Nitrogen*Tillage,
       data=data.frame(sample_data(r_rarefied_mothur_phyloseq)),
       permutations=9999,by="margin")
#Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
#block                   3    0.3321 0.11069  1.1876 0.07103 0.1034483    
#Cover                   2    0.2928 0.14638  1.5706 0.06262 0.0079960 ** 
#Nitrogen                1    0.4125 0.41245  4.4253 0.08823 0.0004998 ***
#Tillage                 1    0.3274 0.32744  3.5132 0.07004 0.0004998 ***
#Cover:Nitrogen          2    0.2413 0.12065  1.2945 0.05162 0.0644678 .  
#Cover:Tillage           2    0.1866 0.09332  1.0013 0.03993 0.4347826    
#Nitrogen:Tillage        1    0.0986 0.09861  1.0580 0.02109 0.3213393    
#Cover:Nitrogen:Tillage  1    0.0807 0.08072  0.8661 0.01727 0.6871564    
#Residuals              29    2.7029 0.09320         0.57817              
#Total                  42    4.6749                 1.00000 

#at phylum level
adonis(t(otu_table(phylum_r_rarefied_mothur_phyloseq))~block+Cover*Nitrogen*Tillage,
       data=data.frame(sample_data(phylum_r_rarefied_mothur_phyloseq)),
       permutations=2000,by="margin")#at least 1000
#Df SumsOfSqs   MeanSqs F.Model      R2   Pr(>F)   
#block                   3  0.010709 0.0035696  0.9282 0.05232 0.522739   
#Cover                   2  0.020507 0.0102536  2.6662 0.10019 0.014493 * 
#Nitrogen                1  0.018905 0.0189054  4.9159 0.09236 0.003998 **
#Tillage                 1  0.018491 0.0184911  4.8082 0.09034 0.003498 **
#Cover:Nitrogen          2  0.010254 0.0051272  1.3332 0.05010 0.225387   
#Cover:Tillage           2  0.009532 0.0047662  1.2393 0.04657 0.264868   
#Nitrogen:Tillage        1  0.002027 0.0020271  0.5271 0.00990 0.718641   
#Cover:Nitrogen:Tillage  1  0.002737 0.0027366  0.7116 0.01337 0.588706   
#Residuals              29  0.111527 0.0038458         0.54486            
#Total                  42  0.204690                   1.00000

#at class level
adonis(t(otu_table(class_r_rarefied_mothur_phyloseq))~block+Cover*Nitrogen*Tillage,data=data.frame(sample_data(class_r_rarefied_mothur_phyloseq)),permutations=2000,by="margin")
#Df SumsOfSqs  MeanSqs F.Model      R2    Pr(>F)    
#block                   3   0.02828 0.009425  1.2839 0.07063 0.2013993    
#Cover                   2   0.03186 0.015930  2.1701 0.07958 0.0214893 *  
#Nitrogen                1   0.05218 0.052179  7.1080 0.13033 0.0004998 ***
#Tillage                 1   0.03141 0.031410  4.2788 0.07846 0.0019990 ** 
#Cover:Nitrogen          2   0.01846 0.009228  1.2571 0.04610 0.2453773    
#Cover:Tillage           2   0.01825 0.009126  1.2432 0.04559 0.2368816    
#Nitrogen:Tillage        1   0.00271 0.002706  0.3687 0.00676 0.9060470    
#Cover:Nitrogen:Tillage  1   0.00432 0.004321  0.5886 0.01079 0.7486257    
#Residuals              29   0.21289 0.007341         0.53176              
#Total                  42   0.40035                  1.00000  

#at class level
adonis(t(otu_table(class_r_rarefied_mothur_phyloseq))~block+Cover*Nitrogen*Tillage,data=data.frame(sample_data(class_r_rarefied_mothur_phyloseq)),permutations=2000,by="margin")
#Df SumsOfSqs  MeanSqs F.Model      R2    Pr(>F)    
#block                   3   0.02828 0.009425  1.2839 0.07063 0.2048976    
#Cover                   2   0.03186 0.015930  2.1701 0.07958 0.0159920 *  
#Nitrogen                1   0.05218 0.052179  7.1080 0.13033 0.0004998 ***
#Tillage                 1   0.03141 0.031410  4.2788 0.07846 0.0009995 ***
#Cover:Nitrogen          2   0.01846 0.009228  1.2571 0.04610 0.2498751    
#Cover:Tillage           2   0.01825 0.009126  1.2432 0.04559 0.2558721    
#Nitrogen:Tillage        1   0.00271 0.002706  0.3687 0.00676 0.9165417    
#Cover:Nitrogen:Tillage  1   0.00432 0.004321  0.5886 0.01079 0.7446277    
#Residuals              29   0.21289 0.007341         0.53176              
#Total                  42   0.40035                  1.00000 

#at order level
adonis(t(otu_table(order_r_rarefied_mothur_phyloseq))~block+Cover*Nitrogen*Tillage,data=data.frame(sample_data(order_r_rarefied_mothur_phyloseq)),permutations=2000,by="margin")
#Df SumsOfSqs  MeanSqs F.Model      R2    Pr(>F)    
#block                   3   0.04283 0.014277  1.2345 0.06754 0.2173913    
#Cover                   2   0.05073 0.025366  2.1933 0.07999 0.0074963 ** 
#Nitrogen                1   0.07707 0.077072  6.6641 0.12153 0.0004998 ***
#Tillage                 1   0.05726 0.057259  4.9510 0.09029 0.0004998 ***
#Cover:Nitrogen          2   0.02789 0.013945  1.2058 0.04398 0.2588706    
#Cover:Tillage           2   0.02916 0.014580  1.2607 0.04598 0.2213893    
#Nitrogen:Tillage        1   0.00631 0.006307  0.5454 0.00995 0.8300850    
#Cover:Nitrogen:Tillage  1   0.00756 0.007557  0.6534 0.01192 0.7426287    
#Residuals              29   0.33539 0.011565         0.52884              
#Total                  42   0.63420                  1.00000 

#at family level
adonis(t(otu_table(family_r_rarefied_mothur_phyloseq))~block+Cover*Nitrogen*Tillage,data=data.frame(sample_data(family_r_rarefied_mothur_phyloseq)),permutations=2000,by="margin")
#Df SumsOfSqs  MeanSqs F.Model      R2    Pr(>F)    
#block                   3   0.05597 0.018658  1.2356 0.06833 0.1859070    
#Cover                   2   0.06057 0.030286  2.0057 0.07395 0.0089955 ** 
#Nitrogen                1   0.09431 0.094312  6.2459 0.11514 0.0004998 ***
#Tillage                 1   0.07671 0.076706  5.0799 0.09364 0.0004998 ***
#Cover:Nitrogen          2   0.03739 0.018697  1.2382 0.04565 0.2048976    
#Cover:Tillage           2   0.03511 0.017556  1.1627 0.04287 0.2768616    
#Nitrogen:Tillage        1   0.01037 0.010367  0.6866 0.01266 0.7346327    
#Cover:Nitrogen:Tillage  1   0.01078 0.010781  0.7140 0.01316 0.6841579    
#Residuals              29   0.43790 0.015100         0.53460              
#Total                  42   0.81912                  1.00000

#at genera level
adonis(t(otu_table(genus_r_rarefied_mothur_phyloseq))~block+Cover*Nitrogen*Tillage,data=data.frame(sample_data(genus_r_rarefied_mothur_phyloseq)),permutations=2000,by="margin")
#Df SumsOfSqs  MeanSqs F.Model      R2    Pr(>F)    
#block                   3   0.07411 0.024704  1.2504 0.06935 0.1454273    
#Cover                   2   0.07683 0.038417  1.9445 0.07190 0.0124938 *  
#Nitrogen                1   0.11933 0.119330  6.0399 0.11167 0.0004998 ***
#Tillage                 1   0.10290 0.102898  5.2082 0.09629 0.0004998 ***
#Cover:Nitrogen          2   0.04885 0.024427  1.2364 0.04572 0.2093953    
#Cover:Tillage           2   0.04228 0.021138  1.0699 0.03956 0.3728136    
#Nitrogen:Tillage        1   0.01637 0.016372  0.8286 0.01532 0.6036982    
#Cover:Nitrogen:Tillage  1   0.01500 0.015001  0.7593 0.01404 0.6986507    
#Residuals              29   0.57296 0.019757         0.53616              
#otal                  42   1.06863                  1.00000 

###Simper analysis--phylum level
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
simper_output_tillage_class <- simper(otutable_class,sampledata_class$Tillage,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
tillage_simper_class <- summary(simper_output_tillage_class,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_tillage_class
write.csv(as.data.frame(tillage_simper_class$Conventional_No_tillage),file = "1stsimper_output_tillage_class.csv")

simper_output_nitrogen_class <- simper(otutable_class,sampledata_class$Nitrogen,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
nitrogen_simper_class <-summary(simper_output_nitrogen_class,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_nitrogen_class
write.csv(as.data.frame(nitrogen_simper_class$N0_N60),file = "1stsimper_output_nitrogen_class.csv")

simper_output_cover_class <- simper(otutable_class,sampledata_class$Cover,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
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

##Simper analysis--family level
sampledata_family <- sample_data(family_r_rarefied_mothur_phyloseq)
sr_family_r_rarefied_mothur_phyloseq <- t(sqrt(otu_table(family_r_rarefied_mothur_phyloseq)))#square root normalization
simper_output_tillage_family <- simper(sr_family_r_rarefied_mothur_phyloseq,sampledata_family$Tillage,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
tillage_simper_family <-summary(simper_output_tillage_family,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_tillage_family
write.csv(as.data.frame(tillage_simper_family$Conventional_No_tillage),file = "1stsimper_output_tillage_family.csv")

simper_output_nitrogen_family <- simper(sr_family_r_rarefied_mothur_phyloseq,sampledata_family$Nitrogen,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
nitrogen_simper_family <- summary(simper_output_nitrogen_family,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_nitrogen_family
write.csv(as.data.frame(nitrogen_simper_family$N0_N60),file = "1stsimper_output_nitrogen_family.csv")

simper_output_cover_family <- simper(sr_family_r_rarefied_mothur_phyloseq,sampledata_family$Cover,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
cover_simper_family <-summary(simper_output_cover_family,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_cover_family
write.csv(as.data.frame(cover_simper_family$No_cover_Vetch),file = "1stsimper_output_cover_family_No_cover_Vetch.csv")
write.csv(as.data.frame(cover_simper_family$No_cover_Wheat),file = "1stsimper_output_cover_family_No_cover_Wheat.csv")
write.csv(as.data.frame(cover_simper_family$Vetch_Wheat),file = "1stsimper_output_cover_family_Vetch_Wheat.csv")

##Simper analysis--genus level
sampledata_genus <- sample_data(genus_r_rarefied_mothur_phyloseq)
sr_genus_r_rarefied_mothur_phyloseq <- t(sqrt(otu_table(genus_r_rarefied_mothur_phyloseq)))#square root normalization
simper_output_tillage_genus <- simper(sr_genus_r_rarefied_mothur_phyloseq,sampledata_genus$Tillage,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
tillage_simper_genus <- summary(simper_output_tillage_genus,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_tillage_genus
write.csv(as.data.frame(tillage_simper_genus$Conventional_No_tillage),file = "1stsimper_output_tillage_genus.csv")

simper_output_nitrogen_genus <- simper(sr_genus_r_rarefied_mothur_phyloseq,sampledata_genus$Nitrogen,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
nitrogen_simper_genus <- summary(simper_output_nitrogen_genus,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_nitrogen_genus
write.csv(as.data.frame(nitrogen_simper_genus$N0_N60),file = "1stsimper_output_nitrogen_genus.csv")

simper_output_cover_genus <- simper(sr_genus_r_rarefied_mothur_phyloseq,sampledata_genus$Cover,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
cover_simper_genus <- summary(simper_output_cover_genus,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_cover_genus
write.csv(as.data.frame(cover_simper_genus$No_cover_Vetch),file = "1stsimper_output_cover_genus_No_cover_Vetch.csv")
write.csv(as.data.frame(cover_simper_genus$No_cover_Wheat),file = "1stsimper_output_cover_genus_No_cover_Wheat.csv")
write.csv(as.data.frame(cover_simper_genus$Vetch_Wheat),file = "1stsimper_output_cover_genus_Vetch_Wheat.csv")

##Simper analysis--otu level
sampledata_otu <- sample_data(r_rarefied_mothur_phyloseq)
sr_otu_r_rarefied_mothur_phyloseq <- t(sqrt(otu_table(r_rarefied_mothur_phyloseq)))#square root normalization
simper_output_tillage_otu <- simper(sr_otu_r_rarefied_mothur_phyloseq,sampledata_otu$Tillage,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
tillage_simper_otu <- summary(simper_output_tillage_otu,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_tillage_otu
write.csv(as.data.frame(tillage_simper_otu$Conventional_No_tillage),file = "1stsimper_output_tillage_otu.csv")

simper_output_nitrogen_otu <- simper(sr_otu_r_rarefied_mothur_phyloseq,sampledata_otu$Nitrogen,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
nitrogen_simper_otu <- summary(simper_output_nitrogen_otu,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_nitrogen_otu
write.csv(as.data.frame(nitrogen_simper_otu$N0_N60),file = "1stsimper_output_nitrogen_otu.csv")

simper_output_cover_otu <- simper(sr_otu_r_rarefied_mothur_phyloseq,sampledata_otu$Cover,permutations =100, trace = FALSE, parallel = getOption("mc.cores"))
cover_simper_otu <- summary(simper_output_cover_otu,ordered = TRUE,digits = max(3,getOption("digits") - 3))
simper_output_cover_otu
write.csv(as.data.frame(cover_simper_otu$No_cover_Vetch),file = "1stsimper_output_cover_otu_No_cover_Vetch.csv")
write.csv(as.data.frame(cover_simper_otu$No_cover_Wheat),file = "1stsimper_output_cover_otu_No_cover_Wheat.csv")
write.csv(as.data.frame(cover_simper_otu$Vetch_Wheat),file = "1stsimper_output_cover_otu_Vetch_Wheat.csv")

write.csv(as.data.frame(tax_table(r_rarefied_mothur_phyloseq)),file = "wholetaxtable.csv")
write.csv(as.data.frame(t(tax_table(phylum_r_rarefied_mothur_phyloseq))),file = "phylum_taxtable.csv")
write.csv(as.data.frame(t(tax_table(class_r_rarefied_mothur_phyloseq))),file = "class_taxtable.csv")
write.csv(as.data.frame(t(tax_table(order_r_rarefied_mothur_phyloseq))),file = "order_taxtable.csv")
write.csv(as.data.frame(t(tax_table(family_r_rarefied_mothur_phyloseq))),file = "family_taxtable.csv")
write.csv(as.data.frame(t(tax_table(genus_r_rarefied_mothur_phyloseq))),file = "genus_taxtable.csv")
