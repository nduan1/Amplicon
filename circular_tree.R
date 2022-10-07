### Read into shared and taxonomy file
library(phyloseq)
sharedfile_tx = "1stning16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.tx.1.subsample.shared"
taxfile_tx = "1stning16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.tx.1.cons.taxonomy"
tree_mothur <- import_mothur(mothur_shared_file = sharedfile, mothur_constaxonomy_file = taxfile)
colnames(tax_table(tree_mothur)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
sort(sample_sums(tree_mothur),decreasing = FALSE)
rms_tree_mothur<-prune_taxa(taxa_sums(tree_mothur)>5, tree_mothur)
sort(sample_sums(rms_tree_mothur),decreasing = FALSE)
#remove the BLANK
sample_names(rms_tree_mothur)
sub_rms_tree_mothur <- prune_samples(sample_names(rms_tree_mothur)[-c(1)],rms_tree_mothur)#remove VCNO1-4 and VNTN60_4
sub_rms_tree_mothur
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 25366 taxa and 44 samples ]
#tax_table()   Taxonomy Table:    [ 25366 taxa by 6 taxonomic ranks ]
rarefied_tree_mothur<-rarefy_even_depth(sub_rms_tree_mothur,sample.size = min(sample_sums(sub_rms_tree_mothur)),rngseed = 1013,replace = TRUE,trimOTUs = TRUE,verbose = TRUE)#normalize the sequences based on the least number of sequences
rarefied_tree_mothur
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 23045 taxa and 44 samples ]
#tax_table()   Taxonomy Table:    [ 23045 taxa by 6 taxonomic ranks ]
genus_rarefied_tree_mothur<-tax_glom(rarefied_tree_mothur,taxrank = rank_names(rarefied_tree_mothur)[6],NArm=TRUE)
tax_table(genus_rarefied_tree_mothur)
taxa_sums(genus_rarefied_tree_mothur)
t2_genus_rarefied_tree_mothur <-  prune_taxa(names(sort(taxa_sums(genus_rarefied_tree_mothur), TRUE)[1:200]), genus_rarefied_tree_mothur)#Prune further, to top 200 most abundant taxa
r_t2_genus_rarefied_tree_mothur<-transform_sample_counts(t2_genus_rarefied_tree_mothur,function(x)x/sum(x))#relative abundance
#merge phyloseq
tree_meta<-read.csv("mothur_meta_tree.csv",sep=",",header=TRUE,row.names = 1)#this file already pruned manually(remove the bad samples)
dim(tree_meta)
sample_names(r_t2_genus_rarefied_tree_mothur)#check the names of sample in both otu table and meta table
tree_moth_merge<-merge_phyloseq(r_t2_genus_rarefied_tree_mothur,sample_data(tree_meta))
merge_tree_moth <- merge_samples(tree_moth_merge,group = "Treat",fun = mean)#as documented, mean function is ignored for the OTU_table,so in fact the function is sum instead of mean.
merge_tree_moth
sample_sums(merge_tree_moth)#equals to sample number of each treatment
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 200 taxa and 12 samples ]
#sample_data() Sample Data:       [ 12 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 200 taxa by 6 taxonomic ranks ]
### Write out the shared file and the taxonomy file for GraPhlAn analysis

genus_graphlan_OTU <- otu_table(merge_tree_moth)
dim(genus_graphlan_OTU)
genus_graphlan_OTU[1:5,1:5]
genus_graphlan_taxonomy<-tax_table(merge_tree_moth)
dim(genus_graphlan_taxonomy)
genus_graphlan_taxonomy[1:5,]
write.csv(genus_graphlan_taxonomy,file="genus_table_for_graphlan.csv")
genus_graphlan_combine<-data.frame(t(genus_graphlan_OTU),genus_graphlan_taxonomy,tax_sum=taxa_sums(genus_graphlan_OTU))
dim(genus_graphlan_combine)
genus_graphlan_combine[1:5,]
write.csv(genus_graphlan_combine,file="otu_tax_combine_for_graphlan.csv")
#!manually revise the genus or class name which contain "uncultured", and then read into r and then do the next command rownames(genus_graphlan_combine)<-genus_graphlan_combine$Genus
genus_graphlan_combine <- read.csv(file = "otu_tax_combine_for_graphlan.csv",sep = ',',row.names = 1)
rownames(genus_graphlan_combine)<-genus_graphlan_combine$Genus#replace the otu rownames as genus name
write.csv(genus_graphlan_combine,file="otu_tax_combine_for_graphlan.csv")
### Generate node size file
#-----Genus------
Genus_node_size<-data.frame(name=tax_table(merge_tree_moth)[,6],size=taxa_sums(merge_tree_moth))
dim(Genus_node_size)
Genus_node_size[1:5,]
write.csv(Genus_node_size,file="Genus_node_size.csv")
#----Family------
family_graphlan_phyloseq<-tax_glom(merge_tree_moth,taxrank = rank_names(merge_tree_moth)[5],NArm=TRUE) 
Family_node_size<-data.frame(name=tax_table(family_graphlan_phyloseq)[,5],size=taxa_sums(family_graphlan_phyloseq))
dim(Family_node_size)
Family_node_size[1:5,]
write.csv(Family_node_size,file="Family_node_size.csv")
#----Order-------
order_graphlan_phyloseq<-tax_glom(merge_tree_moth,taxrank = rank_names(merge_tree_moth)[4],NArm = TRUE)
Order_node_size<-data.frame(name=tax_table(order_graphlan_phyloseq)[,4],size=taxa_sums(order_graphlan_phyloseq))
dim(Order_node_size)
Order_node_size[1:5,]
write.csv(Order_node_size,file="Order_node_size.csv")
#----Class----
class_graphlan_phyloseq<-tax_glom(merge_tree_moth,taxrank = rank_names(merge_tree_moth)[3],NArm = TRUE)
Class_node_size<-data.frame(name=tax_table(class_graphlan_phyloseq)[,3],size=taxa_sums(class_graphlan_phyloseq))
dim(Class_node_size)
Class_node_size[1:5,]
write.csv(Class_node_size,file="Class_node_size.csv")
#----Phylum-----
phylum_graphlan_phyloseq<-tax_glom(merge_tree_moth,taxrank = rank_names(merge_tree_moth)[2],NArm = TRUE)
Phylum_node_size<-data.frame(name=tax_table(phylum_graphlan_phyloseq)[,2],size=taxa_sums(phylum_graphlan_phyloseq))
dim(Phylum_node_size)
Phylum_node_size[1:5,]
write.csv(Phylum_node_size,file="Phylum_node_size.csv")
unique(phylum_graphlan_phyloseq$Phylum)
taxa_names(phylum_graphlan_phyloseq)

