#############################################################
#
# Ref to the ARTICLE
# 
# Alessandro Passera & Davide Bulgarelli
# alessandro.passera@unimi.it
# d.bulgarelli@dundee.ac.uk
# 
# script to reproduce calculations and figures presented in the manuscript
# 
# Disclaimer: this is a temporary file, the script might be subjected to changes derived from the revision process
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#1st time installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")


#required packages 
library("phyloseq")
library("ggplot2")
library("vegan")
library ("ape")


#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#set the working directory, depending on the path in the computer used
setwd("YOUR_PATH_TO_DIRECTORY")

#############################################################
#import the count matrix and the design file
#############################################################


#OTU table generated using QIIME 1.9.0. 
dat_info <- read.delim("JH17_otu_table_SILVA132_97_nc2_nCMC.txt", skip=1, sep = "\t", header=T, row.names=1, blank.lines.skip = FALSE)

#inspect the file 
dim(dat_info)
colnames(dat_info)

#design file
#This file is designated Worksheet_ws1 in Supplementary Information
design <- read.delim("Mapping_File_All.txt", sep = "\t", header=TRUE, row.names=1)
design

#create a dataframe containing taxonomic information
dat_tax <- as.data.frame(dat_info[, 34])
class(dat_tax)
dat_tax
rownames(dat_tax) <- rownames(dat_info)
dat_tax[1:10, ]
dat_info[1:10, 34]
#write.table(dat_tax, file="Mais_tax_ordered.txt", sep ="\t")

#create a file with count data for maize only
dat_count <-dat_info[, rownames(design)] 

#inspect the file 
dim(dat_count)
colnames(dat_count)

#extract the total number of reads clustered at OTU 97% identiy (the number underneath the Hv identifier represents the total number of reads clustered for that sample) 
OTU_97_reads <- sort(colSums(dat_count))
OTU_97_reads

#total reads
OTU_97_reads_sum <- sum(colSums(dat_count))
OTU_97_reads_sum


#############################################################
#Genererate the phyloseq object
#Data required: dat_count_noplants; design, JH09_FC_dat_tax_noPlants_ordered.txt, and 97_otus.tree.gz
#############################################################

#The OTU Table counts
Mais_OTU <- otu_table(dat_count, taxa_are_rows=TRUE)

#The taxonomy information
#it is a tab-delimited file with 8 columns, the column headers are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
#The original taxonomy file has been generated on lines 75-82 and re-named in excel
Mais_taxa_ordered <- read.delim ("Mais_tax_ordered.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = FALSE)
Mais_taxa <- tax_table(as.matrix(Mais_taxa_ordered))
dim(Mais_taxa)

#The mapping file 
Mais_map <- sample_data(design)

#The phylogenetic tree: the OTU table has been generated using a closed reference approach agains the greengenes 13_05 database use the corresponding phylogenetic tree 
SILVA_tree <- read_tree("97_otus.tre")

#check whether the tree is rooted
is.rooted(SILVA_tree)

#merge the files and create the phyloseq object
Mais_data_phyloseq <- merge_phyloseq(Mais_OTU, Mais_taxa, Mais_map,  SILVA_tree)
Mais_data_phyloseq

#remove sample with less than 50 reads
Mais_data_phyloseq_2 <- prune_samples(sample_sums(Mais_data_phyloseq) > 49, Mais_data_phyloseq)
sort(sample_sums(Mais_data_phyloseq_2))

#Transform the count in relative abundance cpm
Mais_data_phyloseq_prop <- transform_sample_counts(Mais_data_phyloseq_2,  function(x) 1e+03 * x/sum(x))

#################
#Alpha diversity
#################
p <- plot_richness(Mais_data_phyloseq_2, x= "Code", measures = c("Observed", "Chao1", "Shannon"), color = "Code")
p <- p + geom_boxplot(data = p$data, aes(x = Code, y = value, color = NULL), alpha = 0.5)
p

##################################
#Beta diversity
#################################

#PCoA weighted unifrac distance
#info Unifrac: https://en.wikipedia.org/wiki/UniFrac
Mais_data_phyloseq_prop_wunifrac <- ordinate(Mais_data_phyloseq_prop, "PCoA", "unifrac", weighted = TRUE)
p=plot_ordination(Mais_data_phyloseq_prop, Mais_data_phyloseq_prop_wunifrac , shape ="Year", color ="Code")
p = p + geom_point(size = 6, alpha = 0.75)
p = p + scale_colour_manual(values = c("#1D91C0", "#67001F", "#CB181D", "#78C679", "#6A51A3", "#000000"))
p + ggtitle("PCoA 16S data, Unifrac distance")

#WU distance adonis
WU <- phyloseq::distance(Mais_data_phyloseq_prop, "unifrac", weighted= TRUE)
adonis(WU ~ Genotype*Year, design[colnames(otu_table(Mais_data_phyloseq_prop)), ], permutations = 10000)
                                                   
#UU distance adonis
UU <- phyloseq::distance(Mais_data_phyloseq_prop, "unifrac", weighted= FALSE)
adonis(UU ~ Genotype*Year, design[colnames(otu_table(Mais_data_phyloseq_prop)), ], permutations = 10000)

#############################
#Histograms
#############################

#Phylum level
Mais_data_Phylum <- tax_glom(Mais_data_phyloseq_2, taxrank = "Phylum")
Mais_data_Phylum_2 <- merge_samples(Mais_data_Phylum, "GxY")
design_glom <- read.delim("Mapping_File_Glom.txt", sep = "\t", header=TRUE, row.names=1)
Mais_data_Phylum_3 <- merge_phyloseq(otu_table(Mais_data_Phylum_2), sample_data(design_glom), tax_table(Mais_data_Phylum_2), phy_tree(Mais_data_Phylum_2))
Mais_data_Phylum_4 <- transform_sample_counts(Mais_data_Phylum_3,  function(x) 1e+02 * x/sum(x))
Mais_data_Phylum_5 <- filter_taxa(Mais_data_Phylum_4, function(x) mean(x) >0.1, TRUE)
Phylum_color <- c("#E5F5F9", "#1D91C0", "#67001F", "#CB181D", "#78C679", "#F46D43", "#A6CEE3", "#FD8D3C", "#A6D854", "#D4B9DA", "#6A51A3")
p <- plot_bar(Mais_data_Phylum_5, fill = "Phylum", x = "Code", facet_grid = "Year")
p = p + scale_fill_manual(values = Phylum_color)
p

#Family level
Mais_data_Fam <- tax_glom(Mais_data_phyloseq_2, taxrank = "Family")
Mais_data_Fam_2 <- merge_samples(Mais_data_Fam, "GxY")
design_glom <- read.delim("Mapping_File_Glom.txt", sep = "\t", header=TRUE, row.names=1)
Mais_data_Fam_3 <- merge_phyloseq(otu_table(Mais_data_Fam_2), sample_data(design_glom), tax_table(Mais_data_Fam_2), phy_tree(Mais_data_Fam_2))
Mais_data_Fam_4 <- transform_sample_counts(Mais_data_Fam_3,  function(x) 1e+02 * x/sum(x))
Top_N_Fams <- names(sort(taxa_sums(Mais_data_Fam_4), TRUE)[1:15])
Mais_data_Fam_Top <- prune_taxa(Top_N_Fams, Mais_data_Fam_4)
Mais_data_Fam_5 <- filter_taxa(Mais_data_Fam_Top, function(x) mean(x) >0.1, TRUE)
Fam_color <- c("#E5F5F9", "#1D91C0", "#67001F", "#CB181D", "#78C679", "#F46D43", "#A6CEE3", "#FD8D3C", "#A6D854", "#D4B9DA", "#6A51A3", "#7F0000", "#D9D9D9", "#FFF7BC", "#000000")
p <- plot_bar(Mais_data_Fam_5, fill = "Family", x = "Code", facet_grid = "Year")
p = p + scale_fill_manual(values = Fam_color)
p
