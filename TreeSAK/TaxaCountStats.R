
################################################################################

# R script for ranking marker proteins for generating concatenated species trees (127 taxa set)
# Dombrowski et al., 2020
# Finalized: February 2020

# modified by Weizhi
# Rscript /Users/songweizhi/PycharmProjects/Sponge_Hologenome/Scripts/TaxaCountStats.R -t treefile_v2.tre -l List_of_trees_2.txt -g mapping_3.txt -x MarkerList.txt -s TaxaCounts_op.txt -r Genes_to_remove.txt -o a.txt

####################################### argument parser ######################################
library(optparse)
option_list = list(
  make_option(c("-t", "--tree"),       type="character", help="combined_contree_file"),
  make_option(c("-l", "--treelist"),   type="character", help="list_of_trees_txt"),
  make_option(c("-g", "--mapping"),    type="character", help="mapping_txt"),
  make_option(c("-x", "--markerlist"), type="character", help="marker_list_txt"),
  make_option(c("-s", "--cstop"),      type="character", help="combined_count_sister_taxa_output"),
  make_option(c("-r", "--removegene"), type="character", help="genes_to_remove_txt"),
  make_option(c("-o", "--output"),     type="character", help="output table")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

combined_contree_file         = opt$tree
list_of_trees_txt             = opt$treelist
mapping_txt                   = opt$mapping
marker_list_txt               = opt$markerlist
combined_count_sister_taxa_op = opt$cstop
genes_to_remove_txt           = opt$removegene
output_table                  = opt$output

################################################################################

#rm(list=ls())
sessionInfo()

################################################################################
#0.1 setting working directory (!adjust wdir accordingly!)
################################################################################

# setting working directory (!adjust wdir accordingly!)
#wdir <- "/Users/songweizhi/Desktop/Anja_paper/Nina/4_151Marker_analyses/127_taxa"
#wdir <- "/Users/songweizhi/Desktop/Input_folder_to_R"
#setwd(wdir)

# Weizhi
# List_of_trees_2.txt                                     id of marker genes (HOGs) (152, mind order)
# treefile_v2.tre                                         tree corresponding to each marker gene (HOG)
# mapping_3.txt                                           taxonomy/cluster/group/color for each genome (Domain is higher than cluster)
# Genes_to_remove.txt                                     basically a list of marker gene ids
# MarkerList.txt                                          basically a list of marker gene ids (152, mind order)
# TaxaCounts_151MarkerGenes_ArcRefv5UAP2_129taxa_v5.txt   concatenated output from count_sister_taxa.py

################################################################################
# load packages
################################################################################

suppressMessages(library("plyr"))
suppressMessages(library("dbplyr"))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gplots))
suppressMessages(library(ape))

################################################################################
#1. read in the treefiles (all concatenated in one large document)
################################################################################
#read in concatenated tree and list of trees
tree_order <- read.table(list_of_trees_txt, sep="\t", header=F, fill=TRUE, quote = "")

trees<-read.tree(combined_contree_file)


#make a table from taxa labels and count how many taxa are in each tree
y <- c()
for(i in 1:length(trees)){
  x <- length(trees[[i]]$tip.label)
  y <- rbind (y,x)
}

Species_in_tree <- as.data.frame(y)
rownames(Species_in_tree) <- as.character(tree_order$V1)
colnames(Species_in_tree) <- "NrSpecies"
#head(Species_in_tree)

################################################################################
#2. read in taxa mapping files and stats from tom's script 
################################################################################

#taxa - taxonomy mapping file
mapping <- read.table(mapping_txt, sep="\t", header=T, fill=TRUE, quote = "", comment.char = "", check.names = FALSE)
#head(mapping)

#cleanup and shorten mapping file
mapping_clean <-unique(mapping[,c("Cluster", "Domain")])
#head(mapping_clean)
#dim(mapping_clean)

#list of genes to remove because in these trees archaea were not monphyletic
genes_to_remove <- read.table(genes_to_remove_txt, sep="\t", header=T, fill=TRUE, quote = "")

#list of total markers used in these analyses
MarkerList <- read.table(marker_list_txt, sep="\t", header=T, fill=TRUE, quote = "")

#read in statistics file from count_sister_taxa.py
SisterCounts <- read.table(combined_count_sister_taxa_op, sep="\t", header=T, fill=TRUE, quote = "")
#colnames(SisterCounts)

################################################################################
#3. transform data
################################################################################
#reduce table on and remove hits with low support (0.1 in this case but can be changed)
SisterCounts_temp0 <- subset(SisterCounts, Normalized2_sum_of_occurances >= 0.1)

#control, whether something is missing in the mapping file, the first setdiff is the relevant one!
List_taxa <- as.character(unique(SisterCounts_temp0$Group_of_interest))
#setdiff(List_taxa, mapping$Cluster)

#add in Group Info (i.e. DPANN, Eury, TACK) for Group of interest (needed to define HGT events)
SisterCounts_temp1 <-  merge(SisterCounts_temp0, mapping_clean, by.x = "Group_of_interest", by.y = "Cluster", all.x = T)
colnames(SisterCounts_temp1) <- c( "Group_of_interest", "MarkerID", "Sister_taxa", "Normalized_sum_of_occurances","splits","Normalized2_sum_of_occurances","Clusters", "Group_of_interest_Group")
#head(SisterCounts_temp1)

#add in Group Info (i.e. DPANN, Eury, TACK) for Sister_taxa (needed to define HGT events)
SisterCounts_temp2 <-  merge(SisterCounts_temp1, mapping_clean, by.x = "Sister_taxa", by.y = "Cluster", all.x = T)
colnames(SisterCounts_temp2) <- c( "Sister_taxa", "Group_of_interest", "MarkerID","Normalized_sum_of_occurances", "splits","Normalized2_sum_of_occurances","Clusters","Group_of_interest_Group", "Sister_taxon_Group")
#head(SisterCounts_temp2)

#resort dataframe for aesthetics
SisterCounts_temp3 <- SisterCounts_temp2[,c("MarkerID","Group_of_interest", "Group_of_interest_Group", "Sister_taxa", "Sister_taxon_Group", "Normalized_sum_of_occurances","splits", "Normalized2_sum_of_occurances","Clusters")]
#head(SisterCounts_temp3)

#count nr of total splits
SisterCounts_temp4 <- cbind(SisterCounts_temp3, count.fields(textConnection(as.character(SisterCounts_temp3$Clusters)), sep = ","))
#head(SisterCounts_temp4)

#make new column and make a remark whether clusters of interest are split or not
#Notice: The column "Clusters" lists if there is a split (i.e. if UAP2 12 then all 12 MAGs are together, if UAP2 has 8,4 then UAP2 is split once with one cluster with 8 and the other with 4 taxa)
SisterCounts_temp4$SplitGroups<- ifelse(grepl(",",SisterCounts_temp3$Clusters), "split", "no")
#head(SisterCounts_temp4)

#rename a column for better readability
names(SisterCounts_temp4)[names(SisterCounts_temp4) == "count.fields(textConnection(as.character(SisterCounts_temp3$Clusters)), "] <- 'NrSplits'
#head(SisterCounts_temp4)

#remove dublicates and keep the Group_of_interest with the best Normalized2_sum_of_occurances value
#this is done to only have one hits per arcog and group of interest to better normalize the data by the total nr of arcogs and do have a consistent link to the total number of phylogenetic clusters
SisterCounts_best <-
  SisterCounts_temp4 %>%
  group_by(MarkerID,Group_of_interest) %>%
  filter(Normalized2_sum_of_occurances == max(Normalized2_sum_of_occurances)) 

#count nr of taxonomic groups (i.e. clusters)  in each tree
Nr_clusters <- Number_of_taxa <- ddply(SisterCounts_best, .(MarkerID), summarize, NrClusters = length(Group_of_interest))
#head(Nr_clusters)

#merge nr of clusters and nr of species with datatable
SisterCounts_best_temp1 <- merge(SisterCounts_best, Nr_clusters, by = "MarkerID")
SisterCounts_best_temp2 <- merge(SisterCounts_best_temp1,Species_in_tree, by.x = "MarkerID", by.y = "row.names" )
#head(SisterCounts_best_temp2)

#print table          
# write.table(SisterCounts_best_temp1, "2_Output/Taxa_Summary_1.txt",  sep = "\t", row.names = F, quote =F)

################################################################################
#4. summarize split events to be able to rank marker genes
################################################################################
#summarize splits/cluster
Split_counts <- ddply(SisterCounts_best_temp1, .(MarkerID,SplitGroups, NrClusters), summarise, quantity = length(NrSplits))
Split_counts_wide <- spread(Split_counts, SplitGroups, quantity)

#make new column to calulate the percentage of split clusters
Split_counts_wide$SplitsPerCluster <- round((Split_counts_wide$split/Split_counts_wide$NrClusters)*100, digits = 1)
#head(Split_counts_wide)

#summarize and coun the total number of splits
Split_Total <- ddply(SisterCounts_best_temp1, .(MarkerID, NrClusters), summarise, TotalSplits = sum(NrSplits))
#head(Split_Total)

#combine the two dataframes generated above
#Summary_temp1 <- merge(Split_counts_wide[,c("MarkerID","NrClusters","SplitsPerCluster")],HGT_Counts, by = "MarkerID" )
Summary_temp2 <- merge(Split_counts_wide[,c("MarkerID","NrClusters","SplitsPerCluster")], Split_Total, by = "MarkerID")
Summary_temp3 <- merge(Summary_temp2, Species_in_tree, by.x = "MarkerID", by.y = "row.names")
Summary_temp3$TotalSplits_to_Species <- round((Summary_temp3$TotalSplits/Summary_temp3$NrSpecies)*100, digits = 1)
#head(Summary_temp3)

#subset to only print relevant info 
Summary_temp4 <- Summary_temp3[,c("MarkerID","NrSpecies", "NrClusters.x", "SplitsPerCluster","TotalSplits", "TotalSplits_to_Species" )]
#head(Summary_temp4)

#if genes were lost already during the tree building step, add that info in
Summary_temp5 <- merge(MarkerList, Summary_temp4, by = "MarkerID", all.x = T)
#head(Summary_temp5)

################################################################################
#5. find highest/lowest 25/50% ranking markers
################################################################################
#make vector of genes that are not good marker genes based on literature and that are not monophyletic
genes_to_remove_vector <- as.character(genes_to_remove$MarkerID) 
#genes_to_remove_vector

#remove genes from dataframe
Stats_temp1A <- Summary_temp4[ ! Summary_temp4$MarkerID %in% genes_to_remove_vector, ]
#dim(Summary_temp4)
#dim(Stats_temp1A)

#define a cutoff to remove gene trees that have less than 50% of the species as we do not want to use these genes for concatenations
cutoff <- mean(Stats_temp1A$NrSpecies)/2
#cutoff

#remove genes that have few species
Stats_temp1B <- subset(Stats_temp1A, Stats_temp1A[ , "NrSpecies"] > cutoff) 
#dim(Stats_temp1A)
#dim(Stats_temp1B)

#rank according to the split clusters in percentage from 1 to xx (lowest/best value = lowest nr) = RankA
Stats_temp2 <- Stats_temp1B %>% mutate(RankA = rank(SplitsPerCluster, ties.method = 'first'))

#rank according to the total splits normalized by the total number of species from 1 to xx (lowest value = lowest nr)
Stats_temp3 <- Stats_temp2 %>% mutate(RankB = rank(TotalSplits_to_Species, ties.method = 'first'))

#combine RankA and RankB to get the best for each method
Stats_temp3$RankA_B <- Stats_temp3$RankA+Stats_temp3$RankB
#dim(Stats_temp3)

#define the concatenated marker sets and create vectors
nr_genes <- length(Stats_temp3$MarkerID)
cutoff_25perc <- round(nr_genes/4, digits = 0)
cutoff_50perc<- round(nr_genes/2, digits = 0)
#cutoff_25perc
#cutoff_50perc

#subset the tables for the different cutoffs
best_50perc <- as.data.frame(Stats_temp3 %>% top_n(-cutoff_50perc, RankA_B))
Stats_temp3$best_50perc <- best_50perc$MarkerID[match(Stats_temp3$MarkerID, best_50perc$MarkerID)]

best_25perc <- Stats_temp3 %>% top_n(-cutoff_25perc, RankA_B)
Stats_temp3$best_25perc <- best_25perc$MarkerID[match(Stats_temp3$MarkerID, best_25perc$MarkerID)]

worst_50perc <- Stats_temp3 %>% top_n(cutoff_50perc, RankA_B)
Stats_temp3$worst_50perc <-  worst_50perc$MarkerID[match(Stats_temp3$MarkerID, worst_50perc$MarkerID)]

worst_25perc <- Stats_temp3 %>% top_n(cutoff_25perc, RankA_B)
Stats_temp3$worst_25perc <-   worst_25perc$MarkerID[match(Stats_temp3$MarkerID, worst_25perc$MarkerID)]

Stats_temp3$FullSet <- Stats_temp3$MarkerID

#merge with original table (to keep the statistics)
Stats_temp4 <- merge(Summary_temp5, Stats_temp3[,c("MarkerID", "RankA", "RankB", "RankA_B", "FullSet", "best_50perc", "best_25perc", "worst_50perc","worst_25perc")], by = "MarkerID", all.x = T)

#print
write.table(Stats_temp4, output_table, sep = "\t", row.names = F, quote =F, na = "")

