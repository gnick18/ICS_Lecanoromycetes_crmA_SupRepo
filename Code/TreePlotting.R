# Scripts for Sigh Crm paper
## Figure 1: The main ics tree, unrooted with gcfs overlaid
#uniform ggtree plotting
library(ggtree)
library(ape)
library(treeio)
library(phytools)
library(hash)
library(ggnewscale)
library(ggtreeExtra)
library(tidyr)
library(ggplot2)
library(tibble)
library(dplyr)
library(svglite)
library(aplot)
library(ggpattern)
library(rgl)
library(treespace)

rm(list=ls())

###########
# PLOT 1: making the full ICS tree with certain clades highlighted
###########
metaPath <- "/Volumes/T7/MiscProjects/Singh/April_FinalWork/ExtractLFFs_BGCs/LFF3_4/ICS_Meta_withTaxonomy.tsv"
meta <- read.table(metaPath, sep="\t", header=TRUE, stringsAsFactors=FALSE)
treePath <- "/Volumes/T7/MiscProjects/Singh/FinalVersion/Filtered_ICS_Protein_Trees/Filtered_ICS_Proteins.trimmed.treefile"
tree <- read.tree(treePath)


# check if it is unrooted, we don't want it to me
is.rooted(tree)

## Step 0: Cleaning up the data
mpos <- match(meta$TreeTip_Name, tree$tip.label)
unmatched <- meta$TreeTip_Name[is.na(mpos)]
# Print them for inspection
if (length(unmatched) > 0) {
  cat("Warning: The following TreeTip_Name entries were not found in the tree:\n")
  print(unmatched)
}
meta <- meta[!is.na(mpos), ]


## Finding all of the tips that are not in the metadata
missing_tips <- setdiff(tree$tip.label, meta$TreeTip_Name)
# Print them for inspection
if (length(missing_tips) > 0) {
  cat("Warning: The following tips are not in the metadata:\n")
  print(missing_tips)
}

# Step 1: removing the duplicate tips from singh and my database
removeThese <- c(
  "Xylographa_vitiligo_MMC15_003146",
  "Xylographa_vitiligo_MMC15_001555",
  "Xylographa_vitiligo_MMC15_000673",
  "Xylographa_trunciseda_MMC11_005954",
  "Xylographa_trunciseda_MMC11_003068",
  "Xylographa_soralifera_MMC17_008632",
  "Xylographa_soralifera_MMC17_003108",
  "Xylographa_soralifera_MMC17_001664",
  "Xylographa_pallens_MMC27_005897",
  "Xylographa_pallens_MMC27_003017",
  "Xylographa_opergraphella_MMC26_003333",
  "Xylographa_opergraphella_MMC26_002878",
  "Xylographa_opergraphella_MMC26_002626",
  "Xylographa_carneopallida_MMC34_006463",
  "Xylographa_carneopallida_MMC34_002691",
  "Xylographa_bjoerkii_MMC18_003552",
  "Xylographa_bjoerkii_MMC18_000573",
  "Variospora_aurantia_LQ341_001612",
  "Variospora_aurantia_LQ341_000328",
  "Usnochroma_carphineum",
  "Trapelia_coarctata_MMC30_006246",
  "Trapelia_coarctata_MMC30_005726",
  "Trapelia_coarctata_MMC30_005530",
  "Sticta_canariensis_MMC29_007091",
  "Sticta_canariensis_MMC29_006734",
  "Seirophora_villosa",
  "Seirophora_lacunosa",
  "Pseudocyphellaria_aurata_FUN_006565",
  "Physcia_stellaris",
  "Lobaria_immixta",
  "Icmadophila_ericetorum",
  "Gyalolechia_ehrenbergii",
  "Caloplaca_aetnensis",
  "Caloplaca_aegaea",
  "Bacidia_gigantensis_KY384_008841"
)


## checking if table and tree line up
metaTree_tips <- meta$TreeTip_Name
# tree
tree_tips <- tree$tip.label
## are there any tree tips not in the meta table?
tips_not_in_meta <- setdiff(tree_tips, metaTree_tips)

# Remove the tips from the tree
tree <- drop.tip(tree, tips_not_in_meta)

#Step 2: Convert the phylo object to a ggtree object
# Rename the TreeTipName column to ID, and making it the first column
meta <- meta %>%
  rename(ID = TreeTip_Name) %>%
  select(ID, everything())


#printing the number of tips and the number of unique accessions
cat("Number of tips in the tree:", length(tree$tip.label), "\n")
cat("Number of unique accessions in the metadata:", length(unique(meta$Accession)), "\n")

meta$tip_offset <- ifelse(meta$Lichen == "Lichen", 2, 0)  # Adjust the value (e.g., 2) to control spacing
# copying the meta table but only keeping the columns we use
meta <- meta %>%
  select(ID, Species, Accession, Lichen, ManualLabel, Protein, tip_offset, Phylum, Fungal_Class) %>%
  mutate(Lichen = factor(Lichen, levels = c("Lichen", "NotLichen")))

#making the heatmap column
meta <- meta %>%
  mutate(LichenStatus = ifelse(Lichen == "Lichen", 1, 0))

#Step 3: Adding the metadata to the tree and converting it to a ggtree
# midpoint root the tree
########
tree <- midpoint.root(tree)
########
icsTree <- ggtree(tree, layout='ape')  %<+% meta + 
  geom_tippoint(aes(
    color = Lichen,
    alpha = ifelse(Lichen == "NotLichen", 1, 1),
    size = ifelse(Lichen == "Lichen", 1.5, 0.5)
  )) +
  scale_color_manual(values = c("Lichen" = "#07f49e", "NotLichen" = "#42047e"),
                     na.value = "transparent") +
  scale_alpha_identity() +
  scale_size_identity() 

icsTree_Cladogram <- ggtree(tree, layout='circular', branch.length = 'none') %<+% meta + 
  geom_tippoint(aes(
    color = Lichen,
    size = ifelse(Lichen == "Lichen", 1, 0.5)
  )) +
  scale_color_manual(values = c("Lichen" = "#07f49e", "NotLichen" = "#42047e"),
                     na.value = "transparent") +
  scale_size_identity() +
  geom_tiplab(size = 0.5, linesize = 0.1, align = TRUE) 

############################################################
##### ADDING IN LFF 3 AND LFF 4
# finding the MRCA of LFF3 based on ones I am handtyping in based on the results in the cladogram svg file
lff3_tips <- c("Emmonsia_crescens_PGH30761.1","Schaereria_dolodes_MCJ1483472.1","Xylaria_flabelliformis_KAI0184952.1",
               "Letrouitia_leprolyta_KAI4189289.1","Xanthoria_sp_2_TBL2021_KAI4283325.1","Xanthoria_aureola_KAI4221412.1",
               "Nannizzia_gypsea_XP_003174548.1")

lff4_tips <- c("Letrouitia_leprolyta_KAI4182038.1", "Usnochroma_carphineum_KAI4132389.1","Xylographa_bjoerkii_MCJ1387245.1",
               "Penicillium_nalgiovense_OQE80469.1", "Penicillium_crustosum_KAF7521114.1","Aspergillus_saccharolyticus_XP_025429098.1",
               "Aspergillus_thermomutatus_RHZ52826.1","Penicillium_subrubescens_OKP03822.1")

lff3_subset2 <- c("Xanthoria_sp_2_TBL2021_KAI4283325.1","Xanthoria_aureola_KAI4221412.1",
                  "Lasallia_pustulata_SLM34376.1","Teloschistes_flavicans_KAI4254493.1")

dit_subset <- c("Exserohilum_turcicum_XP_008027087.1","Zopfia_rhizophila_KAF2195181.1","Neofusicoccum_parvum_XP_007585686.1",
                "Scheffersomyces_stipitis_KAG2735933.1","Colletotrichum_tropicale_KAF4837729.1","Monosporascus_sp_mg162_RYP51665.1",
                "Aspergillus_ibericus_XP_025572841.1")

mrca_lff3 <- findMRCA(tree, tips = lff3_tips, type = "node")
mrca_lff4 <- findMRCA(tree, tips = lff4_tips, type = "node")
mrca_lff3_subset2 <- findMRCA(tree, tips = lff3_subset2, type = "node")
mrca_dit_subset <- findMRCA(tree, tips = dit_subset, type = "node")

## We need to get all of these tips, and annotate them in the ManualLabel column as lichengroup 3 and 4
lff3_all_tips <- unlist(getDescendants(tree, mrca_lff3))
lff4_all_tips <- unlist(getDescendants(tree, mrca_lff4))
lff3_subset2_tips <- unlist(getDescendants(tree, mrca_lff3_subset2))
dit_subset_tips <- unlist(getDescendants(tree, mrca_dit_subset))
lff3_all_tips <- lff3_all_tips[lff3_all_tips <= length(tree$tip.label)]
lff4_all_tips <- lff4_all_tips[lff4_all_tips <= length(tree$tip.label)]
dit_all_tips <- dit_subset_tips[dit_subset_tips <= length(tree$tip.label)]
lff3_subset2_tips <- lff3_subset2_tips[lff3_subset2_tips <= length(tree$tip.label)]
lff3_ids <- tree$tip.label[lff3_all_tips]
lff4_ids <- tree$tip.label[lff4_all_tips]
lff3_subset2_ids <- tree$tip.label[lff3_subset2_tips]
dit_ids <- tree$tip.label[dit_all_tips]

meta$ManualLabel[meta$ID %in% lff3_ids] <- "lichengroup3"
meta$ManualLabel[meta$ID %in% lff4_ids] <- "lichengroup4"
meta$ManualLabel[meta$ID %in% dit_ids] <- "dit_supergroup"
## for the values that are in subset 2, we will call them lichengroup4_subset2
meta$ManualLabel[meta$ID %in% lff3_subset2_ids] <- "lichengroup3_subset2"

#### Save this table
write.table(meta, file="/Volumes/T7/MiscProjects/Singh/April_FinalWork/ExtractLFFs_BGCs/LFF3_4/newVersionDIT_Meta.tsv", sep="\t", quote=FALSE, row.names=FALSE)
######################################################

#Step 5: Adding in the GCFs that were manually annotated
labels <- unique(meta$ManualLabel)
labels <- labels[labels != ""]
# removing the crm_split label

highlightColors <- c("#CC79A7","#E69F00","#009E73",'purple','tan',
                     "#0072B2","#D55E00","#000000",'red','brown',"grey","darkblue", "darkgreen")

#Getting the MRCA nodes for all of the classes we want to highlight
mrcaNodes = hash()
for(i in 1:length(labels)) {
  label<-labels[i]
  idsInGCF_unedited<- meta[meta$ManualLabel == label, "ID"]
  mrca_node = findMRCA(tree, tips=idsInGCF_unedited, type=c("node"))
  mrcaNodes[[label]] <- mrca_node
}

counter = 1
for (label in ls(mrcaNodes)) {    
  nodeToHighlight <- mrcaNodes[[label]]
  print(label)
  # First do the ape tree
  tempTree <- icsTree +
    geom_highlight(node=nodeToHighlight, fill=highlightColors[counter], alpha=0.2)
  icsTree <- tempTree
  # Now do the cladogram
  tempTree <- icsTree_Cladogram +
    geom_highlight(node=nodeToHighlight, fill=highlightColors[counter], alpha=0.8, to.bottom = TRUE) 
  icsTree_Cladogram <- tempTree
  #
  counter <- counter + 1
}

counter = 1
for (label in ls(mrcaNodes)) {    
  nodeToHighlight <- mrcaNodes[[label]]
  # First do the ape tree
  tempTree <- icsTree +
    geom_cladelabel(node=nodeToHighlight, label=label, 
                    fill = "white", geom='label', angle="auto", horizontal=FALSE, hjust=0.5,
                    fontsize = 3, color = highlightColors[counter])
  icsTree <- tempTree
  # Now do the cladogram
  tempTree <- icsTree_Cladogram +
    geom_cladelabel(node=nodeToHighlight, label=label, 
                    fill = "white", geom='label', angle="auto", horizontal=FALSE, hjust=0.5,
                    fontsize = 3, color = highlightColors[counter])
  icsTree_Cladogram <- tempTree
  #
  counter <- counter + 1
}

# Saving this as svg files: Windows version
# saveHere = "E:/MiscProjects/Singh/Plots/ICSTree"
# Saving this as svg files: Mac version
saveHere = "/Volumes/T7/MiscProjects/Singh/Plots/LFF3_VERSIONS/"

ggsave(
  filename = paste0(saveHere, "/ICS_Tree.svg"),
  plot = icsTree,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)
ggsave(
  filename = paste0(saveHere, "/ICS_Tree_Cladogram_lab.svg"),
  plot = icsTree_Cladogram,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)

############ FOR PLOTTING MANUALLY
# make a version of icsTree where the four outgroups are labled with red text 
outgroup_1 <- "Fusarium_tricinctum_KAH7242825.1"
outgroup_2 <- "Didymella_exigua_XP_033446528.1"
outgroup_3 <- "Trichoderma_asperellum_XP_024755881.1"
outgroup_4 <- "Penicillium_expansum_KGO39525.1"
my_outgroups <- c(outgroup_1, outgroup_2, outgroup_3, outgroup_4)
icsTree_outgroupLabels <- ggtree(tree, layout = "ape") + # Ensure layout matches your "unrooted" description
  geom_treescale() +
  geom_tiplab(aes(subset = label %in% my_outgroups), 
              color = "red", 
              size = 3) # Adjust text size as needed

ggsave(
  filename = paste0(saveHere, "/ICS_Tree_OutsLabeled.svg"),
  plot = icsTree_outgroupLabels,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)
  



### Save this metafile and tree: Windows version
# saveHere = "E:/MiscProjects/Singh/Plots/ICSTree"
## Save this metafile and tree: Mac version
saveHere = "/Volumes/T7/MiscProjects/Singh/Plots/ICSTree"
# save the table as a tsv 
write.table(meta, file = paste0(saveHere, "/ICS_Meta.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
# save the phylotree as a newick file
write.tree(tree, file = paste0(saveHere, "/ICS_Tree.tre"), append = FALSE, digits = 10)


icsTree

## TESTING
# Collapsing an internal node (for example, node 5) and drawing the tree
collapsed_tree <- collapse(tree, node = 5)  # Collapse node 5
ggtree(collapsed_tree) +
  geom_treescale() +
  geom_tiplab()


###########
# PLOT 1.5: Creating smaller subsets of key regions we want to highlight
###########
## subset the tree to only include lichenGroup1 with an outgroup
lichenGroup1_tips <- meta[meta$ManualLabel == "lichengroup1", "ID"]
lichenGroup2_tips <- meta[meta$ManualLabel == "lichengroup2", "ID"]
lichenGroup3_tips <- meta[meta$ManualLabel == "lichengroup3_subset2", "ID"]
lichenGroup4_tips <- meta[meta$ManualLabel == "lichengroup4", "ID"]
outgroup_1 <- "Fusarium_tricinctum_KAH7242825.1"
outgroup_2 <- "Fusarium_tricinctum_KAH7242825.1"
outgroup_3 <- "Trichoderma_asperellum_XP_024755881.1"
outgroup_4 <- "Penicillium_expansum_KGO39525.1"
lichenGroup1_tips <- c(lichenGroup1_tips, outgroup_1)
lichenGroup2_tips <- c(lichenGroup2_tips, outgroup_2)
lichenGroup3_tips <- c(lichenGroup3_tips, outgroup_3)
lichenGroup4_tips <- c(lichenGroup4_tips, outgroup_4)

# Subset the tree to only include the tips we want
lichenGroup1_tree <- keep.tip(tree, lichenGroup1_tips)
lichenGroup2_tree <- keep.tip(tree, lichenGroup2_tips)
lichenGroup3_tree <- keep.tip(tree, lichenGroup3_tips)
lichenGroup4_tree <- keep.tip(tree, lichenGroup4_tips)

# root both trees at the outgroup
## unroot all of the trees 
lichenGroup1_tree <- unroot(lichenGroup1_tree)
lichenGroup2_tree <- unroot(lichenGroup2_tree)
lichenGroup3_tree <- unroot(lichenGroup3_tree)
lichenGroup4_tree <- unroot(lichenGroup4_tree)
#
lichenGroup1_tree <- root(lichenGroup1_tree, outgroup = outgroup_1, resolve.root = TRUE)
lichenGroup2_tree <- root(lichenGroup2_tree, outgroup = outgroup_2, resolve.root = TRUE)
lichenGroup3_tree <- root(lichenGroup3_tree, outgroup = outgroup_3, resolve.root = TRUE)
lichenGroup4_tree <- root(lichenGroup4_tree, outgroup = outgroup_4, resolve.root = TRUE)

ggtree_lg1 <- ggtree(lichenGroup1_tree, layout='rectangular') %<+% meta + 
  geom_tippoint(aes(
    color = Lichen,
  )) +
  scale_color_manual(values = c("Lichen" = "#07f49e", "NotLichen" = "#42047e"),
                     na.value = "transparent") +
  scale_size_identity() +
  ## Adding in the bootstrap values as red text
  geom_nodelab(aes(subset=!is.na(label), label = label), size = 2.5, color="red", geom='text', node='internal') +
  geom_tiplab(size=2, align=FALSE, offset=0.05)
  

ggtree_lg2 <- ggtree(lichenGroup2_tree, layout='rectangular') %<+% meta + 
  geom_tippoint(aes(
    color = Lichen,
  )) +
  scale_color_manual(values = c("Lichen" = "#07f49e", "NotLichen" = "#42047e"),
                     na.value = "transparent") +
  scale_size_identity() +
  ## Adding in the bootstrap values as red text
  geom_nodelab(aes(subset=!is.na(label), label = label), size = 2.5, color="red", geom='text', node='internal') +
  geom_tiplab(size=2, align=FALSE, offset=0.05)

ggtree_lg3 <- ggtree(lichenGroup3_tree, layout='rectangular') %<+% meta +
  geom_tippoint(aes(
    color = Lichen,
  )) +
  scale_color_manual(values = c("Lichen" = "#07f49e", "NotLichen" = "#42047e"),
                     na.value = "transparent") +
  scale_size_identity() +
  ## Adding in the bootstrap values as red text
  geom_nodelab(aes(subset=!is.na(label), label = label), size = 2.5, color="red", geom='text', node='internal') +
  geom_tiplab(size=2, align=FALSE, offset=0.05)

ggtree_lg4 <- ggtree(lichenGroup4_tree, layout='rectangular') %<+% meta +
  geom_tippoint(aes(
    color = Lichen,
  )) +
  scale_color_manual(values = c("Lichen" = "#07f49e", "NotLichen" = "#42047e"),
                     na.value = "transparent") +
  scale_size_identity() +
  ## Adding in the bootstrap values as red text
  geom_nodelab(aes(subset=!is.na(label), label = label), size = 2.5, color="red", geom='text', node='internal') +
  geom_tiplab(size=2, align=FALSE, offset=0.05)


## save these plots
cladeFolder = "/Volumes/T7/MiscProjects/Singh/Plots/LFF3_VERSIONS/CladePlots"
ggsave(paste0(cladeFolder, "/LichenGroup1_Tree.svg"), plot=ggtree_lg1, device="svg", width=4.5, height=4)
ggsave(paste0(cladeFolder, "/LichenGroup2_Tree.svg"), plot=ggtree_lg2, device="svg", width=4, height=2)
ggsave(paste0(cladeFolder, "/LichenGroup3_Tree.svg"), plot=ggtree_lg3, device="svg", width=5, height=3)
ggsave(paste0(cladeFolder, "/LichenGroup4_Tree.svg"), plot=ggtree_lg4, device="svg", width=5, height=4)


## Save


##################################################################################################################################
## PLOT 1.55: Pie chart showing the number of species I picked up with crmA before and after this analysis
##################################################################################################################################
library(RColorBrewer)
library(waffle)
color <- brewer.pal(3, "Dark2")

allSpecies = read.csv("E:/MiscProjects/Singh/Paper/SupTable1_BuscoICSCounts.tsv", sep="\t", header=TRUE)
## getting the number of species we would say has crm before vs after
ascoOnly = allSpecies[allSpecies$Phylum == "Ascomycota",]
asco_noSac = ascoOnly[ascoOnly$Fungal_Class != "Saccharomycetes",]
## drop rows with duplicate values in Species
asco_noSac = asco_noSac[!duplicated(asco_noSac$Species),]
allAscos_num = length(asco_noSac$Species)
asco_noSac$crmType = ""
# get all of the species with a fused value 
fused_meta = meta[meta$ManualLabel == "crm_fused",]
split_meta = meta[meta$ManualLabel == "crm_split",]
#for every species in the fused meta, check if it is in the ascoOnly table
for (i in 1:nrow(fused_meta)) {
  # get the species name
  species = fused_meta$Species[i]
  # check if it is in the asco_noSac table
  if (species %in% asco_noSac$Species) {
    # get the index of the row in the asco_noSac table
    index = which(asco_noSac$Species == species)
    # set the crmType to "crm_fused"
    asco_noSac$crmType[index] = "crm_fused"
  }
}
#do the same thing for the split, BUT if it's already set to "crm_fused", set it to "both"
for (i in 1:nrow(split_meta)) {
  # get the species name
  species = split_meta$Species[i]
  # check if it is in the asco_noSac table
  if (species %in% asco_noSac$Species) {
    # get the index of the row in the asco_noSac table
    index = which(asco_noSac$Species == species)
    # set the crmType to "crm_split"
    if (asco_noSac$crmType[index] == "crm_fused") {
      asco_noSac$crmType[index] = "both"
    } else {
      asco_noSac$crmType[index] = "crm_split"
    }
  }
}
# before is both + fused only
beforeNumCounted = length(which(asco_noSac$crmType == "crm_fused" | asco_noSac$crmType == "both"))
plotting_beforeDif = allAscos_num - beforeNumCounted
# after is all three
afterNumCounted = length(which(asco_noSac$crmType == "crm_fused" | asco_noSac$crmType == "both" | asco_noSac$crmType == "crm_split"))
plotting_afterDif = allAscos_num - afterNumCounted

###
nickles2023_table <- read.csv("/Users/gnickles/Desktop/ICS_Fungi_Project/DataAnalysis/PaperStuff/SupTables/SupTable5_BETTER.tsv", sep="\t", header=TRUE)
# to get the species we need to split by _ and rejoin everything with _ BUT the first element (remove the first element)
nickles2023_table$Species <- sapply(strsplit(nickles2023_table$Cluster_Name, "_"), function(x) paste(x[-1], collapse = "_"))
nickles2023_crmGroup <- nickles2023_table[nickles2023_table$FinalGCFs == "crm*",]

nickles2023_crmGroup_species <- unique(nickles2023_crmGroup$Species)
## count how many of these are in the asco_noSac table
nickles2023_count <- 0
for (i in 1:length(nickles2023_crmGroup_species)) {
  species = nickles2023_crmGroup_species[i]
  # check if it is in the asco_noSac table
  if (species %in% asco_noSac$Species) {
    # set the crmType to "crm_fused"
    nickles2023_count = nickles2023_count + 1
  }
}



#saving this using svg light
wBefore <- waffle(c(51,887), title="Before", rows=20, colors = c('#ffba49','#20a39e'))
wAfter <- waffle(c(314,624), title="After", rows=20, colors = c('#ffba49','#20a39e'))

# iron(
#   wBefore,
#   wAfter
# )
#saving the charts
saveHere_before <- "E:/MiscProjects/Singh/Plots/crmA_Counts_Before.svg"
ggsave(saveHere_before, plot=(wBefore + theme(legend.position="none")), device="svg", width=7, height=5)

saveHere_after <- "E:/MiscProjects/Singh/Plots/crmA_Counts_After.svg"
ggsave(saveHere_after, plot=(wAfter + theme(legend.position="none")), device="svg", width=7, height=5)




###########
# PLOT 2: Mapping the GCFs for the lichen and crm GCFs back to a species tree
###########
rm(list=ls())

treePath = "E:/MickeyDomainCollab/ASTRALRun/FungalSpeciesTree_FigTreeExport.tre"
antismashOverlayFile  = "E:/MickeyDomainCollab/EditedOverlay.tsv"
taxOverlayPath = "E:/MickeyDomainCollab/TaxonomyOverlay.tsv"

phyloTree <- read.tree(treePath)
antismashOverlay <- read.csv(antismashOverlayFile, sep="\t")
taxOverlay_unedited <- read.csv(taxOverlayPath, sep="\t")

#doing the edited to the tax file
#honestly a lot of this is trial and error, but the big thing is to make the row that matches the tip labels "ID"
mergeThis <- antismashOverlay[, c("Species","Accession")]
filterAccessions <- mergeThis$Accession
taxOverlay_unedited <- taxOverlay_unedited[taxOverlay_unedited$accession %in% filterAccessions, ] #filtering down the table to include only what is in the tree
taxOverlay_unedited <- taxOverlay_unedited[!duplicated(taxOverlay_unedited$accession), ] #removing any duplicate rows
names(taxOverlay_unedited)[1] <- "Accession"
taxOverlay <- merge(taxOverlay_unedited, mergeThis, by="Accession")
names(taxOverlay)[8] <- "ID"
taxOverlay <- subset(taxOverlay, select = -Accession)
taxOverlay <- taxOverlay %>%  #this just moves the ID column to be first in the dataframe
  select(ID, everything())


## Adding in the lichen status and three columns, crm, lff1, and lff2. 
icsmetaPath <- "E:/MiscProjects/Singh/April_FinalWork/FinalMetaTable_April2025.tsv"
icsmeta <- read.table(icsmetaPath, sep="\t", header=TRUE, stringsAsFactors=FALSE)
icsAccessions <- unique(icsmeta$Accession)
icsAccessions <- icsAccessions[icsAccessions != "None"]
icsSpecies <- unique(icsmeta$Species)
icsSpecies <- icsSpecies[icsSpecies != '']
#removing the IDs in the taxOverlay and the tree that are not in the icsmeta Species column
taxOverlay_trimmed <- taxOverlay[taxOverlay$ID %in% icsSpecies, ]

# Group and summarize
grouped_summarized <- icsmeta %>%
  group_by(Accession,Lichen, Species, ManualLabel) %>%
  summarise(Count = n(), .groups = 'drop')

# Filter out empty species or labels
grouped_summarized <- grouped_summarized %>%
  filter(Species != "", ManualLabel != "")

# Keep only the row with the highest Count for each Species + ManualLabel
grouped_summarized <- grouped_summarized %>%
  group_by(Species, ManualLabel) %>%
  slice_max(order_by = Count, n = 1, with_ties = FALSE) %>%
  ungroup()

#rename species to ID
grouped_summarized$ID <- grouped_summarized$Species
grouped_summarized_wide <- grouped_summarized %>%
  pivot_wider(names_from = ManualLabel, values_from = Count, values_fill = 0)
# keep only the columns we want
grouped_summarized_wide <- grouped_summarized_wide %>%
  select(Accession, Species, ID, Lichen, crm_split, crm_fused, lichengroup1, lichengroup2)

# for all of the grouped_summarized_wide crm_split, crm_fused, lichengroup1, and lichengroup2 columns, if the value is 0, change it to "NotPresent", if it is 1, change it to "Present"
grouped_summarized_wide <- grouped_summarized_wide %>%
  mutate(across(
    c(crm_split, crm_fused, lichengroup1, lichengroup2),
    ~ ifelse(. == 0, "NotPresent", "Present")
  ))

# Step 1: Build full list of species-accession pairs from icsmeta
icsmeta_base <- icsmeta %>%
  filter(Species != "", Accession != "None") %>%
  distinct(Accession, Species, Lichen) %>%
  mutate(ID = Species)

# Step 2: Merge with GCF-wide info, preserving full list of accessions/species
grouped_summarized_full <- left_join(icsmeta_base, grouped_summarized_wide, by = c("Accession", "Species", "ID", "Lichen"))

# Step 3: Replace any NAs in the GCF label columns with "NotPresent"
grouped_summarized_full <- grouped_summarized_full %>%
  mutate(across(c(crm_split, crm_fused, lichengroup1, lichengroup2), ~ replace_na(.x, "NotPresent")))

# matching up these Species to the IDs in the taxOverlay
taxOverlay_trimmed <- merge(taxOverlay_trimmed, grouped_summarized_full, by="ID", all.x=TRUE)
#if the species is NA, then make it the ID column value

# Reshape to long format
heatmap_data_long <- taxOverlay_trimmed %>%
  select(ID, Lichen, crm_split, crm_fused, lichengroup1, lichengroup2) %>%
  pivot_longer(
    cols = c(crm_split, crm_fused, lichengroup1, lichengroup2),
    names_to = "GCF",
    values_to = "Detected"
  )

dup_info <- icsmeta %>%
  filter(Species != "", ManualLabel %in% c("crm_split", "crm_fused")) %>%
  group_by(Species, ManualLabel) %>%
  summarise(GCF_count = n(), .groups = "drop") %>%
  pivot_wider(names_from = ManualLabel, values_from = GCF_count, values_fill = 0) %>%
  rename(ID = Species)  # Match tree tip labels

dup_info <- dup_info %>%
  mutate(border_flag = case_when(
    crm_split >= 2 | crm_fused >= 2 ~ "Duplicate",
    TRUE ~ NA_character_
  ))
# if not a Duplicate, then make it "NotDuplicate"
dup_info$border_flag[is.na(dup_info$border_flag)] <- "NotDuplicate"

taxOverlay_trimmed <- left_join(taxOverlay_trimmed, dup_info[, c("ID", "border_flag")], by = "ID")

# Combine crm_split and crm_fused into a new column
taxOverlay_trimmed$crm_combined <- case_when(
  taxOverlay_trimmed$crm_split == "Present" & taxOverlay_trimmed$crm_fused == "Present" ~ "Both",
  taxOverlay_trimmed$crm_split == "Present" ~ "crm_split",
  taxOverlay_trimmed$crm_fused == "Present" ~ "crm_fused",
  TRUE ~ "None"
)
## 

##########
# Rooting the tree with Rozella allomycis, the only sequenced Crypotmycota which is well documented to be the most ancestral fungi lineage
##########
#removing all of the tips not in the taxOverlay_trimmed table
keepThese <- c(taxOverlay_trimmed$ID, "Rozella_allomycis")
phyloTree_trim <- keep.tip(phyloTree, keepThese)
rootedTree <- root(phyloTree_trim, outgroup = "Rozella_allomycis")

##########
# Converting the tree to a ggtree object and coloring the tips based on the phylum they fall in
##########
fungalTree <- ggtree(rootedTree, layout = "rectangular", size = 0.1) 
taxOverlay_trimmed <- taxOverlay_trimmed %>%
  mutate(LichenStatus = ifelse(Lichen == "Lichen", 1, 0))

fTreeClass <- fungalTree %<+% taxOverlay_trimmed + geom_tippoint(aes(
                    color = Lichen,
                    alpha = ifelse(Lichen == "NotLichen", 1, 1),
                    size = ifelse(Lichen == "Lichen", 1, 0.5)
                  )) +
                    scale_color_manual(values = c("Lichen" = "#07f49e", "NotLichen" = "#42047e"),
                      na.value = "transparent") + # remove the legend 
                    scale_size_continuous(range = c(0.5, 1), guide = "none") +
                    scale_alpha_continuous(range = c(0.5, 1), guide = "none") 

##########
# Adding some hand-labels to the tree based on certain key taxa
##########
#Mislabeled species on NCBI: Keep adding to this list if more are detected
misLabels <- c("Penicillium_sp_occitanis", "Aspergillus_ustus")
#Converting that list to the classes they are in
classesToHighlight = c("Sordariomycetes","Dothideomycetes", "Eurotiomycetes","Saccharomycetes", "Leotiomycetes", "Agaricomycetes", "Lecanoromycetes")
#setting up the color dictionary, this def isn't the easiest way to do it but 
#I have to do it this way unless I want to re-write the code 
colorDic = hash()
colorDic[["Sordariomycetes"]] = "#0072B2"
colorDic[["Dothideomycetes"]] = "#E69F00"
colorDic[["Eurotiomycetes"]] = "#9b1d20"
colorDic[["Saccharomycetes"]] = "black"
colorDic[["Leotiomycetes"]] = "#009E73"
colorDic[["Agaricomycetes"]] = "#F0E442"
colorDic[["Lecanoromycetes"]] = "#ff99c8"

#Getting the MRCA nodes for all of the classes we want to highlight
mrcaNodes = hash()
for(i in 1:length(classesToHighlight)) {         # Head of for-loop
  classss<-classesToHighlight[i]
  idsInClass_unedited<- taxOverlay_trimmed[taxOverlay_trimmed$Fungal_Class == classss, "ID"]
  #removing the id if it is in the mislabel table
  idsInClass <- setdiff(idsInClass_unedited, misLabels)
  mrca_node = findMRCA(rootedTree, tips=idsInClass, type=c("node"))
  mrcaNodes[[classss]] <- mrca_node
}

counter = 1
for (classss in ls(mrcaNodes)) {    
  nodeToHighlight <- mrcaNodes[[classss]]
  fillColor = colorDic[[classss]]
  tempTree <- fTreeClass +
    geom_highlight(node=nodeToHighlight, fill=fillColor, alpha=0.1)
  fTreeClass <- tempTree
  counter <- counter + 1
}
#Adding in the labels!
counter = 1
for (classss in ls(mrcaNodes)) {    
  nodeToHighlight <- mrcaNodes[[classss]]
  fillColor = colorDic[[classss]]
  tempTree <- fTreeClass +
    geom_cladelabel(node=nodeToHighlight, label=classss, 
                    fill = "white", geom='label', angle="auto", horizontal=FALSE,hjust=0.5,
                    fontsize = 3, color =fillColor)
  fTreeClass <- tempTree
  counter <- counter + 1
}

##########
# Adding in two heatmaps, one showing if the tip had a lichen status and the other showing the columns crm_split, crm_fused, lff1, and lff2
##########
# First heatmap (Lichen)
fTree_Heatmap <- fTreeClass +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = ID, fill = Lichen),
    width = 0.5,
    offset = 0.1,
    pwidth = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(
    values = c("Lichen" = "#07f49e", "NotLichen" = "#42047e"),
    na.value = "transparent"
  )

# Second heatmap (crm_split)
fTree_Heatmap <- fTree_Heatmap +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = ID, fill = crm_split),
    width = 0.5,
    offset = 0.1,
    pwidth = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(
    values = c("Present" = "#ffd500", "NotPresent" = "#000000"),
    na.value = "transparent"
  )

##### printing the number of ascomycetes with a crm_split or crm_fused
### first we need to add in ALL of the ascomycetes (go back to OG tax overlay)
fullAsco <- taxOverlay %>%
  filter(Phylum == "Ascomycota") 
## we need to merge taxOverlay and taxOverlay_trimmed by the fullAsco 
columnsToMerge <- c("ID", "Phylum.x","Fungal_Class.x", "crm_split", "crm_fused")
fullAsco <- fullAsco %>%
  left_join(taxOverlay_trimmed, by = "ID") %>%
  select(all_of(columnsToMerge)) 
## rename phylum.x to phylum
fullAsco <- fullAsco %>%
  rename(Phylum = Phylum.x)
## rename Fungal_Class.x to Fungal_Class
fullAsco <- fullAsco %>%
  rename(Fungal_Class = Fungal_Class.x)
## make all na in crm_split and crm_fused to "NotPresent"
fullAsco$crm_split[is.na(fullAsco$crm_split)] <- "NotPresent"
fullAsco$crm_fused[is.na(fullAsco$crm_fused)] <- "NotPresent"
# now we can count the number of ascomycetes with a crm_split or crm_fused
fullAsco$crm_either <- ifelse(fullAsco$crm_split == "Present" | fullAsco$crm_fused == "Present", "Present", "NotPresent")
#print the num present and what percentage that is of all ascomycetes
numPresent <- fullAsco %>%
  filter(crm_either == "Present") %>%
  nrow()
numTotal <- fullAsco %>%
  nrow()
print(paste("Number of Ascomycetes with a crm_split or crm_fused: ", numPresent))
print(paste("Percentage of Ascomycetes with a crm_split or crm_fused: ", round(numPresent/numTotal*100, 2), "%"))

# remove Saccharomyces cerevisiae and recount
numPresentNoSac <- fullAsco %>%
  filter(crm_either == "Present" & Fungal_Class != "Saccharomycetes") %>%
  nrow()
numTotalNoSac <- fullAsco %>%
  filter(Fungal_Class != "Saccharomycetes") %>%
  nrow()
print(paste("Number of Ascomycetes with a crm_split or crm_fused (no Saccharomyces cerevisiae): ", numPresentNoSac))
print(paste("Percentage of Ascomycetes with a crm_split or crm_fused (no Saccharomyces cerevisiae): ", round(numPresentNoSac/numTotalNoSac*100, 2), "%"))



### Save this plot as svg
saveHere <- "E:/MiscProjects/Singh/Plots"
ggsave(
  filename = file.path(saveHere, "FungalTree_Heatmap_Lichen_crm-split.svg"),
  plot = fTree_Heatmap,
  device = "svg",
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

# Third heatmap (crm_fused)
fTree_Heatmap <- fTreeClass +
  new_scale_fill() +  # This is CRUCIAL
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = ID, fill = lichengroup1),
    width = 0.5,
    offset = 0.1,
    pwidth = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(
    values = c("Present" = "#E600D1", "NotPresent" = 'transparent'),
    na.value = "transparent"
  )
# Fourth heatmap (lff1)
fTree_Heatmap <- fTree_Heatmap +
  new_scale_fill() +  # This is CRUCIAL
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = ID, fill = lichengroup2),
    width = 0.5,
    offset = 0.1,
    pwidth = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(
    values = c("Present" = "#FF940C", "NotPresent" = "transparent"),
    na.value = "transparent"
  )
# Saving the plot
saveHere <- "E:/MiscProjects/Singh/Plots"
ggsave(
  filename = file.path(saveHere, "FungalTree_Heatmap_crms_seperated.svg"),
  plot = fTree_Heatmap,
  device = "svg",
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

# Fifth heatmap (lff2)
fTree_Heatmap <- fTreeClass +
  new_scale_fill() +  # This is CRUCIAL
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = ID, fill = lichengroup2),
    width = 0.5,
    offset = 0.1,
    pwidth = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(
    values = c("Present" = "#f8acff", "NotPresent" = "#000000"),
    na.value = "transparent"
  )
# Saving the plot
saveHere <- "E:/MiscProjects/Singh/Plots"
ggsave(
  filename = file.path(saveHere, "FungalTree_Heatmap_lff2.svg"),
  plot = fTree_Heatmap,
  device = "svg",
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

## Saving the table used with the tree
saveHere <- "E:/MiscProjects/Singh/Plots"
write.csv(taxOverlay_trimmed, file.path(saveHere, "FungalTree_Heatmap_table.csv"), row.names = FALSE)

## Doing a version with the crm fused and crm split, where it instead plots it on one heatmap with different colors. If a species has both crm fused and crm split, it will be colored a unique color
fTree_CRMCombined <- fTreeClass +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = ID, fill = crm_combined, color = border_flag),
    width = 0.5,
    offset = 0.1,
    pwidth = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(
    values = c(
      "crm_split" = "#cd0402",
      "crm_fused" = "#00b3ed",
      "Both" = "#ffde20",
      "None" = "transparent"
    ),
    na.value = "transparent"
  ) +
  scale_color_manual(
    values = c("Duplicate" = "white", "NotDuplicate" = "transparent"),
    na.value = "transparent"
  )

# Saving the plot
saveHere <- "E:/MiscProjects/Singh/Plots"
ggsave(
  filename = file.path(saveHere, "FungalTree_Heatmap_crm_combined.svg"),
  plot = fTree_CRMCombined,
  device = "svg",
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

###########
# PLOT 3: Making the phylogeny from the crm ICS domains
###########
rm(list = ls())

crm_treePath <- "E:/MiscProjects/Singh/April_FinalWork/crmTree/DIT1_PvcA_crm_domains_trimmed.faa.treefile"
crmTree <- read.tree(crm_treePath)

crmMetaPath <- "E:/MiscProjects/Singh/April_FinalWork/GetNRPS_crm/NRPS_CRM_FinalTable.tsv"
crmMeta <- read.table(crmMetaPath, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# the proteinsToRemove are already gone in the crmMeta table
proteinsToRemove = c("RFU76422.1",
                    "GKU05918.1", 'GKU06936.1', 'GKU06935.1', 'GKU06937.1',
                    'XP_016263411.1', 'XP_016263409.1'
)
crmTree <- drop.tip(crmTree, proteinsToRemove)

#removing the _ICS from every tree tip label
crmTree$tip.label <- gsub("_ICS", "", crmTree$tip.label)

## Making the column in crmMeta that matches to ID
crmMeta$ID <- paste(crmMeta$Species, crmMeta$crm_protein, sep = "_")
crmMeta <- crmMeta %>%
  select(ID, everything())

# adding the outgroup
crmMetaColumns <- colnames(crmMeta)
newRow <- as.data.frame(matrix(NA, nrow = 1, ncol = length(crmMetaColumns)))
colnames(newRow) <- crmMetaColumns
#
newRow$ID <- "Fusarium_equiseti_CAG7566245.1"
newRow$Species <- "Fusarium_equiseti"
newRow$crm_protein <- "CAG7566245.1"
newRow$Lichen <- "NotLichen"
newRow$Accession <- "GCA_910393935.1"

# Appending the new row
crmMeta <- rbind(crmMeta, newRow)
# change all the NAs to "None"
crmMeta[is.na(crmMeta)] <- "None"

## Adding in a barplot showing the distance between the NRPS and ICS if split
crmMeta$bp_separating[crmMeta$bp_separating == "None"] <- NA
crmMeta$bp_separating <- as.numeric(crmMeta$bp_separating)

## Removing any other tips that are not in the crmMeta table
missingTips <- setdiff(crmTree$tip.label, crmMeta$ID)
if (length(missingTips) > 0) {
  crmTree <- drop.tip(crmTree, missingTips)
}

## Check if the outgroup is in the tree
outgroup <- "Fusarium_equiseti_CAG7566245.1"
if (outgroup %in% crmTree$tip.label) {
  # Outgroup is present in the tree
  cat("Outgroup is present in the tree.\n")
} else {
  # Outgroup is not present in the tree
  cat("Outgroup is NOT present in the tree.\n")
}

## Adding in the column the shows if the gene directions are opposite if they are split
# comparing: crm_direction to nrps_direction
crmMeta$gene_direction <- ifelse(crmMeta$crm_direction == crmMeta$nrps_direction, "Same", "Opposite")
# if the NRPS_attached_ICS is "Yes" then I will set the gene_direction to "None"
crmMeta$gene_direction[crmMeta$NRPS_attached_ICS == "Yes"] <- "None"
crmMeta$gene_direction[crmMeta$NRPS_attached_ICS == "None"] <- "None"

## Adding in a heatmap showing the gene group for the NRPSs
nrpsGroups_path <- "E:/MiscProjects/Singh/April_FinalWork/SynthaserRun/synthaser_domainGrouping.tsv"
nrpsGroups <- read.table(nrpsGroups_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
## remove all rows where there isn't a value in GroupLetter
nrpsGroups <- nrpsGroups %>%
  filter(!is.na(GroupLetter) & GroupLetter != "")

## Change the first column to ID
colnames(nrpsGroups)[1] <- "ID"
## merging this with the crmMeta table by the ID column 
crmMeta <- left_join(crmMeta, nrpsGroups, by = "ID")

### Rooting the tree at Fusarium_equiseti_CAG7566245.1
crmTree_rooted <- root(crmTree, outgroup = "Fusarium_equiseti_CAG7566245.1", resolve.root = TRUE)
colnames(crmMeta)[colnames(crmMeta) == "NRPS_attached_ICS."] <- "NRPS_attached_ICS"

### Checking for bootstrap values in the tree
if (!is.null(crmTree_rooted$node.label)) {
  cat("Bootstrap values are present in the tree!\n")
} else {
  stop("No bootstrap values found in `crmTree$node.label`. Did you run the tree inference with support?")
}
## Converting it to a ggtree
crm_ggtree <-ggtree(crmTree_rooted, layout = "rectangular") +
  geom_nodelab(aes(subset=!is.na(label), label = label), nudge_x = -0.025,nudge_y=3.5, size = 2, color="red", geom='text', node='internal')

crm_ggtree <- crm_ggtree %<+% crmMeta + 
  geom_tippoint(aes(color = Lichen)) +
  scale_color_manual(values = c("Lichen" = "#07f49e", "NotLichen" = "#42047e", "None" = "transparent"),
                     na.value = "transparent") +
  geom_tiplab(size = 0.5, offset = 0.05) 

## The split vs fused is in the NRPS_attached_ICS. column
crm_ggtree_heatmap <- crm_ggtree +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = ID, fill = NRPS_attached_ICS),
    width = 0.1,
    offset = 0.3,
    pwidth = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(
    values = c(
      "Yes" = "#2c6e49",
      "No" = "#70e000",
      "None" = "transparent"
    ),
    na.value = "transparent"
  )

## Adding in a heatmap showing the gene directions of the nrps and ics if they are split
crm_ggtree_heatmap<- crm_ggtree_heatmap + new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y=ID, fill=gene_direction),
    width = 0.05,
    pwidth = 0.1,
    offset = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(values = c("Same" = "#023e8a", "Opposite" = "#00b4d8", "None" = "transparent"),
                    na.value = "transparent")
  
## Adding in the distance barplot with a scale on the X values
crm_ggtree_heatmap <- crm_ggtree_heatmap + 
  geom_fruit(
    geom = geom_col,
    mapping = aes(y = ID, x = bp_separating),
    orientation = "y",
    width = 1,
    offset = 0.05,
    pwidth = 0.5,  
    axis.params = list(
      axis = "x",                     # show x-axis
      text.angle = 90,                # horizontal labels
      text.size = 1,               # small axis text
    )
  )

### Adding in the group letter
crm_ggtree_heatmap<- crm_ggtree_heatmap +
  geom_fruit(
    geom = geom_text,
    mapping = aes(y = ID, label = GroupLetter),
    pwidth = 0.05,
    size = 1,  # Adjust for readability
    hjust = 0.5,  # Center the text
    fontface = "bold",
    offset = 0.1
  )

# saving this as an svg 
ggsave(
  filename = "E:/MiscProjects/Singh/Plots/crmTree_labeled.svg",
  plot = crm_ggtree_heatmap,
  width = 10,
  height = 10,
  units = "in"
)

## Adding in one more layer with colors for the same GroupLetters: letters a through l
crm_ggtree_nrpsGroups<- crm_ggtree +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = ID, fill = GroupLetter),
    width = 0.1,
    offset = 0.1,
    pwidth = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(values = c(
    "a" = "#01befe",
    "b" = "#ffdd00",
    "c" = "#ff7d00",
    "d" = "#ff006d",
    "e" = "#adff02",
    "f" = "#8f00ff",
    "g" = "black",
    "h" = "black",
    "i" = "black",
    "j" = "black",
    "k" = "black"
  ),
  na.value = "transparent")

# saving this as an svg 
ggsave(
  filename = "E:/MiscProjects/Singh/Plots/crmTree_NRPSGroups.svg",
  plot = crm_ggtree_nrpsGroups,
  width = 10,
  height = 10,
  units = "in"
)

## Saving the full crm meta table
write.table(crmMeta, file = "E:/MiscProjects/Singh/Plots/crmMeta.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


###########
# PLOT 4: Aligning the crm tree to the species tree
###########
rm(list = ls())

################# 
# GETTING THE CRM TREE READY
################# 
crm_treePath <- "E:/MiscProjects/Singh/April_FinalWork/crmTree/DIT1_PvcA_crm_domains_trimmed.faa.treefile"
crmTree <- read.tree(crm_treePath)

crmMetaPath <- "E:/MiscProjects/Singh/April_FinalWork/GetNRPS_crm/NRPS_CRM_FinalTable.tsv"
crmMeta <- read.table(crmMetaPath, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# the proteinsToRemove are already gone in the crmMeta table
proteinsToRemove = c("RFU76422.1",
                     "GKU05918.1", 'GKU06936.1', 'GKU06935.1', 'GKU06937.1',
                     'XP_016263411.1', 'XP_016263409.1'
)
crmTree <- drop.tip(crmTree, proteinsToRemove)

#removing the _ICS from every tree tip label
crmTree$tip.label <- gsub("_ICS", "", crmTree$tip.label)

## Making the column in crmMeta that matches to ID
crmMeta$ID <- paste(crmMeta$Species, crmMeta$crm_protein, sep = "_")
crmMeta <- crmMeta %>%
  select(ID, everything())

# adding the outgroup
crmMetaColumns <- colnames(crmMeta)
newRow <- as.data.frame(matrix(NA, nrow = 1, ncol = length(crmMetaColumns)))
colnames(newRow) <- crmMetaColumns
#
newRow$ID <- "Fusarium_equiseti_CAG7566245.1"
newRow$Species <- "Fusarium_equiseti"
newRow$crm_protein <- "CAG7566245.1"
newRow$Lichen <- "NotLichen"
newRow$Accession <- "GCA_910393935.1"

# Appending the new row
crmMeta <- rbind(crmMeta, newRow)
# change all the NAs to "None"
crmMeta[is.na(crmMeta)] <- "None"

## Removing any other tips that are not in the crmMeta table
missingTips <- setdiff(crmTree$tip.label, crmMeta$ID)
if (length(missingTips) > 0) {
  crmTree <- drop.tip(crmTree, missingTips)
}
colnames(crmMeta)[colnames(crmMeta) == "NRPS_attached_ICS."] <- "NRPS_attached_ICS"

### Rooting the crm tree at Fusarium_equiseti_CAG7566245.1 which is the known ICS-TauD outgroup
crmTree_rooted <- root(crmTree, outgroup = "Fusarium_equiseti_CAG7566245.1", resolve.root = TRUE)

################# 
# GETTING THE SPECIES TREE READY
################# 
# treePath = "E:/MickeyDomainCollab/ASTRALRun/FungalSpeciesTree_FigTreeExport.tre"
treePath= "E:/MickeyDomainCollab/ASTRALRun/FungalSpeciesTree.tre"
antismashOverlayFile  = "E:/MickeyDomainCollab/EditedOverlay.tsv"
taxOverlayPath = "E:/MickeyDomainCollab/TaxonomyOverlay.tsv"

phyloTree <- read.tree(treePath)
antismashOverlay <- read.csv(antismashOverlayFile, sep="\t")
taxOverlay_unedited <- read.csv(taxOverlayPath, sep="\t")

#doing the edited to the tax file
#honestly a lot of this is trial and error, but the big thing is to make the row that matches the tip labels "ID"
mergeThis <- antismashOverlay[, c("Species","Accession")]
filterAccessions <- mergeThis$Accession
taxOverlay_unedited <- taxOverlay_unedited[taxOverlay_unedited$accession %in% filterAccessions, ] #filtering down the table to include only what is in the tree
taxOverlay_unedited <- taxOverlay_unedited[!duplicated(taxOverlay_unedited$accession), ] #removing any duplicate rows
names(taxOverlay_unedited)[1] <- "Accession"
taxOverlay <- merge(taxOverlay_unedited, mergeThis, by="Accession")
names(taxOverlay)[8] <- "ID"
taxOverlay <- subset(taxOverlay, select = -Accession)
taxOverlay <- taxOverlay %>%  #this just moves the ID column to be first in the dataframe
  select(ID, everything())

#removing the IDs in the taxOverlay and the tree that are not in the Species column for the crm tree
speciesInCrmTree <- unique(crmMeta$Species)
## Adding in a single Saccharomyces cerevisiae tip as an outgroup for the species tree
speciesInCrmTree <- c(speciesInCrmTree, "Saccharomyces_cerevisiae")
#
taxOverlay_trimmed <- taxOverlay[taxOverlay$ID %in% speciesInCrmTree, ]
keepTheseTips <- taxOverlay_trimmed$ID

phyloTree_trimmed <- keep.tip(phyloTree, keepTheseTips)

### Rooting the species tree at Saccharomyces cerevisiae as it is the outgroup to the non-Saccharomyces Ascomycetes
phyloTree_rooted <- root(phyloTree_trimmed, outgroup = "Saccharomyces_cerevisiae", resolve.root = TRUE)

## there are nan values in phyloTree_rooted$edge.length that almost certainly comes from zero-length branches
# the other tree doesn't have this issue as ASTRAL and IqTree2 treet these null branch lengths differently
# for the purposes of plotting we're going to set all of the NaN edges to 0.001
phyloTree_rooted$edge.length[!is.finite(phyloTree_rooted$edge.length)] <- 0.001

## Adding the Lichen column from the crmMeta table to the taxOverlay table
addThis <- crmMeta[, c("Species", "Lichen")]
addThis <- addThis[!duplicated(addThis), ] #removing any duplicate rows
# rename species to ID and make it first in the dataframe
names(addThis)[1] <- "ID"
addThis <- addThis %>%  
  select(ID, everything())
## merging it with the taxOverlay_trimmed table
taxOverlay_trimmed <- merge(taxOverlay_trimmed, addThis, by="ID", all.x=TRUE)


#### Adding in the taxonomy to the crm from the now edited taxOverlay_trimmed table
# matching columns crmMeta$Species with taxOverlay_trimmed$ID
addTheseToCrmMeta <- taxOverlay_trimmed[, c("ID", "Kingdom", "Phylum", "Fungal_Class", "Order", "Family", "Genus")]
crmMeta <- merge(crmMeta, addTheseToCrmMeta, by.x="Species", by.y="ID", all.x=TRUE)


################# 
# ALIGNING THE TREES TOGETHER
################# 
###
### Update association matrix based on Species matching with the species tree tips
association <- cbind(crmMeta$ID, crmMeta$Species)  # Use Species from crmMeta for association
# association <- cbind(crmTree_rooted$tip.label, crmTree_rooted$tip.label)

### Creating the col to add in colors to the links based on fused, split and lichen
# Assign colors based on NRPS_attached_ICS column
fused_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "Yes"]
split_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "No"]

# Find where those IDs appear in the association matrix
searchresult_fused <- which(association[,1] %in% fused_ids)
searchresult_split <- which(association[,1] %in% split_ids)

# Default color for all lines
col <- rep(make.transparent("#000000", 0.2), nrow(association))

# Color the lines based on the NRPS_attached_ICS values (fused = red, split = green)
col[c(searchresult_fused)] <- make.transparent("#ed1c6e", 0.7)
col[c(searchresult_split)] <- make.transparent("#7ec242", 0.7)

# Now, override the color for any "Lichen" species
lichen_species <- crmMeta$Species[crmMeta$Lichen == "Lichen"]
searchresult_lichen <- which(association[, 2] %in% lichen_species)
col[searchresult_lichen] <- "#4cc9f0"  # Set Lichen-related links to this color

# ## Save the crmMeta table and the taxOverlay table
# write.csv(crmMeta, "E:/MiscProjects/Singh/Plots/crmPlot/crmMeta.csv", row.names = FALSE)
# write.csv(taxOverlay_trimmed, "E:/MiscProjects/Singh/Plots/crmPlot/taxOverlay.csv", row.names = FALSE)

# lichen_species <- crmMeta$Species[crmMeta$Lichen == "Lichen"]
# # Find where those lichen species appear in the association matrix and change their color
# searchresult_lichen <- which(association[,1] %in% lichen_species)
# 
# # Overwrite their colors to the Lichen color
# col[searchresult_lichen] <- "#4cc9f0"  # Lichen color

#### Running the cophylo command
crm_speciesAligned<- cophylo(crmTree_rooted,phyloTree_rooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association)

# colorDic[["Sordariomycetes"]] = "#0072B2"
# colorDic[["Dothideomycetes"]] = "#E69F00"
# colorDic[["Eurotiomycetes"]] = "#9b1d20"
# colorDic[["Saccharomycetes"]] = "black"
# colorDic[["Leotiomycetes"]] = "#009E73"
# colorDic[["Agaricomycetes"]] = "#F0E442"
# colorDic[["Lecanoromycetes"]] = "#ff99c8"

### Adding in colors based on taxonomy for the species tree
class_colors <- c(
  "Dothideomycetes" = "#E69F00",
  "Sordariomycetes" = "#0072B2",
  "Eurotiomycetes" = "#9b1d20",
  "Leotiomycetes" = "#009E73",
  "Lecanoromycetes" = "#ff99c8"
)

### Setting the branch colors for the species tree
right <- rep("#000000", nrow(crm_speciesAligned$trees[[2]]$edge))  # default black

# Loop over each target class
for (target_class in names(class_colors)) {
  # Get tips for this class
  tips_to_color <- taxOverlay_trimmed$ID[taxOverlay_trimmed$Fungal_Class == target_class]
  tips_to_color <- tips_to_color[tips_to_color %in% crm_speciesAligned$trees[[2]]$tip.label]
  
  if (length(tips_to_color) < 2) next  # Need at least 2 tips for MRCA
  
  mrca_node <- getMRCA(crm_speciesAligned$trees[[2]], tips_to_color)
  descendants <- getDescendants(crm_speciesAligned$trees[[2]], mrca_node)
  
  edge_indices <- sapply(descendants, function(x, y) which(y == x), y = crm_speciesAligned$trees[[2]]$edge[,2])
  right[edge_indices] <- class_colors[[target_class]]
}

edge.col <- list(
  left = rep("black", nrow(crm_speciesAligned$trees[[1]]$edge)),  # or color the left too if needed
  right = right
)

## Comment this out if you want it printed instead of saved as an svg
outputFilePath <- file.path("E:/MiscProjects/Singh/Plots/cophylo_crm_species_labeled.svg")
# outputFilePath <- file.path("E:/MiscProjects/Singh/Plots/cophylo_crm_species.svg")
svglite(outputFilePath,
        width = 10, height = 8)

###### PLOTTING
#with tips
plot(crm_speciesAligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, edge.col = edge.col,
     show.tip.label=TRUE, fsize=0.1)
# without tips
# plot(crm_speciesAligned,link.type="curved", link.lty="solid",
#      hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
#      link.col=col, edge.col = edge.col,
#      show.tip.label=FALSE, ftype="off")

### Add class-based tip color on LEFT tree
tip_class_colors <- c(
  "Dothideomycetes" = "#E69F00",
  "Sordariomycetes" = "#0072B2",
  "Eurotiomycetes" = "#9b1d20",
  "Leotiomycetes" = "#009E73",
  "Lecanoromycetes" = "#ff99c8"
)

# Get a named vector of class by species name
tip_classes <- crmMeta$Fungal_Class
names(tip_classes) <- crmMeta$ID

# Initialize empty vector for tip colors
tip_colors_left <- rep("#000000", length(crm_speciesAligned$trees[[1]]$tip.label))

# Match species → color if class is in defined set
for (i in seq_along(tip_colors_left)) {
  tip_name <- crm_speciesAligned$trees[[1]]$tip.label[i]
  fungal_class <- tip_classes[tip_name]
  if (!is.na(fungal_class) && fungal_class %in% names(tip_class_colors)) {
    tip_colors_left[i] <- tip_class_colors[fungal_class]
  }
}
# ✨ Add colored tip points to left tree
tiplabels.cophylo(
  pch = 21,         # shape = circle
  pie = NULL,
  bg = tip_colors_left,  # fill color
  which = "left",
  cex = 1       # adjust size to taste
)

### Add class-based tip color on RIGHT tree
lichenColors <- c(
  "Lichen" = "#452979",
  "NotLichen" = "#69c28d"
)

tip_class_colors_right <- taxOverlay_trimmed$Lichen
names(tip_class_colors_right) <- taxOverlay_trimmed$ID

tip_colors_right <- rep("black", length(crm_speciesAligned$trees[[2]]$tip.label))
for (i in seq_along(crm_speciesAligned$trees[[2]]$tip.label)) {
  tip <- crm_speciesAligned$trees[[2]]$tip.label[i]
  class <- tip_class_colors_right[tip]
  if (!is.na(class) && class %in% names(lichenColors)) {
    tip_colors_right[i] <- lichenColors[[class]]
  }
}
# ✨ Add colored tip points to  the tree on the right
tiplabels.cophylo(
  pch = 21,         # shape = circle
  pie = NULL,
  bg = tip_colors_right,  # fill color
  which = "right",
  cex = 1       # adjust size to taste
)

# ### Adding in pies for the bootstrap values
# # Get node numbers for both trees
# nodes_left <- 1:crm_speciesAligned$trees[[1]]$Nnode + Ntip(crm_speciesAligned$trees[[1]])
# nodes_right <- 1:crm_speciesAligned$trees[[2]]$Nnode + Ntip(crm_speciesAligned$trees[[2]])
# 
# # Convert node labels to numeric (or assign 0 if NA)
# labels_left <- as.numeric(crm_speciesAligned$trees[[1]]$node.label)
# labels_right <- as.numeric(crm_speciesAligned$trees[[2]]$node.label)
# 
# # Add bootstrap pies to left tree
# nodelabels.cophylo(
#   node = nodes_left,
#   pie = cbind(labels_left, 100 - labels_left),
#   piecol = c("#118ab2", "white"),
#   cex = 0.2,
#   which = "left"
# )
# # Add bootstrap pies to right tree
# nodelabels.cophylo(
#   node = nodes_right,
#   pie = cbind(labels_right, 1 - labels_right),
#   piecol = c("#ffd166", "white"),
#   cex = 0.2,
#   which = "right"
# )

dev.off()


###########
# PLOT 5: Aligning the crm tree to the species tree; VERSION WITH ALL FUSED REMOVED
###########
rm(list = ls())

################# 
# GETTING THE CRM TREE READY
################# 
crm_treePath <- "E:/MiscProjects/Singh/April_FinalWork/crmTree/DIT1_PvcA_crm_domains_trimmed.faa.treefile"
crmTree <- read.tree(crm_treePath)

crmMetaPath <- "E:/MiscProjects/Singh/April_FinalWork/GetNRPS_crm/NRPS_CRM_FinalTable.tsv"
crmMeta <- read.table(crmMetaPath, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# the proteinsToRemove are already gone in the crmMeta table
proteinsToRemove = c("RFU76422.1",
                     "GKU05918.1", 'GKU06936.1', 'GKU06935.1', 'GKU06937.1',
                     'XP_016263411.1', 'XP_016263409.1'
)
crmTree <- drop.tip(crmTree, proteinsToRemove)

## removing all fused proteins from the crm tree and table
crmMeta <- crmMeta[(crmMeta$NRPS_attached_ICS. == "Yes"), ]

#removing the _ICS from every tree tip label
crmTree$tip.label <- gsub("_ICS", "", crmTree$tip.label)

## Making the column in crmMeta that matches to ID
crmMeta$ID <- paste(crmMeta$Species, crmMeta$crm_protein, sep = "_")
crmMeta <- crmMeta %>%
  select(ID, everything())

# adding the outgroup
crmMetaColumns <- colnames(crmMeta)
newRow <- as.data.frame(matrix(NA, nrow = 1, ncol = length(crmMetaColumns)))
colnames(newRow) <- crmMetaColumns
#
newRow$ID <- "Fusarium_equiseti_CAG7566245.1"
newRow$Species <- "Fusarium_equiseti"
newRow$crm_protein <- "CAG7566245.1"
newRow$Lichen <- "NotLichen"
newRow$Accession <- "GCA_910393935.1"

# Appending the new row
crmMeta <- rbind(crmMeta, newRow)
# change all the NAs to "None"
crmMeta[is.na(crmMeta)] <- "None"

## Removing any other tips that are not in the crmMeta table
missingTips <- setdiff(crmTree$tip.label, crmMeta$ID)
if (length(missingTips) > 0) {
  crmTree <- drop.tip(crmTree, missingTips)
}
colnames(crmMeta)[colnames(crmMeta) == "NRPS_attached_ICS."] <- "NRPS_attached_ICS"

### Rooting the crm tree at Fusarium_equiseti_CAG7566245.1 which is the known ICS-TauD outgroup
crmTree_rooted <- root(crmTree, outgroup = "Fusarium_equiseti_CAG7566245.1", resolve.root = TRUE)

################# 
# GETTING THE SPECIES TREE READY
################# 
# treePath = "E:/MickeyDomainCollab/ASTRALRun/FungalSpeciesTree_FigTreeExport.tre"
treePath= "E:/MickeyDomainCollab/ASTRALRun/FungalSpeciesTree.tre"
antismashOverlayFile  = "E:/MickeyDomainCollab/EditedOverlay.tsv"
taxOverlayPath = "E:/MickeyDomainCollab/TaxonomyOverlay.tsv"

phyloTree <- read.tree(treePath)
antismashOverlay <- read.csv(antismashOverlayFile, sep="\t")
taxOverlay_unedited <- read.csv(taxOverlayPath, sep="\t")

#doing the edited to the tax file
#honestly a lot of this is trial and error, but the big thing is to make the row that matches the tip labels "ID"
mergeThis <- antismashOverlay[, c("Species","Accession")]
filterAccessions <- mergeThis$Accession
taxOverlay_unedited <- taxOverlay_unedited[taxOverlay_unedited$accession %in% filterAccessions, ] #filtering down the table to include only what is in the tree
taxOverlay_unedited <- taxOverlay_unedited[!duplicated(taxOverlay_unedited$accession), ] #removing any duplicate rows
names(taxOverlay_unedited)[1] <- "Accession"
taxOverlay <- merge(taxOverlay_unedited, mergeThis, by="Accession")
names(taxOverlay)[8] <- "ID"
taxOverlay <- subset(taxOverlay, select = -Accession)
taxOverlay <- taxOverlay %>%  #this just moves the ID column to be first in the dataframe
  select(ID, everything())

#removing the IDs in the taxOverlay and the tree that are not in the Species column for the crm tree
speciesInCrmTree <- unique(crmMeta$Species)
## Adding in a single Saccharomyces cerevisiae tip as an outgroup for the species tree
speciesInCrmTree <- c(speciesInCrmTree, "Saccharomyces_cerevisiae")
#
taxOverlay_trimmed <- taxOverlay[taxOverlay$ID %in% speciesInCrmTree, ]
keepTheseTips <- taxOverlay_trimmed$ID

phyloTree_trimmed <- keep.tip(phyloTree, keepTheseTips)

### Rooting the species tree at Saccharomyces cerevisiae as it is the outgroup to the non-Saccharomyces Ascomycetes
phyloTree_rooted <- root(phyloTree_trimmed, outgroup = "Saccharomyces_cerevisiae", resolve.root = TRUE)

## there are nan values in phyloTree_rooted$edge.length that almost certainly comes from zero-length branches
# the other tree doesn't have this issue as ASTRAL and IqTree2 treet these null branch lengths differently
# for the purposes of plotting we're going to set all of the NaN edges to 0.001
phyloTree_rooted$edge.length[!is.finite(phyloTree_rooted$edge.length)] <- 0.001

## Adding the Lichen column from the crmMeta table to the taxOverlay table
addThis <- crmMeta[, c("Species", "Lichen")]
addThis <- addThis[!duplicated(addThis), ] #removing any duplicate rows
# rename species to ID and make it first in the dataframe
names(addThis)[1] <- "ID"
addThis <- addThis %>%  
  select(ID, everything())
## merging it with the taxOverlay_trimmed table
taxOverlay_trimmed <- merge(taxOverlay_trimmed, addThis, by="ID", all.x=TRUE)


#### Adding in the taxonomy to the crm from the now edited taxOverlay_trimmed table
# matching columns crmMeta$Species with taxOverlay_trimmed$ID
addTheseToCrmMeta <- taxOverlay_trimmed[, c("ID", "Kingdom", "Phylum", "Fungal_Class", "Order", "Family", "Genus")]
crmMeta <- merge(crmMeta, addTheseToCrmMeta, by.x="Species", by.y="ID", all.x=TRUE)


################# 
# ALIGNING THE TREES TOGETHER
################# 
# #### last step, we need to filter any tips with no color in the association, these are ICS proteis with no NRPS (edge case)
# no_color_indices_crmTree <- crmMeta[crmMeta$NRPS_attached_ICS=="None", ]$ID
# no_color_indices_speicesTree <- crmMeta[crmMeta$NRPS_attached_ICS=="None", ]$Species
# 
# ## drop those tips from each tree
# crmTree_rooted <- drop.tip(crmTree_rooted, no_color_indices_crmTree)
# phyloTree_rooted <- drop.tip(phyloTree_rooted, no_color_indices_speicesTree)
# ## drop in the table also
# crmMeta <- crmMeta[!crmMeta$ID %in% no_color_indices_crmTree, ]

#### sanity check, print if there are any tip.labels in each tree that is not in the crmMeta
# crm tree match is ID column
# species tree match is Speceis column
print("Checking if all tips in crmTree_rooted are in crmMeta")
if (all(crmTree_rooted$tip.label %in% crmMeta$ID)) {
  print("All tips in crmTree_rooted are in crmMeta")
} else {
  # print which tips are not in crmMeta
  print(crmTree_rooted$tip.label[!crmTree_rooted$tip.label %in% crmMeta$ID])
}
print("Checking if all tips in phyloTree_rooted are in crmMeta")
if (all(phyloTree_rooted$tip.label %in% crmMeta$Species)) {
  print("All tips in phyloTree_rooted are in crmMeta")
} else {
  # print which tips are not in crmMeta
  print(phyloTree_rooted$tip.label[!phyloTree_rooted$tip.label %in% crmMeta$Species])
}

###
### Update association matrix based on Species matching with the species tree tips
association <- cbind(crmMeta$ID, crmMeta$Species)  # Use Species from crmMeta for association
## drop the Fusarium_equiseti_CAG7566245.1 from the association matrix
association <- association[!association[,1] %in% c("Fusarium_equiseti_CAG7566245.1"), ]
# association <- cbind(crmTree_rooted$tip.label, crmTree_rooted$tip.label)

### Creating the col to add in colors to the links based on fused, split and lichen
# Assign colors based on NRPS_attached_ICS column
fused_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "Yes"]
split_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "No"]

# Find where those IDs appear in the association matrix
searchresult_fused <- which(association[,1] %in% fused_ids)
searchresult_split <- which(association[,1] %in% split_ids)

# Default color for all lines
col <- rep(make.transparent("#000000", 0.2), nrow(association))

# Color the lines based on the NRPS_attached_ICS values (fused = red, split = green)
col[c(searchresult_fused)] <- make.transparent("#ed1c6e", 0.7)
col[c(searchresult_split)] <- make.transparent("#7ec242", 0.7)

# Now, override the color for any "Lichen" species
lichen_species <- crmMeta$Species[crmMeta$Lichen == "Lichen"]
searchresult_lichen <- which(association[, 2] %in% lichen_species)
col[searchresult_lichen] <- "#4cc9f0"  # Set Lichen-related links to this color



#### Running the cophylo command
crm_speciesAligned<- cophylo(crmTree_rooted,phyloTree_rooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association)

### Adding in colors based on taxonomy for the species tree
class_colors <- c(
  "Dothideomycetes" = "#E69F00",
  "Sordariomycetes" = "#0072B2",
  "Eurotiomycetes" = "#9b1d20",
  "Leotiomycetes" = "#009E73",
  "Lecanoromycetes" = "#ff99c8"
)

### Setting the branch colors for the species tree
right <- rep("#000000", nrow(crm_speciesAligned$trees[[2]]$edge))  # default black

# Loop over each target class
for (target_class in names(class_colors)) {
  # Get tips for this class
  tips_to_color <- taxOverlay_trimmed$ID[taxOverlay_trimmed$Fungal_Class == target_class]
  tips_to_color <- tips_to_color[tips_to_color %in% crm_speciesAligned$trees[[2]]$tip.label]
  
  if (length(tips_to_color) < 2) next  # Need at least 2 tips for MRCA
  
  mrca_node <- getMRCA(crm_speciesAligned$trees[[2]], tips_to_color)
  descendants <- getDescendants(crm_speciesAligned$trees[[2]], mrca_node)
  
  edge_indices <- sapply(descendants, function(x, y) which(y == x), y = crm_speciesAligned$trees[[2]]$edge[,2])
  right[edge_indices] <- class_colors[[target_class]]
}

edge.col <- list(
  left = rep("black", nrow(crm_speciesAligned$trees[[1]]$edge)),  # or color the left too if needed
  right = right
)

## Comment this out if you want it printed instead of saved as an svg
# outputFilePath <- file.path("E:/MiscProjects/Singh/Plots/NOFUSED_cophylo_crm_species_labeled.svg")
outputFilePath <- file.path("E:/MiscProjects/Singh/Plots/NOFUSED_cophylo_crm_species.svg")
svglite(outputFilePath,
        width = 10, height = 8)

###### PLOTTING
# #with tips
# plot(crm_speciesAligned,link.type="curved", link.lty="solid",
#      hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
#      link.col=col, edge.col = edge.col,
#      show.tip.label=TRUE, fsize=0.1)
# without tips
plot(crm_speciesAligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, edge.col = edge.col,
     show.tip.label=FALSE, ftype="off")




### Add class-based tip color on LEFT tree
tip_class_colors <- c(
  "Dothideomycetes" = "#E69F00",
  "Sordariomycetes" = "#0072B2",
  "Eurotiomycetes" = "#9b1d20",
  "Leotiomycetes" = "#009E73",
  "Lecanoromycetes" = "#ff99c8"
)

# Get a named vector of class by species name
tip_classes <- crmMeta$Fungal_Class
names(tip_classes) <- crmMeta$ID

# Initialize empty vector for tip colors
tip_colors_left <- rep("#000000", length(crm_speciesAligned$trees[[1]]$tip.label))

# Match species → color if class is in defined set
for (i in seq_along(tip_colors_left)) {
  tip_name <- crm_speciesAligned$trees[[1]]$tip.label[i]
  fungal_class <- tip_classes[tip_name]
  if (!is.na(fungal_class) && fungal_class %in% names(tip_class_colors)) {
    tip_colors_left[i] <- tip_class_colors[fungal_class]
  }
}
# ✨ Add colored tip points to left tree
tiplabels.cophylo(
  pch = 21,         # shape = circle
  pie = NULL,
  bg = tip_colors_left,  # fill color
  which = "left",
  cex = 1       # adjust size to taste
)

### Add class-based tip color on RIGHT tree
lichenColors <- c(
  "Lichen" = "#452979",
  "NotLichen" = "#69c28d"
)

tip_class_colors_right <- taxOverlay_trimmed$Lichen
names(tip_class_colors_right) <- taxOverlay_trimmed$ID

tip_colors_right <- rep("black", length(crm_speciesAligned$trees[[2]]$tip.label))
for (i in seq_along(crm_speciesAligned$trees[[2]]$tip.label)) {
  tip <- crm_speciesAligned$trees[[2]]$tip.label[i]
  class <- tip_class_colors_right[tip]
  if (!is.na(class) && class %in% names(lichenColors)) {
    tip_colors_right[i] <- lichenColors[[class]]
  }
}
# ✨ Add colored tip points to  the tree on the right
tiplabels.cophylo(
  pch = 21,         # shape = circle
  pie = NULL,
  bg = tip_colors_right,  # fill color
  which = "right",
  cex = 1       # adjust size to taste
)

dev.off()

###########
# PLOT 6: Aligning the crm tree to the nrps tree
###########
rm(list = ls())

################# 
# GETTING THE CRM TREE READY
################# 
crm_treePath <- "E:/MiscProjects/Singh/April_FinalWork/crmTree/DIT1_PvcA_crm_domains_trimmed.faa.treefile"
crmTree <- read.tree(crm_treePath)

crmMetaPath <- "E:/MiscProjects/Singh/April_FinalWork/GetNRPS_crm/NRPS_CRM_FinalTable.tsv"
crmMeta <- read.table(crmMetaPath, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# the proteinsToRemove are already gone in the crmMeta table
proteinsToRemove = c("RFU76422.1",
                     "GKU05918.1", 'GKU06936.1', 'GKU06935.1', 'GKU06937.1',
                     'XP_016263411.1', 'XP_016263409.1'
)
crmTree <- drop.tip(crmTree, proteinsToRemove)

#removing the _ICS from every tree tip label
crmTree$tip.label <- gsub("_ICS", "", crmTree$tip.label)

## Drop all the rows where NRPS_attched_ICS is None
crmMeta <- crmMeta[!crmMeta$NRPS_attached_ICS. == "None", ]

## Making the column in crmMeta that matches to ID
crmMeta$ID <- paste(crmMeta$Species, crmMeta$crm_protein, sep = "_")
crmMeta <- crmMeta %>%
  select(ID, everything())

# adding the outgroup
crmMetaColumns <- colnames(crmMeta)
newRow <- as.data.frame(matrix(NA, nrow = 1, ncol = length(crmMetaColumns)))
colnames(newRow) <- crmMetaColumns
#
newRow$ID <- "Fusarium_equiseti_CAG7566245.1"
newRow$Species <- "Fusarium_equiseti"
newRow$crm_protein <- "CAG7566245.1"
newRow$Lichen <- "NotLichen"
newRow$Accession <- "GCA_910393935.1"

# Appending the new row
crmMeta <- rbind(crmMeta, newRow)
# change all the NAs to "None"
crmMeta[is.na(crmMeta)] <- "None"

## Removing any other tips that are not in the crmMeta table
missingTips <- setdiff(crmTree$tip.label, crmMeta$ID)
if (length(missingTips) > 0) {
  crmTree <- drop.tip(crmTree, missingTips)
}
colnames(crmMeta)[colnames(crmMeta) == "NRPS_attached_ICS."] <- "NRPS_attached_ICS"

### Rooting the crm tree at Fusarium_equiseti_CAG7566245.1 which is the known ICS-TauD outgroup
# crmTree_rooted <- root(crmTree, outgroup = "Fusarium_equiseti_CAG7566245.1", resolve.root = TRUE)
## Version where it's not rooted, but we copy over the tree with that variable name to make the code work
crmTree_rooted <- crmTree

################# 
# GETTING THE NRPS TREE READY
################# 
nrpsTree_Path <- "E:/MiscProjects/Singh/April_FinalWork/nrpsTree/nrps_crm_domains_trim.faa.treefile"
nrpsTree <- read.tree(nrpsTree_Path)

nrpsTreeTips <- nrpsTree$tip.label

# Function to handle protein name matching
find_protein_match <- function(tipName, crmMeta) {
  # Split the tip name by "_"
  name_parts <- unlist(strsplit(tipName, "_"))
  # First attempt: Grab the second-to-last element
  proteinName_first_attempt <- name_parts[length(name_parts) - 1]
  # Check if the first attempt matches crm or nrps protein exactly
  matchRow <- crmMeta[(crmMeta$crm_protein == proteinName_first_attempt | crmMeta$nrps_protein == proteinName_first_attempt), ]
  # If there is exactly 1 match, return it
  if (nrow(matchRow) == 1) {
    return(matchRow)
  }
  # If no match or more than one, grab the third and second-to-last parts
  if (length(name_parts) >= 3) {
    proteinName_second_attempt <- paste(name_parts[(length(name_parts) - 2):(length(name_parts) - 1)], collapse = "_")
    # Try to match this combined name in the crm and nrps columns
    matchRow <- crmMeta[(crmMeta$crm_protein == proteinName_second_attempt | crmMeta$nrps_protein == proteinName_second_attempt), ]
    if (nrow(matchRow) == 1) {
      return(matchRow)
    } else {
      # If still no match or multiple matches, print the tipName for manual inspection
      print(paste("Manual inspection needed for tip:", tipName))
    }
  }
  return(NULL)  # Return NULL if no valid match found
}

crmMeta$NRPS_TreeMatch <- NA  # Initialize the NRPS_TreeMatch column
# Now let's loop through the tree tips and perform the matching
for (i in 1:length(nrpsTreeTips)) {
  tipName <- nrpsTreeTips[i]
  
  # Find the matching row in crmMeta based on the logic
  matchRow <- find_protein_match(tipName, crmMeta)
  
  # If we found a valid match (not NULL), assign the tipName to the matching row's NRPS_TreeMatch
  if (!is.null(matchRow)) {
    crmMeta[crmMeta$crm_protein == matchRow$crm_protein | crmMeta$nrps_protein == matchRow$nrps_protein, "NRPS_TreeMatch"] <- tipName
  }
}

## if there are any rows that still don't have a value for NRPS_TreeMatch, 
## get those IDs, remove them from crmMeta, and remove them from the crmTree_rooted
missingMatch <- crmMeta[is.na(crmMeta$NRPS_TreeMatch), ]
if (nrow(missingMatch) > 0) {
  missingIDs <- missingMatch$ID
  # Remove these IDs from crmMeta
  crmMeta <- crmMeta[!crmMeta$ID %in% missingIDs, ]
  # Remove these tips from the tree
  crmTree_rooted <- drop.tip(crmTree_rooted, missingIDs)
}

################# 
# ALIGNING THE TREES TOGETHER
################# 
###
### Update association matrix based on Species matching with the species tree tips
association <- cbind(crmMeta$ID, crmMeta$NRPS_TreeMatch)  # Use Species from crmMeta for association

### Creating the col to add in colors to the links based on fused, split and lichen
# Assign colors based on NRPS_attached_ICS column
fused_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "Yes"]
split_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "No"]

# Find where those IDs appear in the association matrix
searchresult_fused <- which(association[,1] %in% fused_ids)
searchresult_split <- which(association[,1] %in% split_ids)

# Default color for all lines
col <- rep(make.transparent("#000000", 0.2), nrow(association))

# Color the lines based on the NRPS_attached_ICS values (fused = red, split = green)
col[c(searchresult_fused)] <- make.transparent("#ed1c6e", 0.7)
col[c(searchresult_split)] <- make.transparent("#7ec242", 0.7)

# Now, override the color for any "Lichen" species
lichen_species <- crmMeta$ID[crmMeta$Lichen == "Lichen"]
searchresult_lichen <- which(association[, 1] %in% lichen_species)
col[searchresult_lichen] <- "#4cc9f0"  # Set Lichen-related links to this color

#### Running the cophylo command
crm_nrps_Aligned<- cophylo(crmTree_rooted,nrpsTree, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association)


############### Getting the RF distance for the same trees
# Create copies to avoid modifying original trees
crmTree_copy <- crmTree_rooted
nrpsTree_copy <- nrpsTree
## For this we will rename the NRPS tree tips to match the crm tree using the crmMeta$ID and crmMeta$NRPS_TreeMatch columns
for (i in 1:nrow(crmMeta)) {
  # Get the ID and NRPS_TreeMatch for the current row
  id <- crmMeta$ID[i]
  nrps_match <- crmMeta$NRPS_TreeMatch[i]
  
  # Rename the tip in the nrpsTree_copy to match the ID in crmTree_copy
  if (!is.na(nrps_match)) {
    nrpsTree_copy$tip.label[nrpsTree_copy$tip.label == nrps_match] <- id
  }
}
### remove any tips in NRPS tree not in the crm tree, and vice versa
nrpsTree_copy <- drop.tip(nrpsTree_copy, setdiff(nrpsTree_copy$tip.label, crmTree_copy$tip.label))
crmTree_copy <- drop.tip(crmTree_copy, setdiff(crmTree_copy$tip.label, nrpsTree_copy$tip.label))


# --- Step 2: Calculate Robinson-Foulds (RF) Distance ---
trees_for_rf <- list(crmTree_copy, nrpsTree_copy)
class(trees_for_rf) <- "multiPhylo"
rf_matrix <- multiRF(trees_for_rf, quiet = TRUE)
rf_distance <- rf_matrix[1, 2]

num_tips <- length(crmTree_copy$tip.label) # Get number of tips from one tree
max_rf_distance <- 0 # Initialize
if (num_tips >= 3) {
  # Note: This assumes the trees are fully resolved (bifurcating).
  # RF distance is based on differing bipartitions (internal edges).
  # Each tree has n-3 internal edges. Max diff = (n-3) + (n-3) = 2n - 6
  max_rf_distance <- 2 * (num_tips - 3)
}

# Print the RF distance and its maximum possible value
print(paste0("Robinson-Foulds (RF) Distance: ", rf_distance,
             " (Maximum possible: ", max_rf_distance, ")"))


# --- Step 3: Calculate Branch Score (Kuhner-Felsenstein) Distance ---

# Use the trees prepared in the previous block (crmTree_copy, nrpsTree_copy)
# KF.dist calculates the branch score distance, which incorporates branch lengths
library(phangorn)

bs_distance <- tryCatch(
  KF.dist(crmTree_copy, nrpsTree_copy), # Use the trees with matched labels
  error = function(e) {
    warning("Could not calculate Branch Score distance: ", e$message)
    return(NA) # Return NA if calculation fails
  }
)

# Calculate the sum of squared branch lengths for each tree for context
# This represents the distance if the trees shared zero splits.
sum_sq_len_crm <- sum(crmTree_copy$edge.length^2)
sum_sq_len_nrps <- sum(nrpsTree_copy$edge.length^2)
total_sum_sq_len <- sum_sq_len_crm + sum_sq_len_nrps

# Print the Branch Score distance and the total squared length
print(paste0("Branch Score (KF) Distance (Topology & Branch Lengths): ", bs_distance,
             " (Total sum of squared branch lengths in both trees: ", round(total_sum_sq_len, 4), ")"))




## Comment this out if you want it printed instead of saved as an svg
# outputFilePath <- file.path("E:/MiscProjects/Singh/Plots/NOFUSED_cophylo_crm_species_labeled.svg")
outputFilePath <- file.path("E:/MiscProjects/Singh/Plots/NOFUSED_cophylo_crm_nrps.svg")
svglite(outputFilePath,
        width = 10, height = 8)

###### PLOTTING
# #with tips
# plot(crm_speciesAligned,link.type="curved", link.lty="solid",
#      hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
#      link.col=col, show.tip.label=TRUE, fsize=0.1)
# without tips
plot(crm_nrps_Aligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, show.tip.label=FALSE, ftype="off")

 ### Add class-based tip color on LEFT tree
tip_class_colors <- c(
  "Dothideomycetes" = "#E69F00",
  "Sordariomycetes" = "#0072B2",
  "Eurotiomycetes" = "#9b1d20",
  "Leotiomycetes" = "#009E73",
  "Lecanoromycetes" = "#ff99c8"
)

# Get a named vector of class by species name
tip_classes <- crmMeta$Fungal_Class
names(tip_classes) <- crmMeta$ID

# Initialize empty vector for tip colors
tip_colors_left <- rep("#000000", length(crm_speciesAligned$trees[[1]]$tip.label))

# Match species → color if class is in defined set
for (i in seq_along(tip_colors_left)) {
  tip_name <- crm_speciesAligned$trees[[1]]$tip.label[i]
  fungal_class <- tip_classes[tip_name]
  if (!is.na(fungal_class) && fungal_class %in% names(tip_class_colors)) {
    tip_colors_left[i] <- tip_class_colors[fungal_class]
  }
}
# ✨ Add colored tip points to left tree
tiplabels.cophylo(
  pch = 21,         # shape = circle
  pie = NULL,
  bg = tip_colors_left,  # fill color
  which = "left",
  cex = 1       # adjust size to taste
)

### Add class-based tip color on RIGHT tree
lichenColors <- c(
  "Lichen" = "#452979",
  "NotLichen" = "#69c28d"
)

tip_class_colors_right <- taxOverlay_trimmed$Lichen
names(tip_class_colors_right) <- taxOverlay_trimmed$ID

tip_colors_right <- rep("black", length(crm_speciesAligned$trees[[2]]$tip.label))
for (i in seq_along(crm_speciesAligned$trees[[2]]$tip.label)) {
  tip <- crm_speciesAligned$trees[[2]]$tip.label[i]
  class <- tip_class_colors_right[tip]
  if (!is.na(class) && class %in% names(lichenColors)) {
    tip_colors_right[i] <- lichenColors[[class]]
  }
}
# ✨ Add colored tip points to  the tree on the right
tiplabels.cophylo(
  pch = 21,         # shape = circle
  pie = NULL,
  bg = tip_colors_right,  # fill color
  which = "right",
  cex = 1       # adjust size to taste
)

dev.off()


###########
# PLOT 7: Aligning the crm tree to the species tree; VERSION WITH ALL SPLIT REMOVED
###########
rm(list = ls())

################# 
# GETTING THE CRM TREE READY
################# 
crm_treePath <- "E:/MiscProjects/Singh/April_FinalWork/crmTree/DIT1_PvcA_crm_domains_trimmed.faa.treefile"
crmTree <- read.tree(crm_treePath)

crmMetaPath <- "E:/MiscProjects/Singh/April_FinalWork/GetNRPS_crm/NRPS_CRM_FinalTable.tsv"
crmMeta <- read.table(crmMetaPath, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# the proteinsToRemove are already gone in the crmMeta table
proteinsToRemove = c("RFU76422.1",
                     "GKU05918.1", 'GKU06936.1', 'GKU06935.1', 'GKU06937.1',
                     'XP_016263411.1', 'XP_016263409.1'
)
crmTree <- drop.tip(crmTree, proteinsToRemove)

# Keeping only the fused copies
crmMeta <- crmMeta[(crmMeta$NRPS_attached_ICS. == "No"), ]
## remove all the 

#removing the _ICS from every tree tip label
crmTree$tip.label <- gsub("_ICS", "", crmTree$tip.label)

## Making the column in crmMeta that matches to ID
crmMeta$ID <- paste(crmMeta$Species, crmMeta$crm_protein, sep = "_")
crmMeta <- crmMeta %>%
  select(ID, everything())

# adding the outgroup
crmMetaColumns <- colnames(crmMeta)
newRow <- as.data.frame(matrix(NA, nrow = 1, ncol = length(crmMetaColumns)))
colnames(newRow) <- crmMetaColumns
#
newRow$ID <- "Fusarium_equiseti_CAG7566245.1"
newRow$Species <- "Fusarium_equiseti"
newRow$crm_protein <- "CAG7566245.1"
newRow$Lichen <- "NotLichen"
newRow$Accession <- "GCA_910393935.1"

# Appending the new row
crmMeta <- rbind(crmMeta, newRow)
# change all the NAs to "None"
crmMeta[is.na(crmMeta)] <- "None"

## Removing any other tips that are not in the crmMeta table
missingTips <- setdiff(crmTree$tip.label, crmMeta$ID)
if (length(missingTips) > 0) {
  crmTree <- drop.tip(crmTree, missingTips)
}
colnames(crmMeta)[colnames(crmMeta) == "NRPS_attached_ICS."] <- "NRPS_attached_ICS"

### Rooting the crm tree at Fusarium_equiseti_CAG7566245.1 which is the known ICS-TauD outgroup
crmTree_rooted <- root(crmTree, outgroup = "Fusarium_equiseti_CAG7566245.1", resolve.root = TRUE)

################# 
# GETTING THE SPECIES TREE READY
################# 
# treePath = "E:/MickeyDomainCollab/ASTRALRun/FungalSpeciesTree_FigTreeExport.tre"
treePath= "E:/MickeyDomainCollab/ASTRALRun/FungalSpeciesTree.tre"
antismashOverlayFile  = "E:/MickeyDomainCollab/EditedOverlay.tsv"
taxOverlayPath = "E:/MickeyDomainCollab/TaxonomyOverlay.tsv"

phyloTree <- read.tree(treePath)
antismashOverlay <- read.csv(antismashOverlayFile, sep="\t")
taxOverlay_unedited <- read.csv(taxOverlayPath, sep="\t")

#doing the edited to the tax file
#honestly a lot of this is trial and error, but the big thing is to make the row that matches the tip labels "ID"
mergeThis <- antismashOverlay[, c("Species","Accession")]
filterAccessions <- mergeThis$Accession
taxOverlay_unedited <- taxOverlay_unedited[taxOverlay_unedited$accession %in% filterAccessions, ] #filtering down the table to include only what is in the tree
taxOverlay_unedited <- taxOverlay_unedited[!duplicated(taxOverlay_unedited$accession), ] #removing any duplicate rows
names(taxOverlay_unedited)[1] <- "Accession"
taxOverlay <- merge(taxOverlay_unedited, mergeThis, by="Accession")
names(taxOverlay)[8] <- "ID"
taxOverlay <- subset(taxOverlay, select = -Accession)
taxOverlay <- taxOverlay %>%  #this just moves the ID column to be first in the dataframe
  select(ID, everything())

#removing the IDs in the taxOverlay and the tree that are not in the Species column for the crm tree
speciesInCrmTree <- unique(crmMeta$Species)
## Adding in a single Saccharomyces cerevisiae tip as an outgroup for the species tree
speciesInCrmTree <- c(speciesInCrmTree, "Saccharomyces_cerevisiae")
#
taxOverlay_trimmed <- taxOverlay[taxOverlay$ID %in% speciesInCrmTree, ]
keepTheseTips <- taxOverlay_trimmed$ID

phyloTree_trimmed <- keep.tip(phyloTree, keepTheseTips)

### Rooting the species tree at Saccharomyces cerevisiae as it is the outgroup to the non-Saccharomyces Ascomycetes
phyloTree_rooted <- root(phyloTree_trimmed, outgroup = "Saccharomyces_cerevisiae", resolve.root = TRUE)

## there are nan values in phyloTree_rooted$edge.length that almost certainly comes from zero-length branches
# the other tree doesn't have this issue as ASTRAL and IqTree2 treet these null branch lengths differently
# for the purposes of plotting we're going to set all of the NaN edges to 0.001
phyloTree_rooted$edge.length[!is.finite(phyloTree_rooted$edge.length)] <- 0.001

## Adding the Lichen column from the crmMeta table to the taxOverlay table
addThis <- crmMeta[, c("Species", "Lichen")]
addThis <- addThis[!duplicated(addThis), ] #removing any duplicate rows
# rename species to ID and make it first in the dataframe
names(addThis)[1] <- "ID"
addThis <- addThis %>%  
  select(ID, everything())
## merging it with the taxOverlay_trimmed table
taxOverlay_trimmed <- merge(taxOverlay_trimmed, addThis, by="ID", all.x=TRUE)


#### Adding in the taxonomy to the crm from the now edited taxOverlay_trimmed table
# matching columns crmMeta$Species with taxOverlay_trimmed$ID
addTheseToCrmMeta <- taxOverlay_trimmed[, c("ID", "Kingdom", "Phylum", "Fungal_Class", "Order", "Family", "Genus")]
crmMeta <- merge(crmMeta, addTheseToCrmMeta, by.x="Species", by.y="ID", all.x=TRUE)


################# 
# ALIGNING THE TREES TOGETHER
################# 
###
### Update association matrix based on Species matching with the species tree tips
association <- cbind(crmMeta$ID, crmMeta$Species)  # Use Species from crmMeta for association
association <- association[!association[,1] %in% c("Fusarium_equiseti_CAG7566245.1"), ]
#
# association <- cbind(crmTree_rooted$tip.label, crmTree_rooted$tip.label)

### Creating the col to add in colors to the links based on fused, split and lichen
# Assign colors based on NRPS_attached_ICS column
fused_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "Yes"]
split_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "No"]

# Find where those IDs appear in the association matrix
searchresult_fused <- which(association[,1] %in% fused_ids)
searchresult_split <- which(association[,1] %in% split_ids)

# Default color for all lines
col <- rep(make.transparent("#000000", 0.2), nrow(association))

# Color the lines based on the NRPS_attached_ICS values (fused = red, split = green)
col[c(searchresult_fused)] <- make.transparent("#ed1c6e", 0.7)
col[c(searchresult_split)] <- make.transparent("#7ec242", 0.7)

# Now, override the color for any "Lichen" species
lichen_species <- crmMeta$Species[crmMeta$Lichen == "Lichen"]
searchresult_lichen <- which(association[, 2] %in% lichen_species)
col[searchresult_lichen] <- "#4cc9f0"  # Set Lichen-related links to this color

#### Running the cophylo command
crm_speciesAligned<- cophylo(crmTree_rooted,phyloTree_rooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association)

### Adding in colors based on taxonomy for the species tree
class_colors <- c(
  "Dothideomycetes" = "#E69F00",
  "Sordariomycetes" = "#0072B2",
  "Eurotiomycetes" = "#9b1d20",
  "Leotiomycetes" = "#009E73",
  "Lecanoromycetes" = "#ff99c8"
)

### Setting the branch colors for the species tree
right <- rep("#000000", nrow(crm_speciesAligned$trees[[2]]$edge))  # default black

# Loop over each target class
for (target_class in names(class_colors)) {
  # Get tips for this class
  tips_to_color <- taxOverlay_trimmed$ID[taxOverlay_trimmed$Fungal_Class == target_class]
  tips_to_color <- tips_to_color[tips_to_color %in% crm_speciesAligned$trees[[2]]$tip.label]
  
  if (length(tips_to_color) < 2) next  # Need at least 2 tips for MRCA
  
  mrca_node <- getMRCA(crm_speciesAligned$trees[[2]], tips_to_color)
  descendants <- getDescendants(crm_speciesAligned$trees[[2]], mrca_node)
  
  edge_indices <- sapply(descendants, function(x, y) which(y == x), y = crm_speciesAligned$trees[[2]]$edge[,2])
  right[edge_indices] <- class_colors[[target_class]]
}

edge.col <- list(
  left = rep("black", nrow(crm_speciesAligned$trees[[1]]$edge)),  # or color the left too if needed
  right = right
)

## Comment this out if you want it printed instead of saved as an svg
# outputFilePath <- file.path("E:/MiscProjects/Singh/Plots/FUSED-ONLY_cophylo_crm_species_labeled.svg")
outputFilePath <- file.path("E:/MiscProjects/Singh/Plots/FUSED-ONLY_cophylo_crm_species.svg")
svglite(outputFilePath,
        width = 10, height = 8)

###### PLOTTING
# #with tips
# plot(crm_speciesAligned,link.type="curved", link.lty="solid",
#      hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
#      link.col=col, edge.col = edge.col,
#      show.tip.label=TRUE, fsize=0.1)
# without tips
plot(crm_speciesAligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, edge.col = edge.col,
     show.tip.label=FALSE, ftype="off")

### Add class-based tip color on LEFT tree
tip_class_colors <- c(
  "Dothideomycetes" = "#E69F00",
  "Sordariomycetes" = "#0072B2",
  "Eurotiomycetes" = "#9b1d20",
  "Leotiomycetes" = "#009E73",
  "Lecanoromycetes" = "#ff99c8"
)

# Get a named vector of class by species name
tip_classes <- crmMeta$Fungal_Class
names(tip_classes) <- crmMeta$ID

# Initialize empty vector for tip colors
tip_colors_left <- rep("#000000", length(crm_speciesAligned$trees[[1]]$tip.label))

# Match species → color if class is in defined set
for (i in seq_along(tip_colors_left)) {
  tip_name <- crm_speciesAligned$trees[[1]]$tip.label[i]
  fungal_class <- tip_classes[tip_name]
  if (!is.na(fungal_class) && fungal_class %in% names(tip_class_colors)) {
    tip_colors_left[i] <- tip_class_colors[fungal_class]
  }
}
# ✨ Add colored tip points to left tree
tiplabels.cophylo(
  pch = 21,         # shape = circle
  pie = NULL,
  bg = tip_colors_left,  # fill color
  which = "left",
  cex = 1       # adjust size to taste
)

### Add class-based tip color on RIGHT tree
lichenColors <- c(
  "Lichen" = "#452979",
  "NotLichen" = "#69c28d"
)

tip_class_colors_right <- taxOverlay_trimmed$Lichen
names(tip_class_colors_right) <- taxOverlay_trimmed$ID

tip_colors_right <- rep("black", length(crm_speciesAligned$trees[[2]]$tip.label))
for (i in seq_along(crm_speciesAligned$trees[[2]]$tip.label)) {
  tip <- crm_speciesAligned$trees[[2]]$tip.label[i]
  class <- tip_class_colors_right[tip]
  if (!is.na(class) && class %in% names(lichenColors)) {
    tip_colors_right[i] <- lichenColors[[class]]
  }
}
# ✨ Add colored tip points to  the tree on the right
tiplabels.cophylo(
  pch = 21,         # shape = circle
  pie = NULL,
  bg = tip_colors_right,  # fill color
  which = "right",
  cex = 1       # adjust size to taste
)

dev.off()

###########
# Plotting the NRPS tree in isolation
###########
rm(list=ls())

nrpsTree <- read.tree("E:/MiscProjects/Singh/April_FinalWork/nrpsTree/nrps_aln_trimmed.faa.treefile")

crmMetaPath <- "E:/MiscProjects/Singh/April_FinalWork/GetNRPS_crm/NRPS_CRM_FinalTable.tsv"
crmMeta <- read.table(crmMetaPath, sep="\t", header=TRUE, stringsAsFactors=FALSE)
crmMeta$ID <- ""
crmMeta <- crmMeta[, c("ID", setdiff(names(crmMeta), "ID"))]


crmMetaColumns <- colnames(crmMeta)
newRow <- as.data.frame(matrix(NA, nrow = 1, ncol = length(crmMetaColumns)))
colnames(newRow) <- crmMetaColumns
#
newRow$ID <- "XP_748589.2_NRPS-like-enzyme"
newRow$Species <- "Aspergillus_fumigatus"
newRow$nrps_protein <- "XP_748589.2"
newRow$Lichen <- "NotLichen"
newRow$Accession <- "GCF_000002655.1"

# Appending the new row
crmMeta <- rbind(crmMeta, newRow)
# change all the NAs to "None"
crmMeta[is.na(crmMeta)] <- "None"

# root the tree at the outgroup
nrpsTree <- root(nrpsTree, outgroup = "XP_748589.2_NRPS-like-enzyme", resolve.root = TRUE)

## Making the column in crmMeta that matches to ID
nrpsTreeTips <- nrpsTree$tip.label

# Function to handle protein name matching
find_protein_match <- function(tipName, crmMeta) {
  # Split the tip name by "_"
  name_parts <- unlist(strsplit(tipName, "_"))
  # First attempt: Grab the second-to-last element
  proteinName_first_attempt <- name_parts[length(name_parts) - 1]
  # Check if the first attempt matches crm or nrps protein exactly
  matchRow <- crmMeta[(crmMeta$crm_protein == proteinName_first_attempt | crmMeta$nrps_protein == proteinName_first_attempt), ]
  # If there is exactly 1 match, return it
  if (nrow(matchRow) == 1) {
    return(matchRow)
  }
  # If no match or more than one, grab the third and second-to-last parts
  if (length(name_parts) >= 3) {
    proteinName_second_attempt <- paste(name_parts[(length(name_parts) - 2):(length(name_parts) - 1)], collapse = "_")
    # Try to match this combined name in the crm and nrps columns
    matchRow <- crmMeta[(crmMeta$crm_protein == proteinName_second_attempt | crmMeta$nrps_protein == proteinName_second_attempt), ]
    if (nrow(matchRow) == 1) {
      return(matchRow)
    } else {
      # If still no match or multiple matches, print the tipName for manual inspection
      print(paste("Manual inspection needed for tip:", tipName))
    }
  }
  return(NULL)  # Return NULL if no valid match found
}

# Now let's loop through the tree tips and perform the matching
for (i in 1:length(nrpsTreeTips)) {
  tipName <- nrpsTreeTips[i]
  
  # Find the matching row in crmMeta based on the logic
  matchRow <- find_protein_match(tipName, crmMeta)
  
  # If we found a valid match (not NULL), assign the tipName to the matching row's ID
  if (!is.null(matchRow)) {
    crmMeta[crmMeta$crm_protein == matchRow$crm_protein | crmMeta$nrps_protein == matchRow$nrps_protein, "ID"] <- tipName
  }
}

## if there are any rows that still don't have a value for ID
crmMeta <- crmMeta[!is.na(crmMeta$ID), ]

## Adding in a barplot showing the distance between the NRPS and ICS if split
crmMeta$bp_separating[crmMeta$bp_separating == "None"] <- NA
crmMeta$bp_separating <- as.numeric(crmMeta$bp_separating)

## Removing any other tips that are not in the crmMeta table
missingTips <- setdiff(nrpsTree$tip.label, crmMeta$ID)
if (length(missingTips) > 0) {
  nrpsTree <- drop.tip(nrpsTree, missingTips)
}

## Adding in the column the shows if the gene directions are opposite if they are split
# comparing: crm_direction to nrps_direction
crmMeta$gene_direction <- ifelse(crmMeta$crm_direction == crmMeta$nrps_direction, "Same", "Opposite")
# if the NRPS_attached_ICS is "Yes" then I will set the gene_direction to "None"
crmMeta$gene_direction[crmMeta$NRPS_attached_ICS == "Yes"] <- "None"
crmMeta$gene_direction[crmMeta$NRPS_attached_ICS == "None"] <- "None"

## Adding in a heatmap showing the gene group for the NRPSs
nrpsGroups_path <- "E:/MiscProjects/Singh/April_FinalWork/SynthaserRun/synthaser_domainGrouping.tsv"
nrpsGroups <- read.table(nrpsGroups_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
## remove all rows where there isn't a value in GroupLetter
nrpsGroups <- nrpsGroups %>%
  filter(!is.na(GroupLetter) & GroupLetter != "")

## Change the first column to ID
colnames(nrpsGroups)[1] <- "crmID_Match"
## We need to make the crmID_Match column 
crmMeta$crmID_Match <- paste(crmMeta$Species, crmMeta$crm_protein, sep="_")

## merging this with the crmMeta table by the ID column 
crmMeta <- left_join(crmMeta, nrpsGroups, by = "crmID_Match")

### Checking for bootstrap values in the tree
if (!is.null(nrpsTree$node.label)) {
  cat("Bootstrap values are present in the tree!\n")
} else {
  stop("No bootstrap values found in `crmTree$node.label`. Did you run the tree inference with support?")
}
## Converting it to a ggtree
nrps_ggtree <-ggtree(nrpsTree, layout = "rectangular") +
  geom_nodelab(aes(subset=!is.na(label), label = label), nudge_x = -0.025,nudge_y=3.5, size = 2, color="red", geom='text', node='internal')

# rename NRPS_attached_ICS. to NRPS_attached_ICS
colnames(crmMeta)[colnames(crmMeta) == "NRPS_attached_ICS."] <- "NRPS_attached_ICS"

nrps_ggtree <- nrps_ggtree %<+% crmMeta + 
  geom_tippoint(aes(color = Lichen)) +
  scale_color_manual(values = c("Lichen" = "#07f49e", "NotLichen" = "#42047e", "None" = "transparent"),
                     na.value = "transparent") +
  geom_tiplab(size = 0.5, offset = 0.05) 

## The split vs fused is in the NRPS_attached_ICS column
nrps_ggtree_heatmap <- nrps_ggtree +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = ID, fill = NRPS_attached_ICS),
    width = 0.1,
    offset = 0.3,
    pwidth = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(
    values = c(
      "Yes" = "#2c6e49",
      "No" = "#70e000",
      "None" = "transparent"
    ),
    na.value = "transparent"
  )

## Adding in a heatmap showing the gene directions of the nrps and ics if they are split
nrps_ggtree_heatmap <- nrps_ggtree_heatmap + new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y=ID, fill=gene_direction),
    width = 0.05,
    pwidth = 0.1,
    offset = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(values = c("Same" = "#611b0e", "Opposite" = "#f04f6e", "None" = "transparent"),
                    na.value = "transparent")

## Adding in the distance barplot with a scale on the X values
nrps_ggtree_heatmap <- nrps_ggtree_heatmap + 
  geom_fruit(
    geom = geom_col,
    mapping = aes(y = ID, x = bp_separating),
    orientation = "y",
    width = 1,
    offset = 0.05,
    pwidth = 0.5,  
    axis.params = list(
      axis = "x",                     # show x-axis
      text.angle = 90,                # horizontal labels
      text.size = 1,               # small axis text
    )
  )

# saving this as an svg 
ggsave(
  filename = "E:/MiscProjects/Singh/Plots/NRPSTree.svg",
  plot = nrps_ggtree_heatmap,
  width = 10,
  height = 10,
  units = "in"
)


## Adding in one more layer with colors for the same GroupLetters: letters a through l
nrps_ggtree_nrpsGroups<- nrps_ggtree +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = ID, fill = GroupLetter),
    width = 0.1,
    offset = 0.1,
    pwidth = 0.1,
    axis.params = list(axis = "none")
  ) +
  scale_fill_manual(values = c(
    "a" = "#01befe",
    "b" = "#ffdd00",
    "c" = "#ff7d00",
    "d" = "#ff006d",
    "e" = "#adff02",
    "f" = "#8f00ff",
    "g" = "black",
    "h" = "black",
    "i" = "black",
    "j" = "black",
    "k" = "black"
  ),
  na.value = "transparent")

# saving this as an svg 
ggsave(
  filename = "E:/MiscProjects/Singh/Plots/NRPSTree_NRPSGroups.svg",
  plot = nrps_ggtree_nrpsGroups,
  width = 10,
  height = 10,
  units = "in"
)


library(ggplot2)
library(dplyr)
library(readr)
library(forcats) # For reordering factor levels if needed
library(stringr) # For str_detect

#######
# Metabolite Measurement bar plot
#######
rm(list=ls())

# --- 1. Load Data ---
# Define the path to your data file
metaboliteMeasurments_path <- "E:/MiscProjects/Singh/Fusarium data-ICS/20240428_LCMS analysis_fusarium.tsv"

# Read the tab-separated file using base R's read.csv
metaboliteMeasurments <- read.csv(metaboliteMeasurments_path, sep = '\t')

# Optional: Inspect the data structure
print("Data structure:")
str(metaboliteMeasurments)
print("First few rows:")
print(head(metaboliteMeasurments))

# --- 2. Data Processing ---
# Calculate mean intensity and standard deviation for each Strain and Metabolite
metabolite_summary <- metaboliteMeasurments %>%
  group_by(Strain, SearchMetabolite) %>%
  summarise(
    mean_intensity = mean(Intensity, na.rm = TRUE), # Calculate mean, remove NAs if any
    sd_intensity = sd(Intensity, na.rm = TRUE),     # Calculate standard deviation
    n = n(),                                     # Count number of replicates
    se_intensity = sd_intensity / sqrt(n)        # Calculate standard error (optional)
  ) %>%
  ungroup() # Ungroup for further operations

# Define the desired order for the strains on the x-axis
# *** IMPORTANT: Replace the placeholder names below with your actual strain names ***
strain_order <- c("PH-1 (WT)", "KLB34 (_ICS)", "KLB36 (_ICS)", "KLB44 (_NRPS)", "KLB45(_NRPS)") # Example order, replace this!

# Create a new column for Strain Type based on strain name patterns
metabolite_summary <- metabolite_summary %>%
  mutate(
    StrainType = case_when(
      str_detect(Strain, fixed("(WT)")) ~ "WT",        # Detects "(WT)"
      str_detect(Strain, fixed("(_ICS)")) ~ "ICS_Deletion", # Detects "(_ICS)"
      str_detect(Strain, fixed("(_NRPS)")) ~ "NRPS_Deletion", # Detects "(_NRPS)"
      TRUE ~ "Other" # Default category if none of the above match
    ),
    # Convert Strain to an ordered factor using the defined order
    Strain = factor(Strain, levels = strain_order),
    # Convert StrainType to a factor for consistent coloring/legend
    StrainType = factor(StrainType, levels = c("WT", "ICS_Deletion", "NRPS_Deletion", "Other"))
  )

print("Summary statistics (with StrainType and ordered Strain factor):")
print(metabolite_summary)

# --- 3. Create the Plot ---

# Define custom colors for the Strain Types
# *** Adjust colors as desired ***
strain_type_colors <- c(
  "WT" = "grey50",
  "ICS_Deletion" = "steelblue",
  "NRPS_Deletion" = "darkorange",
  "Other" = "lightgrey" # Color for any strains not matching patterns
)

# Create the ggplot object
# Map fill aesthetic to the new StrainType column
metabolite_plot <- ggplot(metabolite_summary, aes(x = Strain, y = mean_intensity, fill = StrainType)) +
  
  # Add bars representing the mean intensity
  geom_bar(stat = "identity", position = position_dodge(), color="black") + # Added black outline to bars
  
  # Add error bars (using standard deviation here, change ymin/ymax for SE if preferred)
  geom_errorbar(
    aes(ymin = mean_intensity - sd_intensity, ymax = mean_intensity + sd_intensity),
    width = 0.2, # Width of the error bar caps
    position = position_dodge(0.9)
  ) +
  
  # Create separate plots for each metabolite, side-by-side
  facet_wrap(~ SearchMetabolite, scales = "free_y", nrow = 1) + # nrow=1 ensures side-by-side
  
  # Apply the custom manual color scale
  scale_fill_manual(values = strain_type_colors) +
  
  # Customize labels and title
  labs(
    title = "Metabolite Intensity by Strain",
    x = "Strain",
    y = "Mean Intensity (+/- SD)", # Update if using SE
    fill = "Strain Type" # Update legend title (though legend is hidden below)
  ) +
  
  # Apply a theme (optional, theme_minimal is clean)
  theme_minimal(base_size = 14) + # Increased base font size
  
  # Further theme customizations (optional)
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # Rotate x-axis labels if needed
    strip.text = element_text(face = "bold"), # Make facet titles bold
    legend.position = "none", # Hide legend (remove this line if you want to see the legend)
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    panel.spacing = unit(1.5, "lines") # Add spacing between facets
  ) +
  
  # Optional: Use scale_fill_brewer or scale_fill_manual for custom colors
  # scale_fill_brewer(palette = "Set2")
  scale_y_continuous(labels = scales::scientific) # Format y-axis in scientific notation

# --- 4. Display the Plot ---
print(metabolite_plot)

# --- 5. Save the Plot (Optional) ---
# ggsave("metabolite_comparison_plot.png", plot = metabolite_plot, width = 12, height = 6, dpi = 300)
saveHere <- "E:/MiscProjects/Singh/Plots/FusariumMetabolites.svg"
ggsave(saveHere, plot = metabolite_plot, width = 12, height = 6) # Save as SVG for scalability

######
######
######
######
###### RESPONSE TO REVIEWERS PLOT: Confirming if every version of the trees agree on split and fused difference
# we also might align them to each other
rm(list=ls())

metaData = read.csv("/Volumes/T7/MiscProjects/Singh/April_FinalWork/ExtractLFFs_BGCs/LFF3_4/ICS_Meta_withTaxonomy.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
crmMetaPath <- "/Volumes/T7/MiscProjects/Singh/April_FinalWork/GetNRPS_crm/NRPS_CRM_FinalTable.tsv"
crmMeta <- read.table(crmMetaPath, sep="\t", header=TRUE, stringsAsFactors=FALSE)

ics_gappyout = read.tree("/Volumes/T7/MiscProjects/Singh/April_FinalWork/crmTree/Tree_ReMake/Gappyout/DIT1_PvcA_crm_domains_gappyout.fas.treefile")
ics_gt0_1 = read.tree("/Volumes/T7/MiscProjects/Singh/April_FinalWork/crmTree/Tree_ReMake/GT0-1/DIT1_PvcA_crm_domains_gt0-1.fas.treefile")
ics_gt0_9 = read.tree("/Volumes/T7/MiscProjects/Singh/April_FinalWork/crmTree/Tree_ReMake/GT0-9/DIT1_PvcA_crm_domains_gt0-9.fas.treefile")

nrps_gappyout = read.tree("/Volumes/T7/MiscProjects/Singh/April_FinalWork/nrpsTree/TreeRemake/Gappyout/nrps_domains_gappyout.fas.treefile")
nrps_gt0_1 = read.tree("/Volumes/T7/MiscProjects/Singh/April_FinalWork/nrpsTree/TreeRemake/GT0-1/nrps_domains_gt0-1.fas.treefile")
nrps_gt0_9 = read.tree("/Volumes/T7/MiscProjects/Singh/April_FinalWork/nrpsTree/TreeRemake/GT0-9/nrps_domains_gt0-9.fas.treefile")

################# 
# GETTING THE CRM TREE READY
################# 
# the proteinsToRemove are already gone in the metaData table
proteinsToRemove = c("RFU76422.1",
                     "GKU05918.1", 'GKU06936.1', 'GKU06935.1', 'GKU06937.1',
                     'XP_016263411.1', 'XP_016263409.1'
)
ics_gappyout <- drop.tip(ics_gappyout, proteinsToRemove)
ics_gt0_1 <- drop.tip(ics_gt0_1, proteinsToRemove)
ics_gt0_9 <- drop.tip(ics_gt0_9, proteinsToRemove)

#removing the _ICS from every tree tip label
ics_gappyout$tip.label <- gsub("_ICS", "", ics_gappyout$tip.label)
ics_gt0_1$tip.label <- gsub("_ICS", "", ics_gt0_1$tip.label)
ics_gt0_9$tip.label <- gsub("_ICS", "", ics_gt0_9$tip.label)

## Drop all the rows where NRPS_attched_ICS is None
crmMeta <- crmMeta[!crmMeta$NRPS_attached_ICS. == "None", ]

## Making the column in crmMeta that matches to ID
crmMeta$ID <- paste(crmMeta$Species, crmMeta$crm_protein, sep = "_")
crmMeta <- crmMeta %>%
  select(ID, everything())

# adding the outgroup
crmMetaColumns <- colnames(crmMeta)
newRow <- as.data.frame(matrix(NA, nrow = 1, ncol = length(crmMetaColumns)))
colnames(newRow) <- crmMetaColumns
#
newRow$ID <- "Fusarium_equiseti_CAG7566245.1"
newRow$Species <- "Fusarium_equiseti"
newRow$crm_protein <- "CAG7566245.1"
newRow$Lichen <- "NotLichen"
newRow$Accession <- "GCA_910393935.1"

# Appending the new row
crmMeta <- rbind(crmMeta, newRow)
# change all the NAs to "None"
crmMeta[is.na(crmMeta)] <- "None"

## Removing any other tips that are not in the crmMeta table
missingTips_gappyout <- setdiff(ics_gappyout$tip.label, crmMeta$ID)
missingTips_gt0_1 <- setdiff(ics_gt0_1$tip.label, crmMeta$ID)
missingTips_gt0_9 <- setdiff(ics_gt0_9$tip.label, crmMeta$ID)

if (length(missingTips_gappyout) > 0) {
  ics_gappyout <- drop.tip(ics_gappyout, missingTips_gappyout)
}
if (length(missingTips_gt0_1) > 0) {
  ics_gt0_1 <- drop.tip(ics_gt0_1, missingTips_gt0_1)
}
if (length(missingTips_gt0_9) > 0) {
  ics_gt0_9 <- drop.tip(ics_gt0_9, missingTips_gt0_9)
}
colnames(crmMeta)[colnames(crmMeta) == "NRPS_attached_ICS."] <- "NRPS_attached_ICS"

### Rooting the crm tree at Fusarium_equiseti_CAG7566245.1 which is the known ICS-TauD outgroup
ics_gappyout_rooted <- root(ics_gappyout, outgroup = "Fusarium_equiseti_CAG7566245.1", resolve.root = TRUE)
ics_gt0_1_rooted <- root(ics_gt0_1, outgroup = "Fusarium_equiseti_CAG7566245.1", resolve.root = TRUE)
ics_gt0_9_rooted <- root(ics_gt0_9, outgroup = "Fusarium_equiseti_CAG7566245.1", resolve.root = TRUE)

################# 
# GETTING THE NRPS TREE READY
################# 

nrps_gappout_labels <- nrps_gappyout$tip.label
nrps_gt0_1_labels <- nrps_gt0_1$tip.label
nrps_gt0_9_labels <- nrps_gt0_9$tip.label

## Check if each list is identical
if (identical(sort(nrps_gappout_labels), sort(nrps_gt0_1_labels))) {
  print("The tip labels of nrps_gappyout and nrps_gt0_1 match.")
}
## They are identical so we only need to do this once

# Function to handle protein name matching
find_protein_match <- function(tipName, crmMeta) {
  # Split the tip name by "_"
  name_parts <- unlist(strsplit(tipName, "_"))
  # First attempt: Grab the second-to-last element
  proteinName_first_attempt <- name_parts[length(name_parts) - 1]
  # Check if the first attempt matches crm or nrps protein exactly
  matchRow <- crmMeta[(crmMeta$crm_protein == proteinName_first_attempt | crmMeta$nrps_protein == proteinName_first_attempt), ]
  # If there is exactly 1 match, return it
  if (nrow(matchRow) == 1) {
    return(matchRow)
  }
  # If no match or more than one, grab the third and second-to-last parts
  if (length(name_parts) >= 3) {
    proteinName_second_attempt <- paste(name_parts[(length(name_parts) - 2):(length(name_parts) - 1)], collapse = "_")
    # Try to match this combined name in the crm and nrps columns
    matchRow <- crmMeta[(crmMeta$crm_protein == proteinName_second_attempt | crmMeta$nrps_protein == proteinName_second_attempt), ]
    if (nrow(matchRow) == 1) {
      return(matchRow)
    } else {
      # If still no match or multiple matches, print the tipName for manual inspection
      print(paste("Manual inspection needed for tip:", tipName))
    }
  }
  return(NULL)  # Return NULL if no valid match found
}

crmMeta$NRPS_TreeMatch <- NA  # Initialize the NRPS_TreeMatch column
# Now let's loop through the tree tips and perform the matching
for (i in 1:length(nrps_gappout_labels)) {
  tipName <- nrps_gappout_labels[i]
  
  # Find the matching row in crmMeta based on the logic
  matchRow <- find_protein_match(tipName, crmMeta)
  
  # If we found a valid match (not NULL), assign the tipName to the matching row's NRPS_TreeMatch
  if (!is.null(matchRow)) {
    crmMeta[crmMeta$crm_protein == matchRow$crm_protein | crmMeta$nrps_protein == matchRow$nrps_protein, "NRPS_TreeMatch"] <- tipName
  }
}
## if there are any rows that still don't have a value for NRPS_TreeMatch, 
## get those IDs, remove them from crmMeta, and remove them from the crmTree_rooted
missingMatch <- crmMeta[is.na(crmMeta$NRPS_TreeMatch), ]
if (nrow(missingMatch) > 0) {
  missingIDs <- missingMatch$ID
  # Remove these IDs from the three ics trees
  crmMeta <- crmMeta[!crmMeta$ID %in% missingIDs, ]
  # Remove these tips from the tree
  ics_gappyout <- drop.tip(ics_gappyout, missingIDs)
  ics_gappyout_rooted <- drop.tip(ics_gappyout_rooted, missingIDs)
  ics_gt0_1 <- drop.tip(ics_gt0_1, missingIDs)
  ics_gt0_1_rooted <- drop.tip(ics_gt0_1_rooted, missingIDs)
  ics_gt0_9 <- drop.tip(ics_gt0_9, missingIDs)
  ics_gt0_9_rooted <- drop.tip(ics_gt0_9_rooted, missingIDs)
}

################# 
# ALIGNING THE TREES TOGETHER ~~~
################# 
# unroot all of the trees
ics_gappyout_unrooted <- unroot(ics_gappyout)
ics_gt0_1_unrooted <- unroot(ics_gt0_1)
ics_gt0_9_unrooted <- unroot(ics_gt0_9)
nrps_gappyout_unrooted <- unroot(nrps_gappyout)
nrps_gt0_1_unrooted <- unroot(nrps_gt0_1)
nrps_gt0_9_unrooted <- unroot(nrps_gt0_9)

############### Update association matrix based on Species matching with the species tree tips
association_ics <- cbind(crmMeta$ID, crmMeta$ID)  # Use Species from crmMeta for association

### Creating the col to add in colors to the links based on fused, split and lichen
# Assign colors based on NRPS_attached_ICS column
fused_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "Yes"]
split_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "No"]

# Find where those IDs appear in the association matrix
searchresult_fused <- which(association_ics[,1] %in% fused_ids)
searchresult_split <- which(association_ics[,1] %in% split_ids)

# Default color for all lines
col <- rep(make.transparent("#000000", 0.2), nrow(association_ics))

# Color the lines based on the NRPS_attached_ICS values (fused = red, split = green)
col[c(searchresult_fused)] <- make.transparent("#ed1c6e", 0.7)
col[c(searchresult_split)] <- make.transparent("#7ec242", 0.7)

# Now, override the color for any "Lichen" species
lichen_species <- crmMeta$ID[crmMeta$Lichen == "Lichen"]
searchresult_lichen <- which(association_ics[, 1] %in% lichen_species)
col[searchresult_lichen] <- "#4cc9f0"  # Set Lichen-related links to this color

###
# First the ICS compared to each other (unrooted only)
###
ics_gap_gt01_Aligned<- cophylo(ics_gappyout_unrooted,ics_gt0_1_unrooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association_ics)
ics_gap_gt09_Aligned<- cophylo(ics_gappyout_unrooted,ics_gt0_9_unrooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association_ics)
ics_gt01_gt09_Aligned<- cophylo(ics_gt0_1_unrooted,ics_gt0_9_unrooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association_ics)
##
# Save the plots as svgs
saveHere <- "/Volumes/T7/MiscProjects/Singh/April_FinalWork/ICS_NRPS_TreeComparisons"
## ICS gap and gt0.1
svglite(paste0(saveHere,"/ics-gap_gt01-aligned.svg"),
        width = 7, height = 6)
plot(ics_gap_gt01_Aligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, show.tip.label=FALSE, ftype="off")
dev.off()
## ICS gap and gt0.9
svglite(paste0(saveHere,"/ics-gap_gt09-aligned.svg"),
        width = 7, height = 6)
plot(ics_gap_gt09_Aligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, show.tip.label=FALSE, ftype="off")
dev.off()
## ICS gt0.1 and gt0.9
svglite(paste0(saveHere,"/ics-gt01_gt09-aligned.svg"),
        width = 7, height = 6)
plot(ics_gt01_gt09_Aligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, show.tip.label=FALSE, ftype="off")
dev.off()
##
############### Update association matrix based on Species matching with the species tree tips
association_nrps <- cbind(crmMeta$ID, crmMeta$NRPS_TreeMatch)  # Use Species from crmMeta for association

### Creating the col to add in colors to the links based on fused, split and lichen
# Assign colors based on NRPS_attached_ICS column
fused_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "Yes"]
split_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "No"]

# Find where those IDs appear in the association matrix
searchresult_fused <- which(association_nrps[,1] %in% fused_ids)
searchresult_split <- which(association_nrps[,1] %in% split_ids)

# Default color for all lines
col <- rep(make.transparent("#000000", 0.2), nrow(association_nrps))

# Color the lines based on the NRPS_attached_ICS values (fused = red, split = green)
col[c(searchresult_fused)] <- make.transparent("#ed1c6e", 0.7)
col[c(searchresult_split)] <- make.transparent("#7ec242", 0.7)

# Now, override the color for any "Lichen" species
lichen_species <- crmMeta$ID[crmMeta$Lichen == "Lichen"]
searchresult_lichen <- which(association_nrps[, 1] %in% lichen_species)
col[searchresult_lichen] <- "#4cc9f0"  # Set Lichen-related links to this color
##
###
# Second the NRPS compared to each other (unrooted only)
###
nrps_gap_gt01_Aligned<- cophylo(nrps_gappyout_unrooted,nrps_gt0_1_unrooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association_nrps)
nrps_gap_gt09_Aligned<- cophylo(nrps_gappyout_unrooted,nrps_gt0_9_unrooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association_nrps)
nrps_gt01_gt09_Aligned<- cophylo(nrps_gt0_1_unrooted,nrps_gt0_9_unrooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association_nrps)
##
# Save the plots as svgs
saveHere <- "/Volumes/T7/MiscProjects/Singh/April_FinalWork/ICS_NRPS_TreeComparisons"
## NRPS gap and gt0.1
svglite(paste0(saveHere,"/nrps-gap_gt01-aligned.svg"),
        width = 7, height = 6)
plot(nrps_gap_gt01_Aligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, show.tip.label=FALSE, ftype="off")
dev.off()
## NRPS gap and gt0.9
svglite(paste0(saveHere,"/nrps-gap_gt09-aligned.svg"),
        width = 7, height = 6)
plot(ics_gap_gt09_Aligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, show.tip.label=FALSE, ftype="off")
dev.off()
## NRPS gt0.1 and gt0.9
svglite(paste0(saveHere,"/nrps-gt01_gt09-aligned.svg"),
        width = 7, height = 6)
plot(ics_gt01_gt09_Aligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, show.tip.label=FALSE, ftype="off")
dev.off()
##
##
############### Update association matrix based on Species matching with the species tree tips
association_ICS_NRPS <- cbind(crmMeta$ID, crmMeta$NRPS_TreeMatch)  # Use Species from crmMeta for association

### Creating the col to add in colors to the links based on fused, split and lichen
# Assign colors based on NRPS_attached_ICS column
fused_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "Yes"]
split_ids <- crmMeta$ID[crmMeta$NRPS_attached_ICS == "No"]

# Find where those IDs appear in the association matrix
searchresult_fused <- which(association_ICS_NRPS[,1] %in% fused_ids)
searchresult_split <- which(association_ICS_NRPS[,1] %in% split_ids)

# Default color for all lines
col <- rep(make.transparent("#000000", 0.2), nrow(association_ICS_NRPS))

# Color the lines based on the NRPS_attached_ICS values (fused = red, split = green)
col[c(searchresult_fused)] <- make.transparent("#ed1c6e", 0.7)
col[c(searchresult_split)] <- make.transparent("#7ec242", 0.7)

# Now, override the color for any "Lichen" species
lichen_species <- crmMeta$ID[crmMeta$Lichen == "Lichen"]
searchresult_lichen <- which(association_ICS_NRPS[, 1] %in% lichen_species)
col[searchresult_lichen] <- "#4cc9f0"  # Set Lichen-related links to this color

#### Running the cophylo commands
ics_nrps_gap_Aligned<- cophylo(ics_gappyout_unrooted,nrps_gappyout_unrooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association_ICS_NRPS)
ics_nrps_gt01_Aligned<- cophylo(ics_gt0_1_unrooted,nrps_gt0_1_unrooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association_ICS_NRPS)
ics_nrps_gt09_Aligned<- cophylo(ics_gt0_9_unrooted,nrps_gt0_9_unrooted, rotate=TRUE, gap =2, space = 30, length.line = 4, assoc=association_ICS_NRPS)
# Save the plots as svgs
saveHere <- "/Volumes/T7/MiscProjects/Singh/April_FinalWork/ICS_NRPS_TreeComparisons"
## ICS and NRPS gappyout
svglite(paste0(saveHere,"/ics-nrps_gappyout-aligned.svg"),
        width = 7, height = 6)
plot(ics_nrps_gap_Aligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, show.tip.label=FALSE, ftype="off")
dev.off()
## ICS and NRPS gt 0.1
svglite(paste0(saveHere,"/ics-nrps_gt01-aligned.svg"),
        width = 7, height = 6)
plot(ics_nrps_gt01_Aligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, show.tip.label=FALSE, ftype="off")
dev.off()
## ICS and NRPS gt0.9
svglite(paste0(saveHere,"/ics-nrps_gt09-aligned.svg"),
        width = 7, height = 6)
plot(ics_nrps_gt09_Aligned,link.type="curved", link.lty="solid",
     hang=-1, print=TRUE, link.lwd=1, lwd=1.5,
     link.col=col, show.tip.label=FALSE, ftype="off")
dev.off()



##### Checking significance of lichen genomes having more ICS
# Create the matrix
# Rows: LFF vs Non-LFF
# Cols: Has ICS vs No ICS
# LFF: 39 Yes, (90-39)=51 No
# Non-LFF: 775 Yes, (3902-775)=3127 No

data <- matrix(c(39, 51, 775, 3127), nrow = 2, byrow = TRUE)
rownames(data) <- c("LFF", "Non-LFF")
colnames(data) <- c("ICS_Present", "ICS_Absent")

# Run Fisher's Exact Test
test <- fisher.test(data)

print(test)

