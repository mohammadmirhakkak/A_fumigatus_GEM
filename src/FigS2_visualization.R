###################
##### Fig. S2 #####
###################

library("ggtree")
library("ape")
library("tidyverse")
library("ggstance")
library("RColorBrewer")

fum_tree <- read.tree("mohammadmirhakkak/A_fumigatus_GEM/dat/core_nucleotide_clonalML.newick")
metadata <- read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/metadata_2021-07-05.csv")
isolate_info <- read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/tree_cluster_metadata_20200826.csv",row.names = 1,header = TRUE)

all_tips = fum_tree$tip.label
tips_252 = isolate_info$Sample
drop_tips = setdiff(all_tips,tips_252)
fum_tree <- drop.tip(fum_tree, drop_tips)

metadata$cluster = as.character(metadata$cluster)

#row.names(metadata) = metadata$id
#metadata = metadata[1:300,]
#metadata = metadata[match(fum_tree$tip.label, metadata$id),]

p2 <- ggtree(fum_tree, branch.length='none') %<+% metadata + geom_tippoint(aes(color = source), size=1.5) + 
  scale_color_manual(values = c("magenta","green"),labels = c("Clinical","Environmental")) + theme_tree2(legend.position = 'top')


acc_reactome = read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/acc_5_genes.csv",row.names = 1,header = TRUE)
#strain_genome$Genome <- factor(strain_genome$Genome)
#strain_genome$Genome = fct_rev(strain_genome$Genome)

#remove AB01 for isolate IDs
for (i in 1:nrow(acc_reactome)){
  acc_reactome$isolates[i] <- str_remove(acc_reactome$isolates[i],"AB01-")
  if (grepl(x = acc_reactome$isolates[i],pattern = ".1",fixed = TRUE)){
    acc_reactome$isolates[i] <- substr(x = acc_reactome$isolates[i],start = 1,stop = (nchar(acc_reactome$isolates[i])-2))
  }
}

acc_reactome <- acc_reactome %>% select(isolates, Accessory, num)


#set color as many as rows for the barplot
colors <- brewer.pal(n=5,"Set3")
new_colors = rep(colors,252)


svg(file = "mohammadmirhakkak/A_fumigatus_GEM/res/FigS2_acc_genes.svg",width = 10,height = 10)
facet_plot(p2,data = acc_reactome,geom = geom_barh,panel = "Accessory",
           mapping = aes(fill=Accessory, x=num), stat = "identity")
dev.off()