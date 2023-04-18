

#################
#### Fig. 2a ####
#################

gem_info = read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/metadata_reformat_supplS3_v2.csv")
gem_info$GEM.information <- factor(gem_info$GEM.information, levels = unique(gem_info$GEM.information))

gem_info <- gem_info[gem_info$GEM.information!="unique genes",]

info_type <- as.factor(c(rep("genes", 252),
                         rep("high-prevalent rxns", 3*252),
                         rep("low-prevalent rxns",504),
                         rep("high-prevalent rxns", 252),
                         rep("low-prevalent rxns",504),
                         rep("metabolites",504)))

gem_info$info_type = info_type

short_info = vector()
for (i in 1:nrow(gem_info)){
  short_info[i] = str_remove(string = gem_info$GEM.information[i], pattern = " genes")
  short_info[i] = str_remove(string = short_info[i], pattern = " reactions")
  short_info[i] = str_remove(string = short_info[i], pattern = " metabolites")
}
gem_info$GEM.information = short_info

p <- ggplot(gem_info, aes(GEM.information, number, color = Sample.origin)) +
  geom_point(position = position_jitterdodge(), alpha = 0.5, size = 0.5) +
  facet_wrap(~ info_type, scales = "free", nrow = 1) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black",angle = 45,hjust = 1,size = 11),
        axis.text.y = element_text(color="black",size = 11),
        strip.text = element_text(size=11),
        panel.grid.major = element_line(colour = "black",linetype="dashed",size=0.1)) +
  scale_color_manual(values = c("magenta","green"),labels = c("Clinical","Environmental")) +
  xlab("")
  #theme(axis.text.x = element_text(angle = 45, hjust=1))

svg(height=3.4, width=8, file="mohammadmirhakkak/A_fumigatus_GEM/res/Fig2a_gem_info.svg")
p
dev.off()


###################
##### Fig. 2b #####
###################
library(RColorBrewer)
library(ggplot2)
pathway <- read.csv('mohammadmirhakkak/A_fumigatus_GEM/res/path_count_v2.csv',row.names = 1)

#substract one from Transport reaction since t_C00140_er_c was deleted and was present in all models. 
pathway[which(pathway$Pathway=='Transport' & pathway$Reactome=='Core'),'num'] = pathway[which(pathway$Pathway=='Transport' & pathway$Reactome=='Core'),'num'] - 1

pathway$Pathway <- factor(pathway$Pathway, levels = pathway$Pathway[1:9])


#vertical version
perc = vector()
for (i in 1:9){
  perc[i] <- as.character(round(pathway[i,1] / (pathway[i,1] + pathway[i+9,1]),digits = 2))
}
perc = c(perc,rep(NA,9))
pathway$perc = perc

svg(height=3.4, width=2.8, file="mohammadmirhakkak/A_fumigatus_GEM/res/Fig2b_path_count.svg")

ggplot(pathway, aes(fill=Reactome, x=num, y=Pathway)) + 
  geom_bar(position="fill", stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black",angle = 90,hjust = 1,size = 8),
        axis.text.y = element_text(color="black",size = 8)) +
  xlab("Number of reactions")+ylab("")+geom_col(color = "black") + 
  geom_text(aes(label = perc),position = position_stack(vjust = 0.5),size = 2.9) +
  scale_fill_brewer(palette = "Spectral") + coord_flip() 
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 15))

dev.off()

###################
##### Fig. 2c #####
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


acc_reactome = read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/acc_5_rxns.csv",row.names = 1,header = TRUE)
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


svg(file = "mohammadmirhakkak/A_fumigatus_GEM/res/Fig2c_acc_rxns.svg",width = 10,height = 10)
facet_plot(p2,data = acc_reactome,geom = geom_barh,panel = "Accessory",
           mapping = aes(fill=Accessory, x=num), stat = "identity")
dev.off()


###################
##### Fig. 3a #####
###################
library(ComplexHeatmap)
library(circlize)
library(viridis)

##import data
df_jac_model <- read.csv('mohammadmirhakkak/A_fumigatus_GEM/res/jac_252.csv',row.names=1,check.names = FALSE,stringsAsFactors = FALSE)
df_labels <- read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/group_jac_252.csv")
#remove first row and first column
#df_jac_model <- df_jac_model[-1,-1]

##row and column annotation
col_annot_model = df_labels$origin

df_jac_model <- as.matrix(df_jac_model)


col_ha_model = HeatmapAnnotation("Sample origin" = col_annot_model,
                                 col = list("Sample origin" = c("Clinical" = "magenta","Environmental"="green")))
                                            



svg("mohammadmirhakkak/A_fumigatus_GEM/res/Fig3a_jaccard.svg",width = 10,height = 7.5)

Heatmap(df_jac_model,
        name = 'Jaccard distance',
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = viridis(5,direction = -1),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        #right_annotation = row_ha_model,
        top_annotation = col_ha_model
)

dev.off()


###################
##### Fig. 3b #####
###################
library(dplyr)
library(tidyr)
library(RANN)
library(ggplot2)
library(ggsignif)
library(vegan)
library(mixOmics)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)


bin_sig_rxns = read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/sig_rxns_252.csv",row.names = 1,check.names = FALSE)
bin_sig_rxns_groups = read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/group_jac_252.csv",row.names = 1)


bin_sig_rxns = as.matrix(bin_sig_rxns)

#column annotation
bin_sig_rxns_groups$names = row.names(bin_sig_rxns_groups)
env = row.names(bin_sig_rxns_groups[bin_sig_rxns_groups$origin=="Environmental",])
cli = row.names(bin_sig_rxns_groups[bin_sig_rxns_groups$origin=="Clinical",])
scatter_annot_final <- cbind(rowSums(bin_sig_rxns[,env])/203 , rowSums(bin_sig_rxns[,cli])/49)
column_ha = HeatmapAnnotation("Presence proportion" = anno_points(ylim = c(0,1),scatter_annot_final, pch = c(16,16), gp = gpar(col = c("green","magenta")),height = unit(1,"cm"),size = unit(0.15,"cm")))

col_fun = colorRamp2(c(0,1),c("gray","black"))

bin_sig_rxns_sorted = cbind(bin_sig_rxns[,env],bin_sig_rxns[,cli])

row_ha_model = HeatmapAnnotation(category = c(rep("Environmental",203),rep("Clinical",49)),
                                 col = list(category = c("Clinical" = "magenta","Environmental"="green")),which = "row",show_legend = FALSE)

GEM_info = read.csv("mohammadmirhakkak/A_fumigatus_GEM/dat/GEM_info.csv",row.names = 2)
row.names(bin_sig_rxns_sorted) = GEM_info[row.names(bin_sig_rxns_sorted),"NAME"]

#heatmap
ht1 = Heatmap(as.matrix(t(bin_sig_rxns_sorted)),
              #column_title_gp = gpar(fontsize = 8),
              column_names_rot = 45,
              row_split = c(rep("",203),rep(" ",49)),
              col = col_fun,
              cluster_row_slices = FALSE,
              show_heatmap_legend = FALSE,
              top_annotation = column_ha,
              name = 'presence/absence plot',
              show_row_names = FALSE,
              right_annotation = row_ha_model)
svg("mohammadmirhakkak/A_fumigatus_GEM/res/Fig3b_sig_rxns.svg",width = 7,height =8)
draw(ht1,padding = unit(c(30, 55, 2, 2), "mm"))
dev.off()



###################
##### Fig. 3c ##### Decision Tree Classification based on metabolite
################### consumption plus presence/absence reactions and vizualization with picharts
library(tibble)
library(stringr)

# Importing the dataset
catabolic_test = read.csv('mohammadmirhakkak/A_fumigatus_GEM/res/catabolic_activity.csv',row.names = 1,check.names = FALSE)
catabolic_test_group = read.csv('mohammadmirhakkak/A_fumigatus_GEM/res/df_groups_catabolic_activity.csv',row.names = 1,check.names = FALSE)
dataset = read.csv('mohammadmirhakkak/A_fumigatus_GEM/res/rxns_252.csv',row.names = 1,check.names = FALSE)
old2new = read.csv("mohammadmirhakkak/A_fumigatus_GEM/dat/old2new.csv")
row_names = row.names(dataset)
row_names_old = vector()
for (i in 1:length(row_names)){
  row_names_old[i] = old2new[old2new$new==row_names[i],'old']
}
row.names(dataset) = row_names_old

#make the dataframe binarized
catabolic_test[catabolic_test>0] = 1


# delete transport and exchanges
removed_rows = vector()
counter = 0
for (i in 1:nrow(dataset)){
  if (substr(x = row.names(dataset)[i],start = 1,stop = 2)=='t_' |
      substr(x = row.names(dataset)[i],start = 1,stop = 6)=='pot_t_' |
      substr(x = row.names(dataset)[i],start = 1,stop = 3)=='EX_'){
    counter = counter + 1
    removed_rows[counter] = i
  }
}
dataset = dataset[-removed_rows,]

#transpose dataset
dataset = as.data.frame(t(dataset))


#unify the row.names
row_names = row.names(dataset)
#for (i in 1:252){
#  row_names[i] = str_replace(string = row_names[i],pattern = ".",replacement = "-")
#}
#row.names(dataset) = row_names



#merge data frames
catabolic_test <- catabolic_test[row.names(dataset),]
dataset = cbind(dataset,catabolic_test)


#environmental -> 0; clinical -> 1
dataset$group = 0
for (i in 1:nrow(dataset)){
  if (catabolic_test_group$niche[i]=='clinical'){
    dataset$group[i] = 1
  }
}

# Encoding the target feature as factor
dataset$group = factor(dataset$group, levels = c(0, 1))


# Make Valid Column Names 
colnames(dataset) <- make.names(colnames(dataset))

library(caTools)
library(rpart)
library(ggplot2)
#control <- rpart.control(minsplit = 20)
classifier_bin_rxn_nutri = rpart(formula = group ~ .,
                                 data = dataset)

#picharts
for (i in 1:nrow(classifier_bin_rxn_nutri$frame)){
  val_env = classifier_bin_rxn_nutri$frame$yval2[i,2]
  val_cln = classifier_bin_rxn_nutri$frame$yval2[i,3]
  df <- data.frame(value = c(val_env, val_cln),
                   Group = c('Environmental','Clinical'))
  
  plot <- ggplot(df, aes(x = "", y = value, fill = Group)) +
    geom_col(color = "black") +
    coord_polar(theta = "y") +
    #geom_text(aes(label = value),
    #          position = position_stack(vjust = 0.5)) +
    theme_void() +
    theme(legend.position = "none") +
    scale_fill_manual(values=c("magenta","green"))
  
  
  dir = paste0("mohammadmirhakkak/A_fumigatus_GEM/res/Fig2e_picharts/pichart",
               as.character(row.names(classifier_bin_rxn_nutri$frame)[i]))
  pdf(paste0(dir,'.pdf'))
  print(plot)
  dev.off()
}