library(ggplot2)


#################
#### Fig. 1b ####
#################

number <- c(1453, 3850, 3051, 2924, 239, 11, 4170, 1957)
gem_info <- c("genes",
              "total rxns",
              "unique rxns",
              "gene-associated rxns",
              "exchange rxns",
              "demand rxns",
              "total mets",
              "unique mets")
information = c(rep("genes",1), rep("reactions",5),rep("metabolites",2))
gem_info = data.frame(number, gem_info, information)

gem_info <- read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/pan_gem_statistics.csv", row.names = 1)

gem_info = gem_info[order(nrow(gem_info):1),]

gem_info <- gem_info[gem_info$gem_info!="unique genes",]
gem_info <- gem_info[gem_info$gem_info!="orphan rxns",]

gem_info$gem_info <- factor(gem_info$gem_info, levels = gem_info$gem_info)

p <- ggplot(gem_info, aes(fill=information, x=number, y=gem_info)) + 
  geom_bar(position="fill", stat = "identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", angle = 45, hjust = 1,size = 10),
        axis.text.y = element_text(color="black",size = 10)) +
  xlab("#")+ylab("")+geom_col(color = "black") + 
  geom_text(aes(label = number),hjust=-.08, size = 4) +
  #scale_fill_brewer(palette = "Spectral") +
  xlim(c(0,6300)) +
  scale_fill_discrete(breaks=c('genes', 'reactions', 'metabolites')) +
  scale_y_discrete(breaks=c("total genes",
                            "total rxns",
                            "unique rxns",
                            "gene-associated rxns",
                            "exchange rxns",
                            "demand rxns",
                            "total mets",
                            "unique mets"),
                   labels=c("total",
                            "total",
                            "unique",
                            "gene-associated",
                            "exchange",
                            "demand",
                            "total",
                            "unique"))

svg("mohammadmirhakkak/A_fumigatus_GEM/res/Fig1b_gem_statistics.svg", width = 3.8, height = 3)
p
dev.off()

#################
#### Fig. 1c ####
#################

#Biomass -> stacked barplot for DNA, RNA, etc. picharts for each stack next to them
biomass = data.frame(bio = rep("Biomass",5),
                     Macromolecule = c("Carbohydrate","Protein","Lipid","Cofactor","DNA and RNA"),
                     value = c(0.428/0.936,0.3/0.936,0.099/0.936,0.066/0.936,0.037/0.936+0.006/0.936))
biomass$Macromolecule <- factor(biomass$Macromolecule, levels = biomass$Macromolecule)

svg("mohammadmirhakkak/A_fumigatus_GEM/res/Fig1c_biomass_component.svg",width = 3,height = 3.5)
ggplot(biomass, aes(fill=Macromolecule, y=value, x=bio)) + 
  geom_bar(position="fill", stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(color="black",size = 12),
        legend.title=element_text(size=12),legend.text=element_text(size=12)) +
  xlab("")+ylab("")+geom_col(color = "black")
dev.off()

#################
#### Fig. 1d ####
#################

library(ggrepel)
pathway = read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/path_count_pan_model.csv",row.names = 1)

# Get the positions
df2 <- pathway %>% 
  mutate(csum = rev(cumsum(rev(num))), 
         pos = num/2 + lead(csum, 1),
         pos = if_else(is.na(pos), num/2, pos))


pathway$pathway = row.names(pathway)
pathway$pathway <- factor(pathway$pathway, levels = pathway$pathway)

svg("mohammadmirhakkak/A_fumigatus_GEM/res/Fig1d_pathway.svg",width = 7,height = 3.8)

ggplot(pathway, aes(x = "", y = num, fill = pathway)) +
  geom_col(color = "black") +
  coord_polar(theta = "y",start = 0,direction = -1) +
  #geom_text(aes(label = num,x = 1.37),
  #          position = position_stack(vjust = 0.5),size = 3) +
  geom_label_repel(data = df2,
                   aes( label = num,y = pos),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  theme_void() +
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12)) +
  scale_fill_manual(values = c("#e6194B",
                               "#f58231",
                               "#ffe119",
                               "#bfef45",
                               "#3cb44b",
                               "#42d4f4",
                               "#4363d8",
                               "#911eb4",
                               "#f032e6",
                               "#fabed4",
                               "#9A6324",
                               "#ffd8b1",
                               "#a9a9a9"))
#theme(legend.position = "none")

dev.off()


#################
#### Fig. 1e ####
#################

compartments = read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/num_compartments.csv",row.names = 1)

compartments$Compartment = row.names(compartments)
compartments$Compartment <- factor(compartments$Compartment, levels = compartments$Compartment)

svg("mohammadmirhakkak/A_fumigatus_GEM/res/Fig1e_compartments.svg",width = 6.5,height = 3.5)

ggplot(compartments, aes(x = Compartment, y = X0, fill = Compartment)) +
  geom_col(color = "black") +
  #coord_polar(theta = "y") +
  geom_text(aes(label = X0),
            position = position_dodge(width = 0.9),vjust = -0.3) +
  theme_void() +  
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12))

dev.off()

#################
#### Fig. 1f ####
#################

df = read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/gr_pred_curated.csv",row.names = 1)

Source <- unique(df$Source)
df_acc <- data.frame(Source,Accuracy = rep(1,length(Source)), SE = rep(0,length(Source)))

filter_df <- df[df$Source=='C',]
accs <- vector()
for (i in 1:nrow(filter_df)){
  accs[i] = (filter_df[i,'TP'] + filter_df[i,'TN']) / sum(filter_df[i,3:6])
}
df_acc[1,'Accuracy'] <- mean(accs)
df_acc[1,'SE'] <- sd(accs) / sqrt(length(accs))

filter_df <- df[df$Source=='N',]
accs <- vector()
for (i in 1:nrow(filter_df)){
  accs[i] = (filter_df[i,'TP'] + filter_df[i,'TN']) / sum(filter_df[i,3:6])
}
df_acc[2,'Accuracy'] <- mean(accs)
df_acc[2,'SE'] <- sd(accs) / sqrt(length(accs))

filter_df <- df[df$Source=='P',]
accs <- vector()
for (i in 1:nrow(filter_df)){
  accs[i] = (filter_df[i,'TP'] + filter_df[i,'TN']) / sum(filter_df[i,3:6])
}
df_acc[3,'Accuracy'] <- mean(accs)
df_acc[3,'SE'] <- sd(accs) / sqrt(length(accs))

filter_df <- df[df$Source=='S',]
accs <- vector()
for (i in 1:nrow(filter_df)){
  accs[i] = (filter_df[i,'TP'] + filter_df[i,'TN']) / sum(filter_df[i,3:6])
}
df_acc[4,'Accuracy'] <- mean(accs)
df_acc[4,'SE'] <- sd(accs) / sqrt(length(accs))

df_acc_curated <- df_acc

df = read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/gr_pred_draft.csv",row.names = 1)

Source <- unique(df$Source)
df_acc <- data.frame(Source,Accuracy = rep(1,length(Source)), SE = rep(0,length(Source)))

filter_df <- df[df$Source=='C',]
accs <- vector()
for (i in 1:nrow(filter_df)){
  accs[i] = (filter_df[i,'TP'] + filter_df[i,'TN']) / sum(filter_df[i,3:6])
}
df_acc[1,'Accuracy'] <- mean(accs)
df_acc[1,'SE'] <- sd(accs) / sqrt(length(accs))

filter_df <- df[df$Source=='N',]
accs <- vector()
for (i in 1:nrow(filter_df)){
  accs[i] = (filter_df[i,'TP'] + filter_df[i,'TN']) / sum(filter_df[i,3:6])
}
df_acc[2,'Accuracy'] <- mean(accs)
df_acc[2,'SE'] <- sd(accs) / sqrt(length(accs))

filter_df <- df[df$Source=='P',]
accs <- vector()
for (i in 1:nrow(filter_df)){
  accs[i] = (filter_df[i,'TP'] + filter_df[i,'TN']) / sum(filter_df[i,3:6])
}
df_acc[3,'Accuracy'] <- mean(accs)
df_acc[3,'SE'] <- sd(accs) / sqrt(length(accs))

filter_df <- df[df$Source=='S',]
accs <- vector()
for (i in 1:nrow(filter_df)){
  accs[i] = (filter_df[i,'TP'] + filter_df[i,'TN']) / sum(filter_df[i,3:6])
}
df_acc[4,'Accuracy'] <- mean(accs)
df_acc[4,'SE'] <- sd(accs) / sqrt(length(accs))

df_acc_draft <- df_acc


#concatenate dataframes
df_gr_acc <- rbind(df_acc_draft,df_acc_curated)
df_gr_acc$GEM = rep(c('draft','curated'),each=4)

df_gr_acc$GEM <- factor(df_gr_acc$GEM, levels = c("draft","curated"))

p<- ggplot(df_gr_acc, aes(x=Source, y=Accuracy, fill=GEM)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=Accuracy, ymax=Accuracy+SE), width=.2,position=position_dodge(.9)) +
  geom_text(aes(label = round(Accuracy,digits = 2)), position = position_dodge(0.9), vjust = 2)
#              position=position_dodge(.9)) 
#print(p)
# Finished bar plot
p <- p+labs("Growth prediction")+
  theme_classic() +
  theme(axis.text = element_text(color="black",size = 13),
        legend.position = "top",plot.title = element_blank(),
        legend.title=element_text(size=13),legend.text=element_text(size=13)) +
  scale_fill_manual(values=c('#999999','#E69F00')) +
  ggtitle("Growth prediction") + ylim(0, 1) + ylab("Growth accuracy")

svg("mohammadmirhakkak/A_fumigatus_GEM/res/Fig1f_growth_prediction_comparison.svg",width = 4, height = 4.5)
p
dev.off()

#################
#### Fig. 1g ####
#################

library(ComplexHeatmap)
library(circlize)

df <- read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/Fig1g_gene_essentiality.csv", row.names = 1)
row.names(df) = c("True", "False")
colnames(df) = c("True", "False")

mat <- as.matrix(df)

col_fun = colorRamp2(c(min(mat),max(mat)), c("white","red"))

svg(height=2, width=2, file="mohammadmirhakkak/A_fumigatus_GEM/res/Fig1g_gene_essentiality.svg")

Heatmap(mat,
        column_names_rot = 45,
        rect_gp = gpar(col = "black", lwd = 2),
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", mat[i, j]), x, y, gp = gpar(fontsize = 20))},
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_heatmap_legend = FALSE,
        row_title = "Essential",
        column_title = "Prediction",
        name = 'Gene essentiality')

dev.off()


#################
#### Fig. 1h ####
#################

df_growth = read.csv("mohammadmirhakkak/A_fumigatus_GEM/res/norm_hypo_growth.csv",row.names = 1)

svg("mohammadmirhakkak/A_fumigatus_GEM/res/Fig1h_norm_hypo.svg",width = 3.6, height = 2.8)

ggplot(df_growth, aes(x=condition, y=growth_rate, fill=type)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge(0.9)) +
  geom_text(aes(label = round(growth_rate,digits = 3)), position = position_dodge(0.9), vjust = 0.5, hjust=1,size = 8) +
  labs("Growth rate") +
  theme_classic() +
  theme(axis.text = element_text(color="black",size = 14),
        axis.title = element_text(size=13),
        legend.position = "top",plot.title = element_blank(),
        legend.title=element_blank(),legend.text=element_text(size=12.5)) +
  scale_fill_manual(values=c('aquamarine2','darkorange')) +
  ggtitle("Growth prediction") + ylab("Growth rate (mmol/grDW/hr)") + xlab(element_blank()) + coord_flip()

dev.off()

