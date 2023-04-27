
library(tidyverse)
library(caret)
library(circlize)
library(ComplexHeatmap)
library(ggvenn)


## select signif for rxns based on min OR max
#Loading data
filter_df <- read_tsv("../dat/rxnIDsSignKoModules_220318.tsv")
pos_max_dataset <- read.csv('../dat/max_fva_pos_20230330.csv')
rownames(pos_max_dataset) <- paste0("max_",pos_max_dataset$Metabolite)
pos_max_dataset$Metabolite<- NULL
pos_min_dataset <- read.csv('../dat/min_fva_pos_20230330.csv')
rownames(pos_min_dataset) <-paste0("min_",pos_min_dataset$Metabolite)
pos_min_dataset$Metabolite<- NULL
before_max_dataset <- read.csv('../dat/max_fva_before_20230330.csv')
rownames(before_max_dataset) <-paste0("max_",before_max_dataset$Metabolite)
before_max_dataset$Metabolite<- NULL
before_min_dataset <- read.csv('../dat/min_fva_before_20230330.csv')
rownames(before_min_dataset) <- paste0("min_",before_min_dataset$Metabolite)
before_min_dataset$Metabolite<- NULL

max_dataset <- cbind(before_max_dataset,pos_max_dataset) %>% t %>% as.data.frame()
min_dataset <- cbind(before_min_dataset,pos_min_dataset) %>% t %>% as.data.frame()

max_dataset_p <- data.frame()
for (metabolite in colnames(max_dataset)){
  wilcox_df <- max_dataset[metabolite]
  wilcox_df[wilcox_df==0] <- 0.00001
  wilcox_v <- wilcox_df[[metabolite]]
  if(var(na.omit(wilcox_v)) == 0){
    max_dataset_p<- rbind(max_dataset_p, data.frame("metabolite"=metabolite, "p"=1,fold_change=1))
    next
  }
  wilcox_x <- wilcox_v[1:49]
  wilcox_y <-wilcox_v[50:98]
  fold_change <- median(na.omit(wilcox_y))/median(na.omit(wilcox_x))
  wilcoxmodel <- wilcox.test(wilcox_x,wilcox_y,paired=TRUE)
  p_value <- wilcoxmodel$p.value
  max_dataset_p<- rbind(max_dataset_p, data.frame("metabolite"=metabolite, "p"=p_value,fold_change=fold_change))
}

max_dataset_p$adjp <- p.adjust(max_dataset_p$p, "fdr")
max_dataset_p$difffc <- abs(max_dataset_p$fold_change-1)

min_dataset_p <- data.frame()
for (metabolite in colnames(min_dataset)){
  wilcox_df <- min_dataset[metabolite]
  wilcox_df[wilcox_df==0] <- 1
  wilcox_v <- wilcox_df[[metabolite]]
  if(var(na.omit(wilcox_v)) == 0){
     min_dataset_p<- rbind(min_dataset_p, data.frame("metabolite"=metabolite, "p"=1,fold_change=1))
     next
  }
  wilcox_x <- wilcox_v[1:49]
  wilcox_y <-wilcox_v[50:98]
  fold_change <-median(na.omit(wilcox_y))/median(na.omit(wilcox_x))
  wilcoxmodel <- wilcox.test(wilcox_x,wilcox_y,paired=TRUE)
  p_value <- wilcoxmodel$p.value
  min_dataset_p<- rbind(min_dataset_p, data.frame("metabolite"=metabolite, "p"=p_value,fold_change=fold_change))
}

min_dataset_p$adjp <- p.adjust(min_dataset_p$p, "fdr")
min_dataset_p$difffc <- abs(min_dataset_p$fold_change-1)

significant_max_dataset_p <- max_dataset_p[max_dataset_p$adjp<0.05&max_dataset_p$difffc>0.1,]
significant_max <- significant_max_dataset_p$metabolite
significant_min_dataset_p <- min_dataset_p[min_dataset_p$adjp<0.05&min_dataset_p$difffc>0.1,]
significant_min <- significant_min_dataset_p$metabolite
significant_overlap <- union(gsub("max_","", significant_max),gsub("min_","", significant_min))
both_overlap <- intersect(significant_overlap, filter_df$ID_new)

sig_reaction <- filter_df %>% filter(ID_new %in% both_overlap)

min_dataset_p_selected <- min_dataset_p %>% filter(metabolite %in% paste0("min_",both_overlap)) %>% arrange(metabolite)
max_dataset_p_selected <-  max_dataset_p %>% filter(metabolite %in% paste0("max_",both_overlap))%>% arrange(metabolite)

s1 <- max_dataset_p_selected[max_dataset_p_selected$adjp < min_dataset_p_selected$adjp,]
s2 <- min_dataset_p_selected[max_dataset_p_selected$adjp >= min_dataset_p_selected$adjp,]

dataset <- cbind(max_dataset[s1$metabolite], min_dataset[s2$metabolite])

preprocessParams <- preProcess(dataset, method=c("range"),rangeBounds=c(-1,1))
dataset <- predict(preprocessParams, dataset)
heatmap_mat<- dataset %>%  t%>% as.matrix

colnames(heatmap_mat) <- paste0("s",c(1:49,1:49))
rownames(heatmap_mat) <- NULL
col_ha = HeatmapAnnotation(Group =  rep(c("A.fumigatus +","A.fumigatus -"),each=49),col = list(Group = c("A.fumigatus -" = "#56B4E9", "A.fumigatus +" = "#E69F00")))
row_ha = rowAnnotation(Bound =  c(rep("Upper bound",dim(s1)[1]),rep("Lower bound",dim(s2)[1])),col = list(Bound = c("Upper bound" = "#a14a80", "Lower bound" = "#92f44c")))
out_sig_reaction <- sig_reaction[match(str_split(min_dataset_p_selected$metabolite, "min_", simplify = T)[,2],sig_reaction$ID_new),]
pdf(file="../res/5C.pdf",width=14, height=10)
ht <- Heatmap(heatmap_mat,colorRamp2(c(-1, 0, 1), c("blue", "#EEEEEE", "red")), top_annotation = col_ha,right_annotation = row_ha,column_title = NULL,  column_split= c(rep("A.fumigatus +",49),rep("A.fumigatus -",49)),row_gap = unit(2, "mm"),height = unit(18, "cm"),heatmap_legend_param = list(title="Normalized value",direction = "horizontal"),row_title_gp = gpar(fontsize = 20), clustering_distance_rows="spearman")
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()