
library(tidyverse)
library(ggplot2)
library(ggpubr)


feat_df <- read_tsv("../dat/results_source_contributions_matrix.tsv")
feat_long_df <- feat_df%>% pivot_longer(!sample) 
feat_long_df$group <- NA
feat_long_df$group[str_detect(feat_long_df$name, "Env")] <- "Environment"
feat_long_df$group[str_detect(feat_long_df$name, "Oral")] <- "Oral"
feat_long_df$group[str_detect(feat_long_df$name, "CF")] <- "Cystic fibrosis"
feat_long_df$group[str_detect(feat_long_df$name, "Unknown")] <- "Unknown"
plot_df <- feat_long_df %>% group_by(sample,group) %>% summarise("Source proportion"=sum(value))

p1 <- ggplot(data=plot_df, mapping = aes(x=group,y=`Source proportion`,fill=group))+ 
  geom_boxplot(lwd=1.1)+
  theme_pubr()+  
  scale_fill_brewer(palette="RdBu")+
theme(plot.title = element_text(size = 30, hjust = 0.5),  text = element_text(size = 30), legend.title = element_blank() , axis.text.x=element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.ticks = element_line(colour = "black", size = 1),axis.ticks.length=unit(.4, "cm"),axis.text=element_text(size=30),
        legend.key.size = unit(1, "cm"),legend.text=element_text(size=18))+
  stat_compare_means(comparisons =list(c("Cystic fibrosis","Environment"),c("Cystic fibrosis", "Oral")),label.x =1.5, label.y=1,size=8,paired=T)
ggsave("../res/S2.pdf",p1, width=10, height = 10)
