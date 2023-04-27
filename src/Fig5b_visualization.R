
library(pacman)
p_load(tidyverse,ggpubr)



gr_df <- read_csv("../dat/growth.csv")
gr_plot_df <- gr_df  %>%dplyr::select(subset=c("growth_rate","cohort"))
colnames(gr_plot_df) <- c("growh_rate","cohort")
p <- wilcox.test(gr_plot_df$growh_rate[1:49],gr_plot_df$growh_rate[50:98], paired = T)$p.value
gr_plot<- ggplot(gr_plot_df, aes(x = cohort, y = `growh_rate`)) + 
  geom_boxplot(aes(fill = cohort)) +
  scale_fill_manual(labels = c(expression(italic("A.fumigatus -")), expression(italic("A.fumigatus +"))),values=c("#56B4E9","#E69F00"))+
  theme_pubr()+
  annotate(geom = "text",label =paste("p = ",formatC(p, 3)),x=1.25,y=18,hjust=0.3,size=13,colour="black")+
  theme(plot.title = element_text(size = 40, hjust = 0.5),  text = element_text(size = 40), legend.title = element_blank() , axis.text.x=element_blank(),axis.line = element_line(colour = 'black', size = 2),axis.ticks = element_line(colour = "black", size = 2),axis.ticks.length=unit(.4, "cm"),axis.text=element_text(size=40),
        legend.key.size = unit(2.3, "cm"),legend.text=element_text(size=24),legend.position = "none") +
  xlab(element_blank()) + ylab('Growth rate of clinical GEMs') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0, binwidth = 0.02)+
  coord_cartesian(ylim = c(14 ,18))

g <- ggarrange(gr_plot,ncol=1,nrow = 1)

ggsave("../res/5B.pdf", width=8, height = 10,gr_plot)
