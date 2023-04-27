library(viridis)
library(ggpubr)
library(vegan)
library(sjmisc)
library(magrittr)
library(tidyverse)

meta_df <- read_csv('../dat/clinical_data.csv')
meta_df <- meta_df%>%arrange(clinic_id)
meta_df <- as.data.frame(meta_df)
meta_df <-  meta_df[c("Ifd.Nummer", "cohort")]
colnames(meta_df) <- c("Sample", "Cohort")
meta_df$Cohort<- factor(meta_df$Cohort, levels = c("before", "pos"))
meta_df$Cohort <- as.factor(meta_df$Cohort)
meta_df$Sample <- paste0("D",meta_df$Sample)
rownames(meta_df) <-meta_df$Sample
meta_df <- meta_df %>% arrange(Cohort)
metaphlan_df <- read_csv("../dat/otu.csv") %>% column_to_rownames(var="...1")

## Alpha Diversity
before_df <- as.data.frame(metaphlan_df)[meta_df[meta_df$Cohort=="before",]$Sample]
pos_df <-  as.data.frame(metaphlan_df)[meta_df[meta_df$Cohort=="pos",]$Sample]
alpha_diversities <- function(df){
  list = list()
  list[["Shannon"]] = vegan::diversity(df,index = "shannon")
  list[["Chao"]] = apply(df, 1, fossil::chao1)
  return(list)
}


Alpha_Div_before = alpha_diversities(t(before_df))
Alpha_Div_pos= alpha_diversities(t(pos_df))

Shannon = data.frame(Shannon_Index = c(Alpha_Div_before$Shannon,
                                       Alpha_Div_pos$Shannon),
                     Group = factor(c(rep("before",length(Alpha_Div_before$Shannon)),
                                      rep("pos",length(Alpha_Div_pos$Shannon))),
                                    levels = c("before","pos")))

Shannon_plot<- ggplot(Shannon, aes(x = Group, y = Shannon_Index)) + 
  geom_boxplot(aes(fill = Group)) +scale_fill_manual(labels =  c(expression(italic("A.fumigatus -")), expression(italic("A.fumigatus +"))),values=c("#56B4E9","#E69F00"))+
  theme_pubr()  +
  theme(plot.title = element_text(size = 25, hjust = 0.5),  text = element_text(size = 28), legend.title = element_blank() , axis.text.x=element_blank(), axis.ticks = element_blank(),
        legend.key.size = unit(2, "cm"),legend.text=element_text(size=24)) + 
  geom_signif(test="wilcox.test", comparisons = list(c("before", "pos")), map_signif_level=function(p) sprintf("p = %.2g", p), test.args=list(paired = T), tip_length = 0,textsize = 10,y_position = 4.7) +
  xlab(element_blank()) + ylab('Shannon index') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0, binwidth = 0.02) +
  coord_cartesian(ylim = c(0 ,5.5))
  #ggtitle("Comparison of Shannon Index ") 


Chao = data.frame(Chao_Index = c(Alpha_Div_before$Chao,
                                       Alpha_Div_pos$Chao),
                     Group = factor(c(rep("before",length(Alpha_Div_before$Chao)),
                                      rep("pos",length(Alpha_Div_pos$Chao))),
                                    levels = c("before","pos")))

Chao_plot<- ggplot(Chao, aes(x = Group, y = Chao_Index)) + 
  geom_boxplot(aes(fill = Group)) +scale_fill_manual(labels =  c(expression(italic("A.fumigatus -")), expression(italic("A.fumigatus +"))),values=c("#56B4E9","#E69F00"))+
  theme_pubr()  +
  theme(plot.title = element_text(size = 25, hjust = 0.5),  text = element_text(size = 28), legend.title = element_blank() , axis.text.x=element_blank(), axis.ticks = element_blank(),
        legend.key.size = unit(2, "cm"),legend.text=element_text(size=24)) + 
  geom_signif(test="wilcox.test", comparisons = list(c("before", "pos")), map_signif_level=function(p) sprintf("p = %.2g", p), test.args=list(paired = T), tip_length = 0,textsize = 10,y_position = 240) +
  xlab(element_blank()) + ylab('Chao index') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0, binwidth = 0.02)+coord_cartesian(ylim = c(0 ,280))
  #ggtitle("Comparison of Chao Index ")


 
g <- ggarrange(plotlist=list(Shannon_plot,Chao_plot),common.legend = TRUE)
ggsave("../res/S3_a.pdf",g, width = 8)


#### Aitchison-PCoA
set.seed(123)
SAM= meta_df
aitch_full_features <- t(metaphlan_df)
aitch_full_features[aitch_full_features == 0] <- 0.00001
dist.aitchison= coda.base::dist(aitch_full_features, method="aitchison", diag=T, upper=T)
ord.aitchison <-  cmdscale(dist.aitchison, k = 2, eig = T)
adonis.aitchison <- adonis2(dist.aitchison ~  Cohort, data = SAM, permutations = 999)

pcoa.aitchison.plotting <- as.data.frame(ord.aitchison$points)
colnames(pcoa.aitchison.plotting) <- c("axis_1", "axis_2")
pcoa.aitchison.plotting$cohort <- SAM$Cohort
pcoa.aitchison.plotting$label <- SAM$Sample
pcoa_ep1 <- round(ord.aitchison$eig[1]/(sum(ord.aitchison$eig))*100,1)
pcoa_ep2 <- round(ord.aitchison$eig[2]/(sum(ord.aitchison$eig))*100,1)
pcoa.aitchison.plotting$id <- c(1:40,1:40)
pcoa.aitchison.plotting$label <- c(1:40,1:40)
pcoa.aitchison.plot <- ggplot(pcoa.aitchison.plotting, aes(x = axis_1, y = axis_2, colour = cohort,label=label)) +
  geom_point(size = 3) +
  geom_line(aes(group=id),col = 'gray',linetype=2)+
  stat_ellipse(type = "norm")+
  scale_color_manual(labels = c(expression(italic("A.fumigatus -")), expression(italic("A.fumigatus +"))),values=c("before"="#56B4E9","pos"="#E69F00")) +
  theme_pubr() +
  geom_text(aes(label=label),hjust=0, vjust=0,,  size=8)+
  xlab(paste("PCoA 1 (",pcoa_ep1,"%)")) +
  ylab(paste("PCoA 2 (",pcoa_ep2,"%)")) +
  annotate(geom = "text",label = paste("p = ",round(adonis.aitchison$`Pr(>F)`[1],3)),x=0,y=max(pcoa.aitchison.plotting$axis_2)+20,hjust=0.3,size=10,colour="black")+
  ggtitle("Aitchison-PCoA")+
  theme(plot.title = element_text(size = 28, hjust = 0.5),  text = element_text(size = 28),legend.key.size = unit(2, "cm"),legend.title = element_blank(),legend.text=element_text(size=28),axis.title.x = element_text(size =30),
  axis.title.y = element_text(size = 30))+ylim(-60,60)+xlim(-60,60)
ggsave("../res/S3_b.pdf",pcoa.aitchison.plot, width = 9,height = 9)
