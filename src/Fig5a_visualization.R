library(pacman)
p_load(tidyverse,phyloseq,reshape2,vegan,ape,grid,gridExtra,ggplot2,magrittr,phylosmith,ggpubr)

meta_df <- read_csv('../dat/clinical_data')
meta_df <- meta_df%>%arrange(clinic_id)
meta_df <- as.data.frame(meta_df)
meta_df <-  meta_df[c("Ifd.Nummer", "cohort")]
colnames(meta_df) <- c("Sample", "Cohort")
meta_df$Age <- round(meta_df$Age,1)
meta_df$Cohort<- factor(meta_df$Cohort, levels = c("before", "pos"))
meta_df$Cohort <- as.factor(meta_df$Cohort)
rownames(meta_df) <-meta_df$Sample

metabolites_df <- read_csv("../dat/metabolites_imputation.csv")
metabolites_df <- as.data.frame(metabolites_df)


metabolites_otu_df <- metabolites_df %>% subset(select=-c(Metabolites)) %>%abs()

metabolites_otu_df <-metabolites_otu_df %>%  t %>% as.data.frame() 
metabolites_meta_df <- meta_df

metabolites_taxa_df <- data.frame(name = colnames(metabolites_otu_df))
rownames(metabolites_taxa_df) <- metabolites_taxa_df$name
metabolites_META = sample_data(metabolites_meta_df)
metabolites_TAX = tax_table(as.matrix(metabolites_taxa_df))
metabolites_OTU = otu_table(t(metabolites_otu_df),taxa_are_rows=TRUE)
metabolites_physeq = phyloseq(metabolites_OTU, metabolites_TAX, metabolites_META)
metabolites_log2_physeq =  transform_sample_counts(metabolites_physeq, function(x) x / log2(x) )
metabolites_otu_meta_df <- merge(metabolites_meta_df["Cohort"], metabolites_otu_df, by=0) %>% subset(select=-c(`Row.names`))

SAM = sample_data(metabolites_physeq) %>% data.frame(.)
ord_bc <- ordinate(metabolites_physeq, method="PCoA", distance="euclidean")
dist_bc <- phyloseq::distance(metabolites_physeq, method = "euclidean")
ord.bc <-  cmdscale(dist_bc, k = 2, eig = T)
set.seed(7)
adonis.bc <- adonis2(dist_bc~  Cohort, data = SAM, permutations = 999)
pcoa.bc.plotting <- as.data.frame(ord.bc$points)
colnames(pcoa.bc.plotting) <- c("axis_1", "axis_2")
pcoa.bc.plotting$cohort <- SAM$Cohort
pcoa.bc.plotting$label <- SAM$Sample
pcoa_ep1 <- round(ord.bc$eig[1]/(sum(ord.bc$eig))*100,1)
pcoa_ep2 <- round(ord.bc$eig[2]/(sum(ord.bc$eig))*100,1)
pcoa.bc.plotting$id <- c(1:40,1:40)
pcoa.bc.plotting$label <- c(1:40,1:40)

plot_ordination(metabolites_physeq, ordinate(metabolites_physeq, method="PCoA", distance="bray"), color = "Cohort") + 
  geom_point(size = 3) 

pcoa.bc.plot <- ggplot(pcoa.bc.plotting, aes(x = axis_1, y = axis_2, colour = cohort,label=label)) +
  geom_point(size = 3) +
  geom_line(aes(group=id),col = 'gray',linetype=2)+
  stat_ellipse(type = "norm",lwd=2)+
  scale_color_manual(labels = c(expression(italic("A.fumigatus -")), expression(italic("A.fumigatus +"))),values=c("before"="#56B4E9","pos"="#E69F00")) +
  theme_pubr() +
  geom_text(aes(label=label),hjust=0, vjust=0,,  size=8)+
  xlab(paste("PCoA 1 (",pcoa_ep1,"%)")) +
  ylab(paste("PCoA 2 (",pcoa_ep2,"%)")) +
  annotate(geom = "text",label = paste("p = ",round(adonis.bc$`Pr(>F)`[1],3)),x=max(pcoa.bc.plotting$axis_1)-100,y=max(pcoa.bc.plotting$axis_2),hjust=0.3,size=14,colour="black")+
  ggtitle("Euclidean-PCoA(Metabolites)")+
  theme(plot.title = element_text(size = 32, hjust = 0.5),  text = element_text(size = 32),legend.key.size = unit(2, "cm"),legend.title = element_blank(),legend.text=element_text(size=28),axis.title.x = element_text(size =36),axis.text.x=element_blank(),axis.line = element_line(colour = 'black', size = 2),axis.ticks = element_line(colour = "black", size = 2),axis.ticks.length=unit(.4, "cm"),
  axis.title.y = element_text(size = 36))

g<-pcoa.bc.plot
ggsave("../res/Fig5a.pdf", width =18, height = 12,g)