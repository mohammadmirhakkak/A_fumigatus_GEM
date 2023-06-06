library(ggvenn)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(ggpubr)
library(vegan)
library(phyloseq)
library(phylosmith)
library(sjmisc)
library(magrittr)
library(tidyverse)

meta_df <- read_csv('../dat/clinical_data.csv')
meta_df <- meta_df%>%arrange(cohort,clinic_id)
meta_df <- as.data.frame(meta_df)
meta_df <-  meta_df[c("Ifd.Nummer", "cohort","clinic_id")]
colnames(meta_df) <- c("Sample", "Cohort","clinic_id")
meta_df$Cohort<- factor(meta_df$Cohort, levels = c("before", "pos"))
meta_df$Cohort <- as.factor(meta_df$Cohort)
meta_df$Sample <- paste0("D",meta_df$Sample)
rownames(meta_df) <-meta_df$Sample
metaphlan_df <- read_csv("../dat/otu.csv") %>% column_to_rownames(var="...1")

taxa<- data.frame("Genus" =str_split(rownames(metaphlan_df), "_", simplify = T)[,1] ,"Species" = rownames(metaphlan_df))
rownames(taxa) <- taxa$Species
META <- sample_data(meta_df)
TAX <- tax_table(as.matrix(taxa))
OTU <- otu_table(metaphlan_df,taxa_are_rows=TRUE)

physeq_relative_abundance <- phyloseq(OTU, TAX, META)
physeq_relative_abundance_genus <- tax_glom(physeq_relative_abundance, taxrank="Genus")
taxa_names(physeq_relative_abundance_genus) <- str_split(taxa_names(physeq_relative_abundance_genus) , "_", simplify = T)[,1]

# Top 10 species
otu_relative_abun_df <- otu_table(physeq_relative_abundance) %>% data.frame
median_otu <- apply(otu_relative_abun_df, 1, median, na.rm=T)
top10_taxa <-names(median_otu)[tail( order(median_otu), 10)]
physeq_top10_relative_abundance=prune_taxa(top10_taxa, physeq_relative_abundance)
top10_out <- otu_table(physeq_top10_relative_abundance) %>% data.frame
top10_out <- top10_out[order(apply(top10_out, 1, median, na.rm=T),decreasing=T),]
SAM_temp <- sample_data(physeq_top10_relative_abundance)
SAM_temp$Sample[SAM_temp$Cohort=="before"] <- paste0("D",seq(1:dim(SAM_temp[SAM_temp$Cohort=="before"])[1]),"a")
SAM_temp$Sample[SAM_temp$Cohort=="pos"] <- paste0("D",seq(1:dim(SAM_temp[SAM_temp$Cohort=="pos"])[1]),"b")
sample_data(physeq_top10_relative_abundance)<- SAM_temp
sample_names(physeq_top10_relative_abundance) <- sample_data(physeq_top10_relative_abundance)$Sample

heatmap_mat<- otu_table(physeq_top10_relative_abundance)%>% data.frame%>% as.matrix
heatmap_mat <- heatmap_mat[order(apply(heatmap_mat, 1, median, na.rm=T),decreasing=T),]
colnames(heatmap_mat) <- c(1:40,1:40)
rownames(heatmap_mat) <- gsub("_"," ",rownames(heatmap_mat) )
ha = HeatmapAnnotation(foo = anno_block(gp =  gpar(fill = c("#56B4E9","#E69F00")), labels = c("A.fumigatus -","A.fumigatus +"),labels_gp=gpar(fontface = "italic",fontsize = 20)))
split = rep(1:2, each = 40)
pdf(file="top10_species.pdf",width=16, height=9)
Heatmap(heatmap_mat,colorRamp2(c(0, 50, 100), c(rev(viridis(3)))),column_split = split, top_annotation = ha, 
    column_title = NULL,row_order=NULL, column_order = NULL,height = unit(8, "cm"),heatmap_legend_param = list(title="Abundance(%)"),column_names_gp = grid::gpar(fontsize = 14),row_names_gp = grid::gpar(fontsize = 13,fontface = "italic"))
dev.off()

# distribution of top 10 species

top <- function(x, n){
    tail( order(x), n )
}
top10_c <- c()

for (i in colnames(otu_relative_abun_df)){
  print(rownames(otu_relative_abun_df)[top(otu_relative_abun_df[i],10)])
  top10_c <- c(top10_c,rownames(otu_relative_abun_df)[top(otu_relative_abun_df[i],10)])
}
summary_table <- table(top10_c)
top10_summary_table <- summary_table[top10_taxa] %>% as.data.frame()
in_all_summary_table <- rowSums(otu_relative_abun_df > 0)[top10_taxa] %>% as.data.frame()
summary_all <- cbind(in_all_summary_table,top10_summary_table)
summary_all[10:1,]

# Top 10 genus
otu_relative_abun_df <- otu_table(physeq_relative_abundance_genus) %>% data.frame
median_otu <- apply(otu_relative_abun_df, 1, median, na.rm=T)
top10_taxa <-names(median_otu)[tail( order(median_otu), 10)]
physeq_top10_relative_abundance=prune_taxa(top10_taxa, physeq_relative_abundance_genus)
SAM_temp <- sample_data(physeq_top10_relative_abundance)
SAM_temp$Sample[SAM_temp$Cohort=="before"] <- paste0("D",seq(1:dim(SAM_temp[SAM_temp$Cohort=="before"])[1]),"a")
SAM_temp$Sample[SAM_temp$Cohort=="pos"] <- paste0("D",seq(1:dim(SAM_temp[SAM_temp$Cohort=="pos"])[1]),"b")
sample_data(physeq_top10_relative_abundance)<- SAM_temp
sample_names(physeq_top10_relative_abundance) <- sample_data(physeq_top10_relative_abundance)$Sample

heatmap_mat<- otu_table(physeq_top10_relative_abundance)%>% data.frame%>% as.matrix
heatmap_mat <- heatmap_mat[order(apply(heatmap_mat, 1, median, na.rm=T),decreasing=T),]
colnames(heatmap_mat) <- c(1:40,1:40)
ha = HeatmapAnnotation(foo = anno_block(gp =  gpar(fill = c("#56B4E9","#E69F00")), labels = c("A.fumigatus -","A.fumigatus +"),labels_gp=gpar(fontface = "italic",fontsize = 20)))
split = rep(1:2, each = 40)
pdf(file="top10_genus.pdf",width=16, height=9)
Heatmap(heatmap_mat,colorRamp2(c(0, 50, 100), c(rev(viridis(3)))),column_split = split, top_annotation = ha, 
    column_title = NULL,row_order=NULL, column_order = NULL,height = unit(8, "cm"),heatmap_legend_param = list(title="Abundance(%)"),column_names_gp = grid::gpar(fontsize = 14),row_names_gp = grid::gpar(fontsize = 13,fontface = "italic"))
dev.off()

# distribution of top 10 genus
top <- function(x, n){
    tail( order(x), n )
}
top10_c <- c()

for (i in colnames(otu_relative_abun_df)){
  print(rownames(otu_relative_abun_df)[top(otu_relative_abun_df[i],10)])
  top10_c <- c(top10_c,rownames(otu_relative_abun_df)[top(otu_relative_abun_df[i],10)])
}
summary_table <- table(top10_c)
top10_summary_table <- summary_table[top10_taxa] %>% as.data.frame()
in_all_summary_table <- rowSums(otu_relative_abun_df > 0)[top10_taxa] %>% as.data.frame()
summary_all <- cbind(in_all_summary_table,top10_summary_table)
summary_all[10:1,]
