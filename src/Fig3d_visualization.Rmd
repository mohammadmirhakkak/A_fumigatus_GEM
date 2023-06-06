
library(pacman)
p_load(tidyverse,ROCR)

roc_df <- read_csv("../dat/95_pos_niche_auc_df.csv")
predictions <- roc_df %>% group_by(fold_id)  %>% dplyr::select(pred) %>% group_split(.keep=F)
predictions <- lapply(predictions,FUN=function(x) x[[1]])
labels <- roc_df %>% group_by(fold_id)  %>% dplyr::select(truth) %>% group_split(.keep=F)
labels <- lapply(labels,FUN=function(x) x[[1]])

pred <- prediction(predictions, labels)
perf <- performance(pred, "prec", "rec")
out_df <- data.frame(preauc = unlist(perf@y.values),fold_id = rep(1:5,each = length(perf@y.values[[1]])) )
write_csv(out_df, "Precision_recal.csv")
prauc <- performance(pred, "aucpr")
mean_pracu <-  round(mean(unlist(prauc@y.values)),2)
pdf("prROC.pdf",width=8,height=8)
plot(perf,
     avg= "threshold",
     colorize=F,
     lwd= 3,col="red",
     ylim=c(0.8,1),
     main= "",xlab="",ylab="",cex.lab=1.5, cex.axis=1.5, yaxis.cex.axis=1.5,xaxis.cex.axis=1.5,cex.main=1.5, cex.sub=1.5)
box(lwd=2)
plot(perf,
     lty=3,
     col="grey78",
     add=TRUE)
text(0.7,0.8,paste0("Mean PR-AUC=",mean_pracu),cex=2,col="black" )
dev.off()