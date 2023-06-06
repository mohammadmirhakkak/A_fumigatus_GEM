library(doParallel)
library(tidyverse)
library(caret)
library(imbalance)
library(MUVR)
library(VSURF)
library(Boruta)
library(mixOmics)

#Loading data
dataset1 <- read_csv('../dat/min_fva_whole_v6_fastcc_95_afm_mm_media_taurine.csv')
dataset2 <- read_csv('../dat/max_fva_whole_v6_fastcc_95_afm_mm_media_taurine.csv')
dataset1$X1 <- paste("min_",dataset1$X1, sep="")
dataset2$X1 <- paste("max_",dataset2$X1, sep="")
dataset <- rbind(dataset1,dataset2)
dataset <- dataset %>% mutate(across(everything(), .fns = ~replace_na(.,0)))
groups <- read_csv('../dat/df_groups_fit_fva_whole_fastcc_v6_95_afm_mm_media_taurine.csv')

groups <- groups %>% filter(!str_detect(X1,"RR"))
dataset <- dataset[c("X1",groups$X1)]
#Filter variables
dataset <- dataset %>% filter(!grepl('^DM',X1)) %>% filter(!grepl('formation',X1)) %>% filter(!grepl('Growth',X1)) %>% drop_na()

#Rename variables(some packages can not recognize [] )
dataset <- t(dataset)
colnames(dataset) <- sub("\\[",'_',dataset['X1',])
colnames(dataset) <- sub("\\]",'_',colnames(dataset) )
dataset <- as_tibble(dataset[-1,])

#Remove zero variance variables
dataset <- dataset[,apply(dataset, 2, var, na.rm=TRUE) != 0]
#nzv <- nearZeroVar(dataset)
#dataset <- dataset[, -nzv]

#Remove correlated variables
descrCor <- cor(data.matrix(dataset))
highlyCorDescr <- findCorrelation(descrCor, cutoff = .8)
dataset <- dataset[,-highlyCorDescr]

#Remove Linear Dependencies
comboInfo <- findLinearCombos(dataset)
if (!is.null(comboInfo$remove)) {dataset <- dataset[, -comboInfo$remove]}

#Scaling
preProcValues<- caret::preProcess(data.matrix(dataset), method = c("center", "scale"))
dataset <- as.data.frame(predict(preProcValues, data.matrix(dataset)))

#Set group
group_niche <- as.factor(groups$niche)
#group_drug <- as.factor(groups$drug)
group <- as_tibble(group_niche)
colnames(group) <- c("group")

#Feature Selection
merged_features = c()
nCore = 60
dataset_with_label <- cbind(dataset, group)
dataset_minor_with_label <- dataset_with_label[dataset_with_label$group== "Clinical",]
dataset_major_with_label <- dataset_with_label[dataset_with_label$group== "Environmental",]
for (i in seq(1,50)){
  set.seed(i)
  sample_dataset_major_with_label <- sample_n(dataset_major_with_label, size =floor(nrow(dataset_major_with_label)/2))
  dat <- rbind(dataset_minor_with_label,sample_dataset_major_with_label)
  set.seed(i)
  dataset_with_oversampling <- imbalance::oversample(dat,method = "ADASYN",classAttr = "group")
  group_dataset_with_oversampling <- dataset_with_oversampling["group"]
  dataset_with_oversampling["group"] <- NULL

  nRep = nCore
  nOuter = 8
  varRatio = 0.8
  method = "RF"
  cl=makeCluster(nCore)
  registerDoParallel(cl)
  set.seed(2021)
  muvrModel = MUVR(X=dataset_with_oversampling, Y=group_dataset_with_oversampling$group, nRep=nRep, nOuter=nOuter, varRatio=varRatio,method=method)
  MUVR_min_VIP = getVIP(muvrModel, model = 'min')$name
  stopCluster(cl)
  ## plot
  #plotMV(muvrModel, model='min')

  ## Boruta method
  set.seed(2021)
  borutaModel <- Boruta(x=dataset_with_oversampling, y=group_dataset_with_oversampling$group)
  boruta_min_VIP <- getSelectedAttributes(borutaModel)

  ## VSURF method
  set.seed(2021)
  vsurfModel <- VSURF(x=dataset_with_oversampling, y=group_dataset_with_oversampling$group,ncores=nCore,parallel=TRUE)
  vsurf_min_VIP <- colnames(dataset_with_oversampling)[vsurfModel$varselect.pred]

  ## sPLS-DA method

  dataset.plsda <- plsda(data.matrix(dataset_with_oversampling), as.factor(group_dataset_with_oversampling$group), ncomp = 10)
  set.seed(2021)
  perf.plsda.dataset <- perf(dataset.plsda, validation = "Mfold", folds = 5, 
                   progressBar = FALSE, auc = TRUE, nrepeat = 10, cpus = nCore) 
  list.keepX <- c(1:100)
  set.seed(2021)
  tune.splsda.dataset<- tune.splsda(data.matrix(dataset_with_oversampling),as.factor(group_dataset_with_oversampling$group), ncomp = perf.plsda.dataset$choice.ncomp[1], validation = 'Mfold', folds = 5, dist = 'max.dist', measure = "BER",test.keepX = list.keepX, nrepeat = 10, cpus = nCore)
  ncomp <- tune.splsda.dataset$choice.ncomp$ncomp
  select.keepX <- tune.splsda.dataset$choice.keepX[1:ncomp]
  splsda.dataset <- splsda(data.matrix(dataset_with_oversampling), as.factor(group_dataset_with_oversampling$group),ncomp = ncomp, keepX = select.keepX) 
  ## only select features of comp 1
  splsda_min_VIP <- selectVar(splsda.dataset,comp=1)$name

  merged_features <- c(merged_features, MUVR_min_VIP,boruta_min_VIP,vsurf_min_VIP,splsda_min_VIP)

}


#Vote for feature selection
merged_features_table <- table(merged_features) %>% as.data.frame() %>% arrange(desc(Freq))
merged_features_table["merged_features"] <- lapply(merged_features_table["merged_features"], as.character)
merged_features_table_freq <- unique(merged_features_table[["Freq"]])

colnames(merged_features_table) = c("Feature", "Freq")
save.image(file = "feature_selection.RData")
write.csv(merged_features_table, "features_table.csv", row.names = FALSE)
write.csv(dataset, "X.csv",row.names = FALSE)
write.csv(group, "y.csv",row.names = FALSE)
stopCluster(cl)