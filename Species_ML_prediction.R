library(caret)
library(randomForest)
library(pROC)

###Load Data & pre-processing###
norm_n_filter <- function(dat){
  dat <- sweep(dat,2,colSums(dat),`/`)
  filter.feat <- names(which(rowSums(dat*100 > 0.1) > .1*ncol(dat)))
  dat <- dat[rownames(dat) %in% filter.feat,]
  dat <- sweep(dat,2,colSums(dat),`/`)
  dat <- as.data.frame(t(dat))
  return(dat)
}
{
  #Meta data
  primary.dat <- read.csv("Data/Primary_data.csv", header = T,row.names = 1)
  colnames(primary.dat) <- gsub("\\.\\..*", "", colnames(primary.dat))
  rownames(primary.dat) <- toupper(rownames(primary.dat))
  primary.dat["R0446","Hypertension.Staging"] <- 1 #Fix Incorrect Label
  primary.dat <- primary.dat[primary.dat$Anti.hypertensive.Agents==0,]
  
  secondary.dat <- read.csv("Data/Secondary_data.csv", header = T,row.names = 1)
  colnames(secondary.dat) <- gsub("\\.\\..*", "", colnames(secondary.dat))
  rownames(secondary.dat) <- toupper(rownames(secondary.dat))
  
  blood.stool.dat <- read.csv("Data/Blood_stool_biomarkers.csv", header = T,row.names = 1)
  colnames(blood.stool.dat) <- gsub("\\.\\..*", "", colnames(blood.stool.dat))
  rownames(blood.stool.dat) <- toupper(rownames(blood.stool.dat))
  colnames(blood.stool.dat)[6] <- "Serum.Vitamin.D.Classification"
  blood.stool.dat[] <- lapply(blood.stool.dat, function(x) as.numeric(as.character(x)))
  #Remove Outliers
  blood.stool.dat[which(blood.stool.dat$Blood.SCFA_Butyric.Acid>20),"Blood.SCFA_Butyric.Acid"] <- NA
  blood.stool.dat[which(blood.stool.dat$Blood.SCFA_Propionic.Acid>750),"Blood.SCFA_Propionic.Acid"] <- NA
  blood.stool.dat[which(blood.stool.dat$Blood.SCFA_Isobutyric.Acid>40),"Blood.SCFA_Isobutyric.Acid"] <- NA
  blood.stool.dat[which(blood.stool.dat$MCP.1>1000),"MCP.1"] <- NA
  blood.stool.dat$Total.Blood.SCFA <- rowSums(blood.stool.dat[,c("Blood.SCFA_Propionic.Acid","Blood.SCFA_Butyric.Acid","Blood.SCFA_Acetic.Acid","Blood.SCFA_Isobutyric.Acid")])
  blood.stool.dat$Total.Stool.SCFA <- rowSums(blood.stool.dat[,c("Stool.Acetic.Acid","Stool.Propionic.Acid","Stool.Butyric.Acid")])
  immune_markers <- colnames(blood.stool.dat)[8:20]
  
  #Dietary Data
  diet.dat <- read.csv("Data/diet_score.csv", header = T, row.names = 1)
  colnames(diet.dat) <- gsub("\\.\\..*", "", colnames(diet.dat))
  rownames(diet.dat) <- toupper(rownames(diet.dat))
  
  #Extra Dietary Data
  food.groups.dat <- read.csv("Data/food_groups.csv", header = T, row.names = 1)
  colnames(food.groups.dat) <- gsub("\\.\\..*", "", colnames(food.groups.dat))
  rownames(food.groups.dat) <- toupper(rownames(food.groups.dat))
  
  #Diet group Data
  diet.groups.dat <- read.csv("Data/diet_groups.csv", header = T, row.names = 1)
  colnames(diet.groups.dat) <- gsub("\\.\\..*", "", colnames(diet.groups.dat))
  rownames(diet.groups.dat) <- toupper(rownames(diet.groups.dat))
  diet.groups.dat$med_original_unhealthy <- ifelse(diet.groups.dat$med_original >= 4, 0, 1)
  diet.groups.dat$DASH_mullen_unhealthy <- ifelse(diet.groups.dat$DASH_mullen >= 4.5, 0, 1)
  
  #Merge all tables
  pat_order <- rownames(primary.dat)
  meta.data <- cbind(primary.dat,blood.stool.dat[pat_order,],secondary.dat[pat_order,],diet.dat[pat_order,],food.groups.dat[pat_order,],diet.groups.dat[pat_order,])
  meta.data <- meta.data[,unlist(lapply(meta.data, is.numeric))]
  
  #Metagenomic data
  read.data <- read.table("Data/metaphlan3_res.tsv", sep = "\t", header = T, row.names = 1)
  read.data <- read.data[grepl("k__Bacteria",rownames(read.data)),]
  read.data <- read.data[,-1]
  colnames(read.data) <- substr(colnames(read.data), 1, 5); colnames(read.data) <- toupper(colnames(read.data))
  common <- intersect(colnames(read.data),rownames(primary.dat))
  read.data <- read.data[,common]
  
  #Species data
  species.data <- read.data[grepl("s__",rownames(read.data)),]
  rownames(species.data) <- gsub(".*s__","",rownames(species.data))
  species.data <- norm_n_filter(species.data)
  species.data <- species.data[rownames(meta.data),]
  #Genus data
  genus.data <- read.data[grepl("g__[^_]+$",rownames(read.data)),]
  rownames(genus.data) <- gsub(".*\\|g__","",rownames(genus.data))
  genus.data <- norm_n_filter(genus.data)
  genus.data <- genus.data[rownames(meta.data),]
  #Phylum data
  phylum.data <- read.data[grepl("p__[^_]+$",rownames(read.data)),]
  rownames(phylum.data) <- gsub(".*\\|p__","",rownames(phylum.data))
  phylum.data <- norm_n_filter(phylum.data)
  phylum.data <- phylum.data[rownames(meta.data),]
  ###Arcsin-sqrt transformation for compositional data###
  transformed.species.data <- asin(sqrt(species.data))
  transformed.genus.data <- asin(sqrt(genus.data))
  transformed.phylum.data <- asin(sqrt(phylum.data))
  
  meta.data$Sex <- ifelse(meta.data$Sex == 0, "Male", 'Female')
  male <- rownames(meta.data[meta.data$Sex=="Male",])
  female <- rownames(meta.data[meta.data$Sex=="Female",])
  meta.data$Hypertension.Staging <- as.factor(meta.data$Hypertension.Staging)
  meta.data$Smoking <- as.factor(meta.data$Smoking)
  meta.data$Menopause_code <- as.factor(meta.data$Menopause_code)
  meta.data$med_original_cat <- as.factor(meta.data$med_original_unhealthy)
  meta.data$DASH_mullen_cat <- as.factor(meta.data$DASH_mullen_unhealthy)
  
  merged.dat <- cbind(meta.data,transformed.species.data[rownames(meta.data),])
}

#################################
###Machine Learning Classifier###
#################################
library(caret)
library(pROC)
library(randomForest)

classification_cv <- function(dat=NULL,outcome,covar = NULL){
  top_feats = 20
  folds <- createMultiFolds(outcome, k = 10, times = 5)
  all_preds <- vector()
  all_obs <- factor()
  all_sig <- vector()
  all_auc <- vector()
  ii=1
  for(i in folds){
    print(paste0("Start fold ", ii))
    train_ind <- i
    if(hasArg(dat)){
      #Feature selection
      merged <- data.frame(y = outcome[train_ind],dat[train_ind,])
      fit.control = trainControl(method = "repeatedcv", number = 10, repeats = 5)
      pls.fit <- train(y ~ .,  data = merged, method = "pls", trControl = fit.control, tuneGrid = data.frame(ncomp = 1:4))
      var_importance <- varImp(pls.fit)$importance
      var_importance <- var_importance[order(var_importance$Overall,decreasing = T),,drop=F]
      sig_feat <- rownames(var_importance)[1:top_feats]
      all_sig <- c(all_sig,sig_feat)
      if(hasArg(covar)){ #Microbiome w/ covar
        input_train <- data.frame(y=as.factor(outcome[train_ind]), dat[train_ind,sig_feat],covar[train_ind,,drop=F])
        input_test <- data.frame(dat[-train_ind,],covar[-train_ind,,drop=F])
      } else { #Only microbiome
        input_train <- data.frame(y=as.factor(outcome[train_ind]), dat[train_ind,sig_feat])
        input_test <- dat[-train_ind,]
      }
    } else { #only covar
      input_train <- data.frame(y=as.factor(outcome[train_ind]), covar[train_ind,,drop=F])
      input_test <- covar[-train_ind,,drop=F]
    }
    mod <- randomForest(y ~ ., data = input_train, mtry = 5, scale = TRUE, importance = T)
    preds <- predict(mod,input_test,type="prob")
    preds <- as.numeric(preds[,2])
    all_preds <- c(all_preds,preds)
    all_obs <- c(all_obs,outcome[-train_ind])
    all_auc <- c(all_auc, auc(outcome[-train_ind],preds))
    ii=ii+1
  }
  return(list(sig_feat = all_sig, all_obs = all_obs, all_preds = all_preds, all_auc = all_auc))
}

###Microbiome + Clinical + SCFA Combinations
set.seed(1)
all_m_res <- classification_cv(dat = transformed.species.data, outcome = meta.data$Hypertension.Staging)
all_c_res <- classification_cv(covar = meta.data[,c("Sex","Age","BMI")],outcome = meta.data$Hypertension.Staging)
all_p_res <- classification_cv(covar = meta.data[,"Blood.SCFA_Propionic.Acid",drop=F], outcome = meta.data$Hypertension.Staging)
all_m_c_res <- classification_cv(dat = transformed.species.data, covar = meta.data[,c("Sex","Age","BMI")], outcome = meta.data$Hypertension.Staging)
all_m_p_res <- classification_cv(dat = transformed.species.data, covar = meta.data[,"Blood.SCFA_Propionic.Acid",drop=F],outcome = meta.data$Hypertension.Staging)
all_p_c_res <- classification_cv(covar = meta.data[,c("Sex","Age","BMI","Blood.SCFA_Propionic.Acid")], outcome = meta.data$Hypertension.Staging)
all_m_p_c_res <- classification_cv(dat = transformed.species.data, covar = meta.data[,c("Sex","Age","BMI","Blood.SCFA_Propionic.Acid")], outcome = meta.data$Hypertension.Staging)

fem_m_res <- classification_cv(dat = transformed.species.data[female,], outcome = meta.data[female,"Hypertension.Staging"])
fem_c_res <- classification_cv(covar = meta.data[female,c("Sex","Age","BMI")],outcome = meta.data[female,"Hypertension.Staging"])
fem_p_res <- classification_cv(covar = meta.data[female,"Blood.SCFA_Propionic.Acid",drop=F], outcome = meta.data[female,"Hypertension.Staging"])
fem_m_c_res <- classification_cv(dat = transformed.species.data[female,], covar = meta.data[female,c("Sex","Age","BMI")], outcome = meta.data[female,"Hypertension.Staging"])
fem_m_p_res <- classification_cv(dat = transformed.species.data[female,], covar = meta.data[female,"Blood.SCFA_Propionic.Acid",drop=F],outcome = meta.data[female,"Hypertension.Staging"])
fem_p_c_res <- classification_cv(covar = meta.data[female,c("Sex","Age","BMI","Blood.SCFA_Propionic.Acid")], outcome = meta.data[female,"Hypertension.Staging"])
fem_m_p_c_res <- classification_cv(dat = transformed.species.data[female,], covar = meta.data[female,c("Sex","Age","BMI","Blood.SCFA_Propionic.Acid")], outcome = meta.data[female,"Hypertension.Staging"])

mal_m_res <- classification_cv(dat = transformed.species.data[male,], outcome = meta.data[male,"Hypertension.Staging"])
mal_c_res <- classification_cv(covar = meta.data[male,c("Sex","Age","BMI")],outcome = meta.data[male,"Hypertension.Staging"])
mal_p_res <- classification_cv(covar = meta.data[male,"Blood.SCFA_Propionic.Acid",drop=F], outcome = meta.data[male,"Hypertension.Staging"])
mal_m_c_res <- classification_cv(dat = transformed.species.data[male,], covar = meta.data[male,c("Sex","Age","BMI")], outcome = meta.data[male,"Hypertension.Staging"])
mal_m_p_res <- classification_cv(dat = transformed.species.data[male,], covar = meta.data[male,"Blood.SCFA_Propionic.Acid",drop=F],outcome = meta.data[male,"Hypertension.Staging"])
mal_p_c_res <- classification_cv(covar = meta.data[male,c("Sex","Age","BMI","Blood.SCFA_Propionic.Acid")], outcome = meta.data[male,"Hypertension.Staging"])
mal_m_p_c_res <- classification_cv(dat = transformed.species.data[male,], covar = meta.data[male,c("Sex","Age","BMI","Blood.SCFA_Propionic.Acid")], outcome = meta.data[male,"Hypertension.Staging"])

plot.roc(all_m_res$all_obs, all_m_res$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.x = 0.2, print.auc.y = 0.6, main = "All cohort",font.lab=2, print.auc.cex = 0.8)
plot.roc(all_c_res$all_obs, all_c_res$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.55, print.auc.cex = 0.8)
plot.roc(all_p_res$all_obs, all_p_res$all_preds, col="purple", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.5, print.auc.cex = 0.8)
plot.roc(all_m_c_res$all_obs, all_m_c_res$all_preds, col="deepskyblue", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.45, print.auc.cex = 0.8)
plot.roc(all_m_p_res$all_obs, all_m_p_res$all_preds, col="firebrick", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.4, print.auc.cex = 0.8)
plot.roc(all_p_c_res$all_obs, all_p_c_res$all_preds, col="limegreen", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.35, print.auc.cex = 0.8)
plot.roc(all_m_p_c_res$all_obs, all_m_p_c_res$all_preds, col="deeppink", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.3, print.auc.cex = 0.8)
legend("topleft", legend=c("M","C","P","M+C","M+P","C+P","M+C+P"), col=c("goldenrod", "salmon","purple","deepskyblue","firebrick","limegreen","deeppink"), lwd=4, cex =0.5, xpd = F, horiz = F,text.font = 2)

plot.roc(fem_m_res$all_obs, fem_m_res$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.x = 0.2, print.auc.y = 0.6, main = "Female cohort",font.lab=2, print.auc.cex = 0.8)
plot.roc(fem_c_res$all_obs, fem_c_res$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.55, print.auc.cex = 0.8)
plot.roc(fem_p_res$all_obs, fem_p_res$all_preds, col="purple", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.5, print.auc.cex = 0.8)
plot.roc(fem_m_c_res$all_obs, fem_m_c_res$all_preds, col="deepskyblue", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.45, print.auc.cex = 0.8)
plot.roc(fem_m_p_res$all_obs, fem_m_p_res$all_preds, col="firebrick", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.4, print.auc.cex = 0.8)
plot.roc(fem_p_c_res$all_obs, fem_p_c_res$all_preds, col="limegreen", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.35, print.auc.cex = 0.8)
plot.roc(fem_m_p_c_res$all_obs, fem_m_p_c_res$all_preds, col="deeppink", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.3, print.auc.cex = 0.8)
legend("topleft", legend=c("M","C","P","M+C","M+P","C+P","M+C+P"), col=c("goldenrod", "salmon","purple","deepskyblue","firebrick","limegreen","deeppink"), lwd=4, cex =0.5, xpd = F, horiz = F,text.font = 2)

plot.roc(mal_m_res$all_obs, mal_m_res$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.x = 0.2, print.auc.y = 0.6, main = "Male cohort",font.lab=2, print.auc.cex = 0.8)
plot.roc(mal_c_res$all_obs, mal_c_res$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.55, print.auc.cex = 0.8)
plot.roc(mal_p_res$all_obs, mal_p_res$all_preds, col="purple", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.5, print.auc.cex = 0.8)
plot.roc(mal_m_c_res$all_obs, mal_m_c_res$all_preds, col="deepskyblue", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.45, print.auc.cex = 0.8)
plot.roc(mal_m_p_res$all_obs, mal_m_p_res$all_preds, col="firebrick", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.4, print.auc.cex = 0.8)
plot.roc(mal_p_c_res$all_obs, mal_p_c_res$all_preds, col="limegreen", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.35, print.auc.cex = 0.8)
plot.roc(mal_m_p_c_res$all_obs, mal_m_p_c_res$all_preds, col="deeppink", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.3, print.auc.cex = 0.8)
legend("topleft", legend=c("M","C","P","M+C","M+P","C+P","M+C+P"), col=c("goldenrod", "salmon","purple","deepskyblue","firebrick","limegreen","deeppink"), lwd=4, cex =0.5, xpd = F, horiz = F,text.font = 2)

t.test(all_c_res$all_auc,all_m_c_res$all_auc)
t.test(female_c_res$all_auc,female_m_c_res$all_auc)
t.test(male_c_res$all_auc,male_m_c_res$all_auc)

parse_freq <- function(vec, main = ""){
  vec = as.data.frame(sort(table(vec), decreasing = T))[1:15,]
  ggplot(vec, aes(reorder(vec,Freq,sum),Freq)) + geom_bar(stat='identity') + coord_flip() + theme(axis.title.y = element_blank()) +ggtitle(main)
}
parse_freq(all_m_res$sig_feat, main = "All Cohort")
parse_freq(female_m_res$sig_feat, main = "Female")
parse_freq(male_m_res$sig_feat, main = "Male")

###With additional Covars
f_met = meta.data[female,]
m_met = meta.data[male,]
m2_vars = c("Sex","Age","BMI")
m3_vars = c("Sex","Age","BMI","Smoking","Sodium.Intake.based.on.spot.urine")
m4_vars = c("Sex","Age","BMI","Smoking","Sodium.Intake.based.on.spot.urine","Fatty.liver.by.CAP.score")

all_m2 <- classification_cv(covar = meta.data[complete.cases(meta.data[,m2_vars]),m2_vars],outcome = meta.data[complete.cases(meta.data[,m2_vars]),'Hypertension.Staging'])
all_m3 <- classification_cv(covar = meta.data[complete.cases(meta.data[,m3_vars]),m3_vars],outcome = meta.data[complete.cases(meta.data[,m3_vars]),'Hypertension.Staging'])
all_m4 <- classification_cv(covar = meta.data[complete.cases(meta.data[,m4_vars]),m4_vars],outcome = meta.data[complete.cases(meta.data[,m4_vars]),'Hypertension.Staging'])
female_m2 <- classification_cv(covar = f_met[complete.cases(f_met[,m2_vars]),m2_vars],outcome = f_met[complete.cases(f_met[,m2_vars]),'Hypertension.Staging'])
female_m3 <- classification_cv(covar = f_met[complete.cases(f_met[,m3_vars]),m3_vars],outcome = f_met[complete.cases(f_met[,m3_vars]),'Hypertension.Staging'])
female_m4 <- classification_cv(covar = f_met[complete.cases(f_met[,m4_vars]),m4_vars],outcome = f_met[complete.cases(f_met[,m4_vars]),'Hypertension.Staging'])
male_m2 <- classification_cv(covar = m_met[complete.cases(m_met[,m2_vars]),m2_vars],outcome = m_met[complete.cases(m_met[,m2_vars]),'Hypertension.Staging'])
male_m3 <- classification_cv(covar = m_met[complete.cases(m_met[,m3_vars]),m3_vars],outcome = m_met[complete.cases(m_met[,m3_vars]),'Hypertension.Staging'])
male_m4 <- classification_cv(covar = m_met[complete.cases(m_met[,m4_vars]),m4_vars],outcome = m_met[complete.cases(m_met[,m4_vars]),'Hypertension.Staging'])

plot.roc(all_m2$all_obs, all_m2$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.x = 0.2, print.auc.y = 0.6, main = "All cohort",font.lab=2, print.auc.cex = 0.8)
plot.roc(all_m3$all_obs, all_m3$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.55, print.auc.cex = 0.8)
plot.roc(all_m4$all_obs, all_m4$all_preds, col="purple", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.5, print.auc.cex = 0.8)
legend("topleft", legend=c("M2","M3","M4"), col=c("goldenrod", "salmon","purple"), lwd=4, cex =0.5, xpd = F, horiz = F,text.font = 2)

plot.roc(female_m2$all_obs, female_m2$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.x = 0.2, print.auc.y = 0.6, main = "Female cohort",font.lab=2, print.auc.cex = 0.8)
plot.roc(female_m3$all_obs, female_m3$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.55, print.auc.cex = 0.8)
plot.roc(female_m4$all_obs, female_m4$all_preds, col="purple", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.5, print.auc.cex = 0.8)
legend("topleft", legend=c("M2","M3","M4"), col=c("goldenrod", "salmon","purple"), lwd=4, cex =0.5, xpd = F, horiz = F,text.font = 2)

plot.roc(male_m2$all_obs, male_m2$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.x = 0.2, print.auc.y = 0.6, main = "Male cohort",font.lab=2, print.auc.cex = 0.8)
plot.roc(male_m3$all_obs, male_m3$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.55, print.auc.cex = 0.8)
plot.roc(male_m4$all_obs, male_m4$all_preds, col="purple", lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.5, print.auc.cex = 0.8)
legend("topleft", legend=c("M2","M3","M4"), col=c("goldenrod", "salmon","purple"), lwd=4, cex =0.5, xpd = F, horiz = F,text.font = 2)

#################################
###Machine Learning Regression###
#################################
regression_cv <- function(dat=NULL,outcome,covar=NULL,cutoff=130){
  top_feats = 20
  folds <- createMultiFolds(outcome, k = 10, times = 5)
  all_preds <- vector()
  all_obs <- factor()
  all_sig <- vector()
  all_auc <- vector()
  ii=1
  for(i in folds){
    print(paste0("Start fold ", ii))
    train_ind <- i
    if(hasArg(dat)){
      #Feature selection
      merged <- data.frame(y = outcome[train_ind],dat[train_ind,])
      fit.control = trainControl(method = "repeatedcv", number = 10, repeats = 5)
      mod.fit <- train(y ~ .,  data = merged, method = "pls", trControl = fit.control, tuneGrid = data.frame(ncomp = 1:4), preProcess = c("scale","center"))
      var_importance <- varImp(mod.fit)$importance
      var_importance <- var_importance[order(var_importance$Overall,decreasing = T),,drop=F]
      sig_feat <- rownames(var_importance)[1:top_feats]
      all_sig <- c(all_sig,sig_feat)
      if(hasArg(covar)){ #Microbiome w/ covar
        input_train <- data.frame(y=outcome[train_ind], dat[train_ind,sig_feat],covar[train_ind,])
        input_test <- data.frame(dat[-train_ind,],covar[-train_ind,])
      } else { #Only microbiome
        input_train <- data.frame(y=outcome[train_ind], dat[train_ind,sig_feat])
        input_test <- dat[-train_ind,]
      }
    } else { #only covar
      input_train <- data.frame(y=outcome[train_ind], covar[train_ind,])
      input_test <- covar[-train_ind,]
    }
    mod <- randomForest(y ~ ., data = input_train, mtry = 5, scale = TRUE)
    preds <- predict(mod,input_test)
    all_preds <- c(all_preds,preds)
    all_obs <- c(all_obs,outcome[-train_ind])
    ii=ii+1
  }
  all_obs <- ifelse(all_obs >= cutoff, 1,0) # 80 for DBP, 130 for SBP
  return(list(sig_feat = all_sig, all_obs = all_obs, all_preds = all_preds))
}

sbp_all_m <- regression_cv(dat = transformed.species.data, outcome = meta.data$X24h.Mean.SBP)
sbp_all_c <- regression_cv(covar = meta.data[,c("Sex","Age","BMI")],outcome = meta.data$X24h.Mean.SBP)
sbp_all_m_c <- regression_cv(dat = transformed.species.data, covar = meta.data[,c("Sex","Age","BMI")], outcome = meta.data$X24h.Mean.SBP)
sbp_female_m <- regression_cv(dat = transformed.species.data[female,], outcome = meta.data[female,"X24h.Mean.SBP"])
sbp_female_c <- regression_cv(covar = meta.data[female,c("Age","BMI")],outcome = meta.data[female,"X24h.Mean.SBP"])
sbp_female_m_c <- regression_cv(dat = transformed.species.data[female,], covar = meta.data[female,c("Age","BMI")], outcome = meta.data[female,"X24h.Mean.SBP"])
sbp_male_m <- regression_cv(dat = transformed.species.data[male,], outcome = meta.data[male,"X24h.Mean.SBP"])
sbp_male_c <- regression_cv(covar = meta.data[male,c("Age","BMI")],outcome = meta.data[male,"X24h.Mean.SBP"])
sbp_male_m_c <- regression_cv(dat = transformed.species.data[male,], covar = meta.data[male,c("Age","BMI")], outcome = meta.data[male,"X24h.Mean.SBP"])

plot.roc(sbp_all_m$all_obs, sbp_all_m$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.y = 0.4, main = "All cohort",font.lab=2)
plot.roc(sbp_all_c$all_obs, sbp_all_c$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.3)
plot.roc(sbp_all_m_c$all_obs, sbp_all_m_c$all_preds, col="limegreen", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.2)
legend("bottom", legend=c("Microbiome","Covariates", "M+C"), col=c("goldenrod", "salmon","limegreen"), lwd=5, cex =0.7, xpd = TRUE, horiz = TRUE,text.font = 2)
plot.roc(sbp_female_m$all_obs, sbp_female_m$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.y = 0.4, main = "Female",font.lab=2)
plot.roc(sbp_female_c$all_obs, sbp_female_c$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.3)
plot.roc(sbp_female_m_c$all_obs, sbp_female_m_c$all_preds, col="limegreen", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.2)
legend("bottom", legend=c("Microbiome","Covariates", "M+C"), col=c("goldenrod", "salmon","limegreen"), lwd=5, cex =0.7, xpd = TRUE, horiz = TRUE,text.font = 2)
plot.roc(sbp_male_m$all_obs, sbp_male_m$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.y = 0.4, main = "Male", font.lab=2)
plot.roc(sbp_male_c$all_obs, sbp_male_c$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.3)
plot.roc(sbp_male_m_c$all_obs, sbp_male_m_c$all_preds, col="limegreen", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.2)
legend("bottom", legend=c("Microbiome","Covariates", "M+C"), col=c("goldenrod", "salmon","limegreen"), lwd=5, cex =0.7, xpd = TRUE, horiz = TRUE,text.font = 2)

dbp_all_m <- regression_cv(dat = transformed.species.data, outcome = meta.data$X24h.Mean.DBP, cutoff = 80)
dbp_all_c <- regression_cv(covar = meta.data[,c("Sex","Age","BMI")],outcome = meta.data$X24h.Mean.DBP, cutoff = 80)
dbp_all_m_c <- regression_cv(dat = transformed.species.data, covar = meta.data[,c("Sex","Age","BMI")], outcome = meta.data$X24h.Mean.DBP, cutoff = 80)
dbp_female_m <- regression_cv(dat = transformed.species.data[female,], outcome = meta.data[female,"X24h.Mean.DBP"], cutoff = 80)
dbp_female_c <- regression_cv(covar = meta.data[female,c("Age","BMI")],outcome = meta.data[female,"X24h.Mean.DBP"], cutoff = 80)
dbp_female_m_c <- regression_cv(dat = transformed.species.data[female,], covar = meta.data[female,c("Age","BMI")], outcome = meta.data[female,"X24h.Mean.DBP"], cutoff = 80)
dbp_male_m <- regression_cv(dat = transformed.species.data[male,], outcome = meta.data[male,"X24h.Mean.DBP"], cutoff = 80)
dbp_male_c <- regression_cv(covar = meta.data[male,c("Age","BMI")],outcome = meta.data[male,"X24h.Mean.DBP"], cutoff = 80)
dbp_male_m_c <- regression_cv(dat = transformed.species.data[male,], covar = meta.data[male,c("Age","BMI")], outcome = meta.data[male,"X24h.Mean.DBP"], cutoff = 80)

plot.roc(dbp_all_m$all_obs, dbp_all_m$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.y = 0.4, main = "All cohort",font.lab=2)
plot.roc(dbp_all_c$all_obs, dbp_all_c$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.3)
plot.roc(dbp_all_m_c$all_obs, dbp_all_m_c$all_preds, col="limegreen", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.2)
legend("bottom", legend=c("Microbiome","Covariates", "M+C"), col=c("goldenrod", "salmon","limegreen"), lwd=5, cex =0.7, xpd = TRUE, horiz = TRUE,text.font = 2)
plot.roc(dbp_female_m$all_obs, dbp_female_m$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.y = 0.4, main = "Female",font.lab=2)
plot.roc(dbp_female_c$all_obs, dbp_female_c$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.3)
plot.roc(dbp_female_m_c$all_obs, dbp_female_m_c$all_preds, col="limegreen", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.2)
legend("bottom", legend=c("Microbiome","Covariates", "M+C"), col=c("goldenrod", "salmon","limegreen"), lwd=5, cex =0.7, xpd = TRUE, horiz = TRUE,text.font = 2)
plot.roc(dbp_male_m$all_obs, dbp_male_m$all_preds, col="goldenrod", lwd=2, print.auc=TRUE, print.auc.y = 0.4, main = "Male", font.lab=2)
plot.roc(dbp_male_c$all_obs, dbp_male_c$all_preds, col="salmon", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.3)
plot.roc(dbp_male_m_c$all_obs, dbp_male_m_c$all_preds, col="limegreen", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y = 0.2)
legend("bottom", legend=c("Microbiome","Covariates", "M+C"), col=c("goldenrod", "salmon","limegreen"), lwd=5, cex =0.7, xpd = TRUE, horiz = TRUE,text.font = 2)

