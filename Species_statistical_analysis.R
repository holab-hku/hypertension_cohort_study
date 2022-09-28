library(ggplot2)
library(ggpubr)
library(vegan)
library(rstatix)
library(pheatmap)

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

###My functions###
stacked_bar <- function(dat,feat,max_feat=7,title=""){
  if(ncol(dat)>=max_feat){
    other_col <- colnames(dat)[which(!colnames(dat) %in% names(sort(colSums(dat),decreasing = T)[1:max_feat-1]))] 
    dat$Other <- rowSums(dat[,other_col])
    dat <- dat[,-which(colnames(dat) %in% other_col)]
  }
  dat.m <- reshape2::melt(as.matrix(dat))
  colnames(dat.m) <- c("Patient","Feat","Abund")
  dat.m$Lab <- as.factor(meta.data[match(dat.m$Patient,rownames(meta.data)),feat])
  dat.m$Sex <- as.factor(meta.data[match(dat.m$Patient,rownames(meta.data)),"Sex"])
  all.facet <- dat.m
  all.facet$Sex <- "All"
  dat.m <- rbind(dat.m,all.facet)
  dat.m <- dat.m[!is.na(dat.m$Lab),]
  print(ggplot(dat.m, aes(fill=Feat, y=Abund, x=Lab)) + 
          geom_bar(position="fill", stat="identity") + xlab(feat) + ylab("Relative Abundance") + scale_fill_brewer(palette="Dark2") + 
          ggtitle(title) + facet_wrap(~Sex))
}
my_boxplot <- function(var,label, ylab="", title="", hide_ns = T){
  dat <- data.frame(y = var, lab = meta.data[,label], Sex = meta.data$Sex)
  all.facet <- dat
  all.facet$Sex <- "All"
  dat <- rbind(dat,all.facet)
  dat <- dat[!is.na(dat$lab),]
  stat.test <- dat %>% group_by(Sex) %>% pairwise_wilcox_test(y ~ lab, p.adjust.method = "BH") %>% add_y_position()
  print(ggboxplot(dat,"lab", "y", color = "lab", add = "jitter") + stat_pvalue_manual(stat.test, label = "p.adj", hide.ns = hide_ns, label.size = 6) + 
          xlab(label) + ylab(ylab) + ggtitle(title) + facet_wrap(~Sex) + 
          theme(axis.title = element_text(size = 14, face="bold"), axis.text = element_text(size = 12, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.position = "none"))
}
lm_models <- function(dat,x,y){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(feature in x){
    #df <- dat[!is.na(dat[,feature])&!is.na(meta.data$Smoking)&!is.na(meta.data$Sodium.Intake.based.on.spot.urine)&!is.na(meta.data$Fatty.liver.by.CAP.score),]
    df <- dat[!is.na(dat[,feature]),]
    males <- rownames(df[df$Sex=="Male",])
    females <- rownames(df[df$Sex=="Female",])
    
    m1.all <- lm(df[,y] ~ df[,feature])
    m1.female <- lm(df[females,y] ~ df[females,feature])
    m1.male <- lm(df[males,y] ~ df[males,feature])
    
    m2.all <- lm(df[,y] ~ df[,feature] + Sex + Age + BMI, df)
    m2.female <- lm(df[females,y] ~ df[females,feature] + df[females,"Age"] + df[females,"BMI"])
    m2.male <- lm(df[males,y] ~ df[males,feature] + df[males,"Age"] + df[males,"BMI"])
    
    m3.all <- lm(df[,y] ~ df[,feature] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine, df)
    m3.female <- lm(df[females,y] ~ df[females,feature] + df[females,"Age"] + df[females,"BMI"] + df[females,"Smoking"] + df[females,"Sodium.Intake.based.on.spot.urine"] + df[females,"Menopause_code"])
    m3.male <- lm(df[males,y] ~ df[males,feature] + df[males,"Age"] + df[males,"BMI"] + df[males,"Smoking"] + df[males,"Sodium.Intake.based.on.spot.urine"])
    
    m4.all <- lm(df[,y] ~ df[,feature] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine + Fatty.liver.by.CAP.score, df)
    m4.female <- lm(df[females,y] ~ df[females,feature] + df[females,"Age"] + df[females,"BMI"] + df[females,"Smoking"] + df[females,"Sodium.Intake.based.on.spot.urine"] + df[females,"Menopause_code"] + df[females,"Fatty.liver.by.CAP.score"])
    m4.male <- lm(df[males,y] ~ df[males,feature] + df[males,"Age"] + df[males,"BMI"] + df[males,"Smoking"] + df[males,"Sodium.Intake.based.on.spot.urine"] + df[males,"Fatty.liver.by.CAP.score"])
      
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1.all)$coefficients[2,4],sign(summary(m1.all)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m1",feature,summary(m1.female)$coefficients[2,4],sign(summary(m1.female)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m1",feature,summary(m1.male)$coefficients[2,4],sign(summary(m1.male)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2.all)$coefficients[2,4],sign(summary(m2.all)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m2",feature,summary(m2.female)$coefficients[2,4],sign(summary(m2.female)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m2",feature,summary(m2.male)$coefficients[2,4],sign(summary(m2.male)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3.all)$coefficients[2,4],sign(summary(m3.all)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m3",feature,summary(m3.female)$coefficients[2,4],sign(summary(m3.female)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m3",feature,summary(m3.male)$coefficients[2,4],sign(summary(m3.male)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4.all)$coefficients[2,4],sign(summary(m4.all)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m4",feature,summary(m4.female)$coefficients[2,4],sign(summary(m4.female)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m4",feature,summary(m4.male)$coefficients[2,4],sign(summary(m4.male)$coefficients[2,1]))
  }
  return(cor.matrix)
}
glm_models <- function(dat,x,y){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(feature in x){
    #df <- dat[!is.na(dat[,feature])&!is.na(meta.data$Smoking)&!is.na(meta.data$Sodium.Intake.based.on.spot.urine)&!is.na(meta.data$Fatty.liver.by.CAP.score),]
    df <- dat[!is.na(dat[,feature]),]
    males <- rownames(df[df$Sex=="Male",])
    females <- rownames(df[df$Sex=="Female",])
    
    m1.all <- glm(df[,y] ~ df[,feature], family = 'binomial')
    m1.female <- glm(df[females,y] ~ df[females,feature], family = 'binomial')
    m1.male <- glm(df[males,y] ~ df[males,feature], family = 'binomial')
    
    m2.all <- glm(df[,y] ~ df[,feature] + Sex + Age + BMI, df, family = 'binomial')
    m2.female <- glm(df[females,y] ~ df[females,feature] + df[females,"Age"] + df[females,"BMI"], family = 'binomial')
    m2.male <- glm(df[males,y] ~ df[males,feature] + df[males,"Age"] + df[males,"BMI"], family = 'binomial')
    
    m3.all <- glm(df[,y] ~ df[,feature] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine, df, family = 'binomial')
    m3.female <- glm(df[females,y] ~ df[females,feature] + df[females,"Age"] + df[females,"BMI"] + df[females,"Smoking"] + df[females,"Sodium.Intake.based.on.spot.urine"] + df[females,"Menopause_code"], family = 'binomial')
    m3.male <- glm(df[males,y] ~ df[males,feature] + df[males,"Age"] + df[males,"BMI"] + df[males,"Smoking"] + df[males,"Sodium.Intake.based.on.spot.urine"], family = 'binomial')
    
    m4.all <- glm(df[,y] ~ df[,feature] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine + Fatty.liver.by.CAP.score, df, family = 'binomial')
    m4.female <- glm(df[females,y] ~ df[females,feature] + df[females,"Age"] + df[females,"BMI"] + df[females,"Smoking"] + df[females,"Sodium.Intake.based.on.spot.urine"] + df[females,"Menopause_code"] + df[females,"Fatty.liver.by.CAP.score"], family = 'binomial')
    m4.male <- glm(df[males,y] ~ df[males,feature] + df[males,"Age"] + df[males,"BMI"] + df[males,"Smoking"] + df[males,"Sodium.Intake.based.on.spot.urine"] + df[males,"Fatty.liver.by.CAP.score"], family = 'binomial')
    
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1.all)$coefficients[2,4],sign(summary(m1.all)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m1",feature,summary(m1.female)$coefficients[2,4],sign(summary(m1.female)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m1",feature,summary(m1.male)$coefficients[2,4],sign(summary(m1.male)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2.all)$coefficients[2,4],sign(summary(m2.all)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m2",feature,summary(m2.female)$coefficients[2,4],sign(summary(m2.female)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m2",feature,summary(m2.male)$coefficients[2,4],sign(summary(m2.male)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3.all)$coefficients[2,4],sign(summary(m3.all)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m3",feature,summary(m3.female)$coefficients[2,4],sign(summary(m3.female)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m3",feature,summary(m3.male)$coefficients[2,4],sign(summary(m3.male)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4.all)$coefficients[2,4],sign(summary(m4.all)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m4",feature,summary(m4.female)$coefficients[2,4],sign(summary(m4.female)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m4",feature,summary(m4.male)$coefficients[2,4],sign(summary(m4.male)$coefficients[2,1]))
  }
  return(cor.matrix)
}
permanova_models <- function(dat,x){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
  for(feature in x){
    met <- dat[!is.na(dat[,feature])&!is.na(dat$Smoking)&!is.na(dat$Sodium.Intake.based.on.spot.urine)&!is.na(dat$Fatty.liver.by.CAP.score),]
    spe.dat <- transformed.species.data[rownames(met),]
    males <- rownames(met[met$Sex=="Male",])
    females <- rownames(met[met$Sex=="Female",])
    
    m1.all <- adonis(spe.dat ~ met[,feature], met, permutations=1000, method = "bray")
    m1.female <- adonis(spe.dat[females,] ~ met[females,feature], permutations=1000, method = "bray")
    m1.male <- adonis(spe.dat[males,] ~ met[males,feature], permutations=1000, method = "bray")
    
    m2.all <- adonis(spe.dat ~ met[,feature] + Sex + Age + BMI, met, permutations=1000, method = "bray")
    m2.female <- adonis(spe.dat[females,] ~ met[females,feature] + met[females,"Age"] + met[females,"BMI"], permutations=1000, method = "bray")
    m2.male <- adonis(spe.dat[males,] ~ met[males,feature] + met[males,"Age"] + met[males,"BMI"], permutations=1000, method = "bray")
    
    m3.all <- adonis(spe.dat ~ met[,feature] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine, met, permutations=1000, method = "bray")
    m3.female <- adonis(spe.dat[females,] ~ met[females,feature] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"], permutations=1000, method = "bray")
    m3.male <- adonis(spe.dat[males,] ~ met[males,feature] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"], permutations=1000, method = "bray")
    
    m4.all <- adonis(spe.dat ~ met[,feature] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine + Fatty.liver.by.CAP.score, met, permutations=1000, method = "bray")
    m4.female <- adonis(spe.dat[females,] ~ met[females,feature] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"] + met[females,"Fatty.liver.by.CAP.score"], permutations=1000, method = "bray")
    m4.male <- adonis(spe.dat[males,] ~ met[males,feature] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"] + met[males,"Fatty.liver.by.CAP.score"], permutations=1000, method = "bray")
    
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,m1.all$aov.tab$`Pr(>F)`[1])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m1",feature,m1.female$aov.tab$`Pr(>F)`[1])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m1",feature,m1.male$aov.tab$`Pr(>F)`[1])
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2.all$aov.tab$`Pr(>F)`[1])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m2",feature,m2.female$aov.tab$`Pr(>F)`[1])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m2",feature,m2.male$aov.tab$`Pr(>F)`[1])
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3.all$aov.tab$`Pr(>F)`[1])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m3",feature,m3.female$aov.tab$`Pr(>F)`[1])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m3",feature,m3.male$aov.tab$`Pr(>F)`[1])
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4.all$aov.tab$`Pr(>F)`[1])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m4",feature,m4.female$aov.tab$`Pr(>F)`[1])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m4",feature,m4.male$aov.tab$`Pr(>F)`[1])
  }
  return(cor.matrix)
}

### ENV Variables
vascular_markers_continuous <- c("X24h.Mean.SBP","X24h.Mean.DBP")
vascular_markers_categorical <- c("Hypertension.Staging")
acids_variables <- grep("Blood.SCFA",colnames(blood.stool.dat),value = T)
diet_variables_continuous <- c(names(which(colSums(is.na(diet.dat))!=nrow(diet.dat))),"med_fish","med_original","DASH_mullen")
diet_variables_categorical <- c("med_original_unhealthy","DASH_mullen_unhealthy")
food_group_variables <- names(which(colSums(is.na(food.groups.dat))!=nrow(food.groups.dat)))
#################
### F/B Ratio ###
#################
meta.data$f.b.ratio <- phylum.data$Firmicutes/phylum.data$Bacteroidetes
###Vascular
cor.mat <- lm_models(meta.data,vascular_markers_continuous,"f.b.ratio")
cor.mat.hyp <- glm_models(meta.data,"f.b.ratio",vascular_markers_categorical)
cor.mat.hyp$feature <- "Hypertension.Staging"
cor.mat <- rbind(cor.mat,cor.mat.hyp)
cor.mat$pval <- as.numeric(cor.mat$pval)
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("F/B Ratio - Vascular variables")
##Blood
cor.mat <- lm_models(meta.data,c(acids_variables,immune_markers),"f.b.ratio")
cor.mat$pval <- as.numeric(cor.mat$pval)
cor.mat$sym <- ifelse(cor.mat$pval>0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("F/B Ratio - Blood variables")
###Diet
cor.mat <- lm_models(meta.data,diet_variables_continuous,"f.b.ratio")
cor.mat$pval <- as.numeric(cor.mat$pval)
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("F/B Ratio - Diet variables")
### All food groups - Diet
cor.mat <- lm_models(meta.data,food_group_variables,"f.b.ratio")
cor.mat$pval <- as.numeric(cor.mat$pval)
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("F/B Ratio - Food group variables")

##############################
###Alpha Diversity analysis###
##############################
meta.data$diversity <- vegan::diversity(transformed.species.data, index = 'shannon')
dat <- data.frame(y = meta.data$diversity, lab = meta.data[,"Hypertension.Staging"], Sex = meta.data$Sex)
all.facet <- dat
all.facet$Sex <- "All"
dat <- rbind(dat,all.facet)
dat <- dat[!is.na(dat$lab),]
dat$lab <- ifelse(dat$lab==1, "Hypertensive", "Normotensive")
stat.test <- dat %>% group_by(Sex) %>% pairwise_wilcox_test(y ~ lab, p.adjust.method = "BH") %>% add_y_position()
ggboxplot(dat,"lab", "y", color = "lab", add = "jitter") + stat_pvalue_manual(stat.test, label = "p.adj", hide.ns = TRUE, label.size = 6, y.position = 4.2) + 
        xlab("") + ylab("Shannon Diversity Index") + facet_wrap(~Sex) + ylim(2.5, 4.4) + 
        theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=16, face="bold"),
              legend.position = "none")

###Vascular
cor.mat <- lm_models(meta.data,vascular_markers_continuous,"diversity")
cor.mat$feature <- gsub("\\."," ",cor.mat$feature)
cor.mat$feature <- gsub("X","",cor.mat$feature)
cor.mat.hyp <- glm_models(meta.data,"diversity",vascular_markers_categorical)
cor.mat.hyp$feature <- "Hypertensive"
cor.mat <- rbind(cor.mat,cor.mat.hyp)
cor.mat$pval <- as.numeric(cor.mat$pval)
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10",midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("Shannon diversity - Vascular variables") + 
        theme(axis.title = element_text(size = 10, face="bold"), axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold"), title = element_text(size = 12, face="bold"))
###Blood
cor.mat <- lm_models(meta.data,c(acids_variables,immune_markers),"diversity")
cor.mat$pval <- as.numeric(cor.mat$pval)
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("Shannon diversity - Blood variables")
###Diet
cor.mat <- lm_models(meta.data,diet_variables_continuous,"diversity")
for(x in diet_variables_categorical){
  cor.mat.cat <- glm_models(meta.data,"diversity",x)
  cor.mat.cat$feature <- x
  cor.mat <- rbind(cor.mat,cor.mat.cat)
}
cor.mat$pval <- as.numeric(cor.mat$pval)
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("Shannon diversity - Diet variables")
### All food groups - Diet
cor.mat <- lm_models(meta.data,food_group_variables,"diversity")
cor.mat$pval <- as.numeric(cor.mat$pval)
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("Shannon diversity - Food group variables")

################################
###Permanova - Beta diversity###
################################=
#Vascular parameters
vasc.perma <- permanova_models(meta.data,c(vascular_markers_categorical,vascular_markers_continuous))
vasc.perma$pval <- as.numeric(vasc.perma$pval)
vasc.perma$sym <- ifelse(vasc.perma$pval>=0.05, "", ifelse(vasc.perma$pval<=0.01,ifelse(vasc.perma$pval<=0.001,"***", "**"), "*"))
vasc.perma$log <- -log10(vasc.perma$pval)
vasc.perma$feature <- gsub("Hypertension.Staging","Hypertensive",vasc.perma$feature)
vasc.perma$feature <- gsub("X","",vasc.perma$feature)
vasc.perma$feature <- gsub("\\."," ",vasc.perma$feature)
ggplot(vasc.perma, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold"))
#Acid variables
acid.immune.perma <- permanova_models(meta.data,c(acids_variables,immune_markers))
acid.immune.perma$pval <- as.numeric(acid.immune.perma$pval)
acid.immune.perma$sym <- ifelse(acid.immune.perma$pval>=0.05, "", ifelse(acid.immune.perma$pval<=0.01,ifelse(acid.immune.perma$pval<=0.001,"***", "**"), "*"))
acid.immune.perma$log <- -log10(acid.immune.perma$pval)
ggplot(acid.immune.perma, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("Acid and Immune variables")
#Dietary variables
diet.perma <- permanova_models(meta.data,c(diet_variables_categorical,diet_variables_continuous))
diet.perma$pval <- as.numeric(diet.perma$pval)
diet.perma$sym <- ifelse(diet.perma$pval>=0.05, "", ifelse(diet.perma$pval<=0.01,ifelse(diet.perma$pval<=0.001,"***", "**"), "*"))
diet.perma$log <- -log10(diet.perma$pval)
ggplot(diet.perma, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
  theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
        legend.text =  element_text(size = 12, face="bold"),legend.title =  element_text(size = 12, face="bold"))
#Dietary variables subset for figure
diet.perma <- permanova_models(meta.data,c('HEI_total','fiber_g','diet_sodium_mg','Energy_kcal','alcohol_g'))
diet.perma$pval <- as.numeric(diet.perma$pval)
diet.perma$sym <- ifelse(diet.perma$pval>=0.05, "", ifelse(diet.perma$pval<=0.01,ifelse(diet.perma$pval<=0.001,"***", "**"), "*"))
diet.perma$log <- -log10(diet.perma$pval)
ggplot(diet.perma, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
  coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
  theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
        legend.text =  element_text(size = 12, face="bold"),legend.title =  element_text(size = 12, face="bold"))
#Food groups variables
food.groups.perma <- permanova_models(meta.data,food_group_variables)
food.groups.perma$pval <- as.numeric(food.groups.perma$pval)
food.groups.perma$sym <- ifelse(food.groups.perma$pval>=0.05, "", ifelse(food.groups.perma$pval<=0.01,ifelse(food.groups.perma$pval<=0.001,"***", "**"), "*"))
food.groups.perma$log <- -log10(food.groups.perma$pval)
ggplot(food.groups.perma, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("Food group variables")

# MDS plots of Beta-diversity
#Hypertension
met <- meta.data[!is.na(meta.data[,"Hypertension.Staging"])&!is.na(meta.data$Sodium.Intake.based.on.spot.urine),]
dat <- transformed.species.data[rownames(met),]
males <- rownames(met[met$Sex=="Male",])
females <- rownames(met[met$Sex=="Female",])
p.nova.all <- adonis(dat ~ met[,"Hypertension.Staging"] + met[,"Sex"] + met[,"Age"] + met[,"BMI"] + met[,"Smoking"] + met[,"Sodium.Intake.based.on.spot.urine"], permutations=5000, method = "bray")
p.nova.male <- adonis(dat[males,] ~ met[males,"Hypertension.Staging"] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"], permutations=5000, method = "bray")
p.nova.female <- adonis(dat[females,] ~ met[females,"Hypertension.Staging"] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"], permutations=5000, method = "bray")
p.nova.all <- p.nova.all$aov.tab[1,"Pr(>F)"]; p.nova.male <- p.nova.male$aov.tab[1,"Pr(>F)"]; p.nova.female <- p.nova.female$aov.tab[1,"Pr(>F)"]
p.labels <- data.frame(label = paste0("p = ",c(signif(p.nova.all,2),signif(p.nova.male,2),signif(p.nova.female,2))), gender= c("All","Male","Female"))

bray.dist.matrix <- vegdist(dat[rownames(met),], method="bray")
mds <- cmdscale(bray.dist.matrix, eig = TRUE, x.ret = TRUE)
mds.var.per <- round(mds$eig/sum(mds$eig)*100, 1)
mds.values <- mds$points
mds.data <- data.frame(Patient = rownames(mds.values), X = mds.values[,1], Y = mds.values[,2], var = meta.data[match(rownames(mds.values),rownames(meta.data)),"Hypertension.Staging"], gender=meta.data[match(rownames(mds.values),rownames(meta.data)),"Sex"])
all.facet <- mds.data
all.facet$gender <- "All"
mds.data <- rbind(mds.data,all.facet)
mds.data$var <- ifelse(mds.data$var==0,"Normotensive","Hypertensive")
ggplot(mds.data, aes(x=X, y=Y, col=var)) + geom_point(size = 4, alpha = 0.8) + theme_classic() + facet_wrap(~gender) +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) +
  geom_text(data = p.labels, mapping = aes(x = -Inf, y = -Inf, label = label), hjust = -0.1, vjust = -1, col = "red", size = 5) + 
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=16, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), legend.position = "bottom")

#Masked hypertension
met <- meta.data[!is.na(meta.data[,"Masked.hypertension"])&!is.na(meta.data$Sodium.Intake.based.on.spot.urine),]
dat <- transformed.species.data[rownames(met),]
males <- rownames(met[met$Sex=="Male",])
females <- rownames(met[met$Sex=="Female",])
p.nova.all <- adonis(dat ~ met[,"Masked.hypertension"] + met[,"Sex"] + met[,"Age"] + met[,"BMI"] + met[,"Smoking"] + met[,"Sodium.Intake.based.on.spot.urine"], permutations=5000, method = "bray")
p.nova.male <- adonis(dat[males,] ~ met[males,"Masked.hypertension"] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"], permutations=5000, method = "bray")
p.nova.female <- adonis(dat[females,] ~ met[females,"Masked.hypertension"] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"], permutations=5000, method = "bray")
p.nova.all <- p.nova.all$aov.tab[1,"Pr(>F)"]; p.nova.male <- p.nova.male$aov.tab[1,"Pr(>F)"]; p.nova.female <- p.nova.female$aov.tab[1,"Pr(>F)"]
p.labels <- data.frame(label = paste0("p = ",c(signif(p.nova.all,2),signif(p.nova.male,2),signif(p.nova.female,2))), gender= c("All","Male","Female"))
bray.dist.matrix <- vegdist(dat[rownames(met),], method="bray")
mds <- cmdscale(bray.dist.matrix, eig = TRUE, x.ret = TRUE)
mds.var.per <- round(mds$eig/sum(mds$eig)*100, 1)
mds.values <- mds$points
mds.data <- data.frame(Patient = rownames(mds.values), X = mds.values[,1], Y = mds.values[,2], var = meta.data[match(rownames(mds.values),rownames(meta.data)),"Masked.hypertension"], gender=meta.data[match(rownames(mds.values),rownames(meta.data)),"Sex"])
all.facet <- mds.data
all.facet$gender <- "All"
mds.data <- rbind(mds.data,all.facet)
mds.data$var <- ifelse(mds.data$var==0,"Normotensive","Masked Hypertension")
ggplot(mds.data, aes(x=X, y=Y, col=var)) + geom_point(size = 4, alpha = 0.8) + theme_classic() + facet_wrap(~gender) +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) +
  geom_text(data = p.labels, mapping = aes(x = -Inf, y = -Inf, label = label), hjust = -0.1, vjust = -1, col = "red", size = 5) + 
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=16, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), legend.position = "bottom")

#########################
###Enterotype analysis###
#########################
dat.species <- transformed.species.data[,order(colSums(transformed.species.data),decreasing = T)[1:20]]
dat.genus <- transformed.genus.data[,order(colSums(transformed.genus.data),decreasing = T)[1:15]]
annot = data.frame(Gender=meta.data$Sex, Hypertensive=meta.data$Hypertension.Staging)
rownames(annot) = rownames(meta.data)
pheatmap(t(dat.species), annotation_col=annot, 
         annotation_names_row=FALSE,
         annotation_names_col=FALSE,
         show_colnames = F,
         fontsize_col=5)
pheatmap(t(dat.genus), annotation_col=annot, 
         annotation_names_row=FALSE,
         annotation_names_col=FALSE,
         show_colnames = F,
         fontsize=12)

bray.dist.matrix <- vegdist(transformed.genus.data, method="bray")
mds <- cmdscale(bray.dist.matrix, eig = TRUE, x.ret = TRUE)
mds.var.per <- round(mds$eig/sum(mds$eig)*100, 1)
mds.values <- mds$points
mds.data <- data.frame(Patient = rownames(mds.values), X = mds.values[,1], Y = mds.values[,2], 
                       Bacteroides = dat.genus$Bacteroides, 
                       Bifidobacterium = dat.genus$Bifidobacterium,
                       Prevotella = dat.genus$Prevotella,
                       Ruminococcus = dat.genus$Ruminococcus)
ggplot(mds.data, aes(x=X, y=Y, col= Bacteroides)) + geom_point(size = 5) + theme_classic() +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) +
  scale_color_gradient(low="blue", high="red") + ggtitle("Bacteroides Abundance") +
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), plot.title = element_text(size = 16, face = "bold"))
ggplot(mds.data, aes(x=X, y=Y, col= Bifidobacterium)) + geom_point(size = 5) + theme_classic() +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) +
  scale_color_gradient(low="blue", high="red") + ggtitle("Bifidobacterium Abundance") +
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), plot.title = element_text(size = 16, face = "bold"))
ggplot(mds.data, aes(x=X, y=Y, col= Prevotella)) + geom_point(size = 5) + theme_classic() +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) + 
  scale_color_gradient(low="blue", high="red") + ggtitle("Prevotella Abundance") +
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), plot.title = element_text(size = 16, face = "bold"))
ggplot(mds.data, aes(x=X, y=Y, col= Ruminococcus)) + geom_point(size = 5) + theme_classic() +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) + 
  scale_color_gradient(low="blue", high="red") + ggtitle("Ruminococcus Abundance") +
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), plot.title = element_text(size = 16, face = "bold"))
kmean_res <- kmeans(transformed.genus.data,3)$cluster
meta.data$enterotype <- kmean_res[rownames(meta.data)]
meta.data$enterotype <- as.factor(meta.data$enterotype)
ggplot(mds.data, aes(x=X, y=Y, col= as.factor(kmean_res))) + geom_point(size = 5) + theme_classic() +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) +
  ggtitle("K Means Clustering") +
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), plot.title = element_text(size = 16, face = "bold"))

#Correlation enterotype with continuous vascular
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in c(vascular_markers_continuous)){
  met <- meta.data[!is.na(meta.data[,feature]),]
  males <- rownames(met[met$Sex=="Male",])
  females <- rownames(met[met$Sex=="Female",])
  
  m1.all <- aov(met[,feature] ~ enterotype, met)
  m1.female <- aov(met[females,feature] ~ met[females,"enterotype"])
  m1.male <- aov(met[males,feature] ~ met[males,"enterotype"])
  
  m2.all <- aov(met[,feature] ~ met[,"enterotype"] + Sex + Age + BMI, met)
  m2.female <- aov(met[females,feature] ~ met[females,"enterotype"] + met[females,"Age"] + met[females,"BMI"])
  m2.male <- aov(met[males,feature] ~ met[males,"enterotype"] + met[males,"Age"] + met[males,"BMI"])
  
  m3.all <- aov(met[,feature] ~ met[,"enterotype"] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine, met)
  m3.female <- aov(met[females,feature] ~ met[females,"enterotype"] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"])
  m3.male <- aov(met[males,feature] ~ met[males,"enterotype"] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"])
  
  m4.all <- aov(met[,feature] ~ met[,"enterotype"] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine + Fatty.liver.by.CAP.score, met)
  m4.female <- aov(met[females,feature] ~ met[females,"enterotype"] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"] + met[females,"Fatty.liver.by.CAP.score"])
  m4.male <- aov(met[males,feature] ~ met[males,"enterotype"] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"] + met[males,"Fatty.liver.by.CAP.score"])
  
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m1",feature,summary(m1.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m1",feature,summary(m1.male)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m2",feature,summary(m2.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m2",feature,summary(m2.male)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m3",feature,summary(m3.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m3",feature,summary(m3.male)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m4",feature,summary(m4.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m4",feature,summary(m4.male)[[1]][1,5])
}
cor.matrix$pval <- as.numeric(cor.matrix$pval)
cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
cor.matrix$log <- -log10(cor.matrix$pval)
ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
  coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "white")+ facet_wrap(~gender) + ggtitle("Enterotype groups")

#Correlation enterotype with acids
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in c(acids_variables,immune_markers)){
  met <- meta.data[!is.na(meta.data[,feature]),]
  males <- rownames(met[met$Sex=="Male",])
  females <- rownames(met[met$Sex=="Female",])
  
  m1.all <- aov(met[,feature] ~ enterotype, met)
  m1.female <- aov(met[females,feature] ~ met[females,"enterotype"])
  m1.male <- aov(met[males,feature] ~ met[males,"enterotype"])
  
  m2.all <- aov(met[,feature] ~ met[,"enterotype"] + Sex + Age + BMI, met)
  m2.female <- aov(met[females,feature] ~ met[females,"enterotype"] + met[females,"Age"] + met[females,"BMI"])
  m2.male <- aov(met[males,feature] ~ met[males,"enterotype"] + met[males,"Age"] + met[males,"BMI"])
  
  m3.all <- aov(met[,feature] ~ met[,"enterotype"] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine, met)
  m3.female <- aov(met[females,feature] ~ met[females,"enterotype"] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"])
  m3.male <- aov(met[males,feature] ~ met[males,"enterotype"] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"])
  
  m4.all <- aov(met[,feature] ~ met[,"enterotype"] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine + Fatty.liver.by.CAP.score, met)
  m4.female <- aov(met[females,feature] ~ met[females,"enterotype"] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"] + met[females,"Fatty.liver.by.CAP.score"])
  m4.male <- aov(met[males,feature] ~ met[males,"enterotype"] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"] + met[males,"Fatty.liver.by.CAP.score"])
  
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m1",feature,summary(m1.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m1",feature,summary(m1.male)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m2",feature,summary(m2.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m2",feature,summary(m2.male)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m3",feature,summary(m3.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m3",feature,summary(m3.male)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m4",feature,summary(m4.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m4",feature,summary(m4.male)[[1]][1,5])
}
cor.matrix$pval <- as.numeric(cor.matrix$pval)
cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
cor.matrix$log <- -log10(cor.matrix$pval)
ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
  coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("Enterotype groups")

#Correlation enterotype with diet
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in c(diet_variables_continuous,food_group_variables)){
  met <- meta.data[!is.na(meta.data[,feature]),]
  males <- rownames(met[met$Sex=="Male",])
  females <- rownames(met[met$Sex=="Female",])
  
  m1.all <- aov(met[,feature] ~ enterotype, met)
  m1.female <- aov(met[females,feature] ~ met[females,"enterotype"])
  m1.male <- aov(met[males,feature] ~ met[males,"enterotype"])
  
  m2.all <- aov(met[,feature] ~ met[,"enterotype"] + Sex + Age + BMI, met)
  m2.female <- aov(met[females,feature] ~ met[females,"enterotype"] + met[females,"Age"] + met[females,"BMI"])
  m2.male <- aov(met[males,feature] ~ met[males,"enterotype"] + met[males,"Age"] + met[males,"BMI"])
  
  m3.all <- aov(met[,feature] ~ met[,"enterotype"] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine, met)
  m3.female <- aov(met[females,feature] ~ met[females,"enterotype"] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"])
  m3.male <- aov(met[males,feature] ~ met[males,"enterotype"] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"])
  
  m4.all <- aov(met[,feature] ~ met[,"enterotype"] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine + Fatty.liver.by.CAP.score, met)
  m4.female <- aov(met[females,feature] ~ met[females,"enterotype"] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"] + met[females,"Fatty.liver.by.CAP.score"])
  m4.male <- aov(met[males,feature] ~ met[males,"enterotype"] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"] + met[males,"Fatty.liver.by.CAP.score"])
  
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m1",feature,summary(m1.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m1",feature,summary(m1.male)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m2",feature,summary(m2.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m2",feature,summary(m2.male)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m3",feature,summary(m3.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m3",feature,summary(m3.male)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4.all)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Female","m4",feature,summary(m4.female)[[1]][1,5])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Male","m4",feature,summary(m4.male)[[1]][1,5])
}
cor.matrix$pval <- as.numeric(cor.matrix$pval)
cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
cor.matrix$log <- -log10(cor.matrix$pval)
ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
  coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("Enterotype groups")

#Barplot with SBP
my_boxplot(meta.data$X24h.Mean.SBP,"enterotype", ylab = "24h Mean SBP", hide_ns = F)
my_boxplot(meta.data$X24h.Mean.DBP,"enterotype", ylab = "24h Mean DBP", hide_ns = F)

###########################
###Correlation analysis####
###########################
#24hr Hypertension Staging vs species
{
cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(bac in colnames(transformed.species.data)){
  cor.mat <- rbind(cor.mat,glm_models(merged.dat,bac,"Hypertension.Staging"))
}
cor.mat$pval <- as.numeric(cor.mat$pval)
filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red") +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold"))
}
#24hr SBP vs species
{
cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(bac in colnames(transformed.species.data)){
  cor.mat <- rbind(cor.mat,lm_models(merged.dat,bac,"X24h.Mean.SBP"))
}
cor.mat$pval <- as.numeric(cor.mat$pval)
filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold"))
}
#24hr DBP vs species
{
cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(bac in colnames(transformed.species.data)){
  cor.mat <- rbind(cor.mat,lm_models(merged.dat,bac,"X24h.Mean.DBP"))
}
cor.mat$pval <- as.numeric(cor.mat$pval)
filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold"))
}

#Blood vs Hypertension
{
cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feat in c(acids_variables,immune_markers)){
  cor.mat <- rbind(cor.mat,glm_models(merged.dat,feat,"Hypertension.Staging"))
}
cor.mat$pval <- as.numeric(cor.mat$pval)
filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)
}
#Blood vs 24hr SBP
{
cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feat in c(acids_variables,immune_markers)){
  cor.mat <- rbind(cor.mat,lm_models(merged.dat,feat,"X24h.Mean.SBP"))
}
cor.mat$pval <- as.numeric(cor.mat$pval)
filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)
}
#Blood vs 24hr DBP
{
cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feat in c(acids_variables,immune_markers)){
  cor.mat <- rbind(cor.mat,lm_models(merged.dat,feat,"X24h.Mean.DBP"))
}
cor.mat$pval <- as.numeric(cor.mat$pval)
filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)
}

#Diet vs Hypertension
{
cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feat in c(diet_variables_continuous,food_group_variables)){
  cor.mat <- rbind(cor.mat,glm_models(merged.dat,feat,"Hypertension.Staging"))
}
cor.mat$pval <- as.numeric(cor.mat$pval)
filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)
}
#Diet vs 24hr SBP
{
cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feat in c(diet_variables_continuous,food_group_variables)){
  cor.mat <- rbind(cor.mat,lm_models(merged.dat,feat,"X24h.Mean.SBP"))
}
cor.mat$pval <- as.numeric(cor.mat$pval)
filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)
}
#Diet vs 24hr DBP
{
cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feat in c(diet_variables_continuous,food_group_variables)){
  cor.mat <- rbind(cor.mat,lm_models(merged.dat,feat,"X24h.Mean.DBP"))
}
cor.mat$pval <- as.numeric(cor.mat$pval)
filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
print(ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender))
}
#Diet vs Species
{
sig_diet_feat <- c('HEI_total','fiber_g','med_fish','total_fruit_g','total_vegetable_g')
for(diet in sig_diet_feat){
  cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(bac in colnames(transformed.species.data)){
    cor.mat <- rbind(cor.mat,lm_models(merged.dat,bac,diet))
  }
  cor.mat$pval <- as.numeric(cor.mat$pval)
  filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
  cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
  cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
  cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
  print(ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
    coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle(diet) +
    theme(axis.text = element_text(size = 12, face="bold"), strip.text = element_text(size=12, face="bold"),
          legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold"), plot.title = element_text(size = 14, face = "bold")))
}
}
#Diet vs Blood
{
for(diet in sig_diet_feat){
  cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(blood in c(acids_variables,immune_markers)){
    cor.mat <- rbind(cor.mat,lm_models(merged.dat,blood,diet))
  }
  cor.mat$pval <- as.numeric(cor.mat$pval)
  filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
  if(identical(filt.feat, character(0))){next}
  cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
  cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
  cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
  print(ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle(diet))
}
}
#Blood SCFAs vs Immune
{
for(blood in acids_variables){
  cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(immune in immune_markers){
    cor.mat <- rbind(cor.mat,lm_models(merged.dat,immune,blood))
  }
  cor.mat$pval <- as.numeric(cor.mat$pval)
  filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
  if(identical(filt.feat, character(0))){next}
  cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
  cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
  cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
  print(ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle(blood))
}
}
#Blood vs Species
{
for(blood in c(acids_variables,immune_markers)){
  cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(bac in colnames(transformed.species.data)){
    cor.mat <- rbind(cor.mat,lm_models(merged.dat,bac,blood))
  }
  cor.mat$pval <- as.numeric(cor.mat$pval)
  filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
  if(identical(filt.feat, character(0))){next}
  cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
  cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
  cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
  print(ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red",name = "-log10") +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle(blood) + 
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")))
}
#For Total blood SCFA Fig
cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(bac in colnames(transformed.species.data)){
  cor.mat <- rbind(cor.mat,lm_models(merged.dat,bac,"Total.Blood.SCFA"))
}
cor.mat$pval <- as.numeric(cor.mat$pval)
filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red",name = "-log10") +
      coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("Total Plasma SCFA") + 
      theme(axis.text = element_text(size = 12, face="bold"), strip.text = element_text(size=12, face="bold"),
            legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold"), plot.title = element_text(size = 14, face = "bold"))
}
#Stool SCFAs vs Diet
{
for(stool in stool.var){
  cor.mat <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(diet in c(diet_variables_continuous,food_group_variables)){
    cor.mat <- rbind(cor.mat,lm_models(merged.dat,diet,stool))
  }
  cor.mat$pval <- as.numeric(cor.mat$pval)
  filt.feat <- cor.mat[cor.mat$pval<=0.01,"feature"]
  if(identical(filt.feat, character(0))){next}
  cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
  cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
  cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
  print(ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle(stool))
}
}

#######################################
### Non-hypertensive Male vs Female ###
#######################################
nonhyp.f <- rownames(meta.data[meta.data$Hypertension.stage==0&meta.data$Sex=="Female",])
nonhyp.m <- rownames(meta.data[meta.data$Hypertension.stage==0&meta.data$Sex=="Male",])

#Stacked barplot
dat <- species.data[c(nonhyp.f,nonhyp.m),]
other_col <- colnames(dat)[which(!colnames(dat) %in% names(sort(colSums(dat),decreasing = T)[1:7]))] 
dat$Other <- rowSums(dat[,other_col])
dat <- dat[,-which(colnames(dat) %in% other_col)]
dat.m <- reshape2::melt(as.matrix(dat))
colnames(dat.m) <- c("Patient","Feat","Abund")
dat.m$Lab <- as.factor(meta.data[match(dat.m$Patient,rownames(meta.data)),"Sex"])
ggplot(dat.m, aes(fill=Feat, y=Abund, x=Lab)) + geom_bar(position="fill", stat="identity") + 
        xlab("") + ylab("Relative Abundance") + scale_fill_brewer(palette="Dark2") + ggtitle("Normotensive individuals") +
        theme(axis.title = element_text(size = 14, face="bold"), axis.text = element_text(size = 12, face="bold"), strip.text = element_text(size=16, face="bold"),
              legend.text = element_text(size = 14, face="bold"), legend.title = element_blank(), plot.title = element_text(size = 14, face="bold"))

#Species sig dif
cor.mat <-  data.frame(model = character(), feature = character(), pval = character(), direction = character())
merged.dat <- cbind(dat,meta.data[rownames(dat),])
merged.dat$Sex <- as.factor(merged.dat$Sex)
for(bac in colnames(dat)){
  mod1 <- glm(Sex ~ merged.dat[,bac], merged.dat, family = 'binomial')
  mod2 <- glm(Sex ~ merged.dat[,bac] + Age + BMI, merged.dat, family = 'binomial')
  mod3 <- glm(Sex ~ merged.dat[,bac] + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine, merged.dat, family = 'binomial')
  mod4 <- glm(Sex ~ merged.dat[,bac] + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine + Fatty.liver.by.CAP.score, merged.dat, family = 'binomial')
  cor.mat[nrow(cor.mat)+1,] <- c("m1",bac,summary(mod1)$coefficients[2,4],sign(summary(mod1)$coefficients[2,1]))
  cor.mat[nrow(cor.mat)+1,] <- c("m2",bac,summary(mod2)$coefficients[2,4],sign(summary(mod2)$coefficients[2,1]))
  cor.mat[nrow(cor.mat)+1,] <- c("m3",bac,summary(mod3)$coefficients[2,4],sign(summary(mod3)$coefficients[2,1]))
  cor.mat[nrow(cor.mat)+1,] <- c("m4",bac,summary(mod4)$coefficients[2,4],sign(summary(mod4)$coefficients[2,1]))
}
cor.mat$pval <- as.numeric(cor.mat$pval)
filt.feat <- cor.mat[cor.mat$pval<=0.05,"feature"]
cor.mat <- cor.mat[cor.mat[,"feature"]%in%filt.feat,]
cor.mat$sym <- ifelse(cor.mat$pval>=0.05, "", ifelse(cor.mat$pval<=0.01,ifelse(cor.mat$pval<=0.001,"***", "**"), "*"))
cor.mat$log <- -log10(cor.mat$pval)*as.numeric(cor.mat$direction)
ggplot(cor.mat, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")

#Beta diversity
p.nova <- adonis2(dat ~ merged.dat[,'Sex'] + merged.dat[,"Age"] + merged.dat[,"BMI"], permutations=999, method = "bray")
p.nova <- p.nova[1,"Pr(>F)"]
bray.dist.matrix <- vegdist(dat, method="bray")
mds <- cmdscale(bray.dist.matrix, eig = TRUE, x.ret = TRUE)
mds.var.per <- round(mds$eig/sum(mds$eig)*100, 1)
mds.values <- mds$points
mds.data <- data.frame(Patient = rownames(mds.values), X = mds.values[,1], Y = mds.values[,2], gender=meta.data[match(rownames(mds.values),rownames(meta.data)),"Sex"])
ggplot(mds.data, aes(x=X, y=Y, col=as.factor(gender))) + geom_point(size = 5, alpha = 0.8) + theme_classic() +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) +
  geom_text(mapping = aes(x = -Inf, y = -Inf, label = paste0("p = ",p.nova)), hjust = -0.1, vjust = -1, col = "red", size = 5) +
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=16, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), legend.position = "bottom")
