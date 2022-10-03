
library(ggplot2)
library(reshape2)
library(vegan)
library(ggpubr)
library(rstatix)
library(mixOmics)
library(gplots)
library(factoextra)

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
  primary.dat["R0446","Hypertension.Staging"] <- 1
  primary.dat <- primary.dat[primary.dat$Anti.hypertensive.Agents==0,]
  
  secondary.dat <- read.csv("Data/Secondary_data.csv", header = T,row.names = 1)
  colnames(secondary.dat) <- gsub("\\.\\..*", "", colnames(secondary.dat))
  rownames(secondary.dat) <- toupper(rownames(secondary.dat))
  
  blood.stool.dat <- read.csv("Data/Blood_stool_biomarkers.csv", header = T,row.names = 1)
  colnames(blood.stool.dat) <- gsub("\\.\\..*", "", colnames(blood.stool.dat))
  rownames(blood.stool.dat) <- toupper(rownames(blood.stool.dat))
  colnames(blood.stool.dat)[6] <- "Serum.Vitamin.D.Classification"
  blood.stool.dat[] <- lapply(blood.stool.dat, function(x) as.numeric(as.character(x)))
  #Removing outliers
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
  
  #Pathway data
  path.data <- read.csv("Data/humann3_ko_unstratified.tsv", sep = "\t", header = T, row.names = 1)
  colnames(path.data) <- substr(colnames(path.data), 1, 5); colnames(path.data) <- toupper(colnames(path.data))
  colnames(path.data) <- gsub("R(\\d*)\\.","R0\\1",colnames(path.data))
  common <- intersect(colnames(path.data),rownames(primary.dat))
  path.data <- path.data[,common]
  
  # Mapping ko ids
  ko_mapping <- read.csv("Data/full_KO_table.csv")
  
  grouped.data <- path.data
  grouped.data$group <- ko_mapping[match(rownames(grouped.data),ko_mapping$KO.ID.D),"KO.Desc.C"]
  grouped.data <- aggregate(. ~ group, data = grouped.data, FUN = sum)
  rownames(grouped.data) <- grouped.data$group; grouped.data <- grouped.data[,-1]
  grouped.data <- norm_n_filter(grouped.data)
  grouped.data <- grouped.data[rownames(grouped.data),]
  transformed.grouped.data <- asin(sqrt(grouped.data))
  
  meta.data$Sex <- ifelse(meta.data$Sex == 0, "Male", 'Female')
  male <- rownames(meta.data[meta.data$Sex=="Male",])
  female <- rownames(meta.data[meta.data$Sex=="Female",])
  meta.data$Hypertension.Staging <- as.factor(meta.data$Hypertension.Staging)
  meta.data$Masked.hypertension <- as.factor(meta.data$Masked.hypertension)
  meta.data$Smoking <- as.factor(meta.data$Smoking)
  meta.data$Menopause_code <- as.factor(meta.data$Menopause_code)
  
  merged.dat <- cbind(meta.data,transformed.grouped.data[rownames(meta.data),])
}

###Plotting functions###
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
  print(ggboxplot(dat,"lab", "y", color = "lab") + stat_pvalue_manual(stat.test, label = "p.adj", hide.ns = hide_ns) + 
          theme(legend.position = "none") + xlab(label) + ylab(ylab) + ggtitle(title) + facet_wrap(~Sex))
}
my_scatterplot <- function(var,label, ylab="", title="", pos="bottom"){
  dat <- data.frame(y = var, lab = meta.data[,label], Sex = meta.data$Sex)
  all.facet <- dat
  all.facet$Sex <- "All"
  dat <- rbind(dat,all.facet)
  dat <- dat[!is.na(dat$lab),]
  print(ggscatter(dat,"lab", "y", add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = T) + facet_wrap(~Sex) + 
          xlab(label) + ylab(ylab) + ggtitle(title) + stat_cor(method = "pearson", label.x.npc = "left", label.y.npc =pos, col = "red"))
}
plot_MDS <- function(dat,label, title = "",cat=T){
  met <- meta.data[!is.na(meta.data[,label])&!is.na(meta.data$Sodium.Intake.based.on.spot.urine),]
  dat <- dat[rownames(met),]
  males <- rownames(met[met$Sex=="Male",])
  females <- rownames(met[met$Sex=="Female",])
  p.nova.all <- adonis(dat ~ met[,label] + met[,"Sex"] + met[,"Age"] + met[,"BMI"] + met[,"Smoking"] + met[,"Sodium.Intake.based.on.spot.urine"], permutations=999, method = "bray")
  p.nova.male <- adonis(dat[males,] ~ met[males,label] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"], permutations=999, method = "bray")
  p.nova.female <- adonis(dat[females,] ~ met[females,label] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"], permutations=999, method = "bray")
  p.nova.all <- p.nova.all$aov.tab[1,"Pr(>F)"]
  p.nova.male <- p.nova.male$aov.tab[1,"Pr(>F)"]
  p.nova.female <- p.nova.female$aov.tab[1,"Pr(>F)"]
  p.labels <- data.frame(label = paste0("Permanova: ",c(p.nova.all,p.nova.male,p.nova.female)), gender= c("All","Male","Female"))
  
  bray.dist.matrix <- vegdist(dat[rownames(met),], method="bray")
  mds <- cmdscale(bray.dist.matrix, eig = TRUE, x.ret = TRUE)
  mds.var.per <- round(mds$eig/sum(mds$eig)*100, 1)
  mds.values <- mds$points
  mds.data <- data.frame(Patient = rownames(mds.values), X = mds.values[,1], Y = mds.values[,2], var = meta.data[match(rownames(mds.values),rownames(meta.data)),label], gender=meta.data[match(rownames(mds.values),rownames(meta.data)),"Sex"])
  all.facet <- mds.data
  all.facet$gender <- "All"
  mds.data <- rbind(mds.data,all.facet)
  if(cat){
    ggplot(mds.data, aes(x=X, y=Y, col=as.factor(var))) + geom_point(size = 5) + theme_classic() + labs(color=label) + facet_wrap(~gender) +
      xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) + ggtitle(title) +
      geom_text(data = p.labels, mapping = aes(x = -Inf, y = -Inf, label = label), hjust = -0.1, vjust = -1, col = "red", size = 5)
  } else {
    ggplot(mds.data, aes(x=X, y=Y, col=var)) + geom_point(size = 5) + theme_classic() + labs(color=label) + facet_wrap(~gender) +
      xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) + ggtitle(title) +
      geom_text(data = p.labels, mapping = aes(x = -Inf, y = -Inf, label = label), hjust = -0.1, vjust = -1, col = "red", size = 5)
  }
}
plot_sig_feat <- function(dat, label, cat=T){
  for(bac in colnames(dat)){
    d <- data.frame(y = dat[,bac], lab = meta.data[,label], Sex = meta.data$Sex)
    d <- d[!is.na(d$lab),]
    all.facet <- d
    all.facet$Sex <- "All"
    d <- rbind(d,all.facet)
    if(cat){
      stat.test <- d %>% group_by(Sex) %>% pairwise_wilcox_test(y ~ lab, p.adjust.method = "BH") %>% add_y_position()
      if(length(unique(stat.test$p.adj.signif))>1|unique(stat.test$p.adj.signif)!="ns"){
        print(ggboxplot(d,"lab", "y", color = "lab") + stat_pvalue_manual(stat.test, label = "p.adj", hide.ns = T) 
              + theme(legend.position = "none") + xlab(label) + ylab(bac) + facet_wrap(~Sex))
      }
    } else {
      all.test <- cor.test(d[d$Sex=="All","y"], d[d$Sex=="All","lab"], method = "pearson")
      m.test <- cor.test(d[d$Sex=="Male","y"], d[d$Sex=="Male","lab"], method = "pearson")
      f.test <- cor.test(d[d$Sex=="Female","y"], d[d$Sex=="Female","lab"], method = "pearson")
      if(m.test$p.value <= 0.05|f.test$p.value <= 0.05|all.test$p.value <= 0.05){
        my_scatterplot(d$y, label, ylab = bac)
      }
    }
  }
}
plot_univar_sig <- function(dat,label,p_val_filt=0.05){
  males <- rownames(meta.data[meta.data$Sex=="Male",])
  females <- rownames(meta.data[meta.data$Sex=="Female",])
  res.table <- data.frame(gender = character(), bacteria = character(), m1 = character(), m2 = character(), m3 = character())
  for(bac in colnames(dat)){
    m1.all <-lm(dat[,bac] ~ meta.data[,label])
    m2.all <-lm(dat[,bac] ~ meta.data[,label] + meta.data[,"Sex"] + meta.data[,"Age"] + meta.data[,"BMI"])
    m3.all <-lm(dat[,bac] ~ meta.data[,label] + meta.data[,"Sex"] + meta.data[,"Age"] + meta.data[,"BMI"] + meta.data[,"Smoking"] + meta.data[,"Anti.hypertensive.Agents"])
    m1.male <-lm(dat[males,bac] ~ meta.data[males,label])
    m2.male <-lm(dat[males,bac] ~ meta.data[males,label] + meta.data[males,"Age"] + meta.data[males,"BMI"])
    m3.male <-lm(dat[males,bac] ~ meta.data[males,label] + meta.data[males,"Age"] + meta.data[males,"BMI"] + meta.data[males,"Smoking"] + meta.data[males,"Anti.hypertensive.Agents"])
    m1.female <-lm(dat[females,bac] ~ meta.data[females,label])
    m2.female <-lm(dat[females,bac] ~ meta.data[females,label] + meta.data[females,"Age"] + meta.data[females,"BMI"])
    m3.female <-lm(dat[females,bac] ~ meta.data[females,label] + meta.data[females,"Age"] + meta.data[females,"BMI"] + meta.data[females,"Smoking"] + meta.data[females,"Anti.hypertensive.Agents"])
    m1.all <- summary(m1.all)$coefficients[2,4]; m2.all <- summary(m2.all)$coefficients[2,4]; m3.all <- summary(m3.all)$coefficients[2,4]
    m1.male <- summary(m1.male)$coefficients[2,4]; m2.male <- summary(m2.male)$coefficients[2,4]; m3.male <- summary(m3.male)$coefficients[2,4]
    m1.female <- summary(m1.female)$coefficients[2,4]; m2.female <- summary(m2.female)$coefficients[2,4]; m3.female <- summary(m3.female)$coefficients[2,4]
    res.table[nrow(res.table)+1,] <- c("All",bac,m1.all,m2.all,m3.all)
    res.table[nrow(res.table)+1,] <- c("Male",bac,m1.male,m2.male,m3.male)
    res.table[nrow(res.table)+1,] <- c("Female",bac,m1.female,m2.female,m3.female)
  }
  filt.feat <- unique(res.table[res.table$m1<p_val_filt|res.table$m2<p_val_filt|res.table$m3<p_val_filt,"bacteria"])
  if(identical(filt.feat, character(0))){return("No significant features")}
  res.table <- res.table[res.table[,"bacteria"]%in%filt.feat,]
  res.table <- tidyr::gather(res.table, "model", "p", 3:5)
  res.table$sig <- ifelse(res.table$p>0.05, "Non-significant", "Significant")
  res.table$sym <- ifelse(res.table$p>0.05, "", ifelse(res.table$p<0.01,"**", "*"))
  print(ggplot(res.table, aes(x=model, y = bacteria, fill=sig)) + geom_tile(color="white", size=0.1) + scale_fill_manual(values = c("blue3", "red3"), name="p value") + 
    coord_equal() + xlab("") + ylab("") + ggtitle(label) + geom_text(aes(label=sym), size = 5, col = "white") + facet_wrap(~gender))
}
splsda_pathway <- function(var,X){
  female <- meta.data$Sex == "Female"
  male <- meta.data$Sex == "Male"
  
  splsda.res.all <- splsda(X, var, ncomp=2, keepX=c(15,10))
  splsda.res.m <- splsda(X[male,], var[male], ncomp=2, keepX=c(15,10))
  splsda.res.f <- splsda(X[female,], var[female], ncomp=2, keepX=c(15,10))
  
  plotIndiv(splsda.res.all, group = var,title = "All cohort", legend = T)
  plotIndiv(splsda.res.m, group = var[male], title = "Male", legend = T)
  plotIndiv(splsda.res.f, group = var[female], title = "Female", legend = T)
  
  plotLoadings(splsda.res.all, contrib = 'max', method = 'mean', comp = 1, title = "Contribution comp1: All",legend = F, layout = c(2,4))
  plotLoadings(splsda.res.m, contrib = 'max', method = 'mean', comp = 1, title = "Contribution comp1: Male",legend = F)
  plotLoadings(splsda.res.f, contrib = 'max', method = 'mean', comp = 1, title = "Contribution comp1: Female")
  plotLoadings(splsda.res.all, contrib = 'max', method = 'mean', comp = 2, title = "Contribution comp2: All",legend = F)
  plotLoadings(splsda.res.m, contrib = 'max', method = 'mean', comp = 2, title = "Contribution comp2: Male",legend = F)
  plotLoadings(splsda.res.f, contrib = 'max', method = 'mean', comp = 2, title = "Contribution comp2: Female")
}
plot_sig_heatmap <- function(dat,title){
  dat$pval <- p.adjust(dat$pval, method='BH')
  dat$sym <- ifelse(dat$pval>0.05, "", "*")
  filt.feat1 <- names(which(tapply(dat$pval,dat$feature1,function(x){ sum(x<0.05)>0})))
  filt.feat2 <- names(which(tapply(dat$pval,dat$feature2,function(x){ sum(x<0.05)>0})))
  dat <- dat[dat$feature1%in%filt.feat1&dat$feature2%in%filt.feat2,]
  dat$log10 <- log10(dat$pval)*-as.numeric(dat$direction)
  ggplot(dat, aes(x=feature1, y = feature2, fill=log10)) + geom_tile(color="white", size=0.1)+scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
    coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 4, col = "black") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(title)
}
plot_sig_network <- function(dat,title, pval = 0.01){
  dat$pval <- p.adjust(dat$pval, method='bonferroni')
  dat$direction <- ifelse(dat$direction==1,"green","red")
  g <- graph_from_edgelist(as.matrix(dat[,1:2]), directed = F)
  E(g)$weight <- dat$pval
  E(g)$color <- dat$direction
  g <- delete.edges(g, E(g)[ weight > pval ])
  g <- delete.vertices(g, which(degree(g)==0))
  plot(g, vertex.label.cex = 0.8, vertex.size=5, vertex.label.color="black", edge.width = 2, layout=layout_with_mds, main = title)
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
    spe.dat <- transformed.grouped.data[rownames(met),]
    males <- rownames(met[met$Sex=="Male",])
    females <- rownames(met[met$Sex=="Female",])
    
    m1.all <- adonis(spe.dat ~ met[,feature], met, permutations=999, method = "bray")
    m1.female <- adonis(spe.dat[females,] ~ met[females,feature], permutations=999, method = "bray")
    m1.male <- adonis(spe.dat[males,] ~ met[males,feature], permutations=999, method = "bray")
    
    m2.all <- adonis(spe.dat ~ met[,feature] + Sex + Age + BMI, met, permutations=999, method = "bray")
    m2.female <- adonis(spe.dat[females,] ~ met[females,feature] + met[females,"Age"] + met[females,"BMI"], permutations=999, method = "bray")
    m2.male <- adonis(spe.dat[males,] ~ met[males,feature] + met[males,"Age"] + met[males,"BMI"], permutations=999, method = "bray")
    
    m3.all <- adonis(spe.dat ~ met[,feature] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine, met, permutations=999, method = "bray")
    m3.female <- adonis(spe.dat[females,] ~ met[females,feature] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"], permutations=999, method = "bray")
    m3.male <- adonis(spe.dat[males,] ~ met[males,feature] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"], permutations=999, method = "bray")
    
    m4.all <- adonis(spe.dat ~ met[,feature] + Sex + Age + BMI + Smoking + Sodium.Intake.based.on.spot.urine + Fatty.liver.by.CAP.score, met, permutations=999, method = "bray")
    m4.female <- adonis(spe.dat[females,] ~ met[females,feature] + met[females,"Age"] + met[females,"BMI"] + met[females,"Smoking"] + met[females,"Sodium.Intake.based.on.spot.urine"] + met[females,"Menopause_code"] + met[females,"Fatty.liver.by.CAP.score"], permutations=999, method = "bray")
    m4.male <- adonis(spe.dat[males,] ~ met[males,feature] + met[males,"Age"] + met[males,"BMI"] + met[males,"Smoking"] + met[males,"Sodium.Intake.based.on.spot.urine"] + met[males,"Fatty.liver.by.CAP.score"], permutations=999, method = "bray")
    
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
diet_variables_continuous <- c(names(which(colSums(is.na(diet.dat))!=nrow(diet.dat))),"med_fish","med_original","DASH_mullen","Selenium_ug")
diet_variables_categorical <- c("med_original_unhealthy","DASH_mullen_unhealthy")
food_group_variables <- names(which(colSums(is.na(food.groups.dat))!=nrow(food.groups.dat)))

################################
###Beta-diversity - PERMANOVA###
################################
#Vascular parameters
vasc.perma <- permanova_models(meta.data,c(vascular_markers_categorical,vascular_markers_continuous))
vasc.perma$pval <- as.numeric(vasc.perma$pval)
vasc.perma$sym <- ifelse(vasc.perma$pval>=0.05, "", ifelse(vasc.perma$pval<=0.01,ifelse(vasc.perma$pval<=0.001,"***", "**"), "*"))
vasc.perma$log <- -log10(vasc.perma$pval)
ggplot(vasc.perma, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("Vascular variables")
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
ggplot(diet.perma, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 12, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold"))
#Dietary variables subset
diet.perma <- permanova_models(meta.data,c('HEI_total','fiber_g','diet_sodium_mg','Energy_kcal','alcohol_g'))
diet.perma$pval <- as.numeric(diet.perma$pval)
diet.perma$sym <- ifelse(diet.perma$pval>=0.05, "", ifelse(diet.perma$pval<=0.01,ifelse(diet.perma$pval<=0.001,"***", "**"), "*"))
diet.perma$log <- -log10(diet.perma$pval)
ggplot(diet.perma, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = '-log10',midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold"))
#Food groups variables
food.groups.perma <- permanova_models(meta.data,food_group_variables)
food.groups.perma$pval <- as.numeric(food.groups.perma$pval)
food.groups.perma$sym <- ifelse(food.groups.perma$pval>=0.05, "", ifelse(food.groups.perma$pval<=0.01,ifelse(food.groups.perma$pval<=0.001,"***", "**"), "*"))
food.groups.perma$log <- -log10(food.groups.perma$pval)
ggplot(food.groups.perma, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + ggtitle("Food group variables")
