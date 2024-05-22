library(limma)
library(Mus.musculus)
library(DESeq2)
library(edgeR)
library(tidyverse)


## Getting PCA from the DEG results not LIMMA

# Expression values
dlNorm <-  read.csv("brain/AMY_counts.csv", row.names = 1)
#remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm)[c(1:67)] <- substr(colnames(dlNorm)[c(1:67)], 7, 13)
#batch, cageno, mouseno.

#Group traits
#Getting metadata ready 
coldata <- read_csv("brain/sample70min_table.csv")
head(coldata)
str(coldata)

#fixing things
coldata$region <- gsub("mPF", "PFC", coldata$region)
coldata$groupEX <- coldata$group


# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)

#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)
# still have 2s and 3s in here, but will filter later.

#just get samples I want
coldata$SampleID <- substr(coldata$SampleName, 7, 13)

# coldata <- coldata %>%  filter(post_idbatch!= "3-4Batch12") for mPFC

coldata <- coldata %>%  dplyr::select(-SampleName) %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != 'ascenders') %>% 
  pivot_wider(
    names_from = region, 
    values_from = region)

row.names <- coldata$SampleID
row.names(coldata) <- row.names #Assigning row names from as sample names  
head(coldata)
rownames(coldata)


#check before normalizing 
dlNorm<- dlNorm[, rownames(coldata)]
all(rownames(coldata) == colnames(dlNorm))

#normalize and filter with all groups 

#dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),] #this would remove any row with at least one NA in it.

d = apply(dlNorm, 2, as.numeric)
dim(d)

#DGEList from limma
d0= DGEList(d, group = coldata$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm) #redo rownames as get lost when making matrix
d0 <- calcNormFactors(d0) #Normalizing step.


## Removing genes from PCA where average count is below cutoff.
cutoff <- 25
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 6571   40 # for counts of 50 
# 1] 13198    40 # for counts of 10 
# [1] 9798   20 # for counts of 25 or more
# Now take out groups that you want
#DOMs first 
dge.dl$samples$group

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("CDOM", "DOM", "DES")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

#make sure coldata matches
coldata %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB") %>%
  dplyr::select(SampleID, condition1) -> var_info  

row.names <- var_info$SampleID

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

#refilter, check that they match.
dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check


##making it a factor
var_info$condition1 %>%
  factor(.,levels = c("CDOM","DOM","DES")) -> group.dl

## Make the design/contrasts
design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

## transform to logcpm counts
v = voom(dge.dl, design.dl, plot = F)


colData <- v$design %>% data.frame()
colData <- var_info %>% rownames_to_column(var = "SampleName") %>% cbind(colData) 

#get expression levels
rv <- rowVars((v$E)) #measure of variation for each gene
select <- order(rv, decreasing = TRUE)[1:100] #top 100 varying genes
pca1 <- prcomp(t((v$E)[select, ])) #do pca

#get groups
condition3.df <- as.data.frame(colData[, 'condition1', 
                                       drop = FALSE])



#put results into df
d1 <- data.frame(PC1 = pca1$x[, 1], PC2 = pca1$x[, 2], PC3 = pca1$x[, 3], PC4 = pca1$x[, 4], 
                 condition3.df, name = colnames(v$E))



head(d1)

# make scree plot
pv1 <- ((pca1$sdev^2) / (sum(pca1$sdev^2)))*100
barplot(pv1, cex.names=1, xlab=paste("Principal component (PC), 1-", length(pca1$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

var1 <- ((pca1$sdev[1]^2) / (sum(pca1$sdev^2)))*100
var2 <- ((pca1$sdev[2]^2) / (sum(pca1$sdev^2)))*100
var3 <- ((pca1$sdev[3]^2) / (sum(pca1$sdev^2)))*100
var4 <- ((pca1$sdev[4]^2) / (sum(pca1$sdev^2)))*100

## actually plots 

xlab=paste("PC1, ", round(pv1[1], 2), "%")
ylab=paste("PC2, ", round(pv1[2], 2), "%")
xlab3=paste("PC3, ", round(pv1[3], 2), "%")
ylab4=paste("PC4, ", round(pv1[4], 2), "%")


p1=ggplot(d1, aes(PC1, PC2, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = ylab)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()

cdom=ggplot(d1, aes(PC1, PC2, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = ylab)+
  labs(color = "Social\nCondition",shape = "Social\nCondition", labels = "Social\nCondition")+
  scale_color_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+
  theme_classic()+
  theme(text =element_text(size = 15))
cdom

ggsave("manuscript/brain/results_figures/CDOM_PCA_top100.png", width =5,height = 4, dpi = 300)

#########
# subs
# Expression values
dlNorm <-  read.csv("brain/AMY_counts.csv", row.names = 1)
#remove zeros - Becca told me to do this since zeros can fuck things up. 
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]

#Lines 20 to 70 is just me getting my group variables together
#trim sample ids
colnames(dlNorm)[c(1:67)] <- substr(colnames(dlNorm)[c(1:67)], 7, 13)

#Group traits
#Getting metadata ready 
coldata <- read_csv("brain/sample70min_table.csv")
head(coldata)
str(coldata)

#fixing things
coldata$region <- gsub("mPF", "PFC", coldata$region)
coldata$groupEX <- coldata$group


# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)

#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)


#just get samples I want
coldata$SampleID <- substr(coldata$SampleName, 7, 13)

coldata <- coldata %>%  dplyr::select(-SampleName) %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != 'ascenders') %>% 
  pivot_wider(
    names_from = region, 
    values_from = region)

row.names <- coldata$SampleID
row.names(coldata) <- row.names #Assigning row names from as sample names  
head(coldata)

#check before normalizing 
dlNorm<- dlNorm[, rownames(coldata)]
all(rownames(coldata) == colnames(dlNorm))

#normalize and filter with all groups -
# I normalize with all groups since I want to compare them after

dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0= DGEList(d, group = coldata$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <-25
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
#6571   40
# [1] 9798   20
# Now take out groups that you want
#DOMs first 
dge.dl$samples$group

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("CSUB", "SUB", "ASC")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

coldata %>% 
  filter(condition1 != "DES") %>%
  filter(condition1 != "DOM") %>%
  filter(condition1 != "CDOM") %>%
  dplyr::select(SampleID, condition1) -> var_info  

row.names <- var_info$SampleID

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code
var_info$condition1 %>%
  factor(.,levels = c("CSUB","SUB","ASC")) -> group.dl


design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)


colData <- v.dl$design %>% data.frame()
colData <- var_info %>% rownames_to_column(var = "SampleName") %>% cbind(colData) 

rv <- rowVars((v.dl$E))
select <- order(rv, decreasing = TRUE)[1:100]
pca1 <- prcomp(t((v.dl$E)[select, ]))

condition3.df <- as.data.frame(colData[, 'condition1', 
                                       drop = FALSE])




d1 <- data.frame(PC1 = pca1$x[, 1], PC2 = pca1$x[, 2], PC3 = pca1$x[, 3], PC4 = pca1$x[, 4], 
                 condition3.df, name = colnames(v.dl$E))

d1$condition1

head(d1)

pv1 <- ((pca1$sdev^2) / (sum(pca1$sdev^2)))*100
barplot(pv1, cex.names=1, xlab=paste("Principal component (PC), 1-", length(pca1$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

var1 <- ((pca1$sdev[1]^2) / (sum(pca1$sdev^2)))*100
var2 <- ((pca1$sdev[2]^2) / (sum(pca1$sdev^2)))*100
var3 <- ((pca1$sdev[3]^2) / (sum(pca1$sdev^2)))*100
var4 <- ((pca1$sdev[4]^2) / (sum(pca1$sdev^2)))*100

## actually plots 

xlab=paste("PC1, ", round(pv1[1], 2), "%")
ylab=paste("PC2, ", round(pv1[2], 2), "%")

xlab3=paste("PC3, ", round(pv1[3], 2), "%")
ylab4=paste("PC4, ", round(pv1[4], 2), "%")
d1$condition1

p1s=ggplot(d1, aes(PC1, PC2, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = ylab)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()

csub=ggplot(d1, aes(PC1, PC2, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = ylab)+
  labs(color = "Social\nCondition",shape = "Social\nCondition", labels = "Social\nCondition")+
  scale_color_manual(values = c("#de7065ff","#f7cb44", "#13306dff"))+
  theme_classic()+
  theme(text =element_text(size = 15))

csub
ggsave("manuscript/brain/results_figures/CSUB_PCA_top100.png", width =5,height = 4, dpi = 300)
