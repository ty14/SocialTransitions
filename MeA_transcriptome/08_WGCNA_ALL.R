# libraries 
library(limma)
library(edgeR)
library(WGCNA)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
library(clusterProfiler)
library(enrichplot)
library(organism, character.only = TRUE)
library(DOSE)
library(tidyverse)
grcm38 # mouse genes


# Expression values
dlNorm <-  read.csv("brain/AMY_counts.csv", row.names = 1)
#remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm)[c(1:67)] <- substr(colnames(dlNorm)[c(1:67)], 7, 13)

#Group traits
#Getting metadata ready 
coldata <- read.csv("manuscript/brain/results_tables/coldata.csv", row.names = 1)

#check before normalizing 
dlNorm<- dlNorm[, rownames(coldata)]
all(rownames(coldata) == colnames(dlNorm))

#normalize and filter with all groups 

dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0= DGEList(d, group = coldata$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

#won used 10 in 90% of samples for brain paper, which is what tutorial suggest. 
#in liver paper she did 50 in 90% samples 
cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
#13198    40 samples = 10 

# Now take out groups that you want
dge.dl$samples$group


coldata %>% 
  dplyr::select(SampleID, condition1) -> var_info  

row.names <- var_info$SampleID

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code
var_info$condition1 %>%
  factor(.,levels = c("CDOM", "DOM", "DES","ASC", "SUB", "CSUB")) -> group.dl


design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)


# Many functions expect the matrix to be transposed
datExpr <- t(v.dl$E)
## check rows/cols
nrow(datExpr)
ncol(datExpr)
rownames(datExpr)

# getting trait data
colnames(coldata)
amy <- coldata[c(2,4,5,10:12,15:19)]


ifelse(amy$condition1 == "CDOM", 3, amy$condition1) -> amy$condition2
ifelse(amy$condition2 == "DOM", 2, amy$condition2) -> amy$condition2
ifelse(amy$condition2 == "DES", 1, amy$condition2) -> amy$condition2
ifelse(amy$condition1 == "ASC", 0, amy$condition2) -> amy$condition2
ifelse(amy$condition2 == "SUB", -1, amy$condition2) -> amy$condition2
ifelse(amy$condition2 == "CSUB", -2, amy$condition2) -> amy$condition2


amy1 <- amy %>% 
  dplyr::select(-condition1)

#make everything numeric 
amy1[1:11] <- lapply(amy1[1:11], as.numeric)

#saving to use for other WGCNA analysis 
write.csv(amy1, "manuscript/brain/results_tables/TraitsforWGCNA.csv", row.names = T)


sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 60, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 60, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExprx = datExpr[keepSamples, ]
nGenes = ncol(datExprx)
nSamples = nrow(datExprx)


# Re-cluster samples with out any samples filtered 
var_info
amy1x <- amy1[c(1:8,10,12:16,18:40),]
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(amy1, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(amy1),
                    main = "AMY dendrogram and trait heatmap(all samples)")

collectGarbage()
#not going to filter because ascenders rreally are not grouping anyways?
saveRDS(datExpr,"manuscript/brain/results_RDS/WGCNA_datExpr_ALL10.RDS")
# saveRDS(datExprx,"manuscript/brain/results_RDS/WGCNA_datExpr_ALL_outlierRemoved.RDS")


# Run WGCNA for each dom, descender, controls 
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, )
# Plot the results:
print(sft)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;


# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main ="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");


# this line corresponds to using an R^2 cut-off of h

abline(h=0.85,col="red")  #changed to 0.8 (originally 0.9 in the tutorial)
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft)



# TOMType = "signed", networkType = "signed",
# TOMType = "unsigned", networkType = "unsigned",
# TOMType = c("signed", networkType = "signed hybrid")
cor <- WGCNA::cor
WGCNA_get_net <- function(my_tissue = "MEA70min",
                          my_power =  4, 
                          my_TOMType ="signed", 
                          my_networkType = "signed hybrid"){
  
  x <- readRDS("manuscript/brain/results_RDS/WGCNA_datExpr_ALL10.RDS") 
  set.seed(312)
  net = blockwiseModules(x, 
                         power = my_power,
                         TOMType = my_TOMType, 
                         networkType = my_networkType,
                         minModuleSize = 50,
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)
  
  
  saveRDS(net, "manuscript/brain/results_RDS/WGCNA_net_ALL_Power4_10.RDS")
  
}

WGCNA_get_net("MEA70min", 4, "signed", "signed hybrid")
# ========================================================================================

net <- readRDS(glue::glue("manuscript/brain/results_RDS/WGCNA_net_ALL_Power4_10.RDS"))

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath

dev.off() # make sure you do this before AND after 
png(file = "manuscript/brain/results_figures/cluster_dendo_ALL_Power4_10.png",
    width=600, height=350)
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "MEA ALL Power4 - cutoff 10")

dev.off()

MEs = net$MEs
##################################
set.seed(312)
datExpr <- readRDS("manuscript/brain/results_RDS/WGCNA_datExpr_ALL10.RDS") 
net <- readRDS("manuscript/brain/results_RDS/WGCNA_net_ALL_Power4_10.RDS")
dim(datExpr)

#PART 1
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
rownames(datExpr)


# getting trait data
amy1 <- read.csv("manuscript/brain/results_tables/TraitsforWGCNA.csv",row.names = 1 )
# Re-cluster samples with out any samples filtered 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(amy1, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(amy1),
                    main = "AMY dendrogram and trait heatmap CSUB")

datTraits <- amy1





#PART 2
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

saveRDS(MEs,"manuscript/brain/results_RDS/WGCNA_ALL_MEs_Power4.RDS")

#make corr matrix for paper
datTraitsx <- datTraits %>% select(SocialCondition = condition2)

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraitsx, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

names(MEs) <- gsub("ME", "", names(MEs))
#PART3
 sizeGrWindow(20,12)

dev.off() # make sure you do this before AND after

png(file = "manuscript/brain/results_figures/module_trait_ALL_Power4_cut10.png",
    width=1000, height=2600, res = 300)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraitsx),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "")
#
## remember - they are correlation, not regression

dev.off() 
#PART 4 =====================================================================================
# Define variable David's score containing the David's score column of datTrait
status = as.data.frame(datTraits$condition2);
names(status) = "trait"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, status, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(status), sep="");
names(GSPvalue) = paste("p.GS.", names(status), sep="");


geneModuleMembership %>% 
  rownames_to_column("ensgene") %>% 
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> gene_MM_TS_all



# to make holistic freaking dataframe ==================================================
gene_MM_TS_all$module <- moduleColors


get_mm <- function(x){
  x$moduleMembership <- x[colnames(x) == paste("MM",x$module,sep = "")] %>%
    unlist %>%
    as.numeric
  xx <- x %>%
    dplyr::select(ensgene,module,moduleMembership,GS.trait)
  return(x)
}


wgcna_whole <- get_mm(gene_MM_TS_all[1,])

for(i in 2:nrow(gene_MM_TS_all)){
  wgcna_whole <- rbind(wgcna_whole,get_mm(gene_MM_TS_all[i,]))
}

wgcna_whole %>%
  rename(GS.status = GS.trait) -> wgcna_whole


wgcna_whole %>%
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> wgcna_all

saveRDS(wgcna_all, "manuscript/brain/results_RDS/WGCNA_ALL_WGCNA_MM_GS_all_Power4.RDS")

##boxplot of modules 

MEs <- readRDS("manuscript/brain/results_RDS/WGCNA_ALL_MEs_Power4.RDS")
df <- read.csv("manuscript/brain/results_tables/coldata.csv", row.names = 1)
ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:27, names_to = "Module") %>% 
  full_join(df) %>% filter(Module != "MEgrey")

head(ME_df)

ME_df$condition1
ME_df <- ME_df %>%
  mutate(status = condition1)
ME_df$Module <- gsub("ME", "", ME_df$Module)

ME_df$status <- factor(ME_df$status, levels = c("DES","DOM","CDOM", "ASC", "SUB", "CSUB"))

# pur <- ME_df %>% filter(Module == "purple")

#boxplot of eigengenes and status 
source("functions/geom_boxjitter.R")


p1<- ME_df %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(6))+
  scale_fill_manual(values = viridis::viridis(6))+
  facet_wrap(~Module, scales = "free_y")+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "MeA")+ theme_classic()+ theme(legend.position = "none")

p1


ggsave(filename = "manuscript/brain/results_figures/ALL_eigengene_boxplotPower4.png",
       p1,
       height = 18, width = 22, dpi = 300)


ME_df %>%
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se = F,  alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::viridis(6))+
  scale_fill_manual(values = viridis::viridis(6))+
  facet_wrap(~Module, scales = "free_y")+
  theme(legend.position = "top")+
  labs(x = "Normalized CORT",
       y = "Module eigengene") + theme_classic() -> p2
p2

ggsave(filename = "manuscript/brain/results_figures/ALL_eigengene_CORTcorrelation_Power4.png",
       p2,
       height = 18, width = 22, dpi = 300)


p1sub <- ME_df %>% filter(status %in% c("CSUB", "ASC", "SUB")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y")+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "MeA")+ theme_classic()+ theme(legend.position = "none")

p1sub


ggsave(filename = "manuscript/brain/results_figures/ALL_CSUB_eigengene_boxplotPower4.png",
       p1sub,
       height = 18, width = 22, dpi = 300)

p1dom <- ME_df %>% filter(status %in% c("CDOM", "DES", "DOM")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y")+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "MeA")+ theme_classic()+ theme(legend.position = "none")

p1dom


ggsave(filename = "manuscript/brain/results_figures/ALL_CDOM_eigengene_boxplotPower4.png",
       p1dom,
       height = 18, width = 22, dpi = 300)


ME_df %>% filter(status %in% c("CSUB", "ASC", "SUB")) %>% 
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se = F,  alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y")+
  theme(legend.position = "top")+
  labs(x = "Normalized CORT",
       y = "Module eigengene") + theme_classic() -> p2sub
p2sub

ggsave(filename = "manuscript/brain/results_figures/ALL_CSUB_eigengene_CORTcorrelation_Power4.png",
       p2sub,
       height = 18, width = 22, dpi = 300)

ME_df %>% filter(status %in% c("CDOM", "DES", "DOM")) %>% 
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se = F,  alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y")+
  theme(legend.position = "top")+
  labs(x = "Normalized CORT",
       y = "Module eigengene") + theme_classic() -> p2dom
p2dom

ggsave(filename = "manuscript/brain/results_figures/ALL_CDOM_eigengene_CORTcorrelation_Power4.png",
       p2dom,
       height = 18, width = 22, dpi = 300)


# Module of interest via significance 
my_trait = "Status"
module = "greenyellow"
module = "pink"
module = "cyan"
module = "green"
module = 'turquoise'
module = 'darkturquoise'
module = 'salmon'
module = 'brown'
module = 'midnightblue'
module = 'darkgreen'
module = 'black'
module = 'blue'
module = "darkgrey"
module = "lightgreen"
module = "royalblue"
module = "lightcyan"
module = "yellow"
module = "darkturquoise"
module = "darkred"
module = "orange"
module_list = c("greenyellow", "pink", "cyan", "green", "darkturquoise", "salmon", "brown", "midnightblue",
                "darkgreen", "black", "blue", "darkgrey", "lightgreen", "royalblue", "lightcyan", "yellow", 
                "darkred", "orange")
hub_gene_list = vector('list', length = length(module_list))
names(hub_gene_list) <- module_list


column = base::match(module, modNames);
moduleGenes = moduleColors==module;
colnames(datExpr)[moduleColors==module] -> module_gene
# 
# 
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = glue::glue("Gene significance for {my_trait}"),
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# datExpr[moduleGenes,]

my_MM_threshold = 0.8
my_GS_threshold = 0.2


for (module in module_list){
  column = base::match(module, modNames);
  moduleGenes = moduleColors==module;
  
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  grcm38 %>% 
    dplyr::select(ensgene, symbol, chr, description) %>% 
    filter(ensgene %in% module_gene) -> module_gene_info
  
  
  
  geneModuleMembership %>% 
    rownames_to_column("ensgene") %>% 
    left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> gene_MM_TS
  
  
  
  gene_MM_TS %>% 
    filter(ensgene %in% module_gene) %>% 
    dplyr::select(ensgene, glue::glue("MM{module}"), GS.trait) %>% 
    filter(abs(GS.trait) >= my_GS_threshold) -> x
  
  x[x[,glue::glue("MM{module}")]>my_MM_threshold,] -> hub_genes
  
  hub_genes %>% 
    left_join(module_gene_info) %>% 
    mutate(moduleName = glue::glue("{module}")) %>% 
    rename(moduleMembership = glue::glue("MM{module}")) -> hub_genes
  hub_gene_list[[module]] <- hub_genes
  
}



hub_gene_list %>% 
  do.call(rbind,.)%>% 
  unique() -> hubgenes_df

hubgenes_df %>% 
  dplyr::arrange(desc(GS.trait)) %>% 
  dplyr::select(-ensgene, -chr) %>% 
  head(15)

write.csv(hubgenes_df, "manuscript/brain/results_tables/WCGNA_ALL_hubgene_list_Power4.csv", row.names = F)


ADJ=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ, moduleColors) 

Alldegrees1 %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>% dplyr::select(ensgene,symbol)) -> kIN_df

xx <- cbind(ensgene = colnames(datExpr), module = moduleColors) %>% 
  as.tibble()


kIN_df %>% 
  left_join(xx) -> kIN_df

saveRDS(kIN_df,"manuscript/brain/results_RDS/kIN/ALL_kIN_dataframe_Power4.RDS")


### GO-Terms==============================================
set.seed(312)
# module_list = c( 'brown', 'pink', 'royalblue', 'purple', 'lightyellow', 'cyan', 'magenta', "red","salmon","lightcyan", "black", "tan")
# datExpr <-datExpr <- readRDS("manuscript/brain/results_RDS/WGCNA_datExpr_ALL.RDS") 
# net <- readRDS("manuscript/brain/results_RDS/WGCNA_net_ALL_Power5.RDS")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)

gettop10GO_WGCNA <- function(module,my_showCategory = 10){
  
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  grcm38 %>% 
    filter(ensgene %in% module_gene) %>% 
    filter(!is.na(entrez)) %>% 
    dplyr::select(entrez) -> go_df_wgcna
  
  
  ggo <- enrichGO(gene = go_df_wgcna$entrez %>% unique(),
                  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = 'BP',
                  readable = T,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.50)
  
  
  fortify(
    ggo,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    dplyr::arrange(desc(GeneRatio)) %>% 
    mutate(module = module) -> temp1
  
  return(rbind(temp1))
  
}

my_ont = "BP"
my_showCategory = 100


module_list %>% unique() -> allcolors

WGCNA_GOs <- vector('list', length(allcolors))

for(i in 1:length(allcolors)){
  gettop10GO_WGCNA(allcolors[i],my_showCategory) -> WGCNA_GOs[[i]]
}

WGCNA_GOs %>% 
  do.call(rbind,.) -> wgcna_all_gos

write.csv(wgcna_all_gos, 
          "manuscript/brain/results_tables/WGCNA_ALL_all_gos_catogeryBP_Power4.csv",
          row.names = F)



#######linearmodel 

# linear model
set.seed(312)
datExpr <- readRDS("manuscript/brain/results_RDS/WGCNA_datExpr_ALL10.RDS") 
net <- readRDS("manuscript/brain/results_RDS/WGCNA_net_ALL_Power4_10.RDS")
MEs <- readRDS("manuscript/brain/results_RDS/WGCNA_ALL_MEs_Power4.RDS")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes

coldata <- read.csv("manuscript/brain/results_tables/coldata.csv", row.names = 1)

coldata$condition1 <- factor(coldata$condition1, levels= c("CDOM","CSUB","DOM", "SUB","DES","ASC"))


MEs00 <- orderMEs(MEs0) %>% 
  rownames_to_column("SampleID")

colnames(MEs00)
#transition modules 
pink <- MEs00[,c(1,17)] %>% as.data.frame()
brown <- MEs00[,c(1,13)] %>% as.data.frame()
#ASC
royalblue <- MEs00[,c(1,18)] %>% as.data.frame()
darkgrey <- MEs00[,c(1,14)] %>% as.data.frame()
#DES/CSUB
darkgreen<- MEs00[,c(1,21)] %>% as.data.frame()
dturq<- MEs00[,c(1,25)] %>% as.data.frame()
#cdom/csub
black <- MEs00[,c(1,10)] %>% as.data.frame()
lightgreen <- MEs00[,c(1,26)] %>% as.data.frame()
salmon <- MEs00[,c(1,19)] %>% as.data.frame()
# same
blue <- MEs00[,c(1,7)] %>% as.data.frame()
green <- MEs00[,c(1,11)] %>% as.data.frame()
turq <- MEs00[,c(1,20)] %>% as.data.frame()
yellow <- MEs00[,c(1,6)] %>% as.data.frame()

grey60 <- MEs00[,c(1,5)] %>% as.data.frame()
red <- MEs00[,c(1,9)] %>% as.data.frame()
magenta <- MEs00[,c(1,15)] %>% as.data.frame()
orange <- MEs00[,c(1,12)] %>% as.data.frame()

library(lme4)
library(lmerTest)

pink%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEpink~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEpink ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEpink~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEpink ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEpink~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)



brown%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEbrown~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEbrown ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEbrown~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEbrown ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEbrown~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)

royalblue%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEroyalblue~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEroyalblue ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEroyalblue~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEroyalblue ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEroyalblue~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)



darkgrey%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEdarkgrey~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEdarkgrey ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEdarkgrey~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEdarkgrey ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEdarkgrey~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)

##################

darkgreen %>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEdarkgreen~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEdarkgreen ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEdarkgreen~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEdarkgreen~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEdarkgreen~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)

dturq%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEdarkturquoise~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEdarkturquoise ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEdarkturquoise~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEdarkturquoise ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEdarkturquoise~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)


################

black%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEblack~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEblack ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEblack~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEblack ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEblack~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)



lightgreen%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MElightgreen~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MElightgreen ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MElightgreen~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MElightgreen ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MElightgreen~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)

salmon%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEsalmon~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEsalmon ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEsalmon~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEsalmon ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEsalmon~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)

###########

blue%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEblue~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEblue ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEblue~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEblue ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEblue~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)



green%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEgreen~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEgreen ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEgreen~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEgreen ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEgreen~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)


yellow%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEyellow~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEyellow ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEyellow~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEyellow ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEyellow~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)
###########

orange%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEorange~ condition1 +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEorange ~ condition2 +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEorange~ condition3 +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEorange ~ condition4 +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEorange~ condition5 +(1|batch) , data = ME_df) -> mod5
summary(mod5)


########### statistics with cort 

# linear model
set.seed(312)
datExpr <- readRDS("manuscript/brain/results_RDS/WGCNA_datExpr_ALL10.RDS") 
net <- readRDS("manuscript/brain/results_RDS/WGCNA_net_ALL_Power4_10.RDS")
MEs <- readRDS("manuscript/brain/results_RDS/WGCNA_ALL_MEs_Power4.RDS")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes

coldata <- read.csv("manuscript/brain/results_tables/coldata.csv", row.names = 1)

MEs00 <- orderMEs(MEs0) %>% 
  rownames_to_column("SampleID")

colnames(MEs00)
mb <- MEs00[,c(1,22)] %>% as.data.frame()

mb%>%  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition2 = factor(condition1,levels= c("CSUB","CDOM","DOM", "SUB","DES","ASC"))) %>% 
  mutate(condition3 = factor(condition1,levels= c("DOM","CSUB","CDOM", "SUB","DES","ASC"))) %>% 
  mutate(condition4 = factor(condition1,levels= c("SUB","CSUB","CDOM", "DOM","DES","ASC"))) %>%
  mutate(condition5 = factor(condition1,levels= c("DES","CSUB","CDOM", "SUB","DOM","ASC"))) %>%
  relocate(condition1,condition2,condition3,condition4,condition5, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

colnames(ME_df)

lmer(MEcyan~ condition1*post_Ncort +(1|batch) , data = ME_df) -> mod1
summary(mod1)
lmer(MEcyan ~ condition2*post_Ncort +(1|batch) , data = ME_df) -> mod2
summary(mod2)
lmer(MEcyan~ condition3*post_Ncort +(1|batch) , data = ME_df) -> mod3
summary(mod3)
lmer(MEdarkgreen ~ condition4*post_Ncort +(1|batch) , data = ME_df) -> mod4
summary(mod4)
lmer(MEdarkgreen~ condition5*post_Ncort +(1|batch) , data = ME_df) -> mod5
summary(mod5)


