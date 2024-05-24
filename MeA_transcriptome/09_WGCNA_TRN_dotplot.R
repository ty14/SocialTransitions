library(WGCNA)
library(tidyverse)

set.seed(312)
datExpr <- readRDS("MeA_transcriptome/results_RDS/WGCNA_datExpr_trans.RDS") 
net <- readRDS("MeA_transcriptome/results_RDS/WGCNA_net_trans_Power4_cutoff100.RDS")
MEs <- readRDS("MeA_transcriptome/results_RDS/WGCNA_trans_MEs_Power4.RDS")

MEs 

MEsx <- MEs[,c(1:20)]

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
rownames(datExpr)
colnames(MEsx)<- gsub("ME", "", colnames(MEsx))


MEDiss = 1-cor(MEsx)
METree = hclust(as.dist(MEDiss), method = "average")

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "")

sizeGrWindow(7, 6)
p <- plot(METree, main = "Clustering of module eigengenes", xlab = "", ylab = "",  
          cex.lab=0.65,lwd = 2, rotation = 90) 
p <- plot(METree, main = "", xlab = "", ylab = "",  
          cex.lab=0.65,lwd = 2, rotation = 90) 

par(mfrow = c(1,2), mar = c(5,6,1,6))

dend <- as.dendrogram(hclust(as.dist(MEDiss), method = 'average'))
dend2 <- dendextend::color_labels(dend)
all_dendtree <- dendextend::plot_horiz.dendrogram(dend, side = F, lwd = 2)



# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes

# get linear model data 
coldata <- read_csv("MeA_transcriptome/sample70min_table.csv")
head(coldata)
str(coldata)

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
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "TRANS", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "TRANS", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)

coldata <- coldata %>%
  filter(condition1 != "CSUB") %>%
  filter(condition1 != "CDOM") %>% 
  filter(Prerank != 3)  %>% filter(condition != 'ascenders') 
#just get samples I want
coldata$SampleID <- substr(coldata$SampleName, 7, 13)

coldata$condition1 <- gsub("TRANS", "TRN", coldata$condition1)
coldata$condition1 <- factor(coldata$condition1, levels= c("TRN", "DOM", "SUB"))
coldata$SampleID


orderMEs(MEs0) %>% 
  rownames_to_column("SampleID") %>%
  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(conditionx = factor(condition1, levels = c("DOM", "TRN", 'SUB'))) %>% 
  relocate(condition1,conditionx, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region)-> ME_df

lm_result_list <- list()

library(lme4)
library(lmerTest)

for(x in 1:length(MEs0)){
  k = x + 8
  ME_df[,c(1:8,k)] -> df
  md <- gsub("ME","",colnames(df)[9])
  colnames(df)[9] <- "module"
  lmer(module ~ condition1 +(1|batch) , data = df) -> mod1
  lmer(module ~ conditionx +(1|batch) , data = df) -> mod2
  summary(mod1)
  summary(mod2)
  df
  rbind(summary(lm(df[,9] ~ df[,1]))$coefficients[2,],
        summary(lm(df[,9] ~ df[,1]))$coefficients[3,],
        summary(lm(df[,9] ~ df[,2]))$coefficients[3,])%>%
    as.data.frame() %>%
    cbind(key = c("TRN-DOM","TRN-SUB","DOM-SUB")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}

lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey") %>% format(., scientific = F) -> d





head(d)
colnames(d)
d$`Pr(>|t|)` <- as.numeric(d$`Pr(>|t|)`)
d %>% 
  ggplot(aes(y = module, x =`Pr(>|t|)`)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3.5)+
  facet_wrap(~key)




moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)
heatmap_df <- d %>% 
  mutate(my_alpha = ifelse(`Pr(>|t|)` < 0.05, 1, 0)) %>% 
  mutate(my_alpha2 = Estimate )

moduleColors %>% 
  table() %>% 
  as.data.frame() %>% arrange(Freq)  -> modnum 
colnames(modnum) <- c("module","count")

heatmap_df %>% as_tibble() %>% 
  left_join(modnum) -> heatmap_dfx



heatmap_df %>% 
  filter(my_alpha > 0.95) %>% 
  .$Estimate %>% as.numeric(.) %>% 
  abs() %>% 
  max() -> my_limit
rev(levels(as.factor(heatmap_df$module))) -> xx
xx[xx != "grey"] -> y_limit
xx[xx != "AggRec70min"] -> y_limit


str(heatmap_df)

heatmap_df$Estimate <- as.numeric(heatmap_df$Estimate)
heatmap_df$`Std. Error` <- as.numeric(heatmap_df$`Std. Error`)
heatmap_df$module <- as.factor(heatmap_df$module)



color_above <- heatmap_df$`Pr(>|t|)` < 0.05
color_below <- heatmap_df$`Pr(>|t|)` > 0.05


heatmap_df$color <- ifelse(heatmap_df$`Pr(>|t|)`< color_above,heatmap_df$`Pr(>|t|)`,.055)
heatmap_df$key <- factor(heatmap_df$key, levels = c("TRN-DOM", "TRN-SUB", "DOM-SUB"))
heatmap_df$key2 <- ifelse(heatmap_df$key == "TRN-DOM", "TRN vs. DOM", heatmap_df$key)
heatmap_df$key2 <- ifelse(heatmap_df$key == "TRN-SUB", "TRN vs. SUB", heatmap_df$key2)
heatmap_df$key2 <- ifelse(heatmap_df$key == "DOM-SUB", "DOM vs. SUB", heatmap_df$key2)
heatmap_df$key2 <- factor(heatmap_df$key2, levels = c("TRN vs. DOM", "TRN vs. SUB", "DOM vs. SUB"))
 heatmap_df$module <- factor(heatmap_df$module, levels = c("yellow", "red", "lightgreen",
"royalblue", "turquoise", "greenyellow","pink","grey60", "green","tan","blue",
"purple","salmon","lightcyan","brown","black","cyan", "midnightblue", "lightyellow", "magenta"))



str(heatmap_df)
trn_dot <- heatmap_df %>% 
  ggplot(aes(y = module, x = Estimate, color =color)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3)+
  facet_wrap(~key2)+
  labs(y="", color = "p-value")+
  xlim(-2,2)+
  scale_color_continuous(low = "red", high = "lightgray", breaks = c(0.01, 0.02,0.03,0.04, 0.05))+
  geom_errorbar(aes(xmin = Estimate-`Std. Error`, xmax = Estimate+`Std. Error`),width = 0.2)+
  scale_y_discrete(limits = rev(levels(heatmap_df$module)))+
  theme_bw()+
  theme(axis.ticks = element_blank(),axis.text.y = element_blank(),text = element_text(size = 15))

trn_dot 

