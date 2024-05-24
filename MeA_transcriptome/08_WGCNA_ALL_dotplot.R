library(WGCNA)
library(tidyverse)

set.seed(312)
datExpr <- readRDS("MeA_transcriptome/results_RDS/WGCNA_datExpr_ALL10.RDS") 
net <- readRDS("MeA_transcriptome/results_RDS/WGCNA_net_ALL_Power4_10.RDS")
MEs <- readRDS("MeA_transcriptome/results_RDS/WGCNA_ALL_MEs_Power4.RDS")

MEsx <- MEs[,c(1:25)]


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
p <- plot(METree, main = "Clustering of module eigengenes", xlab = "", ylab = "",  
          cex.lab=0.65,lwd = 2, rotation = 90) 
p <- plot(METree, main = "", xlab = "", ylab = "",  
          cex.lab=0.65,lwd = 2, rotation = 90) 

par(mfrow = c(1,2), mar = c(5,6,1,6))

dend <- as.dendrogram(hclust(as.dist(MEDiss), method = 'average'))
dend2 <- dendextend::color_labels(dend)
all_dendtree <- dendextend::plot_horiz.dendrogram(dend, side = F, lwd = 2)

# ggsave("manuscript/brain/results_figures/all_dendtree.png", all_dendtree, width = 20, height = 8)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes

coldata <- read_csv("MeA_transcriptome/results_tables/coldata.csv")
# 
# coldata <- coldata %>%
#   filter(condition1 != "ASC") %>%
#   filter(condition1 != "SUB") %>%
#   filter(condition1 != "CSUB")
# 
# coldata$condition1 <- factor(coldata$condition1, levels= c("DES", "DOM", "CDOM"))

 coldata <- coldata %>%
  filter(condition1 != "DES") %>%
   filter(condition1 != "DOM") %>%
   filter(condition1 != "CDOM")

 coldata$condition1 <- factor(coldata$condition1, levels= c("ASC", "SUB", "CSUB"))


orderMEs(MEs0) %>% 
  rownames_to_column("SampleID") %>%
  left_join(coldata) %>% na.omit(.) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
 # mutate(conditionx = factor(condition1, levels = c("DOM", "DES", 'CDOM'))) %>%
  mutate(conditionx = factor(condition1, levels = c("SUB", "ASC", 'CSUB'))) %>%
  relocate(condition1,conditionx, batch, post_Ncort, Postrank, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12)-> ME_df

lm_result_list <- list()

library(lme4)
library(lmerTest)

# lmer(MEbrown ~ condition1 +(1|batch) , data = ME_df) -> mod1
# lmer(MEbrown ~ conditionx +(1|batch) , data = ME_df) -> mod2
# summary(mod1)
# summary(mod2)

for(x in 1:length(MEs0)){
  k = x + 9
  ME_df[,c(1:9,k)] -> d
  md <- gsub("ME","",colnames(df)[10])
  colnames(df)[10] <- "module"
  lmer(module ~ condition1 +(1|batch) , data = df) -> mod1
  lmer(module ~ conditionx +(1|batch) , data = df) -> mod2
  summary(mod1)
  summary(mod2)
  df
  rbind(summary(lm(df[,10] ~ df[,1]))$coefficients[2,],
        summary(lm(df[,10] ~ df[,1]))$coefficients[3,],
        summary(lm(df[,10] ~ df[,2]))$coefficients[3,])%>%
    as.data.frame() %>%
      # cbind(key = c("DES-DOM","DES-CDOM","DOM-CDOM")) %>%
    cbind(key = c("ASC-SUB","ASC-CSUB","SUB-CSUB")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}

lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey") %>% filter(., module != "AggRec70min") -> lm_result_all



d <- lm_result_all

head(d)
colnames(d)

d$module <- as.factor(d$module)
d$key <- as.factor(d$key)
d$`Pr(>|t|)`<- as.numeric(d$`Pr(>|t|)`)
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

# heatmap_dfx %>% 
  # as_tibble() %>% 
  # select(module, count) %>% 
  # distinct() %>% 
  # mutate(module = factor(module)) %>% 
  # arrange(desc(count)) %>% 
  # filter(module!="grey") %>% 
  # filter(module!="AggRec70min") %>%
  # .$module -> my_module_level



str(heatmap_df)

heatmap_df$Estimate <- as.numeric(heatmap_df$Estimate)
heatmap_df$`Std. Error`<- as.numeric(heatmap_df$`Std. Error`)
heatmap_df$module <- as.factor(heatmap_df$module)

color_above <- heatmap_df$`Pr(>|t|)` < 0.05
color_below <- heatmap_df$`Pr(>|t|)` > 0.05
# # 
# heatmap_df$color <- ifelse(heatmap_df$`Pr(>|t|)`< color_above,heatmap_df$`Pr(>|t|)`,.055)
# heatmap_df$key <- factor(heatmap_df$key, levels = c("DES-DOM", "DES-CDOM", "DOM-CDOM"))
# heatmap_df$key2 <- ifelse(heatmap_df$key == "DES-DOM", "DES vs. DOM", heatmap_df$key)
# heatmap_df$key2 <- ifelse(heatmap_df$key == "DES-CDOM", "DES vs. CDOM", heatmap_df$key2)
# heatmap_df$key2 <- ifelse(heatmap_df$key == "DOM-CDOM", "DOM vs. CDOM", heatmap_df$key2)
# heatmap_df$key2 <- factor(heatmap_df$key2, levels = c("DES vs. DOM", "DES vs. CDOM", "DOM vs. CDOM"))

heatmap_df$color <- ifelse(heatmap_df$`Pr(>|t|)`< color_above,heatmap_df$`Pr(>|t|)`,.055)
 heatmap_df$key <- factor(heatmap_df$key, levels = c("ASC-SUB", "ASC-CSUB", "SUB-CSUB"))
 heatmap_df$key2 <- ifelse(heatmap_df$key == "ASC-SUB", "ASC vs. SUB", heatmap_df$key)
 heatmap_df$key2 <- ifelse(heatmap_df$key == "ASC-CSUB", "ASC vs. CSUB", heatmap_df$key2)
heatmap_df$key2 <- ifelse(heatmap_df$key == "SUB-CSUB", "SUB vs. CSUB", heatmap_df$key2)
heatmap_df$key2 <- factor(heatmap_df$key2, levels = c("ASC vs. SUB", "ASC vs. CSUB", "SUB vs. CSUB"))
heatmap_df$module <- factor(heatmap_df$module, levels = c('lightgreen', 'darkturquoise', 'tan', 'lightyellow', 'cyan', 'darkgreen',
                                                           'turquoise', 'salmon', 'royalblue', 'pink', 'darkred', 'magenta', 'darkgrey', 'brown', 
                                                          'orange','green','black', 'red', 'greenyellow', 'blue', 'yellow', 'grey60', 
                                                          'purple', 'midnightblue', 'lightcyan'))



dom <- heatmap_df %>% 
  ggplot(aes(y = module, x = Estimate, color =color)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3)+
  facet_wrap(~key2)+
  labs(y="", color = "p-value")+
  xlim(-2,2)+
  scale_color_continuous(low = "red", high = "lightgray", breaks = c(0.01, 0.02,0.03,0.04, 0.05))+
  geom_errorbar(aes(xmin = Estimate-`Std. Error`, xmax = Estimate+`Std. Error`),width = 0.2)+
  theme_bw()+
  scale_y_discrete(limits = rev(levels(heatmap_df$module)))+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        text = element_text(size = 15))
 dom


