#social transitions with out controls. 

# libraries 
library(limma)
library(edgeR)
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

my_logFC_threshold = 0.2

#descenders
limma_list<- readRDS("manuscript/brain/results_RDS/limma_MEA_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dd <- limma_list$desdom

#Ascenders
limma_list<- readRDS("manuscript/brain/results_RDS/limma_MEA_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

as <- limma_list$asub


#upreg
dd_up <- dd %>% filter(logFC >= 0.2)
as_up <- as %>% filter(logFC >= 0.2)

#downreg
dd_down <- dd %>% filter(logFC <= -0.2)
as_down <- as %>% filter(logFC <= -0.2)


# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(dd_up= dd_up$symbol, as_up=as_up$symbol, dd_down = dd_down$symbol, 
                  as_down = as_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#significant test
mat <- matrix(c(142,8,6,123),ncol=2)

# Perform chi-squared test
chi_square <- chisq.test(mat)

# Calculate phi-coefficient manually
n <- sum(mat)
phi <- sqrt(chi_square$statistic / n)

# phi-coefficient and p-value
phi #0.8920702 
chi_square$p.value #3.270708e-50


# get transition genes, logFC, and pvalues 
trans_up <- dd_up$symbol[dd_up$symbol %in% as_up$symbol] %>% as.data.frame() %>% unique(.)
trans_down <- dd_down$symbol[dd_down$symbol %in% as_down$symbol]%>% as.data.frame() %>% unique(.)

td_up<- dd[dd$symbol %in% trans_up$.,] %>% dplyr::select(symbol, DES_logFC = logFC, DES_pvalue = P.Value)
ta_up <- as[as$symbol %in% trans_up$., ] %>% dplyr::select(symbol, ASC_logFC = logFC, ASC_pvalue = P.Value)
trans_up <- td_up %>% full_join(ta_up) %>% mutate(contrast = "trans vs. same")


td_down<- dd[dd$symbol %in% trans_down$.,] %>% dplyr::select(symbol, DES_logFC = logFC, DES_pvalue = P.Value)
ta_down <- as[as$symbol %in% trans_down$., ] %>% dplyr::select(symbol, ASC_logFC = logFC, ASC_pvalue = P.Value)
trans_down <- td_down %>% full_join(ta_down) %>% mutate(contrast = "trans vs. same")

trans_all <- trans_up %>% rbind(trans_down)
#save in tables 
write.csv(trans_all,"manuscript/brain/results_tables/TransitionGenes_withoutControls.csv", row.names = F)


# genes that go in opposite directions. 
dd_down$symbol[dd_down$symbol %in% as_up$symbol]
# [1] "Bmpr1b"  "Irf5"    "Tmem123" "Tgfa"    "Cog2"    "Gjd4"    "Carmil1" "Agap1"  
dd_up$symbol[dd_up$symbol %in% as_down$symbol]
# [1] "Matn2"  "Fads3"  "Eif3j2" "Kcnk12" "Lins1"  "Zfp850"

#####go terms 
source("functions/gettop10GO.R")
trans_up <- dd_up$symbol[dd_up$symbol %in% as_up$symbol] %>% as.data.frame()
trans_down <- dd_down$symbol[dd_down$symbol %in% as_down$symbol]%>% as.data.frame()
#descenders: 
tdd <- dd[dd$symbol %in% trans_up$.,] 
tddx <- dd[dd$symbol %in% trans_down$.,] 
tdd_go <- tdd %>% rbind(tddx)

gettop10GO(tdd_go, my_showCategory) %>% 
  mutate(comparison = "DES vs. DOM") -> top10go1

#descenders: 
tta <- as[as$symbol %in% trans_up$.,] 
ttax <- as[as$symbol %in% trans_$.down,] 
tas_go <- tta %>% rbind(ttax)

gettop10GO(tdd_go, my_showCategory) %>% 
  mutate(comparison = "ASC vs. SUB") -> top10go2


rbind(top10go1,top10go2) -> top10_GOterms

write.csv(top10_GOterms,"manuscript/brain/results_tables/topBP_GOterms_Trans_withoutControls.csv", row.names = F)

# view
top10_GOterms$Description[1:5]
top10_GOterms$Description[6:10]
top10_GOterms$Description[11:15]
top10_GOterms$Description[16:20]
top10_GOterms$Description[21:25]



limma_list<- readRDS("manuscript/brain/results_RDS/limma_MEA_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~unique(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(.,symbol, entrez))) %>%
  map(~unique(.,)) %>% 
  map(~filter(.,!is.na(entrez)))  


d <- limma_list$desdom


limma_list<- readRDS("manuscript/brain/results_RDS/limma_MEA_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~unique(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(.,symbol, entrez))) %>%
  map(~unique(.,)) %>% 
  map(~filter(.,!is.na(entrez)))  


a <- limma_list$asub





dx <- d %>% select(symbol, DES_logFC = logFC, DES_pv = P.Value)
ax <- a %>% select(symbol, ASC_logFC = logFC, ASC_pv = P.Value)

xx <- dx %>% full_join(ax)
ts <- trans_all$symbol


xx$condition <- ifelse(xx$symbol %in% ts, "TRN", "N.S.")
xx$condition <- ifelse(xx$ASC_logFC < 0.02 & xx$DES_logFC > 0.02, "N.T", xx$condition)
xx$condition <- ifelse(xx$ASC_logFC > 0.02 & xx$DES_logFC < 0.02, "N.T", xx$condition)
table(xx$condition)


xx0 <- xx %>% filter(condition=="N.S.")
x1 <- xx %>% filter(condition=="N.T") %>% filter(ASC_pv < 0.05 & DES_pv < 0.05)
xx2 <- xx %>% filter(condition=="TRN")

x <- xx2 %>% arrange(-DES_logFC) %>% head(.,3)
y <- xx2 %>% arrange(DES_logFC) %>% head(.,3)
for_label <- x %>% rbind(y)
tg <- ggplot()+
  # geom_point(data=xx, aes(ASC_logFC, DES_logFC)) +
  geom_point(data=xx0, aes(ASC_logFC, DES_logFC), color='lightgrey', alpha = 0.1,shape = 21, size=3) +
  geom_point(data=xx2, aes(ASC_logFC, DES_logFC), color='blue',fill='blue', alpha = 0.1,shape = 21, size=4)+
  geom_point(data=x1, aes(ASC_logFC, DES_logFC), color='red',fill='red', alpha = 0.1,shape = 21, size=4)+
  ylim(-3, 3)+
  xlim(-3,3)+
  geom_vline(xintercept=c(0,0), linetype="dotted")+
  geom_hline(yintercept=c(0,0), linetype="dotted")+
  annotate("text",label = paste0("Upregulated", '\n', "DES + ASC"), x = 2.5, y = 2.5, size = 5)+
  annotate("text",label = paste0("Downregulated", '\n', "DES + ASC"), x = -2, y = -2.5, size = 5)+
  xlab("ASC LogFC")+
  ylab("DES LogFC")+
  theme_classic()+ theme(text=element_text(size = 10))
tg+ geom_text_repel(aes(for_label$symbol))

# ggsave(filename = "manuscript/brain/results_figures/transition_scatterplot2.png",
      # tg,
       # height = 5, width = 5, dpi = 600)



tg2 <- ggplot(xx, aes(ASC_logFC, DES_logFC, group = condition))+
  geom_point(data=xx0, aes(ASC_logFC, DES_logFC), color='lightgrey', alpha = 0.1,shape = 21, size=3) +
  geom_point(data=xx2, aes(ASC_logFC, DES_logFC), color='blue',fill='blue', alpha = 0.1,shape = 21, size=4)+
  geom_point(data=x1, aes(ASC_logFC, DES_logFC), color='red',fill='red', alpha = 0.1,shape = 21, size=4)+
  geom_text_repel(data= for_label, aes(label = symbol), size =3)+
  ylim(-3, 3)+
  xlim(-3,3)+
  geom_vline(xintercept=c(0,0), linetype="dotted")+
  geom_hline(yintercept=c(0,0), linetype="dotted")+
  # annotate("text",label = paste0("Upregulated", '\n', "DES + ASC"), x = 2.2, y = 2.3, size = 3.5)+
  # annotate("text",label = paste0("Downregulated", '\n', "DES + ASC"), x = -2.2, y = -2.3, size = 3.5)+
  xlab("ASC LogFC")+
  ylab("DES LogFC")+
  theme_classic()+ theme(text=element_text(size = 10))
tg2
ggsave(filename = "manuscript/brain/results_figures/transition_scatterplot3.png",
tg2,
height = 4, width = 4, dpi = 600)
