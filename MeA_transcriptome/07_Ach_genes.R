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

#Descenders 
my_logFC_threshold = 0.2

limma_list<- readRDS("MeA_transcriptome/results_RDS/limma_MEA_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dd <- limma_list$desdom
cdom <- limma_list$cdom
cdes <- limma_list$cdes

cdom[cdom$symbol %in%chol,]

limma_list<- readRDS("MeA_transcriptome/results_RDS/limma_MEA_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

as <- limma_list$asub
cs <- limma_list$asub
ca <- limma_list$asub

ca[ca$symbol %in%chol,]



#chol genes: 
chol <- c('Chrm2', 'Chrm4','Gnai1','Gnai2','Chrm1','Chrm3','Chrm5','Gna11', 'Gna14','Gnaq','Ache',
          'Chat','Grk2', 'Grk5', 'Rgs2', 'Rgs4', 'Rgs6', 'Slc18a3', 'Slc5a7', 'Nat1', 'Lhx8', 'Slc10a4', 
          'Gbx1','Chrna2', 'Chrna3', 'Chrna6', 'Chrna7', 'Chrnb4', 'Chrnb3', 'Agrn', 'Chrna1', 'Chrna10', 
          'Chrna4', 'Chrna5', 'Chrna9', 'Chrnb1', 'Chrnb2', 'Chrnd', 'Chrne', 'Chrng', 'Dok7', 'Lrp4', 'Musk',
          'Rapsn')
x<- dd  %>%   filter(symbol %in% chol) %>% arrange(logFC)
y <- cdom %>% filter(symbol %in% chol) %>% arrange(logFC)
x$symbol[x$symbol %in% y$symbol]

cdes %>% filter(symbol %in% chol) %>% arrange(logFC)

dd  %>% arrange(logFC) %>%  head(.,20) %>% 
  filter(symbol %in% chol) %>% arrange(logFC)

as %>%  filter(symbol %in% chol) %>% arrange(logFC)

chol[chol %in% rsx]

col2 <- c("Lhx8", "Slc18a3", "Slc10a4", "Gbx1", "Chat", "Chrm2", "Slc5a7", "Ache")

dd  %>% filter(symbol %in% col2) %>% arrange(logFC)
as %>%  filter(symbol %in% col2) %>% arrange(logFC)

ex <- readRDS("MeA_transcriptome/results_RDS/limma_vdl_MeA_CDOM.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)


x %>% 
  filter(symbol %in% col2) -> xex



xex2 <- xex %>% pivot_longer(cols = 2:21, names_to = "ids")

p <- xex2 %>% full_join(id)

p$group <- factor(p$group, levels = c("DES", "DOM", "CDOM"))

source('functions/geom_boxjitter.R')
library(viridis)

p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+
  facet_wrap(~symbol, scales = 'free', ncol =8)+
  scale_y_continuous(expand = c(0, 2))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 15))
p1
