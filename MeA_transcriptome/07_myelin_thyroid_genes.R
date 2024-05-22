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

limma_list<- readRDS("manuscript/brain/results_RDS/limma_MEA_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dd <- limma_list$desdom

cdes <- limma_list$cdes

cdom <- limma_list$cdom

#stable data 
limma_list<- readRDS("brain/results/Won_MeA_data/limma_MEA.RDS")%>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 


as <- limma_list$alphasub

ab<- limma_list$alphasubdom

bs <- limma_list$subdomsub



my <- c("Mog","Mal","Mobp","Mbp","Tspan2", "Nkx6-2","Cntn2","Cnp","Arhgef10", "Lpar1")

my <- c('Cnp', 'Mobp', 'Mbp', 'Myrf', 'Mal', 'Bcas1', 'Mog', 'Mag', 'Lpar1',  'Plp1', 'Tspan2', 'Cntn2', 'Opalin', 'Arhgef10')

thy <- c("Dio2", "Thra", "Thrb", "Trh", "Trhr")

dd_down <- dd %>% filter(logFC <=-.2) %>% arrange(logFC) %>%  head(., 50)

dd_down%>% filter(symbol %in% my)


dd %>% filter(symbol %in% my)
cdom %>% filter(symbol %in% my)
cdes %>% filter(symbol %in% my)



dd %>% filter(symbol %in% ar)
cdom %>% filter(symbol %in% ar)
cdes %>% filter(symbol %in% ar)


as %>% filter(symbol %in% my)
ab %>% filter(symbol %in% my)
bs %>% filter(symbol %in% my)



ex <- readRDS("brain/results/Won_MeA_data/limma_vdl_MeA.RDS")
head(ex)
x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 


ht2 <- x %>% filter(symbol %in% c("Mog", "Tspan2", "Cnp"))
id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

colnames(ht2)
ht2 <- ht2 %>% pivot_longer(cols = 2:34, names_to = "ids")

p <- ht2 %>% full_join(id)

# p$group <- factor(p$group, levels = c("Alpha", "Subdominant", "Subdorinate"))

source('functions/geom_boxjitter.R')
library(viridis)

p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~symbol,scales="free_y", ncol =7)+
  ylab("Normalized Expression") +
  theme_bw()+
  theme(legend.position = "none")
p1


ex <- readRDS("manuscript/brain/results_RDS/limma_vdl_MeA_CDOM.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)


x %>% 
  filter(symbol %in% my) -> xex



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
  facet_wrap(~symbol, scales = 'free', ncol =7)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 15))
p1
ggsave("manuscript/brain/results_figures/DES_myelin2_boxplots.png", width =18 , height = 5, dpi = 600)
