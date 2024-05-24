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

y1a <- limma_list$desdom

y2a <- limma_list$cdes

y3a <- limma_list$cdom


dd <- y1a %>%  as_tibble() %>% .$symbol

cdes <-  y2a %>% as_tibble() %>% .$symbol

cdom <- y3a %>%   as_tibble() %>% .$symbol

d_a <- c(dd, cdes,cdom) 
d_ax <- unique(d_a)

source("functions/gettop10GO.R")

gettop10GO(y1a, my_showCategory) %>% 
  mutate(comparison = "DES -DOM") -> top10go1

gettop10GO(y2a, my_showCategory ) %>% 
  mutate(comparison = "DES-CDOM") -> top10go2

gettop10GO(y3a, my_showCategory ) %>% 
  mutate(comparison = "DOM-CDOM") -> top10go3


rbind(top10go1,top10go2,top10go3) -> top10_GOterms


#Ascenders 

limma_list<- readRDS("manuscript/brain/results_RDS/limma_MEA_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

y1a <- limma_list$csub

y2a <- limma_list$casc

y3a <- limma_list$asub

cs <- y1a %>%  as_tibble() %>% .$symbol

ca <-  y2a %>% as_tibble() %>% .$symbol

as <- y3a %>%   as_tibble() %>% .$symbol

a_a <- c(cs, ca,as) 
a_ax <- unique(a_a)

all <- c(a_a, d_a)
allx <- unique(all)
gettop10GO(y1a, my_showCategory) %>% 
  mutate(comparison = "SUB-CSUB") -> top10go1

gettop10GO(y2a, my_showCategory ) %>% 
  mutate(comparison = "ASC-CSUB") -> top10go2

gettop10GO(y3a, my_showCategory ) %>% 
  mutate(comparison = "ASC-SUB") -> top10go3


rbind(top10go1,top10go2,top10go3) -> top10_GOterms


