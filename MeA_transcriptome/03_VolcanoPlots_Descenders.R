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
library(EnhancedVolcano)
library(tidyverse)
grcm38 # mouse genes


my_logFC_threshold = .2

limma_list<- readRDS("MeA_transcriptome/results_RDS/limma_MEA_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~unique(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(.,symbol, entrez))) %>%
  map(~unique(.,)) %>% 
  map(~filter(.,!is.na(entrez)))  


y1a <- limma_list$cdom
y2a <- limma_list$cdes
y3a <- limma_list$desdom

###volano plots
head(y1a)
dc <- y1a %>% mutate(contrast = "DOM vs. CDOM") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 1.75)%>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1.6) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_cdom <- ggplot(data = dc, 
       aes(x = logFC, 
           y = log10, 
           colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange", "grey","purple4"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-2.5,2.5),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-.75,.75),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black",vjust =1, hjust =.45)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 1, vjust = -.55)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))+
  theme_bw() +
  annotate(geom="text", x=2.5, y=.5, label=paste0(" ", "\n", "in DOM"),
           color="black", size = 5)+
  annotate(geom="text", x=-2.5, y=.5, label=paste0(" ", "\n", "in CDOM"),
           color="black",, size = 5)+
  scale_x_continuous(limits = c(-3.5,3.5),breaks = c(-3,-2,-1,0,1,2,3))+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 20)
        
  )

vp_cdom




### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #802 
dc %>% filter(between(logFC, .2, .75))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 676
dc %>% filter(between(logFC, .75, 1.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #107 
dc %>% filter(between(logFC, 1.5, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #19 

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 832 
dc %>% filter(between(logFC, -0.75, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #719
dc %>% filter(between(logFC, -1.5, -0.75)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 104
dc %>% filter(between(logFC, -3, -1.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 9 

#top 10 genes up in dom
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes down in dom 
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)
################
head(y2a)

dc <- y2a %>% mutate(contrast = "DES vs. CDOM") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 1.75) %>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1.5) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_cdes <- ggplot(data = dc, 
       aes(x = logFC, 
           y = log10, 
           colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange", "grey","purple4"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-.75,.75),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black",hjust =1, vjust =.9)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", vjust = -1)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))  +
  theme_bw() +
  annotate(geom="text", x=2.5, y=.5, label=paste0(" ", "\n", "in DES"),
           color="black", size = 5)+
  annotate(geom="text", x=-2.5, y=.5, label=paste0(" ", "\n", "in CDOM"),
           color="black",, size = 5)+
  scale_x_continuous(limits = c(-3.5,3.5),breaks = c(-3,-2,-1,0,1,2,3))+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 20)
        
  )

vp_cdes


### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #575
dc %>% filter(between(logFC, .2, .75))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 463
dc %>% filter(between(logFC, .75, 1.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #100
dc %>% filter(between(logFC, 1.5, 4))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #11 

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 526
dc %>% filter(between(logFC, -0.75, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #472
dc %>% filter(between(logFC, -1.5, -0.75)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 51
dc %>% filter(between(logFC, -3, -1.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 3


#top 10 genes up in des
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes down in des 
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)
################


#################
head(y3a)

dc <- y3a %>% mutate(contrast = "DES vs. DOM") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 2) %>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -2) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_dd <- ggplot(data = dc, 
       aes(x = logFC, 
           y = log10, 
           colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange", "grey","purple4"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-.75,.75),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = .7)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust =.1, vjust = .4)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR))) +
  theme_bw() +
  annotate(geom="text", x=2.5, y=.5, label=paste0(" ", "\n", "in DES"),
           color="black", size = 5)+
  annotate(geom="text", x=-2.5, y=.5, label=paste0(" ", "\n", "in DOM"),
           color="black",, size = 5)+
  scale_x_continuous(limits = c(-3.5,3.5),breaks = c(-3,-2,-1,0,1,2,3))+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 20)
        
  )


vp_dd



### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #1187
dc %>% filter(between(logFC, .2, .75))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 990
dc %>% filter(between(logFC, .75, 1.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #176
dc %>% filter(between(logFC, 1.5, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #21 

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 1130
dc %>% filter(between(logFC, -0.75, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #888
dc %>% filter(between(logFC, -1.5, -0.75)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 212
dc %>% filter(between(logFC, -4, -1.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 28

#top 10 genes up in dom
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes down in dom 
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)
