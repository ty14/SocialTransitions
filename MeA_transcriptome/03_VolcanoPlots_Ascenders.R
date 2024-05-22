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


my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results_RDS/limma_MEA_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~unique(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(.,symbol, entrez))) %>%
  map(~unique(.,)) %>% 
  map(~filter(.,!is.na(entrez)))  


y1a <- limma_list$csub
y2a <- limma_list$casc
y3a <- limma_list$asub

###volano plots
head(y1a)
dc <- y1a %>% mutate(contrast = "SUB vs. CSUB") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 1.7)%>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1.6) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_csub <- ggplot(data = dc, 
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
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black",  hjust = -1,vjust =.8)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 1.5)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))+
  theme_bw() +
  annotate(geom="text", x=2.5, y=.5, label=paste0(" ", "\n", "in SUB"),
           color="black", size = 5)+
  annotate(geom="text", x=-2.5, y=.5, label=paste0(" ", "\n", "in CSUB"),
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

vp_csub


 ggsave("manuscript/brain/results_figures/vplot_csub.png",vp_csub,height =5.5, width =6.5, dpi=600)

### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #309
dc %>% filter(between(logFC, .2, .75))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 276
dc %>% filter(between(logFC, .75, 1.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #36 
dc %>% filter(between(logFC, 1.5, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #6

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 322 
dc %>% filter(between(logFC, -0.75, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #248
dc %>% filter(between(logFC, -1.5, -0.75)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 62
dc %>% filter(between(logFC, -3, -1.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 12

#top 10 genes up in sub
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes down in sub
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)
################
head(y2a)

dc <- y2a %>% mutate(contrast = "ASC vs. CSUB") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 1.5) %>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1.5) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_casc <- ggplot(data = dc, 
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
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", vjust = 1.1, hjust = 1.1)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = .25, vjust = -.7)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))  +
  theme_bw() +
  annotate(geom="text", x=2.5, y=.5, label=paste0(" ", "\n", "in ASC"),
           color="black", size = 5)+
  annotate(geom="text", x=-2.5, y=.5, label=paste0(" ", "\n", "in CSUB"),
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

vp_casc

ggsave("manuscript/brain/results_figures/vplot_casc.png",vp_casc,height =5.5, width =6.5, dpi=600)

### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #280
dc %>% filter(between(logFC, .2, .75))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 250
dc %>% filter(between(logFC, .75, 1.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #26
dc %>% filter(between(logFC, 1.5, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #4

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 225
dc %>% filter(between(logFC, -0.75, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #185
dc %>% filter(between(logFC, -1.5, -0.75)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 37
dc %>% filter(between(logFC, -3, -1.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #3


#top 10 genes up in asc
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes down in asc
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)
################


#################
head(y3a)

dc <- y3a %>% mutate(contrast = "ASC vs. CSUB") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 1.55) %>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1.7) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_as <- ggplot(data = dc, 
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
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = -.3,vjust =.7)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust =1.5, vjust =.3)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR))) +
  theme_bw() +
  annotate(geom="text", x=2.5, y=.5, label=paste0(" ", "\n", "in ASC"),
           color="black", size = 5)+
  annotate(geom="text", x=-2.5, y=.5, label=paste0(" ", "\n", "in SUB"),
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


vp_as

ggsave("manuscript/brain/results_figures/vplot_as.png",vp_as,height =5.5, width =6.5, dpi=600)

### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #658
dc %>% filter(between(logFC, .2, .75))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 576
dc %>% filter(between(logFC, .75, 1.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #77
dc %>% filter(between(logFC, 1.5, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #5

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 727
dc %>% filter(between(logFC, -0.75, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #656
dc %>% filter(between(logFC, -1.5, -0.75)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 61
dc %>% filter(between(logFC, -3, -1.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 10

#top 10 genes up in dom
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes down in dom 
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)
