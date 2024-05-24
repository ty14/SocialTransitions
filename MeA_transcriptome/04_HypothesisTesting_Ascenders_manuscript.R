library(EnhancedVolcano)
library(tidyverse)

my_logFC_threshold = 0.2

limma_list<- readRDS("MeA_transcriptome/results_RDS/limma_MEA_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~unique(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(.,symbol, entrez))) %>%
  map(~unique(.,)) %>% 
  map(~filter(.,!is.na(entrez)))  

y1a <- limma_list$ca
y2a <- limma_list$cs
y3a <- limma_list$asub


 # ARGs from Tyssowski 2018 - (c)
pr <- read_csv("MeA_transcriptome/gene_sets/PrimaryResponseGenes_fromTyssowski.csv")
pr %>%
  as_tibble() %>%
  .$`gene name` -> prx

prx <- unique(prx)
prx_names <- c("Primary Response Genes")


###1
head(y3a)

y3a  %>% 
  filter(symbol %in% prx) %>%
  dplyr::select(symbol, logFC, P.Value) %>% 
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"ASC genes","SUB genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","SUB genes", "N.S."))) %>% 
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
  unique() -> df

keyvals <- ifelse(
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
         'grey'))

table(keyvals) %>% 
  as.data.frame() %>% 
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>% 
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts

if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}

keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'


df %>% 
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  top_n(20) %>% 
  # top_frac(.,1) %>%
  # filter(P.Value <0.01 | abs(logFC) >0.5) %>%
  .$symbol -> for_label

EnhancedVolcano(df,
                lab = df$symbol,
                selectLab= for_label,
                x = 'logFC',
                y = 'P.Value',
                title = "",
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                # drawConnectors = TRUE, # just for a few plots
                # widthConnectors = 0.1, # just for a few plots
                vline = c(-0.15, 0.15),
                vlineCol = c('grey90'),
                vlineType = c( 'dashed'),
                vlineWidth = c(0.3),
                hline = c(0.05),
                hlineCol = c('grey90'),
                hlineType = c( 'dashed'),
                hlineWidth = c(0.5),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2.5,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x = 1.25, y = 4.4,color = "purple4" , size = 3,
           label =glue::glue("Relatively Upregulated in ASC", "\n", "{my_texts[3,4]}"),)+
  annotate("text", x = -1.25, y = 4.4,color = "orange" , size = 3,
           label = glue::glue("Relatively Upregulated in SUB","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-2,2),breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 7),
        plot.subtitle = element_blank()) -> a
print(a)


# 2. 
y3a  %>% 
  filter(symbol %in% prx) %>%
  dplyr::select(symbol, logFC, P.Value) %>% 
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"ASC genes","SUB genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","SUB genes", "N.S."))) %>% 
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
  unique() -> df

df_a <- df %>%  filter(Sig == "ASC genes")

keyvals <- ifelse(
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
         'grey'))

table(keyvals) %>% 
  as.data.frame() %>% 
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>% 
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts

if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}

keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'


df %>% 
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  top_n(20) %>% 
  # top_frac(.,1) %>%
  # filter(P.Value <0.01 | abs(logFC) >0.5) %>%
  .$symbol -> for_label

EnhancedVolcano(df,
                lab = df$symbol,
                selectLab= for_label,
                x = 'logFC',
                y = 'P.Value',
                title = "",
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                # drawConnectors = TRUE, # just for a few plots
                # widthConnectors = 0.1, # just for a few plots
                vline = c(-0.15, 0.15),
                vlineCol = c('grey90'),
                vlineType = c( 'dashed'),
                vlineWidth = c(0.3),
                hline = c(0.05),
                hlineCol = c('grey90'),
                hlineType = c( 'dashed'),
                hlineWidth = c(0.5),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2.5,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x = 2, y = 4.4,color = "purple4" , size = 3,
           label =glue::glue("Relatively Upregulated in ASC", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -2, y = 4.4,color = "orange" , size = 3,
           label = glue::glue("Relatively Upregulated in SUB","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-3,3),breaks = c(-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3.0))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 7),
        plot.subtitle = element_blank()) -> b
print(b)


# c
y3a  %>% 
  filter(symbol %in% prx) %>%
  dplyr::select(symbol, logFC, P.Value) %>% 
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"ASC genes","SUB genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","SUB genes", "N.S."))) %>% 
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
  unique() -> df

apr <- df %>% filter(Sig == "SUB genes")

keyvals <- ifelse(
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
         'grey'))

table(keyvals) %>% 
  as.data.frame() %>% 
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>% 
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts

if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}

keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'


df %>% 
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  top_n(20) %>% 
  # top_frac(.,1) %>%
  # filter(P.Value <0.01 | abs(logFC) >0.5) %>%
  .$symbol -> for_label

EnhancedVolcano(df,
                lab = df$symbol,
                selectLab= for_label,
                x = 'logFC',
                y = 'P.Value',
                title = "",
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                # drawConnectors = TRUE, # just for a few plots
                # widthConnectors = 0.1, # just for a few plots
                vline = c(-0.15, 0.15),
                vlineCol = c('grey90'),
                vlineType = c( 'dashed'),
                vlineWidth = c(0.3),
                hline = c(0.05),
                hlineCol = c('grey90'),
                hlineType = c( 'dashed'),
                hlineWidth = c(0.5),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2.5,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x = 1.25, y = 4.4,color = "purple4" , size = 3,
           label =glue::glue("Relatively Upregulated in ASC", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -1.25, y = 4.4,color = "orange" , size = 3,
           label = glue::glue("Relatively Upregulated in SUB","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-2,2),breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 7),
        plot.subtitle = element_blank()) -> c
print(c)

# 
# ##Save plots 
# invisible(dev.off())
# rs_plot <- gridExtra::grid.arrange(c,b,a, ncol = 3)
# ggsave("manuscript/brain/results_figures/HypothesisTesting_ASC.png",rs_plot,height =5, width =15, dpi=600)


