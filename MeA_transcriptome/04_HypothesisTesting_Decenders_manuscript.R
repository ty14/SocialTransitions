library(EnhancedVolcano)
library(tidyverse)

my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results_RDS/limma_MEA_CDOM.RDS") %>% 
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
# ARGs from Tyssowski 2018 - (c)
pr <- read_csv("manuscript/brain/gene_sets/PrimaryResponseGenes_fromTyssowski.csv")
pr %>%
  as_tibble() %>%
  .$`gene name` -> prx

prx <- unique(prx)
prx_names <- c("Primary Response Genes")

# 4. aggressive genes

ar <- c("Maoa", "Erbb4", "Gria3", 'Mecp2', 'Prnp', 'Avpra1', 'Chrmp2b', "En2", "Fgf14", "Hdac4", "kcnj18", "Lrrc7", "Ache", "Cadm1", 
        "Crhr1", "Dnajb5","Ecm1","Eef1a2", "Ehmt1","Gad2", "Gdl1",
        "Grid1", "Grn", "Gsk3a", "Hsf1", "Lama2", "Mapk15", "Mme", 
        "Nfkb1", "Npy1r", "Osmr", "Pnoc", "Rbfox1", "Spast","Syn1", "Wdr62")

ar %>%
  as.data.frame() -> arx



prx <- unique(prx)


###1
head(y3a)

y3a  %>% 
  filter(symbol %in% ssx) %>%
  dplyr::select(symbol, logFC, P.Value) %>% 
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"DES genes","DOM genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","DOM genes", "N.S."))) %>% 
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
                title = ssx_names,
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
           label =glue::glue("Relatively Upregulated in DES", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -1.25, y = 4.4,color = "orange" , size = 3,
           label = glue::glue("Relatively Upregulated in DOM","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-2,2),breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust =0.5),
        plot.subtitle = element_blank()) -> a
print(a)


# 2. 
y3a  %>% 
  filter(symbol %in% rsx) %>%
  dplyr::select(symbol, logFC, P.Value) %>% 
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"DES genes","DOM genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","DOM genes", "N.S."))) %>% 
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
for_label <- c(for_label, "Mybpc1")

EnhancedVolcano(df,
                lab = df$symbol,
                selectLab= for_label,
                x = 'logFC',
                y = 'P.Value',
                title = rsx_names,
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
           label =glue::glue("Relatively Upregulated in DES", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -2, y = 4.4,color = "orange" , size = 3,
           label = glue::glue("Relatively Upregulated in DOM","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-3,3),breaks = c(-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3.0))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust =0.5),
        plot.subtitle = element_blank()) -> b
print(b)


# c
y3a  %>% 
  filter(symbol %in% prx) %>%
  dplyr::select(symbol, logFC, P.Value) %>% 
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"DES genes","DOM genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","DOM genes", "N.S."))) %>% 
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
  unique() -> df

dpr <- df %>% filter(Sig == "DES genes")
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
                title = prx_names,
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
           label =glue::glue("Relatively Upregulated in DES", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -1.25, y = 4.4,color = "orange" , size = 3,
           label = glue::glue("Relatively Upregulated in DOM","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-2,2),breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust =0.5),
        plot.subtitle = element_blank()) -> c
print(c)


##Save plots 
invisible(dev.off())
rs_plot <- gridExtra::grid.arrange(c,b,a, ncol = 3)
ggsave("manuscript/brain/results_figures/HypothesisTesting_DES.png",rs_plot,height =5, width =15, dpi=600)




#4. 
y3a  %>% 
  filter(symbol %in% ar) %>%
  dplyr::select(symbol, logFC, P.Value) %>% 
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"DES genes","DOM genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","DOM genes", "N.S."))) %>% 
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
for_label <- c(for_label, "Mybpc1")

EnhancedVolcano(df,
                lab = df$symbol,
                selectLab= for_label,
                x = 'logFC',
                y = 'P.Value',
                title = "Aggressive Genes",
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
           label =glue::glue("Relatively Upregulated in DES", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -2, y = 4.4,color = "orange" , size = 3,
           label = glue::glue("Relatively Upregulated in DOM","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-3,3),breaks = c(-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3.0))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust =0.5),
        plot.subtitle = element_blank()) -> b
print(b)

