#boxplots for 8 transition genes 

source("MeA_transcriptome/05_TransitionGenes_withControls.R")

trans_down <- ctrans_down$.
trans_up <- ctrans_up$.
trans <- c(trans_up, trans_down)

#DES 
ex <- readRDS("MeA_transcriptome/results_RDS/limma_vdl_MeA_CDOM.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)


x %>% 
  filter(symbol %in% trans) -> xex



xex2 <- xex %>% pivot_longer(cols = 2:21, names_to = "ids")

p <- xex2 %>% full_join(id)

p$group <- factor(p$group, levels = c("DES", "DOM", "CDOM"))
p$symbols <- factor(p$symbol, levels = c("Unc13a", 'Car12', 'Tesc', 'Lpl', 'Pitpnm2', 'Rai2', "Kif21a", "Nap1l5"))
source('functions/geom_boxjitter.R')
library(viridis)

p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+
  facet_wrap(factor(symbol,levels = c("Unc13a", 'Car12', 'Tesc', 'Lpl', 'Pitpnm2', 'Rai2', "Kif21a", "Nap1l5")) ~ ., scales = 'free', ncol =8)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))
p1

#DES 
ex <- readRDS("MeA_transcriptome/results_RDS/limma_vdl_MeA_CSUB.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)


x %>% 
  filter(symbol %in% trans) -> xex

xex2 <- xex %>% pivot_longer(cols = 2:21, names_to = "ids")

p <- xex2 %>% full_join(id)

p$group <- factor(p$group, levels = c("ASC", "SUB", "CSUB"))
p$symbols <- factor(p$symbol, levels = c("Unc13a", 'Car12', 'Tesc', 'Lpl', 'Pitpnm2', 'Rai2', "Kif21a", "Nap1l5"))

p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#de7065ff","#f7cb44", "#13306dff"))+
  scale_fill_manual(values = c("#de7065ff", "#f7cb44", "#13306dff"))+
  facet_wrap(factor(symbol,levels = c("Unc13a", 'Car12', 'Tesc', 'Lpl', 'Pitpnm2', 'Rai2', "Kif21a", "Nap1l5")) ~ ., scales = 'free', ncol =8)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))

p1
