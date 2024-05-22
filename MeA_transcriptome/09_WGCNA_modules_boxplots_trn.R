# brown, black, cyan, magenta


#transition modules: brown, pink, blue, turquoise, darkorange, green, paleturquoise
#Group traits
#Getting metadata ready 
coldata <- read_csv("brain/sample70min_table.csv")
head(coldata)
str(coldata)

#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "TRANS", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "TRANS", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)


#just get samples I want
coldata$SampleID <- substr(coldata$SampleName, 7, 13)

coldata <- coldata %>%  dplyr::select(-SampleName) %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != 'ascenders') %>% 
  filter(condition1 != "CDOM") %>% 
  filter(condition1 != "CSUB") %>% 
  filter(region == "AMY")

row.names <- coldata$SampleID
row.names(coldata) <- row.names #Assigning row names from as sample names  
head(coldata)

MEs <- readRDS("manuscript/brain/results_RDS/WGCNA_trans_MEs_Power4.RDS")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:22, names_to = "Module") %>% 
  full_join(coldata) %>% filter(Module != "MEgrey") %>% 
  mutate(status = condition1)

head(ME_df)

ME_df$status <- gsub("TRANS", "TRN", ME_df$status)
ME_df$status <- factor(ME_df$status, levels = c("TRN", "DOM", "SUB"))
ME_df$Module <- gsub("ME", "", ME_df$Module)
#up
trans_dfx <- ME_df %>% filter(Module %in% c( 'brown', 'black', 'cyan', 'magenta'))
#down
trans_dfx1  <- ME_df %>% filter(Module%in% c('yellow', 'red', 'lightgreen','grey60', 'lightyellow'))
#dom and sub
trans_dfx <- ME_df %>% filter(Module %in% c("midnightblue", "royalblue"))

red <- ME_df %>% filter(Module == "red")
yellow <- ME_df %>% filter(Module == "yellow")
brown <- ME_df %>% filter(Module == "brown")
cyan <- ME_df %>% filter(Module == "cyan")



#boxplot of eigengenes and status 
source("functions/geom_boxjitter.R")
source("functions/newggtheme.R")

rp <-red%>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # facet_wrap(~Module, ncol =5)+
  scale_color_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+
  scale_fill_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+          
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Red: 574 genes")+ newggtheme_with_legends+ 
  theme(legend.position = "none", text =element_text(size =15))
rp 
ggsave(filename = "manuscript/brain/results_figures/WGCNA_TRN_red.png",
       rp,
       height = 4, width = 4, dpi = 600)


yp <-yellow%>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # facet_wrap(~Module, ncol =5)+
  scale_color_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+
  scale_fill_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+          
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Yellow: 686 genes")+ newggtheme_with_legends+ 
  theme(legend.position = "none", text =element_text(size =15))
yp 
ggsave(filename = "manuscript/brain/results_figures/WGCNA_TRN_yellow.png",
       yp,
       height = 4, width = 4, dpi = 600)

bp <-brown%>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # facet_wrap(~Module, ncol =5)+
  scale_color_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+
  scale_fill_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+          
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Brown: 963 genes")+ newggtheme_with_legends+ 
  theme(legend.position = "none", text =element_text(size =15))
bp
ggsave(filename = "manuscript/brain/results_figures/WGCNA_TRN_brown.png",
       bp,
       height = 4, width = 4, dpi = 600)
cp <-cyan%>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # facet_wrap(~Module, ncol =5)+
  scale_color_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+
  scale_fill_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+          
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Cyan: 392 genes")+ newggtheme_with_legends+ 
  theme(legend.position = "none", text =element_text(size =15))
cp

ggsave(filename = "manuscript/brain/results_figures/WGCNA_TRN_cyan.png",
       cp,
       height = 4, width = 4, dpi = 600)
