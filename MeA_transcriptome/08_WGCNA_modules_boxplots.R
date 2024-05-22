library(WGCNA)
library(tidyverse)


d_lm <- readRDS("manuscript/brain/results_RDS/WGCNA_CDOM_lm_result_all_outliersRemoved.RDS")
d_lm

#DES: yellow, red, green
#cdom: pink, black 

df <- read_csv("manuscript/brain/results_tables/coldata.csv")
MEs <- readRDS("manuscript/brain/results_RDS/WGCNA_CDOM_MEs_outliersRemoved.RDS")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:51, names_to = "Module") %>% 
  full_join(df) %>% filter(Module != "MEgrey")

head(ME_df)

ME_df$condition1
ME_df <- ME_df %>%
  mutate(status = condition1) %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB")
ME_df$Module <- gsub("ME", "", ME_df$Module)

ME_df$status <- factor(ME_df$status, levels = c("DES", "DOM", "CDOM"))

#DES 
des_df <- ME_df %>% filter(Module %in% c("yellow", "green"))
des_dfx <- ME_df %>% filter(Module %in% c("yellow","red", "green", "darkorange"))
#CDOM

cdom_df <- ME_df %>% filter(Module %in% c("black"))
yellow <- ME_df %>% filter(Module == "yellow")
red <- ME_df %>% filter(Module == "red")
green <- ME_df %>% filter(Module == "green")
black <- ME_df %>% filter(Module == "black")


#boxplot of eigengenes and status 
source("functions/geom_boxjitter.R")
source("functions/newggtheme.R")

y <- yellow %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Yellow Module: 617 genes")+ newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))


ggsave(filename = "manuscript/brain/results_figures/WGCNA_DES_yellowmodule.png",
       y,
       height = 4, width = 4, dpi = 600)

  
r <- red %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Red Module: 590 genes")+ newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))


ggsave(filename = "manuscript/brain/results_figures/WGCNA_DES_redmodule.png",
       r,
       height = 4, width = 4, dpi = 600)

b <- black %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Black Module: 548 genes")+ newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))


ggsave(filename = "manuscript/brain/results_figures/WGCNA_DES_blackmodule.png",
       b,
       height = 4, width = 4, dpi = 600)                                                               
       
g <- green %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Green Module: 593 genes")+ newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))


ggsave(filename = "manuscript/brain/results_figures/WGCNA_DES_greenmodule.png",
       g,
       height = 4, width = 4, dpi = 600)                                                               

   
              

pxx <- cdom_df %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  facet_wrap(~Module, scales = "free_y", ncol = 1)+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Social Disruption in Dominants")+ theme_classic()+ theme(legend.position = "none")


ggsave(filename = "manuscript/brain/results_figures/WGCNA_CDOM_black_modules.png",
       pxx,
       height = 4, width = 4, dpi = 600)


######################

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

MEs <- readRDS("manuscript/brain/results_RDS/WGCNA_trans_MEs_Power4_outliersRemoved.RDS")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:38, names_to = "Module") %>% 
  full_join(coldata) %>% filter(Module != "MEgrey") %>% 
  mutate(status = condition1)

head(ME_df)

ME_df$status <- gsub("TRANS", "TRN", ME_df$status)
ME_df$status <- factor(ME_df$status, levels = c("TRN", "DOM", "SUB"))
ME_df$Module <- gsub("ME", "", ME_df$Module)

trans_dfx <- ME_df %>% filter(Module %in% c( 'brown', 'pink', 'blue', 'turquoise', 'darkorange', 'green', 'paleturquoise'))

tur <- ME_df %>% filter(Module == "turquoise")
blue <- ME_df %>% filter(Module == "blue")
brown <- ME_df %>% filter(Module == "brown")
pink <- ME_df %>% filter(Module == "pink")
green <- ME_df %>% filter(Module == "green")
darkorange <- ME_df %>% filter(Module == "darkorange")
#boxplot of eigengenes and status 
source("functions/geom_boxjitter.R")
t <- tur %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
          scale_color_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+
          scale_fill_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+          
          labs(x = "Social condition",
          y = "Module eigengene",
          title = "Turquoise Module: 2715 genes")+ newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))



ggsave(filename = "manuscript/brain/results_figures/WGCNA_TRANs_tru_module.png",
       t,
       height = 4, width = 4, dpi = 600)

bl <- blue %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+
  scale_fill_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+          
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Blue Module: 1218 genes")+ newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))



ggsave(filename = "manuscript/brain/results_figures/WGCNA_TRANs_blue_module.png",
       bl,
       height = 4, width = 4, dpi = 600)

br <- brown %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+
  scale_fill_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+          
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Brown Module: 962 genes")+ newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))



ggsave(filename = "manuscript/brain/results_figures/WGCNA_TRANs_brown_module.png",
       br,
       height = 4, width = 4, dpi = 600)


p <- pink %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+
  scale_fill_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+          
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Pink Module: 382 genes")+ newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))



ggsave(filename = "manuscript/brain/results_figures/WGCNA_TRANs_pink_module.png",
       p,
       height = 4, width = 4, dpi = 600)


g <- green %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+
  scale_fill_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+          
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Green Module: 470 genes")+ newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))



ggsave(filename = "manuscript/brain/results_figures/WGCNA_TRANs_green_module.png",
       g,
       height = 4, width = 4, dpi = 600)

do <- darkorange %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+
  scale_fill_manual(values = c("#932667", "#29AF7FFF","#f7cb44"))+          
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Darkorange Module: 143 genes")+ newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))



ggsave(filename = "manuscript/brain/results_figures/WGCNA_TRANs_do_module.png",
      do,
       height = 4, width = 4, dpi = 600)


## Ascenders
a_lm <- readRDS("manuscript/brain/results_RDS/WGCNA_CSUB_lm_result_Power5_OutliersRemoved.RDS")
a_lm


df <- read_csv("manuscript/brain/results_tables/coldata.csv")
MEs <- readRDS("manuscript/brain/results_RDS/WGCNA_CSUB_MEs_Power5_OutliersRemoved.RDS")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:36, names_to = "Module") %>% 
  full_join(df) %>% filter(Module != "MEgrey")

head(ME_df)

ME_df$condition1
ME_df <- ME_df %>%
  mutate(status = condition1)
ME_df$Module <- gsub("ME", "", ME_df$Module)

ME_df$status <- factor(ME_df$status, levels = c("ASC", "SUB", "CSUB"))


asc_dfx <- ME_df %>% filter(Module %in% c( 'brown',  'blue', 'turquoise', 'lightgreen', 'darkmagenta'))

blue <- ME_df %>% filter(Module == "blue")
lightg <- ME_df %>% filter(Module == "lightgreen")

#boxplot of eigengenes and status 
source("functions/geom_boxjitter.R")
b <-  blue %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#de7065ff","#f7cb44", "#13306dff"))+
  scale_fill_manual(values = c("#de7065ff", "#f7cb44", "#13306dff"))+      

  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Blue Module: 1780 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))


ggsave(filename = "manuscript/brain/results_figures/WGCNA_ascension_bluemodule.png",
       b,
       height = 4, width = 4, dpi = 600)

lg <-  lightg%>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#de7065ff","#f7cb44", "#13306dff"))+
  scale_fill_manual(values = c("#de7065ff", "#f7cb44", "#13306dff"))+      
  
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Lightgreen Module: 131 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

lg
ggsave(filename = "manuscript/brain/results_figures/WGCNA_ascension_lgmodule.png",
       lg,
       height = 4, width = 4, dpi = 600)




##social distruption 
#Group traits
#Getting metadata ready 
coldata <- read_csv("brain/sample70min_table.csv")
head(coldata)
str(coldata)

#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "SD", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SD", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "SD", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "SD", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CON", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CON", coldata$condition1)


#just get samples I want
coldata$SampleID <- substr(coldata$SampleName, 7, 13)

coldata <- coldata %>%  dplyr::select(-SampleName) %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != 'ascenders') %>% 
  filter(region == "AMY")

row.names <- coldata$SampleID
row.names(coldata) <- row.names #Assigning row names from as sample names  
head(coldata)

MEs <- readRDS("manuscript/brain/results_RDS/WGCNA_socialdistruption_MEs_Power4_outliersRemoved.RDS")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:28, names_to = "Module") %>% 
  full_join(coldata) %>% filter(Module != "MEgrey") %>% 
  mutate(status = condition1)

head(ME_df)


ME_df$status <- factor(ME_df$status, levels = c("SD", "CON"))
ME_df$Module <- gsub("ME", "", ME_df$Module)

sd_dfx <- ME_df %>% filter(Module %in% c( 'red', 'pink',"lightyellow", "lightcyan" ))

red <- sd_dfx %>% filter(Module == "red")
pink <- sd_dfx %>% filter(Module == "pink")
ly <- sd_dfx %>% filter(Module == "lightyellow")
lc <- sd_dfx %>% filter(Module == "lightcyan")

#boxplot of eigengenes and status 
source("functions/geom_boxjitter.R")
source("functions/newggtheme.R")
r <-  red %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#21908CFF","#404688FF"))+
  scale_fill_manual(values = c("#21908CFF","#404688FF"))+      
  
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Red Module: 502 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

r
ggsave(filename = "manuscript/brain/results_figures/WGCNA_sd_redmodule.png",
       r,
       height = 5, width = 5, dpi = 300)


p <-  pink %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#21908CFF","#404688FF"))+
  scale_fill_manual(values = c("#21908CFF","#404688FF"))+      
  
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Pink Module: 338 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

p
ggsave(filename = "manuscript/brain/results_figures/WGCNA_sd_pinkmodule.png",
       p,
       height = 5, width = 5, dpi = 300)

lyp <- ly %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#21908CFF","#404688FF"))+
  scale_fill_manual(values = c("#21908CFF","#404688FF"))+      
  
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Lightyellow Module: 152 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

lyp
ggsave(filename = "manuscript/brain/results_figures/WGCNA_sd_lightyellow_module.png",
       lyp,
       height = 5, width = 5, dpi = 300)

lcp <- lc %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#21908CFF","#404688FF"))+
  scale_fill_manual(values = c("#21908CFF","#404688FF"))+      
  
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Lightcyan Module: 174 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

lcp
ggsave(filename = "manuscript/brain/results_figures/WGCNA_sd_lightcyan_module.png",
       lcp,
       height = 5, width = 5, dpi = 300)




MEs <- readRDS("manuscript/brain/results_RDS/WGCNA_ALL_MEs_Power4.RDS")
df <- read.csv("manuscript/brain/results_tables/coldata.csv", row.names = 1)
ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:27, names_to = "Module") %>% 
  full_join(df) %>% filter(Module != "MEgrey")

head(ME_df)

ME_df$condition1
ME_df <- ME_df %>%
  mutate(status = condition1)
ME_df$Module <- gsub("ME", "", ME_df$Module)

ME_df$status <- factor(ME_df$status, levels = c("DES","DOM","CDOM", "ASC", "SUB", "CSUB"))

# Transition all groups
pink <- ME_df %>% filter(Module == "pink")
brown <- ME_df %>% filter(Module == "brown")
green <- ME_df %>% filter(Module == "green")%>% filter(status %in% c("DES", "DOM", "CDOM"))


#ASC 
royalblue <- ME_df %>% filter(Module == "royalblue")
darkgrey <- ME_df %>% filter(module == "darkgrey")
  
#DES <- 
darkgreen <- ME_df %>% filter(Module == "darkgreen")
yellow <- ME_df %>% filter(Module == "yellow") %>% filter(status %in% c("DES", "DOM", "CDOM"))

#ASC w/cort   
lightcyan <- ME_df %>% filter(Module == "lightcyan")
midblue <- ME_df %>% filter(Module == "midnightblue")
black <- ME_df %>% filter(Module == "black")
darkred <- ME_df %>% filter(Module == "darkred")
#DOM w/cort 
darkgreen <- ME_df %>% filter(Module == "darkgreen")
greenyellow <-   ME_df %>% filter(Module == "greenyellow")
# cdom
black <- ME_df %>% filter(Module == "black")
salmon <- ME_df %>% filter(Module == "salmon")

#



source("functions/geom_boxjitter.R")
source("functions/newggtheme.R")

p<-  pink %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::plasma(6)) +
  scale_fill_manual(values = viridis::plasma(6))+    
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Pink: 415 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))
p

ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_pinkmodule.png",
       p,
       height = 5, width = 6, dpi = 600)


pb<-  brown %>% filter(status %in% c("DES", "DOM", "CDOM")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+     
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Brown: 1222 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

pb
ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALLdes_brownmodule.png",
       pb,
       height = 4, width = 4, dpi = 600)

pg<- green %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+ 
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Green: 521 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

pg
ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALLdes_greenmodule.png",
       pg,
       height = 4, width = 4, dpi = 600)


prb<-  royalblue %>% filter(status %in% c("ASC", "SUB", "CSUB")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#de7065ff","#f7cb44", "#13306dff"))+
  scale_fill_manual(values = c("#de7065ff", "#f7cb44", "#13306dff"))+       
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Royalblue: 131 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

prb
ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALLasc_royalbluemodule.png",
       prb,
       height = 4, width = 4, dpi = 600)
darkturq <- ME_df %>% filter(Module == "darkturquoise")
pdt<-  darkturq %>% filter(status %in% c("DES", "DOM", "CDOM")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+   
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Darkturquoise: 64 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))
pdt

ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_desdarkturqmodule.png",
       pdt,
       height = 4, width = 4, dpi = 600)

gy <- ME_df %>% filter(Module == "greenyellow")
dg<-  gy %>% filter(status %in% c("ASC", "SUB", "CSUB")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#de7065ff","#f7cb44", "#13306dff"))+
  scale_fill_manual(values = c("#de7065ff", "#f7cb44", "#13306dff"))+      
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Greenyellow: 295 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

dg
ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_ascgreenyellowmodule.png",
       dg,
       height = 4, width = 4, dpi = 600)


mag <- ME_df %>% filter(Module == "magenta")
magx<-  mag %>% filter(status %in% c("ASC", "SUB", "CSUB")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#de7065ff","#f7cb44", "#13306dff"))+
  scale_fill_manual(values = c("#de7065ff", "#f7cb44", "#13306dff"))+      
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Magenta: 361 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

magx
ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALLasc_magentamodule.png",
       magx,
       height = 4, width = 4, dpi = 600)


g <- ME_df %>% filter(Module == "orange")
gx<-  g %>% filter(status %in% c("ASC", "SUB", "CSUB")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#de7065ff","#f7cb44", "#13306dff"))+
  scale_fill_manual(values = c("#de7065ff", "#f7cb44", "#13306dff"))+      
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Orange: 60 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

gx
ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALLasc_orangemodule.png",
       gx,
       height = 4, width = 4, dpi = 600)




y<-  yellow %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+   
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Yellow: 1178 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

y
ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALLdes_yellowmodule.png",
       y,
       height = 4, width = 4, dpi = 600)

b<-  black %>% filter(status %in% c("DES", "DOM", "CDOM")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+    
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Black: 474 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

b
ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALLdes_blackmodule.png",
       b,
       height = 4, width = 4, dpi = 600)

s<-  salmon %>% filter(status %in% c("DES", "DOM", "CDOM")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+   
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Salmon: 238 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))

s
ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALLdes_salmonmodule.png",
       s,
       height = 4, width = 4, dpi = 600)

#ASC w/cort   
lightcyan <- ME_df %>% filter(Module == "lightcyan")
midblue <- ME_df %>% filter(Module == "midnightblue")
black <- ME_df %>% filter(Module == "black")
darkred <- ME_df %>% filter(Module == "darkred")
#DOM w/cort 
darkgreen <- ME_df %>% filter(Module == "darkgreen")
greenyellow <-   ME_df %>% filter(Module == "greenyellow")


lightcyan %>%
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se = F,  alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::plasma(6)) +
  scale_fill_manual(values = viridis::plasma(6))+    
  labs(x = "Normalized Corticosterone",
       y = "Module eigengene",
  title = "Lighcyan: 204 genes")+  newggtheme_with_legends + 
  theme( text =element_text(size =15)) -> p2
p2

ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_cort_lightcyanmodule.png",
       p2,
       height = 6.5, width = 6.5, dpi = 600)

midblue%>%
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se = F,  alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::plasma(6)) +
  scale_fill_manual(values = viridis::plasma(6))+    
  labs(x = "Normalized Corticosterone",
       y = "Module eigengene",
       title = "Midnightblue: 212 genes")+  newggtheme_with_legends + 
  theme( text =element_text(size =15)) -> p2x


ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_cort_midbluemodule.png",
       p2x,
       height = 6.5, width = 6.5, dpi = 600)

black %>%
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se = F,  alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::plasma(6)) +
  scale_fill_manual(values = viridis::plasma(6))+    
  labs(x = "Normalized Corticosterone",
       y = "Module eigengene",
       title = "Black: 474 genes")+  newggtheme_with_legends + 
  theme( text =element_text(size =15)) -> p2b


ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_cort_blackmodule.png",
       p2b,
       height = 6.5, width = 6.5, dpi = 600)

darkred%>%
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se = F,  alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::plasma(6)) +
  scale_fill_manual(values = viridis::plasma(6))+    
  labs(x = "Normalized Corticosterone",
       y = "Module eigengene",
       title = "Darkred: 118 genes")+  newggtheme_with_legends + 
  theme( text =element_text(size =15)) -> p2dr


ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_cort_darkredmodule.png",
       p2dr,
       height = 6.5, width = 6.5, dpi = 600)

darkgreen %>%
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se = F,  alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::plasma(6)) +
  scale_fill_manual(values = viridis::plasma(6))+    
  labs(x = "Normalized Corticosterone",
       y = "Module eigengene",
       title = "Darkgreen: 89 genes")+  newggtheme_with_legends + 
  theme( text =element_text(size =15)) -> p2dg


ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_cort_darkgreenmodule.png",
       p2dg,
       height = 6.5, width = 6.5, dpi = 600)

greenyellow%>%
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se = F,  alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::plasma(6)) +
  scale_fill_manual(values = viridis::plasma(6))+    
  labs(x = "Normalized Corticosterone",
       y = "Module eigengene",
       title = "Greenyellow: 295 genes")+  newggtheme_with_legends + 
  theme( text =element_text(size =15)) -> p2y


ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_cort_greenyellow.png",
       p2y,
       height = 6.5, width = 6.5, dpi = 600)





darkturq <- ME_df %>% filter(Module == "darkturquoise")
pdt<-  darkturq %>% filter(status %in% c("DES", "DOM", "CDOM")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+   
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Darkturquoise: 64 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))
pdt

ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_desdarkturqmodule.png",
       pdt,
       height = 4, width = 4, dpi = 600)



dg <- ME_df %>% filter(Module == "darkgreen")
pdt<-  dg%>% filter(status %in% c("DES", "DOM", "CDOM")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+   
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Darkgreen: 89 genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))
pdt

ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_desdarkgreenmodule.png",
       pdt,
       height = 4, width = 4, dpi = 600)






dg <- ME_df %>% filter(Module == "pink")
pdt<-  dg %>% filter(status %in% c("DES", "DOM", "CDOM")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#440154FF",  "#29AF7FFF", "#2D708EFF"))+
  scale_fill_manual(values = c("#440154FF","#29AF7FFF", "#2D708EFF"))+   
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Pink: 415  genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))
pdt

ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_desPINK.png",
       pdt,
       height = 4, width = 4, dpi = 600)




dg <- ME_df %>% filter(Module == "pink")
pdt<-  dg %>% filter(status %in% c("ASC", "SUB", "CSUB")) %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#de7065ff","#f7cb44", "#13306dff"))+
  scale_fill_manual(values = c("#de7065ff", "#f7cb44", "#13306dff"))+    
  labs(x = "Social condition",
       y = "Module eigengene",
       title = "Pink: 415  genes")+  newggtheme_with_legends+ theme(legend.position = "none", text =element_text(size =15))
pdt

ggsave(filename = "manuscript/brain/results_figures/WGCNA_ALL_subPINK.png",
       pdt,
       height = 4, width = 4, dpi = 600)



