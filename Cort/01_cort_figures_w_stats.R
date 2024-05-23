# cort analysis 
library(viridis)
library(glue)
library(lme4)
library(MASS)
library(car)
library(tidyverse)

#source
df <- read_csv("Cort/FullCort.csv") 
df<- df %>% dplyr::select(-Well, -final_conc_ng_ul) %>% unique(.)


#changing times
df$time <- ifelse(df$time == "1 hr", "70 min", df$time)
df$time <- ifelse(df$time == "25hr", "25 hr", df$time)


##Changing period
df$period <- ifelse(df$period == "Pre", 0, df$period)
df$period <- ifelse(df$period == "Post" & df$time == "70 min", 1, df$period)
df$period <- ifelse(df$period == "Post" & df$time == "25 hr", 2, df$period)

df <- df %>% 
  dplyr::select( -time)

##adding domgroup 
df$domgroup <- ifelse(df$Prerank == 3|df$Prerank==4, "Subordinate", "Dominant")  
## post domgroup 
df$domgroupPost <- ifelse(df$Postrank == 3|df$Postrank==4, "Subordinate", "Dominant")  

## Getting rid of duplicates 
df <- unique(df) #no duplicates


## making sure everything is a factor 
df$period <- as.factor(df$period)
df$Prerank <- factor(df$Prerank, levels = c("1","3", "4"))
df$Postrank <- factor(df$Postrank, levels = c("1", "3","4"))
df$batch <- as.factor(df$batch)
df$plate <- as.factor(df$plate)
df$pre_idbatchcage <- as.factor(df$pre_idbatchcage)
df$condition<- factor(df$condition, levels = c("control", "same", "ascenders", "descenders"))
df$post_idbatch <- as.factor(df$post_idbatch)
df$group <- as.factor(df$group)

#period 0 = pre, 1 = 70min, 2 = 25hr

dfx <- df %>% dplyr::select(plate, mean_con_ng_ul, batch, period, precage, post_idbatch, group, Prerank, Postrank, condition, domgroup, domgroupPost)
head(dfx)


source('functions/geom_boxjitter.R')

#### No differences between controls and reorganized group pre reorganization 

g1 <- df %>%filter(period == "0")
g1$group <- ifelse(g1$group == "reorganized", "Reorganized","Control")
g1$group <- factor(g1$group, levels = c("Control", "Reorganized"))

g1$precage <- gsub("Cage", "", g1$precage)
#stats
hist(g1$mean_con_ng_ul) #not normal 

# I think this should be cage not id. 
pre.glm <-glmer(mean_con_ng_ul~group +(1|precage)+(1|batch)+(1|plate), data =g1,family = Gamma(link = "log"))

# pre.glm <-glmer(mean_con_ng_ul~group +(1|pre_id)+(1|batch)+(1|plate), data =g1,family = Gamma(link = "log"))


summary(pre.glm)
AIC(pre.glm)

#checks 
acf(resid(pre.glm))
qqPlot(resid(pre.glm))
hist(resid(pre.glm))
plot(pre.glm)

durbinWatsonTest(resid(pre.glm))


#figure 
p1 <- ggplot(g1, aes(group, mean_con_ng_ul,color =group, fill = group))+
  geom_boxjitter(outlier.color = NA, 
                 jitter.shape = 21,
                 alpha = 0.4,
                 jitter.height = 0.02, 
                 jitter.width = 0.030, 
                 errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "",
       y = paste0("Corticosterone (ng/ml)", "\n", "Before Reorganzation"))+
  scale_color_manual(values = c("#22A884FF", "#414487FF")) +
  scale_fill_manual(values = c("#22A884FF", "#414487FF")) +theme_bw() +
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),    
        axis.title.x = element_text(hjust = 1),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        text = element_text(size = 20)
        
  )

p1

#ggsave("manuscript/cort/results_figures/CORTB4_suppl.png", p1,height=4.5, width = 4.5, dpi = 600)



#### no difference between ranks pre reorganization 
g1$prerank2 <- factor(g1$Prerank, levels = c(3,4,1))

# compared to rank 1
prerank.glm <-glmer(mean_con_ng_ul~domgroup +(1|pre_id)+(1|batch)+(1|plate), data =g1,family = Gamma(link = "log"))
summary(prerank.glm)
AIC(prerank.glm)

#checks 
acf(resid(prerank.glm))
qqPlot(resid(prerank.glm))
hist(resid(prerank.glm))
plot(prerank.glm)

durbinWatsonTest(resid(prerank.glm))


# domgroup figure 
p2<- ggplot(g1, aes(domgroup, mean_con_ng_ul,color =domgroup, fill = domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "",
       y = paste0("Corticosterone (ng/ml)", "\n", "Before Reorganzation"))+
  scale_color_manual(values = viridis::viridis(2)) +
  scale_fill_manual(values = viridis::viridis(2)) +theme_bw() +
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),    
        axis.title.x = element_text(hjust = 1),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        text = element_text(size = 20)
        
  )

p2

# ggsave("manuscript/cort/results_figures/CORTB4_domgroup_suppl.png", p2,height=4.5, width = 4.5, dpi = 600)



# CORT suppl figure together
pp <- gridExtra::grid.arrange(p1,p2, nrow=1)
pp

#ggsave("manuscript/cort/results_figures/groupDOMtogether_suppl.png", pp,height=4, width = 10, dpi = 600)


#### Differences in post - pre cort for controls and reorganized groups ####
df1 <- read_csv("cort/FullCort.csv") 
head(df1)
df1<- df1 %>% dplyr::select(-Well, -final_conc_ng_ul) %>%unique(.)
df1$time <- ifelse(df1$time == "1 hr", "Reorganized 70 min", df1$time)
df1$time <- ifelse(df1$time == "25hr", "Reorganized 25 hr", df1$time)
## Time for the control data 
df1$time <- ifelse(df1$time == "Reorganized 70 min" & df1$group == "control","Control 70 min", df1$time)
df1$time <- ifelse(df1$time == "Reorganized 25 hr" & df1$group == "control","Control 25 hr", df1$time)
df1$pre_id <- gsub("Batch", "-", df1$pre_idbatch)

##adding domgroup 
df1$domgroup <- ifelse(df1$Prerank == 3|df1$Prerank==4, "Subordinate", "Dominant")  

##adding final condition groups
df1$condition1<- ifelse(df1$condition == "ascenders",  "ASC", df1$condition)
df1$condition1 <-ifelse(df1$condition == "descenders",  "DES", df1$condition1)
df1$condition1 <- ifelse(df1$condition == "same" & df1$domgroup== "Dominant", "DOM", df1$condition1)
df1$condition1 <- ifelse(df1$condition == "same" & df1$domgroup== "Subordinate", "SUB", df1$condition1)
df1$condition1 <- ifelse(df1$condition == "control" & df1$time == "Control 70 min", "Control 70 min", df1$condition1) 
df1$condition1 <- ifelse(df1$condition == "control" & df1$time == "Control 25 hr", "Control 25 hr", df1$condition1)
df1$condition1 <- factor(df1$condition1, levels = c("DES","DOM","SUB","ASC", "Control 70 min", "Control 25 hr"))


dp <- df1 %>% dplyr::select(mean_con_ng_ul, plate, pre_id,batch, period, time,pre_idbatchcage,post_idbatch, Prerank, Postrank, condition1)
head(dp)

dp <-dp %>% 
  group_by(pre_idbatchcage) %>% 
  pivot_wider(values_from = mean_con_ng_ul, names_from = period) %>% 
  mutate(diff = Post - Pre) %>% ungroup()
head(dp)

dp$time2 <- ifelse(grepl("70 min", dp$time), "70 min", "25 hr")
dp$time2 <- factor(dp$time2, levels = c("70 min", "25 hr"))

dp$condition1 <- gsub("25 hr", "", dp$condition1)
dp$condition1 <- gsub("70 min", "", dp$condition1)

dp$condition1 <- factor(dp$condition1, levels =c( "DOM", "DES","SUB", "ASC","Control "))
  
# Z-Score Standardization: (X – μ) / σ
dm <- mean(dp$diff, na.rm =T)
ds <- sd(dp$diff, na.rm =T)
dp$zdiff <- (dp$diff -dm)/ds
hist(dp$zdiff)

#Min-Max Normalization: (X – min(X)) / (max(X) – min(X))
dp$ndiff <- (dp$diff - min(dp$diff, na.rm =T))/ (max(dp$diff, na.rm = T) - min(dp$diff, na.rm =T)) 
hist(dp$ndiff)

#Min-Max Normalization: (X – min(X)) / (max(X) – min(X))
dp$nPost <- (dp$Post - min(dp$Post, na.rm =T))/ (max(dp$Post, na.rm = T) - min(dp$Post, na.rm =T)) 
hist(dp$nPost)




diff <- ggplot(dp, aes(condition1, diff,color =condition1, fill = condition1))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::plasma(6)) +
  scale_fill_manual(values = viridis::plasma(6))+
  labs(y = paste0("Difference in","\n", "Corticosterone (ng/ml)"), 
       x = "")+
  facet_wrap(~time2) +
  ylim(-400, 600)+
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),    
        axis.title.x = element_text(hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 20)
        
  ) 


fig2D <- diff + geom_hline(yintercept=0, lty=2, col='gray68')
# ggsave("manuscript/cort/results_figures/final_CORT_fig2D.png", fig2D,height=6, width = 11, dpi = 600)



## Stats 70 min

head(dp)
dp70 <- dp %>% filter(time2 == "70 min")
hist(dp70$diff)
dp70$condition1 <- gsub(" ", "", dp70$condition1)

dp70$condition1 <- factor(dp70$condition1, levels = c('Control',"ASC", "SUB", "DES", "DOM"))
dp70$condition2 <- factor(dp70$condition1, levels = c("ASC", "SUB", "DES", "DOM", "Control"))
dp70$condition3<- factor(dp70$condition1, levels = c("SUB", "DES", "DOM", "Control", "ASC"))
dp70$condition4<- factor(dp70$condition1, levels = c("DES", "DOM", "Control", "ASC", "SUB"))


## significantly different from 0 ?

dp70 %>%
  split(.$condition1) %>%
  map(~t.test(.x$diff))



#stats that compare groups to controls

c70.lmx <-lmer(diff~condition1, data =dp70)

dp70$pre_idbatchcage
dp70$post_idbatch
as.data.frame(dp70)

#need to make prebatchcage & postbatchcage (batch + Prerank)
dp70$postbatchcageID<-paste0(unlist(regmatches(dp70$post_idbatch, gregexpr("Batch\\S*", dp70$post_idbatch))),"-",dp70$Prerank)

dp70$prebatchcageID <- sub("^.", "", dp70$pre_idbatchcage)



c70.lmx  <-lmer(diff~condition1 +(1|prebatchcageID)+(1|postbatchcageID), data =dp70)
summary(c70.lmx)
AIC(c70.lmx)  

#checks 
acf(resid(c70.lmx))
qqPlot(resid(c70.lmx))
hist(resid(c70.lmx))
plot(c70.lmx)
durbinWatsonTest(resid(c70.lmx))



#stats that compare groups to ascenders

c70.lm2  <-lmer(diff~condition2 +(1|prebatchcageID)+(1|postbatchcageID), data =dp70)

# c70.lm2 <-lm(diff~condition2, data =dp70)
summary(c70.lm2)
#checks
AIC(c70.lm2)  
acf(resid(c70.lm2))
qqPlot(resid(c70.lm2))
hist(resid(c70.lm2))
plot(c70.lm2)
durbinWatsonTest(resid(c70.lm2))


#stats that compare groups to subordinates
c70.lm3  <-lmer(diff~condition3 +(1|prebatchcageID)+(1|postbatchcageID), data =dp70)

# c70.lm3 <-lm(diff~condition3, data =dp70)
summary(c70.lm3)
#checks
AIC(c70.lm3)  
acf(resid(c70.lm3))
qqPlot(resid(c70.lm3))
hist(resid(c70.lm3))
plot(c70.lm3)
durbinWatsonTest(resid(c70.lm3))


#stats that compare groups to descenders
c70.lm4  <-lmer(diff~condition4 +(1|prebatchcageID)+(1|postbatchcageID), data =dp70)

# c70.lm4 <-lm(diff~condition4, data =dp70)
summary(c70.lm4)
AIC(c70.lm4)
#checks
acf(resid(c70.lm4))
qqPlot(resid(c70.lm4))
hist(resid(c70.lm4))
plot(c70.lm4)
durbinWatsonTest(resid(c70.lm4))






#25 hr stats 
dp25 <- dp %>% filter(time2 == "25 hr")
hist(dp25$diff)
dp25$condition1 <- gsub(" ", "", dp25$condition1)
dp25$condition1 <- factor(dp25$condition1, levels = c('Control',"ASC", "SUB", "DES", "DOM"))
dp25$condition2 <- factor(dp25$condition1, levels = c("ASC", "SUB", "DES", "DOM", "Control"))
dp25$condition3<- factor(dp25$condition1, levels = c("SUB", "DES", "DOM", "Control", "ASC"))
dp25$condition4<- factor(dp25$condition1, levels = c("DES", "DOM", "Control", "ASC", "SUB"))



## significantly different from 0 ?

dp25 %>%
  split(.$condition1) %>%
  map(~t.test(.x$diff))


#need to make prebatchcage & postbatchcage (batch + Prerank)
dp25$postbatchcageID<-paste0(unlist(regmatches(dp25$post_idbatch, gregexpr("Batch\\S*", dp25$post_idbatch))),"-",dp25$Prerank)
dp25$prebatchcageID <- sub("^.", "", dp25$pre_idbatchcage)



#stats that compare groups to controls
c25.lmx  <-lmer(diff~condition1 +(1|prebatchcageID)+(1|postbatchcageID), data =dp25)
#c25.lmx <-lm(diff~condition1, data =dp25)
summary(c25.lmx)
AIC(c25.lmx)  
#checks
acf(resid(c25.lmx))
qqPlot(resid(c25.lmx))
hist(resid(c25.lmx))
plot(c25.lmx)
durbinWatsonTest(resid(c25.lmx))

#stats that compare groups to ascenders
c25.lm2  <-lmer(diff~condition2 +(1|prebatchcageID)+(1|postbatchcageID), data =dp25)
#c25.lm2 <-lm(diff~condition2, data =dp25)
summary(c25.lm2)
AIC(c25.lm2)  


