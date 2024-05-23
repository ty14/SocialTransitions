
#### David's score graph for manuscript code

# load packages and data
library(tidyverse)
library(compete)

df <- read_csv("Behavior/Pre_WL.csv")

head(df)

# how many cages in study
sum(table(df$batch, df$Cage)>0)


## 1. Create aggregated WL matrices for each cage.
# these need to be ordered by David's Scores

# Examining the number of rows of data for each cage/batch
table(df$batch, df$Cage)
table(df$pre_batchcage)

# average fights 
table(df$batch, df$Cage) %>% as.data.frame() %>%  summarise(sum = sum(Freq), mean = sum/65)
# create list of dataframes for each individual cage
l <- split(df, df$pre_batchcage)

# just keep winner loser columns
l.wl <- l %>% map(~ select(., winner,loser) %>% as.data.frame) 

# convert list to win-loss matrices organized by David's Scores
l.mat <- l.wl %>% map(~ compete::get_wl_matrix(.)) %>% map(~ compete::org_matrix(., method="ds"))


## test directional consistency
l.mat.dc <- lapply(l.mat, dc_test)
sapply(l.mat.dc, function(x) x[[1]]) #all p<.001
dcs <- sapply(l.mat.dc, function(x) x[[7]]) #dc values
median(dcs)
quantile(dcs,.25)
quantile(dcs,.75)



# Which are control and which are reorganized
df.groups <- df %>% select(pre_batchcage, group, time) %>% unique

df.groups

# DC descriptives for  control/reorganized
dcs.df <- data.frame(dc = as.numeric(dcs), pre_batchcage = names(dcs))
dcs.df$group <- df.groups$group[match(df.groups$pre_batchcage,dcs.df$pre_batchcage)]
aggregate(dc ~ group, data = dcs.df, FUN = median)
aggregate(dc ~ group, data = dcs.df, FUN = quantile, .25)
aggregate(dc ~ group, data = dcs.df, FUN = quantile, .75)


# split l.mat into 'control' and 'reorganized'
ids <- df.groups$group[match(names(l.mat), df.groups$pre_batchcage)]
l.mat.control <- l.mat[which(ids=="control")]
l.mat.reorg <- l.mat[which(ids=="reorganized")]


# Get average David's Scores for each rank. - for reorg vs control.
ds.reorg <- l.mat.reorg %>% map(~ compete::ds(.)) %>% map(~ matrix(., ncol=4)) %>% do.call("rbind",.)
ds.control <- l.mat.control %>% map(~ compete::ds(.)) %>% map(~ matrix(., ncol=4)) %>% do.call("rbind",.)
colnames(ds.control)<-colnames(ds.reorg)<-c("Rank1", "Rank2", "Rank3", "Rank4")


## Make boxplot of DS

# Reorganize Data
ds.reorg <- as.data.frame(ds.reorg)
ds.control <- as.data.frame(ds.control)
ds.reorg$group <- "Reorganized" 
ds.control$group <- "Control" 
ds.reorg$time <- rep(c("25hr", "70min", "25hr","70min"), times= c(4,20,24,4))
ds.control$time <- rep(c("25hr", "70min"), times= c(7,6))
ds.reorg.long <- ds.reorg %>% pivot_longer(1:4, names_to = "Rank")
ds.reorg.long$id <- rep(1:(nrow(ds.reorg.long)/4), each=4)
ds.reorg.long$Rank <- gsub("Rank", "Rank ", ds.reorg.long$Rank)
ds.control.long <- ds.control %>% pivot_longer(1:4, names_to = "Rank")
ds.control.long$id <- rep(1:(nrow(ds.control.long)/4), each=4)
ds.control.long$Rank <- gsub("Rank", "Rank ", ds.control.long$Rank)

## All together.
ds.all.long <- rbind(ds.control.long, ds.reorg.long)

ds.all.long %>% 
  group_by(Rank,group) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  )

# looking at David Scores by group 
ds.all.rank <- split(ds.all.long, ds.all.long$Rank)
t.test(value ~ group, data = ds.all.rank[[1]]) #alpha ds 
t.test(value ~ group, data = ds.all.rank[[2]]) #beta ds 
t.test(value ~ group, data = ds.all.rank[[3]]) #gamma ds
t.test(value ~ group, data = ds.all.rank[[4]]) #delta ds

# Descriptives of all
ds.all.long %>% 
  group_by(Rank) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  )

ds.all.long %>% filter(Rank=='Rank 2') %>% filter(value<0)
ds.all.long %>% filter(Rank=='Rank 3') %>% filter(value>0)


#get summary stats
ds.reorg.long %>% 
  group_by(Rank,time) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  ) -> ds.reorg.summary


ds.control.long %>% 
  group_by(Rank,time) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  ) -> ds.control.summary


#############################################################################################################

### Body weight.

# import body weight data
bw <- read_csv("Behavior/groups.csv")
bw$pre_idbatchcage <- paste0("Batch",bw$batch,"Cage",bw$precage)
bw <- bw %>% select(pre_idbatchcage, preid,wt_d4,wt_d8,wt_12)

#dataframe of DS
bw.ds <- do.call(rbind, lapply(l.mat %>% map(~ compete::ds(.)), function(x) data.frame(preid = names(x), ds = unname(x))))
bw.ds$Prerank<-as.numeric(gsub("^.*?(\\d)$", "\\1", rownames(bw.ds)))
bw.ds$pre_idbatchcage <- gsub("\\s+", "", sub("\\..*", "", rownames(bw.ds)))
bw.ds$preid <- as.numeric(bw.ds$preid)
bw <- full_join(bw,bw.ds)
bw <- bw[!is.na(bw$ds),]  # remove non observed mice
bw$batch <- sub("C.*", "", bw$pre_idbatchcage)
bw

bw %>% arrange(-wt_d8)  # presume this is 31.8  ???
bw %>% arrange(wt_12)  # presume this is 42.7   ???

bw[which(bw$wt_12==4.7), 5] <- 42.7
bw[which(bw$wt_d8==318), 4] <- 31.8


ggplot(bw, aes(x=factor(Prerank), y=wt_d4, fill=factor(Prerank))) + geom_boxplot() + theme_classic()
ggplot(bw, aes(x=factor(Prerank), y=wt_d8, fill=factor(Prerank))) + geom_boxplot() + theme_classic()
ggplot(bw, aes(x=factor(Prerank), y=wt_12, fill=factor(Prerank))) + geom_boxplot() + theme_classic()


library(lme4)
library(lmerTest)

# Fit mixed effects model
model1 <- lmer(ds ~ wt_d4 + (1|batch/pre_idbatchcage), data = bw)
summary(model1)

model2 <- lmer(ds ~ wt_d8 + (1|batch/pre_idbatchcage), data = bw)
summary(model2)

model3 <- lmer(ds ~ wt_12 + (1|batch/pre_idbatchcage), data = bw)
summary(model3)

#



# make supplemental figure.

library(viridis)
source("functions/geom_boxjitter.R")

range(bw$wt_d4)
range(bw$wt_d8)

pW1<-ggplot(bw,aes(x=factor(Prerank),y=wt_d4,color=factor(Prerank),fill= factor(Prerank)))+
  geom_boxjitter(outlier.color = NA, 
                 jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.025,  
                 errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  ylab("Body Weight g - Day 4") +
  xlab("Rank") +
  scale_y_continuous(limits = c(30, 46), breaks = c(30, 35, 40, 45)) +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  theme_classic() +
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        text = element_text(size = 20)
        
  )

pW2<-ggplot(bw,aes(x=factor(Prerank),y=wt_d8,color=factor(Prerank),fill= factor(Prerank)))+
  geom_boxjitter(outlier.color = NA, 
                 jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.025,  
                 errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  ylab("Body Weight g - Day 8") +
  xlab("Rank") +
  scale_y_continuous(limits = c(30, 46), breaks = c(30, 35, 40, 45)) +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  theme_classic() +
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        text = element_text(size = 20)
        
  )

library(gridExtra)
grid.arrange(pW1,pW2,nrow=1)




##########################################################################################################




### Make Figures of David's Scores

#combine
ds.all.long <- rbind(ds.reorg.long,ds.control.long) %>% mutate(id = paste0(group,id))
ds.all.long %>% group_by(Rank) %>% count()

ds.reorg.long$time <- factor(ds.reorg.long$time, levels = c("70min", "25hr"))
ds.reorg.summary$time <- factor(ds.reorg.summary$time, levels = c("70min", "25hr"))
ds.reorg.long70<- ds.reorg.long %>% filter(time == "70min")
ds.reorg.long25 <-ds.reorg.long %>% filter(time == "25hr")
ds.reorg.summary70<- ds.reorg.summary%>% filter(time == "70min")
ds.reorg.summary25<- ds.reorg.summary%>% filter(time == "25hr")



p1 <- ggplot() + 
  geom_line(data=ds.reorg.long, aes(x=Rank, y=value, group=id), alpha=.2, color="gray57") +
  theme_classic() +
  xlab("") +
  ylab("") +
  geom_hline(yintercept=0, lty=2, color="red", alpha=.5) +
  geom_errorbar(data=ds.reorg.summary, 
                aes(x=Rank, ymin=lqr, ymax=uqr), width=0.0, size=1, color="firebrick") +
  geom_point(data=ds.reorg.summary, 
             aes(x=Rank, y=median), size=3, shape=21, fill="white") + 
  geom_line(data= ds.reorg.summary,aes(x=Rank,y=median,group=time), size=1.25, color="firebrick")+
  ggtitle("Pre Reorganization") +
  facet_wrap(~time)+
  ylim(-6,6)+
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


pre <- p1

ds.control.long$time <- factor(ds.control.long$time, levels = c("70min", "25hr"))
ds.control.summary$time <- factor(ds.control.summary$time, levels = c("70min", "25hr"))

pre_con <- ggplot() + 
  geom_line(data=ds.control.long, aes(x=Rank, y=value, group=id), alpha=.2, color="gray57") +
  theme_classic() +
  xlab("") +
  ylab("") +
  geom_hline(yintercept=0, lty=2, color="red", alpha=.5) +
  geom_errorbar(data=ds.control.summary, 
                aes(x=Rank, ymin=lqr, ymax=uqr), width=0.0, size=1, color="firebrick") +
  geom_point(data=ds.control.summary, 
             aes(x=Rank, y=median), size=3, shape=21, fill="white") +
  geom_line(data= ds.control.summary,aes(x=Rank,y=median,group=time), size=1.25, color="firebrick")+
  ggtitle("Control") +
  facet_wrap(~time)+
  ylim(-6,6)+
  theme(axis.text.x = element_text(vjust = 1, size = 20),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),    
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text = element_text(size = 20)
        
  )

pre_con

