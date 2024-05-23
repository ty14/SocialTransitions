
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


### Body weight.






#combine
ds.all.long <- rbind(ds.reorg.long,ds.control.long) %>% mutate(id = paste0(group,id))

ds.all.long %>% 
  group_by(Rank,time) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  ) -> ds.all.summary

ds.all.summary

ds.all.long %>% group_by(Rank) %>% count()

ds.all.long %>% filter(Rank == "Rank2", value>0)
ds.all.long %>% filter(Rank == "Rank2", value<=0)

ds.all.long %>% filter(Rank == "Rank3", value>0)
ds.all.long %>% filter(Rank == "Rank3", value<=0)


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
pre



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













########################### get post data 
library(tidyverse)

df <- read_csv("behavior/Post_WL.csv")
head(df)


fix10 <- df %>% filter(post_batchcage == "Batch 10 Cage 3")
fix10$loser
fix10["loser"][fix10["loser"]=="Cage 3-3"]<- "Cage 3-1"

df <- df %>% filter(post_batchcage != "Batch 10 Cage 3")
xf <- df %>% rbind(fix10)


## Just assuming that lowest ztime on first day is start of observations:

xf %>% filter(time=="1 hour") %>% group_by(batch) %>% summarise(diff = max(ztime)-min(ztime))

# get data within 70 minutes (4200 seconds) of first day of reorganization.

xf %>% 
  group_by(batch) %>% 
  filter(date == min(date)) %>% 
  filter(ztime < (4200+min(ztime))) -> xf1
table(xf1$group)
table(xf1$batch)

xf %>% 
  group_by(batch) %>% 
  filter(batch %in% c(1:7)) -> xf2
table(xf2$group)
table(xf2$batch)


# reorganized only
xf1r <- xf1 %>% filter(group=="reorganized")
xf2r <- xf2 %>% filter(group=="reorganized")


# create list of dataframes for each individual cage
l <- split(xf1r, xf1r$post_batchcage)
lx <- split(xf2r, xf2r$post_batchcage)

# just keep winner loser columns
l.wl <- l %>% map(~ ungroup(.)) %>% map(~ select(., winner,loser) %>% as.data.frame)
l.wlx <- lx %>% map(~ ungroup(.)) %>% map(~ select(., winner,loser) %>% as.data.frame)

# convert list to win-loss matrices organized by David's Scores
l.mat <- l.wl %>% map(~ compete::get_wl_matrix(.)) %>% map(~ compete::org_matrix(., method="ds"))
l.matx <- l.wlx %>% map(~ compete::get_wl_matrix(.)) %>% map(~ compete::org_matrix(., method="ds"))

# Get average David's Scores for each rank.
l.ds70 <- l.mat %>% map(~ compete::ds(.)) 
l.ds25 <- l.matx %>% map(~ compete::ds(.))

## we have to reinsert "0" DS for two animals that don't appear in any wins/losses
#code below helps id and fix this:
l.ds70  
which(names(l.ds70)=="Batch 4 Cage 4") #[32]
which(names(l.ds70)=="Batch 2 Cage 3") #[27]

# Batch 4 Cage 4  2nd position  [32]  animal id = 2-4
# Batch 2 Cage 3  2nd position  [27]  animal id = 1-4

# manually inputting
x = l.ds70[32]
x1 <- c(x[[1]][1],0,x[[1]][2:3])
names(x1)[2]<-"Cage 2-4"
l.ds70[[32]] <- x1

x = l.ds70[27]
x1 <- c(x[[1]][1],0,x[[1]][2:3])
names(x1)[2]<-"Cage 1-4"
l.ds70[[27]] <- x1


###
ds.df <- l.ds70  %>% map(~ matrix(., ncol=4)) %>% do.call("rbind",.)
colnames(ds.df) <-c("Rank 1", "Rank 2", "Rank 3", "Rank 4")
ds.df



ds.df25 <- l.ds25  %>% map(~ matrix(., ncol=4)) %>% do.call("rbind",.)
colnames(ds.df25) <-c("Rank 1", "Rank 2", "Rank 3", "Rank 4")
ds.df25



# Reorganize Data

ds.long <- ds.df %>% as.data.frame %>% pivot_longer(1:4, names_to = "Rank")
ds.long$id <- rep(1:(nrow(ds.long)/4), each=4)

ds.long25 <- ds.df25 %>% as.data.frame %>% pivot_longer(1:4, names_to = "Rank")
ds.long25$id <- rep(1:(nrow(ds.long25)/4), each=4)


ds.long$group <- "70min" 
ds.long25$group <- "25hr" 

#combine
ds.all.long <- rbind(ds.long,ds.long25) %>% mutate(id = paste0(group,id))

ds.all.long %>% 
  group_by(Rank,group) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  ) -> ds.all.summary

ds.all.summary

ds.all.long %>% group_by(Rank,group) %>% count()

ds.all.long$group <- factor(ds.all.long$group, levels = c("70min", "25hr"))
ds.all.summary$group <- factor(ds.all.summary$group, levels = c("70min", "25hr"))

post_reg <- ggplot() + 
  geom_line(data=ds.all.long, aes(x=Rank, y=value, group=id), alpha=.2, color="gray57") +
  theme_classic() +
  xlab("") +
  ylab("") +
  geom_hline(yintercept=0, lty=2, color="red", alpha=.5) +
  geom_errorbar(data=ds.all.summary, 
                aes(x=Rank, ymin=lqr, ymax=uqr), width=0.0, size=1, color="firebrick") +
  geom_point(data=ds.all.summary, 
             aes(x=Rank, y=median), size=3, shape=21, fill="white") +
  geom_line(data= ds.all.summary,aes(x=Rank,y=median,group=group), size=1.25, color="firebrick")+
  ggtitle("Post Reorganization") +
  facet_wrap(~group)+
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
  

bottom <- grid::textGrob("Mouse Rank", gp = grid::gpar(fontsize = 20))
left <- grid::textGrob("David's Score", gp = grid::gpar(fontsize = 20), rot = 90)

library(gridExtra)
ds <- grid.arrange(pre,post_reg,pre_con, ncol=1,bottom = bottom, left =left)
# ggsave("manuscript/behavior/results_figures/pre_ds_manuscript.png",ds,height =12, width =10, dpi=600)
