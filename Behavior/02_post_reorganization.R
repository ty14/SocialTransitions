
#### David's score graph for manuscript code

source("manuscript/behavior/01_pre_organization.R")


##### Get Post Data 70 min reorg

df <- read_csv("behavior/Post_WL.csv")
head(df)

# Fix typo
fix10 <- df %>% filter(post_batchcage == "Batch 10 Cage 3")
fix10$loser
fix10["loser"][fix10["loser"]=="Cage 3-3"]<- "Cage 3-1"


df <- df %>% filter(post_batchcage != "Batch 10 Cage 3")

#xf <- df %>% rbind(fix10)


# get data within 70 minutes (4200 seconds) of first day of reorganization.
# 
# xf %>% 
#   group_by(batch) %>% 
#   filter(date == min(date)) %>% 
#   filter(ztime < (4200+min(ztime))) -> xf1

## use starttime
zt <- read_csv("behavior/Post_starttimes.csv")
zts <- zt[zt$times=="startTime",]

df <- df %>% rbind(fix10) #needs to be added

xf <- full_join(df, zts %>% dplyr::select(batch, starttime = ztime))

xf %>% 
  group_by(batch) %>% 
  filter(date == min(date)) %>% 
  filter(ztime < (4200+starttime)) -> xf1



# 25 hour data batches [includes 70 min data]
xf %>% 
  group_by(batch) %>% 
  filter(batch %in% c(1:7)) -> xf2

# reorganized only
xf1r <- xf1 %>% filter(group=="reorganized") #all 70mins
xf2r <- xf2 %>% filter(group=="reorganized") #25hrs - all data

# create list of dataframes for each individual cage
l <- split(xf1r, xf1r$post_batchcage) #all 70mins
lx <- split(xf2r, xf2r$post_batchcage) #25hrs -all data

# just keep winner loser columns
l.wl <- l %>% map(~ ungroup(.)) %>% map(~ select(., winner,loser) %>% as.data.frame)   #all 70 mins
l.wlx <- lx %>% map(~ ungroup(.)) %>% map(~ select(., winner,loser) %>% as.data.frame) #25hrs all data

# convert list to win-loss matrices organized by David's Scores
l.mat <- l.wl %>% map(~ compete::get_wl_matrix(.)) %>% map(~ compete::org_matrix(., method="ds"))   #all 70mins
l.matx <- l.wlx %>% map(~ compete::get_wl_matrix(.)) %>% map(~ compete::org_matrix(., method="ds")) #25hrs all data


# Directional consistency test - 70 mins

l.mat.dc <- lapply(l.mat, dc_test)
sapply(l.mat.dc, function(x) x[[1]]) #all p<.001
dcs <- sapply(l.mat.dc, function(x) x[[7]]) #dc values
median(dcs)
quantile(dcs,.25)
quantile(dcs,.75)



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


### David's Scores after 70 mins
ds.df <- l.ds70  %>% map(~ matrix(., ncol=4)) %>% do.call("rbind",.)
colnames(ds.df) <-c("Rank 1", "Rank 2", "Rank 3", "Rank 4")
ds.df


# David's Scores after 25hrs
ds.df25 <- l.ds25  %>% map(~ matrix(., ncol=4)) %>% do.call("rbind",.)
colnames(ds.df25) <-c("Rank 1", "Rank 2", "Rank 3", "Rank 4")
ds.df25


## Put in cage type
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


## 70 minute DS summary:
# Check...  these DS for betas seem higher than Fig 2A - check.

ds.dfdf <-cbind(
  as.data.frame(ds.df),
 group = c("Alphas", "Betas", "Gammas", "Deltas")[as.numeric(substrRight(names(l),1))]
)

ds.dfdf %>%
  pivot_longer(cols=1:4) %>%
  group_by(group, name) %>%
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  ) -> ds.dfdf.sum

ds.dfdf.sum

## 70 minute DS Supplemental Figure X
library(viridis)
source("functions/geom_boxjitter.R")

ds.dfdf %>%
  pivot_longer(cols=1:4) -> ds.dfdf.x


#make sure in correct order of levels
ds.dfdf.x$group <- factor(ds.dfdf.x$group, levels=c("Alphas","Betas","Gammas","Deltas"))

fig2C <- ggplot(ds.dfdf.x,aes(x=name,y=value,color=group,fill= group))+
  geom_boxjitter(outlier.color = NA, 
                 jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.025,  
                 errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  ylab("David's Score (70 min)") +
  xlab("New Rank") +
  scale_color_manual(name=paste0("Previous",'\n', "Rank"),values = viridis::viridis(6)) +
  scale_fill_manual(name=paste0("Previous",'\n', "Rank"),values = viridis::viridis(6))+
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),
        legend.position = c(.9,.8),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        text = element_text(size = 20)) + 
  geom_hline(yintercept = 0, lty=2, lwd=1, color="gray44")
ggsave("manuscript/behavior/results_figures/final_fig2C_ds.png",fig2C,height =6, width =7, dpi=600)


### Fig 2A - David's Scores Across Time.

#check DS
# lapply(l.wl, function(x) compete::ds(compete::get_wl_matrix(x)))  #70min
# lapply(l.wlx, function(x) compete::ds(compete::get_wl_matrix(x))) #all


### This should be done on 'xf' - all data, not the extracted 70min or 25hr cages....


## function to get DS across time for each cage
ds_time <- function(x){
  results <- NULL
  for(i in 1:nrow(x)){
    results[[i]] <- compete::ds(compete::get_wl_matrix(x[1:i,1:2])) # rows 1 thru i
  }
  ids <- unlist(lapply(results,names))
  vals <- unlist(results)
  dx <- data.frame(ids,vals, time = rep(1:nrow(x), times =     unlist(lapply(results, length)) ) )
  return(dx)     
}

#ds_time(l.wl[[1]])
#ds_time(l.wlx[[1]])


## Create dataframe of DS for all cages - both 70 min and 25hrs.

xfr <- xf %>% filter(group=="reorganized") #all - reorganized only
l.all <- split(xfr, xfr$post_batchcage) #all - create list of dataframes for each individual cage
l.wl.all <- l.all %>% map(~ ungroup(.)) %>% map(~ select(., winner,loser) %>% as.data.frame)   #all - # just keep winner loser columns
l.mat.all <- l.wl.all %>% map(~ compete::get_wl_matrix(.)) %>% map(~ compete::org_matrix(., method="ds"))   #all -# convert list to win-loss matrices organized by David's Scores

#ds_time(l.wl.all[[1]])

## do for all cages...
results.all <- NULL
for(i in 1:length(l.wl.all)){
  results.all[[i]] <-  ds_time(l.wl.all[[i]])
}


# ## do for all cages...
# results1 <- NULL
# for(i in 1:length(l.wl)){
#   results1[[i]] <-  ds_time(l.wl[[i]])
# }
# 
# results1
# 
# results25 <- NULL
# for(i in 1:length(l.wlx)){
#   results25[[i]] <-  ds_time(l.wlx[[i]])
# }
# 
# results25



### add in ranks... 
rank_ds <- function(x){
  x %>% filter(time==max(time))  %>% mutate(rank = rank(-vals, ties.method = 'first')) %>%
    select(ids,rank) %>% full_join(x)
}



results.all2 <- results.all %>% map(rank_ds)

# check have all 4 ranks for each cage.
results.all2  %>% map(~table(.$rank))
results.all2  %>% map(~table(.$rank)) %>% map(~length(.)) %>% unlist



# results2 <- results1 %>% map(rank_ds)
# 
# results25x <- results25 %>% map(rank_ds) 

# but remember for two cages, animals are rank 3/4...

# 
# # Batch 4 Cage 4  2nd position  [32]  animal id = 2-4
# # Batch 2 Cage 3  2nd position  [27]  animal id = 1-4
# results2 %>% map(~table(.$rank)) #
# results2[[27]]$rank <- ifelse(results2[[27]]$rank==1, 1, results2[[27]]$rank +1)
# results2[[32]]$rank <- ifelse(results2[[32]]$rank==1, 1, results2[[32]]$rank +1)
# 
# results25x %>% map(~table(.$rank)) #25 hrs look okay. 



## add in batch/cage id, and if alpha/beta/gamma/delta cage...

# substrRight <- function(x, n){
#   substr(x, nchar(x)-n+1, nchar(x))
# }

results.all2 <- Map(cbind, results.all2, type = c("Alphas", "Betas", "Gammas", "Deltas")[as.numeric(substrRight(names(l.wl.all),1))])
results.all2 <- Map(cbind, results.all2, cage=names(l.wl.all))

results.all2


# results2 <- Map(cbind, results2, type = c("Alphas", "Betas", "Gammas", "Deltas")[as.numeric(substrRight(names(l.wl),1))])
# results2 <- Map(cbind, results2, cage=names(l.wl))
# 
# results2
# 
# results25x <- Map(cbind, results25x, type = c("Alphas", "Betas", "Gammas", "Deltas")[as.numeric(substrRight(names(l.wlx),1))])
# results25x <- Map(cbind, results25x, cage=names(l.wlx))
# 
# results25x



## Make one dataframe

# DF <- data.table::rbindlist(results2)
# DF25 <- data.table::rbindlist(results25x)


DFall <- data.table::rbindlist(results.all2)  #all data



## get CIs over time....

# DF %>% 
#   group_by(rank,time,type) %>%
#   summarise(median = median(vals),
#             lqr = quantile(vals,.25),
#             uqr = quantile(vals,.75)
#   ) -> DF.summary
# 
# 
# DF.summary

DFall %>% 
  group_by(rank,time,type) %>%
  summarise(median = median(vals),
            lqr = quantile(vals,.25),
            uqr = quantile(vals,.75)
  ) -> DFall.summary


DFall.summary


DFall.summary$type <- factor(DFall.summary$type, levels = c("Alphas","Betas","Gammas","Deltas"))

fig2A <- ggplot(DFall.summary, aes(x=time,y=median,color=factor(rank))) +
  geom_line() +
  facet_wrap(~type) +
  geom_ribbon(aes(ymin = lqr, ymax = uqr, fill = factor(rank)),alpha=.1, color = NA) +
  theme_classic() +
  ylab("Median David's Scores") +
  xlab("Contest") +
  scale_fill_manual(values = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"),
                    name= paste0("New",'\n', "Rank"), labels=c("1","2","3","4")) +
  scale_color_manual(values = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"),
                     name=paste0("New",'\n', "Rank"), labels=c("1","2","3","4")) +
   theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # axis.line = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        text = element_text(size = 20))
# +
# ggtitle("Reemergence of Ranks after Social Reorganization") 


ggsave("manuscript/behavior/results_figures/final_fig2A.png",fig2A,height =7, width =9, dpi=600)




### 70 mins - are DS sig different within each cage type (a/b/g/d) ?

ds.dfdf # 70 minute data


ds.dfdf.long <- ds.dfdf %>% pivot_longer(cols=1:4)


#rewrite for LSD


library(DescTools)
anovaDS <- function(xx){
# Perform one-way ANOVA with LSD posthocs
anova_results <- aov(formula = value ~ name, data = xx)
pairwise <- PostHocTest(anova_results, method="lsd")
return(pairwise)
}

split(ds.dfdf.long, ds.dfdf.long$group) %>%
  map(~anovaDS(.))



### Differences Across Cage Types in David's Scores

anovaDS2 <- function(xx){
  # Perform one-way ANOVA with LSD posthocs
  anova_results <- aov(formula = value ~ group, data = xx)
  pairwise <- PostHocTest(anova_results, method="lsd")
  return(pairwise)
}

split(ds.dfdf.long, ds.dfdf.long$name) %>%
  map(~anovaDS2(.))





### Latencies  (only include 70 minute data)


# l - 70 min data

l %>% 
  map(~ ungroup(.)) %>% 
  map(~ dplyr::select(., post_idbatchcage,ztime,starttime )) %>% 
  map(~ mutate(., ttime = ztime-starttime, contest = row_number())) -> ll


ll <- Map(cbind, ll, type = c("Alphas", "Betas", "Gammas", "Deltas")[as.numeric(substrRight(names(l),1))])
ll.df <- data.table::rbindlist(ll)

#summary stats for latencies of each contest
ll.df %>% group_by(type, contest) %>%
  summarise(median = median(ttime),
            lqr = quantile(ttime,.25),
            uqr = quantile(ttime,.75)
  ) -> ll.df.sum

range(ll.df.sum$contest)
ll.df.sum$type <- factor(ll.df.sum$type, c("Alphas","Betas", "Gammas", "Deltas"))

# Fig 2C
fig2B <- ggplot(ll.df.sum, aes(contest,median, color = type))+
  geom_step(size=1.25)+
  ylab("Average Latency (s)") +
  xlab("Behavioral Event") +
  labs(color = paste0("Previous", "\n", "Rank"))+
  scale_color_manual(values = viridis::viridis(6)) +
  scale_fill_manual(values = viridis::viridis(6))+
  theme_classic() +
  # xlim(c(1,60))+
#  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55))+
   scale_x_continuous(breaks = seq(0,60,5))+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        legend.position = c(.8,.3),
        text = element_text(size = 20)
  )
ggsave("manuscript/behavior/results_figures/final_fig2B_lat.png",fig2B,height =7, width =10, dpi=600)


### Stats for whether latencies are different between previous social ranks.

# get individual cage level data
ll.df %>% 
  separate(post_idbatchcage, sep=" ", into = c("b","batch","c","cageid")) %>%
  mutate(batch_cageid = paste(batch,cageid,sep="-")) ->latx
latx

#set up comparisons
latx$type2 <- factor(latx$type, c("Deltas", "Alphas", "Betas", "Gammas"))
latx$type3 <- factor(latx$type, c("Gammas", "Alphas", "Betas", "Deltas"))
latx$type4 <- factor(latx$type, c("Betas", "Alphas", "Deltas", "Gammas"))



# first fight
latxx1 <- latx %>% filter(contest == 1)
head(latxx1)

library(lme4)

mixed.lmer1 <- glmer(ttime ~ type + (1|batch), data = latxx1, family = Gamma(link ="log"))
summary(mixed.lmer1)

mixed.lmer1 <- glmer(ttime ~ type2 + (1|batch), data = latxx1, family = Gamma(link ="log"))
summary(mixed.lmer1)

mixed.lmer1 <- glmer(ttime ~ type3 + (1|batch), data = latxx1, family = Gamma(link ="log"))
summary(mixed.lmer1)

mixed.lmer1 <- glmer(ttime ~ type4 + (1|batch), data = latxx1, family = Gamma(link ="log"))
summary(mixed.lmer1)


### Check significant differences over N contests...

#alphas...

lat.alpha <- vector('list',max(latx$contest))

for(i in 1:max(latx$contest)){
  tryCatch({
    
mixed.lmer.out <- glmer(ttime ~ type + (1|batch), 
                     data = latx %>% filter(contest == i), 
                     family = Gamma(link ="log"))


lat.alpha[[i]] <- summary(mixed.lmer.out)}
, error = function(e) {
  # if an error occurs, store NULL in the list
  lat.alpha[[i]] <- NULL
})
}

lat.alpha[1:10] 
lat.alpha[11:20] 
#alpha/gamma NS at 7.
#alpha/delta NS at 9.



#betas...

lat.beta <- vector('list',max(latx$contest))

for(i in 1:max(latx$contest)){
  tryCatch({
    
    mixed.lmer.out <- glmer(ttime ~ type4 + (1|batch), 
                            data = latx %>% filter(contest == i), 
                            family = Gamma(link ="log"))
    
    
    lat.beta[[i]] <- summary(mixed.lmer.out)}
    , error = function(e) {
      # if an error occurs, store NULL in the list
      lat.beta[[i]] <- NULL
    })
}

lat.beta[1:10] 
lat.beta[11:20] 
#betas/gamma NS at 8.
#betas/delta NS at 9.



#deltas..

lat.delta <- vector('list',max(latx$contest))

for(i in 1:max(latx$contest)){
  tryCatch({
    
    mixed.lmer.out <- glmer(ttime ~ type2 + (1|batch), 
                            data = latx %>% filter(contest == i), 
                            family = Gamma(link ="log"))
    
    
    lat.delta[[i]] <- summary(mixed.lmer.out)}
    , error = function(e) {
      # if an error occurs, store NULL in the list
      lat.delta[[i]] <- NULL
    })
}

lat.delta[1:10] 
lat.delta[11:20] 
lat.delta[21:30] 
lat.delta[30:40] 









### Total Fights Per Cage in 70 minutes.

latx  #summary of 70 minute data

latx %>%
  group_by(batch,cageid,type,batch_cageid,type2,type3,type4) %>%
  filter(contest==max(contest)) %>%
  ungroup() %>%
  group_by(type) %>%
  summarise(median = median(contest),
            lqr = quantile(contest,.25),
            uqr = quantile(contest,.75)
  ) 

## Supplemental Figure X - total fights

latx %>%
  group_by(batch,cageid,type,batch_cageid,type2,type3,type4) %>%
  filter(contest==max(contest)) -> latx.fights

#make sure in correct order of levels
latx.fights$type <- factor(latx.fights$type, levels=c("Alphas","Betas","Gammas","Deltas"))

tf <- ggplot(latx.fights,aes(x=type,y=contest,color=type,fill= type))+
  geom_boxjitter(outlier.color = NA, 
                 jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.025,  
                 errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  ylab("Total fights in 70 minutes") +
  xlab("Previous Social Rank") +
  scale_color_manual(values = viridis::viridis(6)) +
  scale_fill_manual(values = viridis::viridis(6))+
  theme_classic() +
  ylim(0,65)+
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
ggsave("manuscript/behavior/results_figures/supp_totalfights.png",tf,height =5, width =5.5, dpi=600)


# Model

mixed.fight <- glmer(contest ~ type + (1|batch), 
                        data = latx.fights, 
                        family = poisson)

summary(mixed.fight)


mixed.fight4 <- glmer(contest ~ type4 + (1|batch), 
                     data = latx.fights, 
                     family = poisson)

summary(mixed.fight4)

###










# ggplot(DF.summary, aes(x=time,y=median,color=factor(rank))) +
#   geom_line() +
#   facet_wrap(~type) +
#   geom_ribbon(aes(ymin = lqr, ymax = uqr, fill = factor(rank)),alpha=.1, color = NA) +
#   theme_classic() +
#   ylab("Median David's Scores") +
#   xlab("Behavioral Event") +
#   scale_fill_manual(values = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"),
#                     name="Rank", labels=c("1","2","3","4")) +
#   scale_color_manual(values = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"),
#                      name="Rank", labels=c("1","2","3","4")) +
#   theme(axis.text.x = element_text(vjust = 1),
#         axis.text.y = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         panel.background = element_blank(),
#         plot.background = element_blank(),
#         axis.text = element_text(color="#3C3C3C", size=rel(0.8)),
#         strip.background = element_blank() 
#   ) +
#   ggtitle("Reemergence of Ranks after Social Reorganization")









## Make Fig 1C

#70 min
ds.long <- ds.df %>% as.data.frame %>% pivot_longer(1:4, names_to = "Rank")
ds.long$id <- rep(1:(nrow(ds.long)/4), each=4)

#25hr
ds.long25 <- ds.df25 %>% as.data.frame %>% pivot_longer(1:4, names_to = "Rank")
ds.long25$id <- rep(1:(nrow(ds.long25)/4), each=4)

ds.long$group <- "70min"
ds.long25$group <- "25hr"

ds.all.longX <- rbind(ds.long,ds.long25)

#combine
ds.all.longX %>% 
  group_by(Rank,group) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  ) -> ds.all.summary

ds.all.summary

ds.all.longX$group <- factor(ds.all.longX$group, levels = c("70min", "25hr"))
ds.all.summary$group <- factor(ds.all.summary$group, levels = c("70min", "25hr"))

# ds.all.long <- rbind(ds.long,ds.long25) %>% mutate(id = paste0(group,id))

# ds.all.long %>% group_by(Rank,time) %>% count()
# 
# ds.all.long$group <- factor(ds.all.long$group, levels = c("70min", "25hr"))
# ds.all.summary$group <- factor(ds.all.summary$group, levels = c("70min", "25hr"))

post_reg <- ggplot() + 
  geom_line(data=ds.all.longX, aes(x=Rank, y=value, group=id), alpha=.2, color="gray57") +
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
  

# arrange plot

 bottom <- grid::textGrob("Mouse Rank", gp = grid::gpar(fontsize = 20))
 left <- grid::textGrob("David's Score", gp = grid::gpar(fontsize = 20), rot = 90)
# 
 library(gridExtra)
 dds <- grid.arrange(pre,post_reg,pre_con, ncol=1,bottom = bottom, left =left)
 dds
 # 
# #ggsave("manuscript/behavior/results_figures/pre_ds_manuscript.png",dds,height =12, width =10, dpi=600)
