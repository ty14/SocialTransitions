
### What predicts hierarchy position post-reorganization ?


# do not need to run these, if already run all code in 01/02.

source("Behavior/01_pre_organization.R")
source("Behavior/02_post_reorganization.R")


## CAGE OF ORIGIN:

## need list of final ranks by pre batch cage.

l.ds70 # david scores of all 52 cages

l.hds <- lapply(l.ds70, function(x) data.frame(name = names(x), value = unname(x))) # get into dataframe
l.hds <- data.table::rbindlist(Map(cbind, l.hds, postbatchcage=names(l.ds70))) # add in info about post cage
l.hds$post_cage <- substr(l.hds$postbatchcage, nchar(l.hds$postbatchcage), nchar(l.hds$postbatchcage)) #post cage
l.hds$post_batch <- as.numeric(stringr::str_extract(l.hds$postbatchcage, "\\d+")) #post batch
l.hds$preid <- substr(l.hds$name, nchar(l.hds$name), nchar(l.hds$name)) #pre ID
l.hds$pre_cage <- as.numeric(stringr::str_extract(l.hds$name, "\\d+")) #pre cage

# add post rank
l.hds <- l.hds %>%
  group_by(postbatchcage) %>%
  mutate(post_rank = dense_rank(-value))

# function to get differences between ranks
get_difs <- function(X) {
p <- X %>% arrange(pre_cage) %>% .$post_rank
groups <- split(p, ceiling(seq_along(p)/4))
group.difs <- lapply(groups, function(x) as.numeric(dist(x))) # calculate difs over all 4 'cages'
return(list(total_sum = sum(unlist(group.difs)),
            total_zeros = sum(unlist(group.difs)==0)))
}


## Observered Results
dif_res <- lapply(split(l.hds, l.hds$post_batch), get_difs)
obs.total.difs <- sum(sapply(dif_res, function(x) x[[1]]))
obs.total.zeros <- sum(sapply(dif_res, function(x) x[[2]]))


####  Generate 10,000 samples of 13 batches to generate random distribution.
set.seed(1002)
nperms <- 130000
val <- NULL
for(i in 1:nperms){val[[i]] <- sample(rep(1:4,4))}
val

#function to get sum of differences and 0s
get_vals <- function(q){
  g <- split(q, ceiling(seq_along(q)/4))
  gd <- lapply(g, function(x) as.numeric(dist(x)))
  gd.sum <- sum(unlist(gd)) # total sum
  gd.0 <- sum(unlist(gd)==0) # total zeros
  return(list(gd.sum,gd.0))
}

results <- lapply(val, get_vals)
gd.sums <- unlist(lapply(results, `[[`, 1))
gd.0s <- unlist(lapply(results, `[[`, 2))


dist.sums<-
unlist(
  lapply(
  split(gd.sums, rep(seq_along(gd.sums), each = 13, length.out = length(gd.sums))), 
  sum)
)

hist(dist.sums)
range(dist.sums)
obs.total.difs


dist.zeros<-
  unlist(
    lapply(
      split(gd.0s, rep(seq_along(gd.0s), each = 13, length.out = length(gd.0s))), 
      sum)
  )
hist(dist.zeros)
range(dist.zeros)
obs.total.zeros


## Make Supplemental Figure: ggplot versions of figures:
dist.sums <- as.data.frame(dist.sums)

p1 <- ggplot() +
  geom_histogram(data=dist.sums, aes(x=dist.sums), color='#009933',fill="#ccffdd", binwidth = 4) +
  geom_vline(xintercept = obs.total.difs, lty=2, color='#000066') +
  theme_classic() +
  xlab("Summed Differences in Ranks from Random Samples") +
  ylab("Frequency")+
  theme(text = element_text(size = 15))

# ggsave("manuscript/behavior/results_figures/SumRankDiff_supp6B.png",p1,height =5, width =6, dpi=600)


## BODY WEIGHT / PREVIOUS ds:

#do alpha,beta,gamma,delta separately
l.hds$pre_idbatchcage <-  paste0("Batch",l.hds$post_batch,"Cage",l.hds$pre_cage)
l.hds<-as.data.frame(l.hds)
l.hds$preid<-as.numeric(l.hds$preid)
post.df <- left_join(l.hds,bw)
l.post <- split(post.df, post.df$Prerank) #put into list by rank

#graph
dw <- post.df %>% select(name, wt_d4,wt_d8,Prerank)
dw2 <- dw %>% pivot_longer(cols = 2:3, names_to = "wt")
dw2$wt <- ifelse(dw2$wt == "wt_d4", "Day 4", "Day 8")
dw2$Prerank <- as.factor(dw2$Prerank)
dwp <- ggplot(dw2,aes(x=Prerank,y=value,color=Prerank,fill= Prerank))+
  geom_boxjitter(outlier.color = NA, 
                 jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.025,  
                 errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  ylab("Body Mass (g)") +
  xlab("Rank") +
  facet_wrap(~wt)+
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

dwp
# ggsave("manuscript/behavior/results_figures/supp_weight.png",dwp,height =5, width =8, dpi=600)



# Models
library(lme4)
library(lmerTest)

# Alphas
model1 <- lmer(value ~ wt_d8 + ds + (1 | postbatchcage), data = l.post[[1]])
summary(model1)
anova(model1)

# Betas
model2 <- lmer(value ~ wt_d8 + ds + (1 | postbatchcage), data = l.post[[2]])
summary(model2)
anova(model2)

# Gammas
model3 <- lmer(value ~ wt_d8 + ds + (1 | postbatchcage), data = l.post[[3]])
summary(model3)
anova(model3)

# Deltas
model4 <- lmer(value ~ wt_d8 + ds + (1 | postbatchcage), data = l.post[[4]])
summary(model4)
anova(model4)

## day 12

# Alphas
model1 <- lmer(value ~ wt12 + ds + (1 | postbatchcage), data = l.post[[1]])
summary(model1)
anova(model1)

# Betas
model2 <- lmer(value ~ wt12 + ds + (1 | postbatchcage), data = l.post[[2]])
summary(model2)
anova(model2)

# Gammas
model3 <- lmer(value ~ wt12 + ds + (1 | postbatchcage), data = l.post[[3]])
summary(model3)
anova(model3)

# Deltas
model4 <- lmer(value ~ wt12 + ds + (1 | postbatchcage), data = l.post[[4]])
summary(model4)
anova(model4)


### total aggression

# import start end times of pre-reorganization behavior.

sends <- read.csv('Behavior/sends.csv')

sends <- sends %>% split(.$batch)


## Function to get differences in starts/ends
get_time_dif<-function(x){
  # Convert timestamp to proper date-time format
  timestamp <- as.POSIXct(x$Timestamp, format="%m/%d/%Y %H:%M:%S")
  # Calculate number of seconds since midnight
  seconds_since_midnight <- difftime(timestamp, trunc(timestamp, "day"), units="secs")
  # Compute differences between consecutive elements
  diffs <- diff(seconds_since_midnight)
  # Extract every second difference starting from the first
  even_diffs <- diffs[seq(1, length(diffs), by=2)]
  return(even_diffs)
}


# from pre_startend_times
obs.secs <- lapply(lapply(sends,get_time_dif),function(x) as.numeric(sum(x)))


# get preorganization behavior data
df1 <- read_csv("Behavior/Pre_WL.csv")

# create list of dataframes for each individual cage
l1 <- split(df1, df1$pre_batchcage)

# just keep winner loser columns
l1.wl <- l1 %>% map(~ select(., winner,loser) %>% as.data.frame) 

# convert list to win-loss matrices organized by David's Scores
l1.mat <- l1.wl %>% map(~ compete::get_wl_matrix(.)) %>% map(~ compete::org_matrix(., method="ds"))

# melt matrices.

library(reshape2)


# Define function to melt a matrix and calculate sum of row values and column values
melt_matrix <- function(mat) {
  mat_melted <- melt(mat)   # Melt matrix to long format
  row_sum <- aggregate(value ~ winner, data = mat_melted, sum)   # Sum row values
  col_sum <- aggregate(value ~ loser, data = mat_melted, sum)   # Sum column values
  sum_df <- merge(row_sum, col_sum, by.x = "winner", by.y = "loser", all = TRUE)   # Combine row sum and column sum
  names(sum_df) <- c("preid", "given", "received")  # Rename columns
  return(sum_df)
}

l.aggr <- lapply(l1.mat, melt_matrix)
l.aggr.df <- data.table::rbindlist(Map(cbind, l.aggr, name = names(l.aggr)))
l.aggr.df$batch <- stringr::str_extract(l.aggr.df$name, "\\d+")
l.aggr.df$pre_cage <- stringr::str_extract(l.aggr.df$name, "\\d+$")
l.aggr.df <- l.aggr.df %>% left_join(
data.frame(time=unlist(obs.secs), batch = stringr::str_extract(names(obs.secs), "\\d+$"))
)
l.aggr.df$pre_idbatchcage <- gsub("\\s+", "", l.aggr.df$name)

# generate per hour value 
l.aggr.df$given1 <- 3600*(l.aggr.df$given / l.aggr.df$time)
l.aggr.df$received1 <- 3600*(l.aggr.df$received / l.aggr.df$time)

l.aggr.df$batch <- NULL
l.aggr.df$name <- NULL
l.aggr.df$pre_cage <- as.numeric(l.aggr.df$pre_cage)


#join per hour to pre_org data predicting post aggr  [value = post DS;  ds = pre DS]
post.df <- post.df %>% left_join(l.aggr.df)

l.post1<- split(post.df, post.df$Prerank) #put into list by rank


# Alphas
model1 <- lmer(value ~ given1 + received1 + (1 | postbatchcage), data = l.post1[[1]])
summary(model1)
anova(model1)

# Betas
model2 <- lmer(value ~ given1 + received1 + (1 | postbatchcage), data = l.post1[[2]])
summary(model2)
anova(model2)

model2 <- lmer(value ~ received1 + (1 | postbatchcage), data = l.post1[[2]])
summary(model2)
anova(model2)

# Gammas
model3 <- lmer(value ~ given1 + received1 + (1 | postbatchcage), data = l.post1[[3]])
summary(model3)
anova(model3)

model3 <- lmer(value ~ received1 + (1 | postbatchcage), data = l.post1[[3]])
summary(model3)
anova(model3)


# Deltas
model4 <- lmer(value ~ given1 + received1 + (1 | postbatchcage), data = l.post1[[4]])
summary(model4)
anova(model4)

model4 <- lmer(value ~  received1 + (1 | postbatchcage), data = l.post1[[4]])
summary(model4)
anova(model4)

pp1<-ggplot(post.df, aes(x=given1, y=value)) + 
  geom_point() +
  facet_wrap(~Prerank, nrow=1) +
  xlab("Aggression Given (per hour)") +
  ylab("Post Reorganization David's Score") +
  stat_smooth(method='lm',se=F)+theme(text = element_text(size = 12))
  
pp2<-ggplot(post.df, aes(x=received1, y=value)) + 
  geom_point() +
  facet_wrap(~Prerank, nrow=1) +
  xlab("Aggression Received (per hour)") +
  ylab("Post Reorganization David's Score") +
  stat_smooth(method='lm',se=F)+theme(text = element_text(size = 12))

library(gridExtra)
p2 <- grid.arrange(pp1,pp2,nrow=2)
# ggsave("manuscript/behavior/results_figures/supp_aggpredict_sup6A.png",p2,height =6, width =12, dpi=600)


### total aggression in cage predict final hierarchy position ?

post.df %>%
  group_by(pre_idbatchcage) %>%
  mutate(totalaggr = (sum(given)+sum(received))/2) %>%
  mutate(rate_totalaggr = 3600*(totalaggr/time)) %>%
  data.frame() -> post.df

l.post2<- split(post.df, post.df$Prerank) #put into list by rank


#plots
p3 <- ggplot(post.df, aes(x=rate_totalaggr, y=value)) + 
  geom_point() +
  facet_wrap(~Prerank, nrow=1) +
  xlab("Aggression per Cage (per hour)") +
  ylab("Post Reorganization David's Score") +
  stat_smooth(method='lm',se=F) + theme(text = element_text(size = 12))
# ggsave("manuscript/behavior/results_figures/AggCage_supp6C.png",p3,height =3, width =12, dpi=600)


# Alphas
model1 <- lmer(value ~ rate_totalaggr + (1 | postbatchcage), data = l.post2[[1]])
summary(model1)
anova(model1)

# Betas

model2 <- lmer(value ~ rate_totalaggr + (1 | postbatchcage), data = l.post2[[2]])
summary(model2)
anova(model2)


# Gammas
model3 <- lmer(value ~ rate_totalaggr + (1 | postbatchcage), data = l.post2[[3]])
summary(model3)
anova(model3)

# Deltas
model4 <- lmer(value ~ rate_totalaggr + (1 | postbatchcage), data = l.post2[[4]])
summary(model4)
anova(model4)

