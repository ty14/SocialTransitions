### Start/End Times for Pre-Aggression

vals <-c("Start","End")


## Import PRE Files:

# Get the full path to the data_raw subfolder
data_raw_path <- file.path(getwd(), "data_raw")

# Get a list of all csv files in the data_raw subfolder
files <- list.files(path = data_raw_path, pattern = "Pre.*\\.csv", full.names = TRUE)

# Initialize empty list to store csv data
csv_list <- list()

# Loop over csv files and import them
for (file in files) {
  # Import csv data
  csv_data <- read.csv(file)
  
  # Add csv data to list
  csv_list[[file]] <- csv_data
}

## List of files
sends <- lapply(csv_list, function(x) x[x$Start.End %in% vals,c(1,8)])

names(sends) <- stringr::str_extract(names(sends), "Batch\\d+")


sends

## Fixing Start/Ends

sends[[1]]<-rbind(sends[[1]][1:10,], 
         data.frame(Timestamp="9/30/2019 15:30:12",Start.End="Start"),
         sends[[1]][11:19,])

sends[[2]]

sends[[3]]

sends[[4]]<-rbind(sends[[4]][1:16,], 
                  data.frame(Timestamp="1/15/2020 10:53:21",Start.End="Start"),
                  sends[[4]][17:19,])


sends[[5]]<-rbind(sends[[5]][1:14,], 
                  data.frame(Timestamp="1/15/2020 15:58:58",Start.End="Start"),
                  sends[[5]][15:17,])

sends[[6]]<-rbind(sends[[6]][1:4,], 
                  data.frame(Timestamp="1/26/2020 13:23:16",Start.End="Start"),
                  sends[[6]][6:18,])


sends[[7]]<-rbind(sends[[7]][1:5,], 
                  data.frame(Timestamp="10/19/2019 12:28:10",Start.End="End"),
                  sends[[7]][6:21,],
                  data.frame(Timestamp="10/23/2019 12:16:42",Start.End="Start"),
                  sends[[7]][23,]
)


sends[[8]]<-rbind(sends[[8]][1:13,], 
                  data.frame(Timestamp="10/29/2019 14:20:12",Start.End="End"),
                  sends[[8]][14:15,])


sends[[9]]<-rbind(sends[[9]][1:10,], 
                  data.frame(Timestamp="11/4/2019 12:47:52",Start.End="Start"),
                  sends[[9]][11:21,])


sends[[10]]<-rbind(sends[[10]][1:13,], 
                  data.frame(Timestamp="11/11/2019 15:36:45",Start.End="End"),
                  sends[[10]][14:21,]) #has erroneous end 2 days after end of observations. Removed.


sends[[11]]<-rbind(sends[[11]][1:18,], 
                  data.frame(Timestamp="11/20/2019 12:30:21",Start.End="Start"),
                  sends[[11]][19:21,])


sends[[12]]<-rbind(sends[[12]][1:2,], 
                   data.frame(Timestamp="11/23/2019 12:03:41",Start.End="Start"),
                   sends[[12]][4:8,],
                   data.frame(Timestamp="11/25/2019 13:15:21",Start.End="Start"),
                   sends[[12]][10:14,],
                   data.frame(Timestamp="11/26/2019 11:59:37",Start.End="Start"),
                   sends[[12]][16:20,]
)


sends[[13]]<-rbind(sends[[13]][1:4,], 
                   sends[[13]][6:11,],
                   sends[[13]][13:18,],
                   sends[[13]][20:22,],
                   data.frame(Timestamp="12/4/2019 11:24:03",Start.End="End"),
                   sends[[13]][23,],
                   sends[[13]][25:27,]
)


sends[[14]]<-rbind(sends[[14]][1:2,], 
                   data.frame(Timestamp="12/1/2019 14:31:04",Start.End="Start"),
                   sends[[14]][5:8,],
                   data.frame(Timestamp="12/2/2019 13:02:38",Start.End="End"),
                   sends[[14]][9:15,],
                   data.frame(Timestamp="12/4/2019 12:31:33",Start.End="End"),
                   sends[[14]][16,],
                   data.frame(Timestamp="12/4/2019 15:02:45",Start.End="End"),
                   sends[[14]][17:20,]
)


#
#
sends


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

#Apply to all batches
lapply(sends,get_time_dif)

lapply(lapply(sends,get_time_dif),function(x) as.numeric(sum(x)))

ts <- lapply(lapply(sends,get_time_dif),function(x) as.numeric(sum(x))) 
tsx <- do.call(rbind, ts) %>% as.data.frame()
tsx$V1 <- as.numeric(tsx$V1)
mean(tsx$V1) #36720.21 secs 
sec <- 36720.21/10
min <- sec/60
min/60
