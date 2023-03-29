#####################################################
# ~ predPharma: parameters from middur ~
#
# reproduces parameter calculations from SCRAP.m
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# rationale is: predPharma Shiny app will be difficult to use if user has to import .mat file from Jason's MATLAB scripts
# best is for user to import RAWs.csv so it is more coherent with FramebyFrame
# then we will convert to middur and calculate parameters internally

# no need to re-write all of SCRAP.m
# go directly to calculations of summary parameters
# will then compare carefully we have the same values


# function diff2(...) -----------------------------------------------------
# I want to 'correct' R's diff function
# e.g. diff(c(FALSE, FALSE, TRUE, TRUE)) gives 0 1 0
# so if look for first TRUE using diff, will give #2
# when actually it is #3
# I think it should return NA 0 1 0, then positions will match original vector
diff2 <- function(vector) {
  
  return(c(NA, diff(vector)))
  
}



# one_fish functions ------------------------------------------------------

# one function per parameter
# each function works on middur data for one fish
# then we will loop through columns

# mc is one column of middur datapoints

# about sleep detection: SCRAP.m uses a threshold at middur 0.1
# if >, not asleep
# if <=, asleep

#### averageActivity ####
averageActivity_onefish <- function(mc) {
  
  return( mean(mc) )
  
}


#### averageWaking ####
# definition is average activity when not asleep
# for some reason, when detecting sleep, uses threshold at 0.1 middur
# but when excluding "asleep" minutes for averageWaking, it uses threshold at 0
averageWaking_onefish <- function(mc) {
  
  # only keep active datapoints
  mca <- mc[which(mc>0)]
  
  # calculate & return mean
  mean(mca)
  
}


#### sleep ####
# total time asleep in minutes
# which here is simply counting number of values <= 0.1
sleep_onefish <- function(mc) {
  
  return( length(which(mc <= 0.1)) )
  
}


#### sleepBout ####
# number of sleep bouts
# here, we will count number of sleep bout starts
sleepBout_onefish <- function(mc) {
  
  # first convert to TRUE = asleep frame & FALSE = not asleep frame
  mcb <- (mc <= 0.1) # b for booleans
  
  # then convert to transitions, which we can do by diff the booleans
  # FALSE is 0, TRUE is 1
  # diff works as (n+1) - n
  # e.g. FALSE, FALSE, TRUE, TRUE, FALSE
  # becomes 0, 1, 0, -1
  # so any 1 marks the start of a sleep bout
  # (as it marks FALSE followed by TRUE)
  mct <- diff2(mcb) # t for transitions
  
  # now count the number of 1s
  # this is the number of sleep bouts that started during that window
  return( length(which(mct==1)) )
  
}


#### sleepLatency ####
# how many minutes before the first sleep bout?
# so this is simply the row index of the first datapoint <= 0.1
# ! SCRAP.m actually reports last not asleep value, not exactly the index of the first sleep bout
# so add minus 1 here
# for some reason there is another minus 1, I do not know where from, but just will copy it here
# (it is consistent for 4 larvae I tried, every time SCRAP.m finds same - 1)
sleepLatency_onefish <- function(mc) {
  
  return( which(mc < 0.1)[1] -1 -1)
  
}


#### sleepLength ####
# average sleep bout length, in minutes
# starts with similar procedure as sleepBout
# then get indices of start and stop for each sleep bout
sleepLength_onefish <- function(mc) {
  
  # first convert to TRUE = asleep frame & FALSE = not asleep frame
  mcb <- (mc <= 0.1) # b for booleans
  
  # then convert to transitions, which we can do by diff the booleans
  # FALSE is 0, TRUE is 1
  # diff works as (n+1) - n
  # e.g. FALSE, FALSE, TRUE, TRUE, FALSE
  # becomes 0, 1, 0, -1
  # so any 1 marks the start of a sleep bout
  # (as it marks FALSE followed by TRUE)
  mct <- diff2(mcb) # t for transitions
  
  # make sure we only have complete sleep bouts
  # by complete, I mean which starts and stops during that window
  
  # start the window at the first sleep bout START
  # stop the window at the last sleep bout STOP
  mctc <- mct[ which(mct==1)[1] : which(mct==-1)[length(which(mct==-1))] ] # c for cropped
  
  # record all the START, i.e. the minute index at which each sleep bout started
  sta <- which(mctc==1)
  
  # record all the STOP, i.e. the minute index at which each sleep bout stopped
  sto <- which(mctc==-1)
  
  # for each START, there should be a STOP
  if(length(sta) != length(sto)) stop('\t \t \t \t >>> Error sleepLength_onefish: different number of sleep bout STARTs and sleep bout STOPs \n')
  
  # now we can simply do the delta STOP - START to find how long each sleep bout lasted
  # and we return the mean of those durations
  return( mean(sto - sta) )
  
}