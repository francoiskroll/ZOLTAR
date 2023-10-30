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


# packages ----------------------------------------------------------------

library(dplyr)
library(tibble)
library(tidyr)

inaThr <- 0.1 # this is a threshold which was added by Jason
# below middur 0.1 (i.e. active for 0.1 second or less within that minute), minute is considered inactive



# function legacyFingerprintMid(...) --------------------------------------

# function legacyFingerprint was calculating fingerprint from .mat file
# functions for this are in legacyFingerprint.R

# now, takes _middur.csv as input, hence 'Mid'

# this will essentially run the entire pipeline of paramsFromMid
# i.e. takes middur file, calculate parameters and assign genotypes, calculate fingerprint

# note, for a given parameter, I think a better approach to calculating Z-scores is to calculate Z-score for each larva separately
# then do mean +/- SEM of those Z-scores, which allows to have an error bar
# I checked in Clustering 2.m: Z-scores were not calculated this way
# so I am assuming Z-scores in the drug database were also calculated like this
# so will stick to this approach here

# note, original fingerprint also have averageWaking daymean & averageWaking nightmean
# averageWaking daymean is the mean of averageWaking day0 (incomplete), day1, day2, day3 (incomplete)
# averageWaking nightmean is the mean of averageWaking night0, night1, night2
# I do not think these parameters are a good idea because they depend on when the experiment is started/stopped (for day)
# I also do not know what they add that is not captured by averageWaking day1/day2 and night1/night2
# Also, why averageWaking only? Why not the same big average for sleep etc?
# I decided not to calculate them here
# will need to delete them from the drug fingerprint too

# for all the other parameters, I get almost exactly the same Z-scores as SCRAP.m, checked for experiment 220531_14
# I also compared fingerprint plots (gglegacyFingerprint from middur or from mat), I cannot see a slightest point moving
# (except from the day/nightmean_averageWaking disappearing)
legacyFingerprintMid <- function(mid,
                                 genopath,
                                 treGrp,
                                 conGrp,
                                 nights=c('night1','night2'),
                                 days=c('day1','day2'),
                                 suns=NA) {
  
  ## calculate parameters
  # we get a list, where each slot is one parameter table
  paral <- calculateParameters(mid=mid, genopath=genopath, suns=suns)
  
  # check that treGrp & conGrp are correct
  if(! treGrp %in% unique(paral[[1]]$grp)) stop('\t \t \t \t Error legacyFingerprintMid: treGrp', treGrp, 'was not found in data. \n')
  if(! conGrp %in% unique(paral[[1]]$grp)) stop('\t \t \t \t Error legacyFingerprintMid: conGrp', conGrp, 'was not found in data. \n')
  
  ## for each parameter, calculate
  # mean of controls
  # sd of controls
  # mean of treated
  # Z-score
  # see note above, would have been better to calculate individual Z-score for each larva then mean ?? SEM, but copying SCRAP.m
  zsl <- lapply(paral, function(pa) {
    
    # loop through the windows we need to calculate
    sapply(c(nights, days), function(win) {
      # win is e.g. night1 or day2
      
      # for each win,
      # calculate the mean of control datapoints
      meanCon <- mean( pa[(which(pa$grp==conGrp)), win] , na.rm=TRUE)
      
      # calculate the sd of control datapoints
      sdCon <- sd( pa[(which(pa$grp==conGrp)), win] , na.rm=TRUE)
      
      # calculate the mean of treatment datapoints
      meanTre <- mean( pa[(which(pa$grp==treGrp)), win] , na.rm=TRUE)
      
      # calculate Z-score
      zsco <- (meanTre - meanCon) / sdCon
      
    })
  })
  # now we have a list zsl
  # each slot is a small row: Z-score for night1, Z-score for night2, Z-score for day1, Z-score for day2
  # turn this into a clean dataframe
  zs <- as.data.frame(do.call(rbind, zsl))
  zs <- cbind(parameter=row.names(zs), zs)
  row.names(zs) <- NULL
  
  # pivot longer so that we have just one Z-score column
  
  # we will call the Z-score column YYMMDD_BX_treGrp
  # get the YYMMDD_BX from the genopath
  print(basename(genopath))
  ybt <- paste(substr(basename(genopath), 1, 9), treGrp, sep='_')
  
  zs <- zs %>%
    pivot_longer(-parameter,
                 names_to='win',
                 values_to=ybt)
  
  ### here we need to add parameters total day waking and total night waking
  # (I do not know why this parameter is there, I do not think that doing a quick operation on two existing parameters counts as a new parameter)
  # (but it is in the original fingerprint, so here it is)
  # from what I understand, it is, for each fish, the sum day1 waking activity + day2 waking activity
  # we did not calculate it in calculateParameters because it is not defined for *one* day or night
  # EDIT: decided to delete those parameters, see comments for rationale
  
  # now just add column uparam
  zs <- zs %>%
    mutate(uparam=paste(win, parameter, sep='_'), .before=1)
  zs <- as.data.frame(zs)
  
  # and we essentially have the same format as the original legacyFingerprint (from .mat file)
  return(zs)
}


# function calculateParameters(...) ---------------------------------------

# calculates parameter for each well and each time window

calculateParameters <- function(mid,
                                genopath,
                                suns) {
  
  # parameters are
  allparameters <- c('sleep',
                     'sleepBout',
                     'sleepLength',
                     'sleepLatency',
                     'averageActivity',
                     'averageWaking')
  
  
  ### split by day/night
  # if suns is not given, we do the automatic way
  if(is.na(suns[1])) {
    # NOTE, as of 30/10/2023
    # ZOLTAR never calls splitMidbyDayNight
    # it calls splitMidbyWoi with defaults set in app.R
    # (default settings will give the same output)
    dn <- splitMidbyDayNight(mid)
  } else {
    dn <- splitMidbyWoi(mid,
                        suns)
  }
  # then the rest can proceed as normal either way
  
  # how many time columns do we have?
  timecols <- which(grepl("^f+[[:digit:]]", colnames(dn[[1]])))[1] - 1
  # above finds first column which is called fX (usually f1 or f97)
  
  ### for each parameter, for each window, for each well, calculate the parameter
  
  # loop through parameters
  paral <- lapply(allparameters, function(param) {
    
    # which is the function we should be using for this parameter?
    # it will be called param_onefish, e.g. sleep_onefish
    fun2use <- paste0(param, '_onefish')
    # we now have the name, here is to link it to the actual function
    onefishFun <- match.fun(fun2use)
    
    # loop through windows
    pal <- lapply(1:length(dn), function(win) {
      
      # loop through wells
      sapply( (timecols+1):ncol(dn[[win]]), function(w) { # w is for well
        
        # make sure the middur values are just a plain vector
        mc <- as.numeric(dn[[win]][,w])
        # run the _onefish function on it
        onefishFun(mc=mc)
        
      } )
      
    })
    # here, after lapply window, we have a list for that parameter
    # where each slot is parameter datapoints (typically 96 of them) for that window
    # create a cleaner dataframe
    # put day/night as columns, i.e. night0 / day1 / ...
    pa <- as.data.frame(do.call(cbind, pal))
    # add day/night as column names
    colnames(pa) <- names(dn)
    # add parameter name
    # add fish number
    pa <- pa %>%
      mutate(fish=colnames(dn[[1]])[(timecols+1):ncol(dn[[1]])], .before=1) %>%
      mutate(parameter=param, .before=1)
    
    ### import genotype file and add grp column
    # import genotype file
    geno <- importGenotype(genopath=genopath)
    # make sure every well is mentioned in genotype file, adding an 'excluded' column if necessary
    # need to give box number here (1 or 2)
    # can guess from the fish IDs
    # assume boxnum=1, unless first fish ID is f97
    boxnum <- 1
    if(pa$fish[1]=='f97') {
      boxnum <- 2
    }
    geno <- fillGenotype(geno=geno, boxnum=boxnum)
    # add genotype column, using function assignGenotype
    pa <- pa %>%
      mutate(grp=assignGenotype(fs=pa$fish, geno=geno), .after='fish')
  })
  
  # here, after lapply parameters, we have a list (paral) where each slot is a parameter table
  # each parameter table, is columns = day/night; rows = fish
  # can add names to this list
  names(paral) <- allparameters
  
  # return this list
  return(paral)
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
  
  return( length(which(mc <= inaThr)) )
  
}
# therefore, units is minutes
# (not hours like in FramebyFrame)


#### sleepBout ####
# number of sleep bouts
# here, we will count number of sleep bout starts
sleepBout_onefish <- function(mc) {
  
  # first convert to TRUE = asleep frame & FALSE = not asleep frame
  mcb <- (mc <= inaThr) # b for booleans
  
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

# Note, I think it would be more accurate to report the first sleep bout, *after* the fish was active once
# so that if the fish is inactive during the transition, we do not report 0
# but this is not how SCRAP.m is coded, and goal here is to replicate exactly
sleepLatency_onefish <- function(mc) {

  lat <- which(mc <= inaThr)[1] -1 -1
  
  # exception: if first datapoint is inactive, it seems that SCRAP.m returns 1
  # in fact, it seems it never returns 0, always 1
  # copy this
  if(!is.na(lat) & lat<=0) { # lat may be NA if fish was never asleep during the window (we return NA in that case)
    lat <- 1
  }
  
  # exception: if fish was never asleep, SCRAP.m returns the last minute of the window
  # this is objectively inaccurate (the fish never slept, so by definition we cannot report a value)
  # but well, can only copy what it does here
  if(is.na(lat)) {
    lat <- length(mc)
  }
  
  return( lat )
  
}


#### sleepLength ####
# average sleep bout length, in minutes
# starts with similar procedure as sleepBout
# then get indices of start and stop for each sleep bout
sleepLength_onefish <- function(mc) {
  
  # first convert to TRUE = asleep frame & FALSE = not asleep frame
  mcb <- (mc <= inaThr) # b for booleans
  
  # then convert to transitions, which we can do by diff the booleans
  # FALSE is 0, TRUE is 1
  # diff works as (n+1) - n
  # e.g. FALSE, FALSE, TRUE, TRUE, FALSE
  # becomes 0, 1, 0, -1
  # so any 1 marks the start of a sleep bout
  # (as it marks FALSE followed by TRUE)
  mct <- diff2(mcb) # t for transitions
  
  # ! if no sleep transition whatsoever (can happen in empty wells)
  # just return NA
  if( length(which(mct==1))==0 | length(which(mct==-1))==0 ) {
    return(NA)
  }
  
  # make sure we only have complete sleep bouts
  # by complete, I mean which starts and stops during that window
  
  # start the window at the first sleep bout START
  # stop the window at the last sleep bout STOP
  mctc <- mct[ which(mct==1)[1] : which(mct==-1)[length(which(mct==-1))] ] # c for cropped
  
  # record all the START, i.e. the minute index at which each sleep bout started
  sta <- which(mctc==1)
  
  # record all the STOP, i.e. the minute index at which each sleep bout stopped
  sto <- which(mctc==-1)
  
  # we should have at least one full sleep bout
  # otherwise there is no point calculating the parameter
  if(length(sta)==0 | length(sto)==0) {
    return(NA)
  }
  
  # for each START, there should be a STOP
  if(length(sta) != length(sto)) stop('\t \t \t \t >>> Error sleepLength_onefish: different number of sleep bout STARTs and sleep bout STOPs \n')
  
  # now we can simply do the delta STOP - START to find how long each sleep bout lasted
  # and we return the mean of those durations
  return( mean(sto - sta) )
  
}


# function splitMidbyDayNight(...) ----------------------------------------

# NOTE, as of 30/10/2023
# this function is not called by the ZOLTAR app anymore
# it calls splitMidbyWoi with default as
# 24, 38
# 38, 48
# 48, 62
# 62, 72
# which gives the same output
# will leave it as other code uses it

# this function splits middur data as a big table, where
# columns are fish (after a few time columns)
# rows are minutes
# into a list where each element is middur data for one day or night
# analogous to splitFramesbyDayNight in FramebyFrame package

# here, we can rely on zth column (number of hours since 9AM day0)
# assuming dayduration is 14hrs
# transitions will happen at 14, 24, 38, etc.

# will assume experiment starts during a day
splitMidbyDayNight <- function(mid) {
  
  dayduration <- 14
  
  # sunsets and sunrises times (in number of hours since day0 sunrise at 9AM) for 10 days
  suns <- c(rbind(seq(dayduration, 720, 24), seq(24, 720, 24)))
  # typically dayduration = 14
  # so will go 14, 24, 38, 48 etc.
  
  # how many transitions do actually have in the data?
  # only keep transitions that happened within the experiment
  # (or in other words, that happened before the end of the middur data)
  suns <- suns[ which(suns <= max(mid$zhrs)) ]
  
  # for each light transition in zth (suns), detect the closest minute
  # we want the minute just before the transition
  # e.g. for 24 hrs, we want ~ 23.99
  # tramin = transition minutes
  tramin <- sapply(suns, function(sun) {
    findLightTransitionMinute(Zeitgeberdurations=mid$zhrs,
                              transitionHour=sun)
  })
  
  # preallocate list dn for day/night
  # each slot will be middur data for one day or night
  dn <- vector(mode='list', length=length(tramin)-1)
  # ! this will skip any day or night that is not full
  # put appropriate names
  # will assume started experiment during the day (i.e. between 9AM and 11PM)
  # so we can simply go day0, night0, day1, night1, etc.
  dn_nms <- c(rbind(sprintf('night%i', 0:30), sprintf('day%i', 1:31))) # names for a one-month long experiment
  # ! assumes starts during day; finishes during day
  # and therefore that day at the start is incomplete and day at the end is incomplete
  # dn_nms above goes night0, day1, night1, day2, etc.; so easy to exclude night0
  # typically gives day1 / day2 + night1 / night2 for the standard Rihel lab experimental design
  
  # take the first few names according to how many we need
  names(dn) <- dn_nms[1:length(dn)]
  
  # fill the list
  # for each window, we take from row just after transition until row just before transition
  # e.g. night0 will be something like: 14.01 hr until 23.99 hr
  # as we are dealing with one-minute binned data, there are not too many rows, so convert data.table to dataframe which I find easier
  for(w in 1:length(dn)) {
    dn[[w]] <- as.data.frame(mid[ (tramin[w] + 1) : (tramin[w+1]) , ])
  }
  
  # return the list
  return(dn)
  
}

#### function findLightTransitionMinute ####
# small function to help splitMidbyDayNight above
# copied from FramebyFrame package findLightTransitionFrame

# Zeitgeberdurations = Zeitgeber durations (i.e. number of hours since day0 9AM), typically a column of middur data
# transitionHour = transition to look for, usually e.g. 14 or 24 or 38 or 48, etc.

# returns the minute (index) just before the light transition
# e.g. if looking for first sunset = 14 zth hours
# may find frame zth = 13.99999 (Case1)
# in that case it returns that frame
# or may find frame zth = 14.00001 or frame with exactly zth = 14 (Case2)
# in that case, assumes this frame is the frame just after the transition and it returns the frame just before
findLightTransitionMinute <- function(Zeitgeberdurations,
                                      transitionHour) {
  
  tfra <- which.min(abs(Zeitgeberdurations - transitionHour)) # transition frame
  
  if (Zeitgeberdurations[tfra] >= transitionHour) { # this is Case2
    
    tfra <- tfra - 1
    
    # check we are now before the transition
    if (Zeitgeberdurations[tfra] >= transitionHour)
      stop('\t \t \t \t >>> Something unexpected: the frame preceding the closest frame to the transition is still after the transition \n')
    
    return(tfra)
    
  } else {
    return(tfra)
  }
  
}



# function splitMidbyWoi(...) ------------------------------------------

# alternative to splitMidbyDayNight is to split by custom times given by the user
# ! woi should still represent, as best as possible, days and nights as we are matching to Rihel et al. 2010.

# this function splits middur data as a big table, where
# columns are fish (after a few time columns)
# rows are minutes
# into a list where each element is middur data for one woi (day/night defined by the user)

# will get a vector start1, end1, start2, end2, etc.
# will assume alternating day/night/...

# suns are sunrises/sunsets given by user

splitMidbyWoi <- function(mid,
                          suns) {
  
  # for each light transition in zth (suns), detect the closest minute
  # we want the minute just before the transition
  # e.g. for 24 hrs, we want ~ 23.99
  # tramin = transition minutes
  tramin <- sapply(suns, function(sun) {
    findLightTransitionMinute(Zeitgeberdurations=mid$zhrs,
                              transitionHour=sun)
  })
  # make tramin in a small table
  # where each row is one day/night
  # and two columns: start / end
  # will simplify below
  # column `start` are all uneven indices; column `stop` are all even indices
  # total length should be even, so we can do:
  tramin <- data.frame(start=tramin[seq(1, length(tramin)-1, 2)],
                       stop=tramin[seq(2, length(tramin), 2)])
  
  
  
  # preallocate list dn for day/night
  # each slot will be middur data for one day or night
  # each day/night is defined by one start and one stop, so just divide by 2 to know number of days/nights
  dn <- vector(mode='list', length=nrow(tramin))
  # ! this will skip any day or night that is not full
  # put appropriate names
  # for now, will force user in app to simply give day1, night1, day2, night2
  # so will set it this way
  # but could be more flexible in the future
  names(dn) <- c('day1', 'night1', 'day2', 'night2')
  
  # take the first few names according to how many we need
  # names(dn) <- dn_nms[1:length(dn)]
  
  # fill the list
  # for each window, we take from row just after transition until row just before transition
  # e.g. night0 will be something like: 14.01 hr until 23.99 hr
  # as we are dealing with one-minute binned data, there are not too many rows, so convert data.table to dataframe which I find easier
  
  # loop through small table tramin
  # for each row, we take start and stop
  # and that slice that chunk of data from mid
  for(w in 1:nrow(tramin)) {
    dn[[w]] <- as.data.frame(mid[ (tramin[w, 'start'] + 1) : (tramin[w, 'stop']) , ])
  }
  
  # return the list
  return(dn)
  
}


# genotype utilities ------------------------------------------------------

# calculateParameters prepares a list
# where each slot is one parameter table
# we need to add the grp column to this list

# some of these functions are copied from package FramebyFrame
# but best to keep projects somewhat independent


### function importGenotype(...) ###
# small command to import .txt genotype file
importGenotype <- function(genopath) {
  if( substrEnding(genopath, 3) != 'txt') stop('\t \t \t \t >>> Genotype file does not finish by .txt \n')
  
  if(!file.exists(genopath)) stop('\t \t \t \t Error: file ', genopath, ' does not exist! \n')
  
  geno <- read.delim(genopath, header=TRUE, skip=1, na.strings=c("","NA"))
  return(geno)
}


### function fillGenotype(...) ###
# makes sure genotype file contains all fish, adding an 'excluded' grp/column if necessary
# Note, boxnum is required here as genotype always given 1>>96 regardless of box1 or box2
fillGenotype <- function(geno,
                         boxnum=1) {
  
  genall <- as.vector(unlist(geno)) # pool all the fish in genotype file
  genall <- sort(genall[!is.na(genall)]) # remove any NA & sort
  
  # check we have a sensical number of fish
  if (length(genall)==0) stop('\t \t \t \t >>> There is no fish ID in this genotype file. \n')
  if (length(genall) > 96) stop('\t \t \t \t >>> More than 96 fish IDs in this genotype file, is this right? \n')
  cat('\t \t \t \t >>> Found N = ', length(genall), 'fish IDs in genotype file \n')
  
  # check fish IDs are given as 1 to 96
  if ( ! identical(unique(genall %in% 1:96), TRUE) ) stop('\t \t \t \t >>> Some fish IDs are not within 1--96. Maybe you gave fish IDs of second box as 97--192? Please change them to 1--96. \n')
  
  # check no duplicates in genotype file
  if(sum(duplicated(genall))!=0) stop('\t \t \t \t >>> Some wells are duplicated in genotype file. \n')
  
  # any fish missing from genotype file? if yes -- add them to excluded column in genotype file
  exclu <- (1:96) [which(! 1:96 %in% genall)] # i.e. which of 1, 2, 3, ..., 96 is not in genotype file
  
  if (length(exclu) != 0) {
    # preallocate excluded column as all NA
    geno$excluded <- NA
    # add any excluded fish
    geno[1:length(exclu), 'excluded'] <- exclu
    
    # now all fish should be present in genotype file
    genall <- as.vector(unlist(geno)) # pool all the wells in genotype file
    genall <- sort(genall[!is.na(genall)]) # remove any NA & sort
    
    # check looks good
    if( length((1:96) [which(! 1:96 %in% genall)]) != 0 )
      stop('\t \t \t \t >>> Still some wells missing from genotype file even after adding excluded column \n')
  }
  
  # if second box, add 96 to the all the fish IDs
  if (boxnum==2) {
    geno <- geno + 96
  }
  
  # now ready to return
  return(geno)
}


### function assignGenotype(...) ###
# assignGenotype takes a vector of fish IDs as input (in the form f1, f2, f3, ...)
# for each, it looks up its genotype in the genotype file (geno) and returns its genotype
# eg. f1, f1, f2 >>> scr, scr, ko
assignGenotype <- function(fs, geno) {
  sapply(fs,
         function(f){ # wn = well name
           colnames(geno)[which(geno==substr(f, 2, 99), arr.ind=TRUE)[2]]
         })
}
# details
# substr(): strips the f from column names, f1, f2, ... >> 1, 2, ..., so can match vs genotype file
# which(arr.ind=TRUE): find that well number in the genotype file (arr.ind=TRUE gives row and column coordinate, then [2] takes column number)
# colnames(geno): get the column name of that element
# eg. f12 >> look for 12 in genotype file >> it is in column 2 >> column 2 name is SCR
# any fish which has now grp = NA was not mentioned in genotype file, so it must have been empty or excluded



# misc functions ----------------------------------------------------------

### function diff2(...) ###
# I want to 'correct' R's diff function
# e.g. diff(c(FALSE, FALSE, TRUE, TRUE)) gives 0 1 0
# so if look for first TRUE using diff, will give #2
# when actually it is #3
# I think it should return NA 0 1 0, then positions will match original vector
diff2 <- function(vector) {
  
  return(c(NA, diff(vector)))
  
}

### function substrEnding ###
# take n last characters of a string
# e.g. substrEnding('210907_12_RAWs.csv', 3) >>> csv
substrEnding <- function(x, n){ # x = string (or vector of strings); n = last n characters
  substr(x, nchar(x)-n+1, nchar(x))
}

