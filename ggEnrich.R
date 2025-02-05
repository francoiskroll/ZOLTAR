#####################################################
# ~ ZFAD: predictive pharmacology, plotting functions ~
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################



# packages ----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)


# function ggDraws(...) ---------------------------------------------------

# ! make sure to match as best as possible the same settings as drugEnrichment(...)
# so that plot actually illustrates the analysis you did

# example
# ggDraws(vdbr=read.csv(here('predPharma', 'scn1labF0/scn1labRanked_f0both.csv')),
#         namesPath=here('annotateDrugDb', 'drugAnnotations_TTD.xlsx'),
#         annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
#         annotation='zebrafishSTITCH',
#         testAnnotation='ENSDARP00000010895',
#         ndraws=100000,
#         minScore=700,
#         whichRank='rankeq',
#         exportPath=here('predPharma', 'scn1labF0both_enrichscn1lab.pdf'),
#         width=100,
#         height=100)

ggDraws <- function(vdbr,
                    namesPath,
                    annotationPath,
                    annotation,
                    testAnnotation,
                    ndraws,
                    minScore=NA,
                    whichRank='rank0',
                    exportPath,
                    width=100,
                    height=100) {
  
  # note, in sampleEnrich(...) we do not record the results from the draws to save (a lot) of time
  # but it means we need to perform the random draws again here
  # so function repeats some steps from drugEnrich(...) and sampleEnrich(...)
  
  
  ### we assume vdbr is list of targets ranked from most positive cosine to most negative cosine
  # confirm this:
  if(vdbr[1,'cos'] < 0)
    stop('\t \t \t \t >>> Error sampleEnrich: top drug does not have a positive cosine.
         sampleEnrich expects list of drugs ranked from most positive cosine at the top to most negative cosine at the bottom \n')
  if(vdbr[nrow(vdbr),'cos'] > 0)
    stop('\t \t \t \t >>> Error sampleEnrich: last drug does not have a negative cosine.
         sampleEnrich expects list of drugs ranked from most positive cosine at the top to most negative cosine at the bottom \n')
  
  
  ### swap each drug for its annotations
  ### swap each drug for its annotations
  if(annotation=='indications') {
    annotationCol <- 'indication'
  } else if(annotation=='TTDtargets') {
    annotationCol <- 'TARGETID'
  } else if(annotation=='KEGG') {
    annotationCol <- 'keggid'
  } else if(annotation=='zebrafishSTITCH') {
    annotationCol <- 'ENSP'
  } else if(annotation=='humanSTITCH') {
    annotationCol <- 'ENSP'
  } else if(annotation=='zebrafishGO') {
    annotationCol <- 'go_id'
  } else if(annotation=='humanGO') {
    annotationCol <- 'go_id'
  }
  
  radf <- swapDrugsforAnnotations(vdbr=vdbr,
                                  namesPath=namesPath,
                                  annotationPath=annotationPath,
                                  annotation=annotation,
                                  annotationCol=annotationCol,
                                  minScore=minScore)
  
  
  ### how many instances of the annotation of interest (testAnnotation) are there?
  nano <- length(which(radf[,annotationCol]==testAnnotation))
  # if 0 instances, stop, there must be a typo or else
  if(nano==0) stop('\t \t \t \t >>> Error ggDraws: annotation "', testAnnotation, '" does not exist. Check spelling. \n')
  
  
  ### calculate sum of ranks for given annotation
  # this gives a measure of "how distant from 0, on either side"
  # which type of rank we use depends on setting whichRank
  sr <- sum(radf[which(radf[,annotationCol]==testAnnotation), whichRank]) # sr for sum of ranks
  
  cat('\t \t \t \t >>> Sum of', whichRank, '=', sr, '\n')
  
  
  ### what is the best possible sum of ranks?
  # say there are 20 instances of a given annotation
  # it is simply the top20 whichRank, i.e. most distant from 0
  # this could represent an extreme left enrichment or extreme right enrichment or extreme binomial
  
  # EDIT: I think simply taking the top e.g. 10 overestimates slightly the best possible sum of ranks
  # e.g. aspirin is top correlating drug, it has rankeq = 1000 and it has 10 targets
  # all of aspirin's targets will be assigned rankeq = 1000, so best possible sum of ranks will be = 1000 + 1000 + 1000 + ...
  # but in practice, one compound (usually*) does not have twice the same target
  # [* note this does happen as of 01/02/2023 because we were quite conservative when removing duplicates in STITCH]
  # [so we left e.g. same compound interacting with same protein, but once inhibition and once activation]
  # rather, best possible sum of ranks should reflect a situation where top 10 compounds all interact with target of interest
  # to simulate this, simply take the top 10 *unique* whichRank
  # e.g. 1000 + 998.5 + ... (even if 1000 is actually repeated multiple times in radf)
  # so take *unique* whichRank, sort them, sum the top `number of instances`
  bsr <- sum ( rev(sort(unique(radf[,whichRank]))) [1:nano] )
  
  
  ### SAMPLING PROCEDURE ###
  
  # preallocate a list to store the results
  rdra <- vector(mode='list', length=ndraws)
  
  for(nd in 1:ndraws) {
    
    # draw new `number of instances` ranks
    # and sum them
    dsr <- sum(sample(radf[,whichRank], size=nano)) # draw sum of ranks
    
    # store result for this draw in a one-row dataframe
    # which we put in list rdra
    rdra[[nd]] <- data.frame(nexamples=nano, ndra=nd, RS=dsr)
    
  }
  
  # currently a list where each element is results from one draw
  # turn it into a dataframe where each row is one draw
  rdra <- as.data.frame(rbindlist(rdra))
  
  # rank the draws from lowest sum of ranks to highest sum of ranks
  # (will matter below)
  rdra <- rdra[order(rdra$RS),]
  
  
  ### calculate p-value
  # how many simulations gave a more extreme sum of ranks than the one we observed?
  # remember we do not care about small ranks (i.e. around cos ~ 0), so we only need to test if sum of ranks is larger than what we observed
  nsim <- length(which(rdra$RS > sr))
  # approximate p-value is
  pv <- nsim / ndraws
  
  # note, we should not be creating any p-value = 0
  # instead, the minimum p-value depends on the number of draws
  # e.g. if we do N=1000 draws, then the minimum p-value is 1/1000 = 0.001
  if(nsim==0) {
    pv <- 1/ndraws
  }
  
  
  ### plot an 'enrichment plot'
  # which is simply the distribution of sum of ranks across the random draws
  
  # where is p=0.05?
  # i.e. 95th percentile, sum of ranks value where 95% of the random draws are to the left of it, i.e. have a lower sum of ranks
  p95 <- rdra[round(0.95*nrow(rdra)), 'RS']
  
  # where is p=0.01?
  # i.e. 99th percentile, sum of ranks value where 99% of the random draws are to the left of it, i.e. have a lower sum of ranks
  p99 <- rdra[round(0.99*nrow(rdra)), 'RS']
  # (we ranked the draws from lowest sum of ranks to highest, that is why we can simply take first 1% of all rows, etc.)
  
  # prepare the data to plot histogram
  drahist <- freqInBins(df=rdra,
                        datcol='RS',
                        binwidth=100)
  
  # plot histogram
  # I find histograms more intuitive than geom_density()
  ggES <- ggplot(drahist, aes(x=middle, y=freq)) +
    geom_bar(stat='identity') +
    # add observed mark
    geom_vline(xintercept=sr, linewidth=0.75, linetype=1, colour='#cb2a20') +
    # add 95th percentile
    geom_vline(xintercept=p95, linewidth=0.5, linetype=2, colour='#aeb3b4') +
    # add 99th percentile
    geom_vline(xintercept=p99, linewidth=0.5, linetype=2, colour='#585e60') +
    theme_minimal() +
    theme(
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      axis.title.x=element_text(size=9, margin=margin(t=-2, r=0, b=0, l=0)),
      axis.title.y=element_text(size=9),
      axis.text.x=element_text(size=7, margin=margin(t=-2, r=0, b=0, l=0)),
      axis.text.y=element_text(size=7, margin=margin(t=0, r=-2, b=0, l=0))
    ) +
    {if(whichRank=='abscos') xlab('sum of cosines')} +
    {if(whichRank %in% c('rank', 'rankeq', 'rank0', 'ranks')) xlab('sum of ranks')} +
    ylab('frequency')
  
  
  # print some statistics for user
  cat('\t \t \t \t >>>', nsim, 'draws out of', ndraws, 'gave a more extreme sum of ranks.\n
      \t \t \t Estimated p-value = ', pv, '\n')
  
  # export plot
  ggsave(exportPath, ggES, width=width, height=height, units='mm', device=cairo_pdf)
  
  return(ggES)
  
}



# ggBarcode(...) ----------------------------------------------------------

ggBarcode <- function(vdbr,
                      namesPath,
                      annotationPath,
                      annotation,
                      testAnnotation,
                      minScore=NA,
                      barwidth1=2,
                      barwidth2=25,
                      exportPath,
                      width=100,
                      height=100) {
  
  
  ### we assume vdbr is list of annotations ranked from most positive cosine to most negative cosine
  # confirm this:
  if(vdbr[1,'cos'] < 0)
    stop('\t \t \t \t >>> Error sampleEnrich: top drug does not have a positive cosine.
         sampleEnrich expects list of drugs ranked from most positive cosine at the top to most negative cosine at the bottom \n')
  if(vdbr[nrow(vdbr),'cos'] > 0)
    stop('\t \t \t \t >>> Error sampleEnrich: last drug does not have a negative cosine.
         sampleEnrich expects list of drugs ranked from most positive cosine at the top to most negative cosine at the bottom \n')
  
  
  ### swap each drug for its annotations
  if(annotation=='indications') {
    annotationCol <- 'indication'
  } else if(annotation=='TTDtargets') {
    annotationCol <- 'TARGETID'
  } else if(annotation=='KEGG') {
    annotationCol <- 'keggid'
  } else if(annotation=='zebrafishSTITCH') {
    annotationCol <- 'ENSP'
  } else if(annotation=='humanSTITCH') {
    annotationCol <- 'ENSP'
  } else if(annotation=='zebrafishGO') {
    annotationCol <- 'go_id'
  } else if(annotation=='humanGO') {
    annotationCol <- 'go_id'
  }
  
  radf <- swapDrugsforAnnotations(vdbr=vdbr,
                                  namesPath=namesPath,
                                  annotationPath=annotationPath,
                                  annotation=annotation,
                                  annotationCol=annotationCol,
                                  minScore=minScore)
  
  
  ### which rows have the annotation of interest?
  radf <- radf %>%
    mutate(isSet=(radf[,annotationCol]==testAnnotation))
  
  
  # print the rows which are isSet
  # sometimes useful to label plot manually
  cat('\t \t \t \t >>> Fingerprints with given annotation:\n')
  radfSet <- radf[which(radf$isSet),]
  print(radfSet)
  
  # some summary metrics
  
  if(annotation=='indications') {
    cat('\t \t \t \t >>> total', nrow(radf),'drug-indication pairs (', length(unique(radf$cid)), 'unique compounds )\n')
  } else if(annotation=='TTDtargets') {
    cat('\t \t \t \t >>> total', nrow(radf),'drug-target pairs (', length(unique(radf$cid)), 'unique compounds )\n')
  } else if(annotation=='KEGG') {
    cat('\t \t \t \t >>> total', nrow(radf),'drug > target > KEGG pathway interactions (', length(unique(radf$cid)), 'unique compounds )\n')
  } else if(annotation=='humanSTITCH') {
    cat('\t \t \t \t >>> total', nrow(radf),'drug-indication pairs (', length(unique(radf$cid)), 'unique compounds )\n')
  } else if(annotation=='zebrafishSTITCH') {
    cat('\t \t \t \t >>> total', nrow(radf),'drug-indication pairs (', length(unique(radf$cid)), 'unique compounds )\n')
  }
  
  cat('\t \t \t \t \t >>> of which,', nrow(radfSet),'fingerprints;', length(unique(radfSet$cid)),'unique compounds have given annotation.\n')
  cat('\t \t \t \t \t >>>', length(which(radfSet$cos>=0.2)), '/', nrow(radfSet),'fingerprints have a cosine >= 0.2.\n')
  cat('\t \t \t \t \t >>>', length(which(radfSet$cos<=-0.2)), '/', nrow(radfSet),'fingerprints have a cosine <= -0.2.\n')
  
  # we will label the position with maximum cos, the position with minimum cos
  # and cos ~ 0, i.e. the position at which cos is closest to 0
  # rows to label:
  cos0pos <- which.min(abs(radf$cos - 0))
  # so rows to label are:
  rw2l <- c(1, cos0pos, nrow(radf))
  
  ggBco <- ggplot()
  # add lines to label a few cos
  for (i in 1:length(rw2l)) {
    ggBco <- ggBco + geom_segment(aes(y=1, yend=0.45, x=radf[rw2l, 'rank'], xend=radf[rw2l, 'rank']),
                                  linetype=2, linewidth=0.75, colour='#ebebeb')
  }
  
  # add the labels
  ggBco <- ggBco +
    geom_text(data=radf[rw2l,], aes(x=rank, y=0.40, label=round(cos, 2)), colour='#4d4d4d', size=2.5) # note ***
  # *** for cos at row cos0pos to not be exactly = 0, the smallest cos needs to be > 0.009
  # it seems like usually it is smaller than that, but we calculate it instead of assuming it is = 0 so we can detect cases when it is not
  
  # add all the bars *not* from the set
  ggBco <- ggBco +
    geom_tile(data=filter(radf, !isSet), aes(x=rank, y=1), fill='#aeb3b4', width=barwidth1)
  
  
  # add the bars from the set
  ggBco <- ggBco +
    geom_tile(data=filter(radf, isSet), aes(x=rank, y=1), fill='#fcb505', colour='white', width=barwidth2, linewidth=0.001) +
    # add the arrows
    geom_segment(data=filter(radf, isSet), aes(x=rank, y=1.52, xend=rank, yend=1.51),
                 colour='#fcb505',
                 lineend='butt',
                 linejoin='mitre', 
                 arrow=arrow(length=unit(2, 'mm'),
                             type='closed')) +
    # and polish
    theme_minimal() +
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_blank(), 
      axis.text.y=element_blank(),
      panel.grid=element_blank(),
      legend.position='none'
    )
  
  ### export barcode plot, if required
  if(!is.na(exportPath)) {
    ggsave(exportPath, ggBco, width=width, height=height, units='mm') 
  }
  
  ### return so displays in RStudio
  return(ggBco)
  
  
}



# function freqInBins(...) ------------------------------------------------

# counts & frequencies within bins of a continuous variable (each bin = one range of values) per group
freqInBins <- function(df, datcol, grpcol=NA, binwidth) {
  
  ###
  # if no grouping
  if (is.na(grpcol)) {
    dat <- as.numeric(df[, datcol]) # data we need to bin
    hist <- hist(dat, breaks=binwidth, include.lowest=TRUE, plot=FALSE)
    lowbound <- hist$breaks # lower boundary of each bin
    # we delete last one which is the the upper boundary of the last bin
    lowbound <- lowbound[-length(lowbound)]
    mids <- hist$mids # centre of each bin
    cou <- hist$counts # count in each bin
    fre <- cou / sum(cou) # frequency in each bin, i.e. proportion of the entire dataset
    gdf <- data.frame(lowbound=lowbound, middle=mids, counts=cou, freq=fre) # prepare small dataframe, two columns: counts and frequencies
    
    return(gdf)
    
  }
  
  ###
  # below taken from PRNP x NANOPORE project, but should update to match above
  # # if grouping
  # 
  # else {
  #   
  #   ngrps <- length(unique(df[,grpcol])) # number of groups
  #   grpnms <- unique(df[,grpcol]) # group names
  #   
  #   cfpg <- vector(mode='list', length=ngrps) # preallocate list counts + frequencies per group
  #   
  #   for(i in 1:length(grpnms)) { # loop through groups
  #     dat <- as.numeric(df[ which(df[,grpcol] == grpnms[i]) , datcol]) # data for this group (the data we need to bin)
  #     cou <- hist(dat, breaks=binwidth, include.lowest=TRUE, plot=FALSE)$counts # counts
  #     fre <- cou / sum(cou) # frequencies
  #     
  #     gdf <- data.frame(grp=grpnms[i], bin=binwidth[-1], counts=cou, freq=fre) # small dataframe for this group, two columns counts and frequencies
  #     
  #     # group column should be called the same as grpcol
  #     colnames(gdf)[1] <- grpcol
  #     
  #     # add it to the list
  #     cfpg[[i]] <- gdf
  #   }
  #   
  #   # stick the list together and return it
  #   fbing <- as.data.frame(rbindlist(cfpg)) # frequencies in each bin by group
  #   return(fbing)
  # }
  
}
