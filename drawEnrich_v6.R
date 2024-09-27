###################################################
# ~~~ predictive pharmacology functions ~~~

# function sampleEnrich
# input is a ranked list of annotations,
# from maximum (closest to 1.0) to minimum (closest to -1.0) correlation with fingerprint
# function calculates enrichment for one annotation, namely:
# positive enrichment, i.e. is annotation enriched in the positive correlations
# negative enrichment, i.e. is annotation enriched in the negative correlations

# Francois Kroll 2022
# francois@kroll.be
###################################################

# up to v3 was called sampleEnrich.R but I did not like the name anymore

### v2
# no more splitting the ranked list in positive and negative cosines
# this approach could miss a situation where e.g. there are 100 Dopamine D2 receptor drugs; they are ALL in negative cosines but not in extreme ranks
# instead, have ranks represent distance from 0 (but all positive) and do sum of ranks as before

### v3
# idea to speed up the draws
# take 10,000 draws
# as soon as 500 draws have a more extreme sum of ranks, we know the final p-value will be > 0.05
# so we could keep a count and report the minimum possible p-value as soon as there is no point trying further

### v4
# common function swapDrugsforAnnotations instead of unique function for each annotation

### v5
# some edits to make it work within Shiny app

### v6
# testing if works fine using drugdbSUM, which is version where each compound is present as only one fingerprint
# prepared by summariseDrugDb.R
# if one CID had multiple fingerprints, we averaged them

# using cos directly as rank
# I think it makes more sense as it takes more directly how well the fingerprints correlate
# example to think about it: a situation where the max cos is 0.9 and many fingerprints with the given annotation are there
# should give a more significant (lower p-value) result than a situation where the max cos is 0.4 and many fingerprints are there
# but using ranks does not make distinction between these situations
# I think using cosines should be a bit more quantitative in that sense
# also, average cosine of fingerprints with that annotation makes sense I think
# is imbalance in number of fingerprints on positive & negative sides an issue?
# imagine there are 10 fingerprints on negative side / 100 fingerprints on positive side
# and the fingerprints w/ annotation are the top 3 negative, cos = 0.6
# and there are many fingerprints on positive side after cos = 0.6
# random draws will almost always pick fingerprints on positive side, so likely to go above 3 * 0.6
# I think that is by design? It is saying "it looks like a surprising result, but it is in fact easy to go above 3 * cos 0.6"


# packages ----------------------------------------------------------------

library(ggpubr)
library(lsa)
library(openxlsx)
library(data.table)
library(stringr)
library(tibble)
library(dplyr)


# function testAllAnnotations(...) ----------------------------------------

testAllAnnotations <- function(vdbr,
                               namesPath=NA,
                               annotationDir,
                               whichRank='rank0',
                               minScore=NA,
                               minNex,
                               ndraws,
                               alphaThr=0.05,
                               maxPval=NA,
                               statsDir,
                               exportPrefix) {
  # 1-- indications
  drugEnrichment(vdbr=vdbr,
                 namesPath=namesPath,
                 annotationPath=paste0(annotationDir, 'TTDindications.csv'),
                 annotation='indications',
                 whichRank=whichRank,
                 minNex=minNex,
                 ndraws=ndraws,
                 alphaThr=alphaThr,
                 statsExport=paste0(statsDir, exportPrefix, '_indications.csv'))
  
  # 2-- TTD targets
  drugEnrichment(vdbr=vdbr,
                 namesPath=namesPath,
                 annotationPath=paste0(annotationDir, 'TTDtargets.csv'),
                 annotation='TTDtargets',
                 whichRank=whichRank,
                 minNex=minNex,
                 ndraws=ndraws,
                 alphaThr=alphaThr,
                 statsExport=paste0(statsDir, exportPrefix, '_TTDtargets.csv'))
  
  # 3-- KEGG pathways
  drugEnrichment(vdbr=vdbr,
                 namesPath=namesPath,
                 annotationPath=paste0(annotationDir, 'TTDkegg.csv'),
                 annotation='KEGG',
                 whichRank=whichRank,
                 minNex=minNex,
                 ndraws=ndraws,
                 alphaThr=alphaThr,
                 statsExport=paste0(statsDir, exportPrefix, '_KEGG.csv'))
  
  # 4-- zebrafish STITCH targets
  drugEnrichment(vdbr=vdbr,
                 namesPath=namesPath,
                 annotationPath=paste0(annotationDir, 'zebrafishSTITCH.csv'),
                 annotation='zebrafishSTITCH',
                 whichRank=whichRank,
                 minScore=minScore,
                 minNex=minNex,
                 ndraws=ndraws,
                 alphaThr=alphaThr,
                 maxPval=maxPval,
                 statsExport=paste0(statsDir, exportPrefix, '_zSTITCH.csv'))
  
  # 5-- human STITCH targets
  drugEnrichment(vdbr=vdbr,
                 namesPath=namesPath,
                 annotationPath=paste0(annotationDir, 'humanSTITCH.csv'),
                 annotation='humanSTITCH',
                 whichRank=whichRank,
                 minScore=minScore,
                 minNex=minNex,
                 ndraws=ndraws,
                 alphaThr=alphaThr,
                 maxPval=maxPval,
                 statsExport=paste0(statsDir, exportPrefix, '_hSTITCH.csv'))
  
  # 6-- zebrafish GO
  drugEnrichment(vdbr=vdbr,
                 namesPath=namesPath,
                 annotationPath=paste0(annotationDir, 'zebrafishGO.csv'),
                 annotation='zebrafishGO',
                 whichRank=whichRank,
                 minScore=minScore,
                 minNex=minNex,
                 ndraws=ndraws,
                 alphaThr=alphaThr,
                 maxPval=maxPval,
                 statsExport=paste0(statsDir, exportPrefix, '_zGO.csv'))
  
  # 7-- human GO
  drugEnrichment(vdbr=vdbr,
                 namesPath=namesPath,
                 annotationPath=paste0(annotationDir, 'humanGO.csv'),
                 annotation='humanGO',
                 whichRank=whichRank,
                 minScore=minScore,
                 minNex=minNex,
                 ndraws=ndraws,
                 alphaThr=alphaThr,
                 maxPval=maxPval,
                 statsExport=paste0(statsDir, exportPrefix, '_hGO.csv'))
  
}


# function drugEnrichment(...) --------------------------------------------

# overarching function to test enrichment of annotations
# vdb is ranked drugDb

drugEnrichment <- function(vdbr,
                           namesPath=NA,
                           annotationPath,
                           annotation,
                           whichRank='rank0',
                           minScore=NA,
                           minNex,
                           ndraws,
                           alphaThr=0.05,
                           maxPval=NA,
                           statsExport=NA) {
  
  # make sure line below is active when using in ZOLTAR app
  source('drawEnrich_v6.R')
  
  
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
  
  
  ### test enrichments with random draws
  # loop through unique annotations and test each
  uan <- sort(unique(radf[,which(colnames(radf)==annotationCol)]))
  
  anr <- lapply(1:length(uan), function(a) {
    
    ua <- uan[a]
    
    # ## loading bar
    # # how far are we in %?
    # perpro <- round((a/length(uan)) * 100)
    # # if 0, replace by 1 to keep below simple
    # if(perpro==0) {
    #   perpro <- 1
    # }
    # # that is how many _ with replace by # in the loading bar
    # lobar <- paste0( paste(rep('#', perpro), collapse='') , paste(rep('_', 100-perpro), collapse='') , collapse='')
    # # cat('\t \t \t \t', lobar, '\n')
    # # cat('\t \t \t \t >>> Testing', annotation,':', ua,'\n')
    
    sampleEnrich(radf=radf,
                 annotationCol=annotationCol,
                 testAnnotation=ua,
                 whichRank=whichRank,
                 minNex=minNex,
                 ndraws=ndraws,
                 alphaThr=alphaThr,
                 maxPval=maxPval)
    
  })
  
  anr <- as.data.frame(rbindlist(anr))
  
  cat('\n \t \t \t \t Correction of random draws p-values \n')
  
  ### Bonferonni correction of pvals
  siBon <- signiBonferonni(pvs=anr$pval,
                           alphaThr=0.05)
  anr <- anr %>%
    add_column(bonferonni=siBon, .after='pval')
  
  ### Benjamini-Hochberg correction of pvals
  siBH <- signiBenjaminiHochberg(pvs=anr$pval,
                                 alphaThr=0.05)
  
  anr <- anr %>%
    add_column(benhoch=siBH, .after='bonferonni')
  
  #######
  
  cat('\n \t \t \t \t Correction of KS test p-values \n')
  
  ### Bonferonni correction of KS pvals
  kssiBon <- signiBonferonni(pvs=anr$kspval,
                             alphaThr=0.05)
  anr <- anr %>%
    mutate(KSbonferonni=kssiBon, .after='kspval')
  
  
  ### Benjamini-Hochberg correction of pvals
  kssiBH <- signiBenjaminiHochberg(pvs=anr$kspval,
                                   alphaThr=0.05)
  anr <- anr %>%
    add_column(KSbenhoch=kssiBH, .after='KSbonferonni')
  
  
  ## add back some information about the annotations
  if(annotation %in% c('zebrafishSTITCH', 'humanSTITCH')) {
    # import zebrafish STITCH
    sti <- fread(here(annotationPath), na.strings='')
    # keep only the protein annotations
    # we only need one example for each protein (cf. duplicated below)
    sti <- sti[!duplicated(sti$ENSP), c('species', 'ENSP', 'ENSG', 'gene_name', 'gene_symbol', 'gene_biotype')]
    # change column name so we can merge with statistics results
    colnames(sti)[which(colnames(sti)==annotationCol)] <- 'annotation'
    # merge
    # we want to keep all rows in anr (and not necessarily take all rows of sti)
    nrowbefore <- nrow(anr)
    anr <- right_join(sti, anr, by='annotation')
    nrowafter <- nrow(anr)
    if(nrowbefore!=nrowafter) stop('\t \t \t \t >>> Error drugEnrichment: unexpected result when adding information to the statistics report \n')
    
  } else if(annotation=='TTDtargets') {
    # import TTD targets
    tar <- read.csv(annotationPath)
    # keep only the protein annotations
    # we only need one example for each protein/target (cf. duplicated below)
    tar <- tar[!duplicated(tar$TARGETID), c('TARGETID', 'UNIPROID', 'TARGNAME', 'GENENAME', 'TARGTYPE', 'BIOCLASS')]
    # change column name so we can merge with statistics results
    colnames(tar)[which(colnames(tar)==annotationCol)] <- 'annotation'
    # merge
    # we want to keep all rows in anr (and not necessarily take all rows of tar)
    nrowbefore <- nrow(anr)
    anr <- right_join(tar, anr, by='annotation')
    nrowafter <- nrow(anr)
    if(nrowbefore!=nrowafter) stop('\t \t \t \t >>> Error drugEnrichment: unexpected result when adding information to the statistics report \n')
    
  } else if(annotation=='KEGG') {
    # import TTD targets
    keg <- read.csv(annotationPath)
    # keep only the KEGG names
    # we only need one example for each KEGG ID (cf. duplicated below)
    keg <- keg[!duplicated(keg$keggid), c('keggid', 'keggname')] # we only keep keggid and keggname, all we need to do is add the names here
    # change column name so we can merge with statistics results
    colnames(keg)[which(colnames(keg)==annotationCol)] <- 'annotation'
    # merge
    # we want to keep all rows in anr (and not necessarily take all rows of tar)
    nrowbefore <- nrow(anr)
    anr <- right_join(keg, anr, by='annotation')
    nrowafter <- nrow(anr)
    if(nrowbefore!=nrowafter) stop('\t \t \t \t >>> Error drugEnrichment: unexpected result when adding information to the statistics report \n')
  } else if(annotation %in% c('zebrafishGO', 'humanGO')) {
    # import GO annotations
    go <- fread(here(annotationPath), na.strings='')
    # keep go_term and go_linkagetype
    # we only need one example for each GO term (cf. duplicated below)
    go <- go[!duplicated(go$go_id), c('go_id', 'go_term', 'go_linkagetype')]
    # change column name so we can merge with statistics results
    colnames(go)[which(colnames(go)==annotationCol)] <- 'annotation'
    # merge
    # we want to keep all rows in anr (and not necessarily take all rows of sti)
    nrowbefore <- nrow(anr)
    anr <- right_join(go, anr, by='annotation')
    nrowafter <- nrow(anr)
    if(nrowbefore!=nrowafter) stop('\t \t \t \t >>> Error drugEnrichment: unexpected result when adding information to the statistics report \n')
  }
  
  
  ## order with lowest pval on top
  anr <- anr[order(anr$pval),]
  
  
  ## write the statistics report to drive
  if(!is.na(statsExport)) {
    write.csv(anr, statsExport, row.names=FALSE)
  }
  
  
  ## also return it
  return(anr)
  
}


# function rankDrugDb(...) ------------------------------------------------

# takes legacy fingerprint as input and ranks the drug database
# from drug generating closest phenotype on top to drug generating most opposite phenotype at bottom

rankDrugDb <- function(legacyFgp,
                       dbPath,
                       metric) {
  
  ### import drugDb
  ddb <- read.csv(dbPath)
  
  # 30/03/2023, decided to delete parameters daymean_averageWaking & nightmean_averageWaking
  # can read rationale in paramsFromMid.R
  # will delete them here so we do not alter the data
  ddb$daymean_averageWaking <- NULL
  ddb$nightmean_averageWaking <- NULL
  
  # detect column where Z-scores start
  zcol <- min(which(startsWith(colnames(ddb), c('day', 'night'))))
  
  # detect column where Z-scores end
  zcoll <- max(which(startsWith(colnames(ddb), c('day', 'night'))))
  
  # (this way if we add columns it does not break the function)
  
  ## ! makes it so legacyFgp gives uparam in same order as ddb
  # first, check both give the same parameters
  if(!identical( sort(colnames(ddb)[zcol:zcoll]) , sort(legacyFgp$uparam)))
             stop('\t \t \t \t >>> Error rankDrugDb: drugDb and legacyFgp do not give the same parameters. \n')
  # next, ensure legacyFgp gives them in same order
  legacyFgp <- legacyFgp[match(colnames(ddb)[zcol:zcoll], legacyFgp$uparam),]
  
  # check this (it is crucial for below)
  if(!identical( colnames(ddb)[zcol:zcoll] , legacyFgp$uparam))
    stop('\t \t \t \t >>> Error rankDrugDb: drugDb and legacyFgp do not give parameters in the same order, even after re-ordering. \n')
  
  ### calculate correlation vs legacyFgp
  dcor <- sapply(1:nrow(ddb), function(d) {
    # drug fingerprint is:
    dfp <- as.numeric(ddb[d, zcol:ncol(ddb)])
    # correlation ( drug fingerprint , input fingerprint )
    cor( dfp , legacyFgp$zsco )
  })
  
  ### calculate cosine similarity vs legacyFgp
  dcos <- sapply(1:nrow(ddb), function(d) {
    # drug fingerprint is:
    dfp <- as.numeric(ddb[d, zcol:ncol(ddb)])
    # correlation ( drug fingerprint , input fingerprint )
    cosine( dfp , legacyFgp$zsco )
  })
  
  # add columns
  vdb <- ddb %>%
    add_column(cor=dcor, .after='name') %>%
    add_column(cos=dcos, .after='name')
  # v for drugs VS fingerprint
  
  
  ### order drug db according to metric chosen by user
  if (metric=='correlation') {
    simCol <- which(colnames(vdb)=='cor')
  } else if (metric=='cosine') {
    simCol <- which(colnames(vdb)=='cos')
  } else {
    stop('\t \t \t \t >>> Error rankDrugDb: unknown metric \n')
  }
  # order drug db
  vdbr <- vdb[rev(order(vdb[,simCol])),]
  
  
  ### add rank columns
  
  # "rank" is simply 1 to last row
  vdbr <- vdbr %>%
    add_column(rank=1:nrow(vdbr), .after=1)
  
  # "ranks" is rank/distance from 0, using negative ranks for negative cosines
  # s for signed
  # i.e. 1, 2, 3, ..., -1, -2, -3, ...
  # preallocate the column
  vdbr <- vdbr %>%
    add_column(ranks=NA, .after='rank')
  # set all the ranks of the positive cosines
  vdbr[which(vdbr$cos>0), 'ranks'] <- length(which(vdbr$cos>0)) : 1 # i.e. we go e.g. 2000, 1999, 1998, etc. so that rank 1 is something like cos ~ +0.001
  # set all the ranks of the negative cosines
  vdbr[which(vdbr$cos<0), 'ranks'] <- -(1 : length(which(vdbr$cos<0))) # i.e. we go e.g. 1, 2, 3, etc. so that rank 1 is something like cos ~ -0.001
  
  # "rank0" is absolute rank from 0, i.e. not using negative ranks
  # we simply copy ranks but remove the sign
  vdbr <- vdbr %>%
    mutate(rank0=abs(ranks), .after='rank')
  
  # "rankeq", e stands for "equally distant"
  # an issue arises when there is a strong imbalance
  # in number of annotations with negative cos vs number of annotations with positive cos
  # say there is only 10 annotations on negative cos side
  # and 100 annotations on positive cos side
  # and annotationTest (e.g. protein X) was top3 negatives, i.e. rank 10, 9, 8, sum of ranks = 27
  # random draws will almost always draw ranks in the positives, simply because there are more of them
  # and even mediocre positive ranks will beat 27, e.g. 37 + 40 + 48
  # summing the cos instead of the ranks may marginally help but is likely not a complete solution
  # say most of the annotations are around cos = 0.8, the random draws will always find a bigger sum than the real one,
  # even if the real one is the maximum negative it can be, say 0.7 + 0.7 + 0.7
  # attempt is to essentially weigh the ranks based on the number of annotations in positives or negatives
  # imagine the ranked list is a ruler of 2 cm, going from -1.0 to 1.0, we put cos ~ 0 at 0 cm
  # then we spread out all the annotations on each side,
  # putting as much distance as possible between them to fill the 1-cm of space that we have
  # the distance of each annotation from 0 is now the rank
  # in example above, each annotation in negative (N = 10) would count for 1 cm/10 = 0.1 cm
  # so starting at 0 and moving left, ranks will be -0.1, -0.2, -0.3, etc.
  # each annotation in positive (N = 100) would count for 1 cm/100 = 0.01 cm
  # so starting at 0 and moving right, ranks wil be +0.01, +0.02, +0.03, etc.
  # this does not solve the imbalance per se (there are still more annotations in positive cosines)
  # but the sum of ranks will be down-weighted, i.e. each annotation drawn from the positive cos is worth less
  # I think this makes intuitive sense because in case there is no imbalance, annotations on either side will count for the same
  
  # preallocate the column
  vdbr <- vdbr %>%
    add_column(rankeq=NA, .after='rank')
  # below: does not necessarily need to be 1000, we could pick any number
  # e.g. 1/ ~ 3000 will make a tiny number, which when summed may look like a scale from 0 to 1 (like a correlation or p-value)
  # using e.g. 1000 makes it so so we usually have some number above 1 when summing a few ranks
  posstep <- 1000/length(which(vdbr$cos>0))
  negstep <- 1000/length(which(vdbr$cos<0))
  vdbr[which(vdbr$cos>0), 'rankeq'] <- rev(seq(posstep, posstep*length(which(vdbr$cos>0)), posstep))
  # from top of the list going down up to cos ~ 0, rankeq goes something like 1000, 999.5, 999, etc.
  vdbr[which(vdbr$cos<0), 'rankeq'] <- seq(negstep, negstep*length(which(vdbr$cos<0)), negstep)
  # from cos ~ 0 going down, rankeq goes something like 0.5, 1.0, 1.5, etc.
  # check no more NA:
  if(sum(is.na(vdbr$rankeq))>0) stop('\t \t \t \t >>> Error: some NA in rankeq column. \n')
  
  # note, ranks or rank0 or rankeq = 0 does not exist
  # right in the middle of the list you have something like:
  # ranks = 1 / rank0 = 1 / rankeq = 0.5 / cos = +0.001
  # ranks = -1 / rank0 = 1 / rankeq = 0.25 / cos = -0.001
  
  # v6 20/08/2024: we also add absolute cosine column
  vdbr <- vdbr %>%
    mutate(abscos=abs(cos), .after='cos')
  
  return(vdbr)
  
}


# function swapDrugsforAnnotations(...) -----------------------------------
# this is a common function to switch the compounds for their annotations
# previously, was a bunch of functions such as rankIndications and rankSTITCHzebrafish
# but they were mostly redundant, best to perform the task with one common function

swapDrugsforAnnotations <- function(vdbr,
                                    namesPath,
                                    annotationPath,
                                    annotation,
                                    annotationCol,
                                    minScore=NA) {
  
  ### check we can deal with this annotation
  if(! annotation %in% c('indications', 'TTDtargets', 'KEGG', 'zebrafishSTITCH', 'humanSTITCH', 'zebrafishGO', 'humanGO'))
    stop('\t \t \t \t >>> Error drugEnrichment: does not currently support this annotation.
         Possible annotations are: "indications", "TTDtargets", "KEGG", "zebrafishSTITCH", "humanSTITCH", "zebrafishGO", "humanGO"')
  
  ### import names
  dnms <- read.csv(namesPath)
  
  
  ### import annotations
  # method depends slightly on the annotation
  if(annotation %in% c('zebrafishSTITCH', 'humanSTITCH', 'zebrafishGO', 'humanGO')) {
    ano <- fread(annotationPath, na.strings='')
  } else if (annotation %in% c('indications', 'TTDtargets', 'KEGG')) {
    ano <- read.csv(annotationPath)
  }
  
  # data.table format is causing some trouble later
  ano <- as.data.frame(ano)
  
  # if we are looking at STITCH targets and minScore is given, trim some interactions
  if(annotation %in% c('zebrafishSTITCH', 'humanSTITCH', 'zebrafishGO', 'humanGO') & !is.na(minScore)) {
    cat('\t \t \t \t >>> minScore threshold =', minScore, ': keeping', nrow(subset(ano, score >= minScore)), 'interactions out of', nrow(ano), '\n')
    ano <- subset(ano, score >= minScore)
    cat('\t \t \t \t >>>' , length(unique(ano[, annotationCol])),'unique proteins or GO terms left \n')
  }
  
  # now to order annotations from top correlating to worst
  # cdb should basically be used as entries to database ano
  # for every CID, we return all annotations and repeat the cosine
  # if that CID occurs again, we return the same annotations again
  # 30/01/2023: we now assign same rank to every target of a given compound
  # e.g. aspirin is rank 10 and has 5 targets, each of these targets gets rank 10
  # previously, we were assigning ranks on the annotations
  # I think this is more accurate, because e.g. aspirin's targets could have been rank 56, 57, 58, 59, 60
  # while the difference between rank 56 and rank 60 is meaningless (they all have the same cos, so are at the same distance from 0)
  # graphically, it is like these 5 targets should all be on top of each other at the same distance from 0
  stl <- lapply(1:nrow(vdbr), function(cr) {
    
    # indices with annotations for this drug (look by CID)
    it <- which(ano$cid==vdbr[cr, 'cid'])
    
    # if no targets
    if (length(it)==0) {
      NULL
      
      # if some targets
    } else {
      # return the rows
      return ( cbind(ano[it,],
                     cos=vdbr[cr, 'cos'],
                     abscos=vdbr[cr, 'abscos'],
                     rank=vdbr[cr, 'rank'],
                     rank0=vdbr[cr, 'rank0'],
                     ranks=vdbr[cr, 'ranks'],
                     rankeq=vdbr[cr, 'rankeq']))
    }
  })
  stl <- rbindlist(stl)
  
  # about replicates (same drug/different experiments):
  # say one drug has 3 replicates, so 3 cos values
  # and this drug has 15 STITCH targets
  # above will record all 15 STITCH targets for replicate1 / those will be assign all same cos & all same rank, from replicate1
  # then will record all 15 STITCH targets again for replicate2 / those will be assign all same cos & all same rank, from replicate2
  # etc.
  # which I think is what we want
  # alternative would be to average the 3 cos but I think better to keep all the data
  
  # as we added rows, we can order again
  # (but should be useless as still ordered by cos)
  stl <- stl[rev(order(stl$cos)),]
  # remember, one drug x target interaction can be at multiple positions
  # as each replicate of this drug has its own cos
  
  return(as.data.frame(stl))
  
  
}

# function sampleEnrich(...) ----------------------------------------------

# are these scores surprising?
# we can tell by drawing some rows at random and calculating the sum of ranks, many times
# small function to do this
sampleEnrich <- function(radf,
                         annotationCol,
                         testAnnotation,
                         whichRank,
                         minNex,
                         ndraws,
                         alphaThr=0.05,
                         maxPval=NA) {
  # ! we assume radf is ranked list of annotations
  # from max cos to min cos in comparison to a given fingerprint
  
  # remember radf is ranked from cos ~ 1 to cos ~ -1
  # confirm this:
  if(radf[1,'cos'] < 0)
    stop('\t \t \t \t >>> Error sampleEnrich: top drug does not have a positive cosine.
         sampleEnrich expects list of drugs ranked from most positive cosine at the top to most negative cosine at the bottom \n')
  if(radf[nrow(radf),'cos'] > 0)
    stop('\t \t \t \t >>> Error sampleEnrich: last drug does not have a negative cosine.
         sampleEnrich expects list of drugs ranked from most positive cosine at the top to most negative cosine at the bottom \n')
  
  
  ### do we have enough examples of testAnnotation to calculate statistics?
  # if not, do not waste time and stop now
  # number of instances of the annotation to be tested?
  nano <- length(which(radf[,annotationCol]==testAnnotation))
  # is that enough?
  if (nano < minNex) {
    
    # incredibly, cat() below was the source of a very challenging issue when using callr::r_bg to run long tasks in background
    # it would somehow cause process to go on forever!
    # cat('\t \t \t \t >>> Fewer than', minNex, 'examples, skip \n')
    
    statres <- data.frame(annotation=testAnnotation, nexamples=nano,
                          sumRanks=NA, sumRanksnorm=NA, bestPos=NA, sumRanksFracPos=NA,
                          sumRanksDir=NA, sumRanksDirnorm=NA,
                          ndraws=NA, nhigher=NA, pval=NA,
                          ksD=NA, kspval=NA )
    return(statres)
  }
  # if enough examples, will continue below
  
  ### calculate sum of ranks for given annotation
  # this gives a measure of "how distant from 0, on either side"
  # which type of rank we use depends on setting whichRank
  sr <- sum(radf[which(radf[,annotationCol]==testAnnotation), whichRank]) # sr for sum of ranks
  
  # ! say Schizophrenia has 50 drugs and Autism has 5 drugs
  # then we expect bigger scores for Schizophrenia regardless of whether it is actually more enriched
  # we can normalise the scores by dividing by the number of times that annotation was present
  srn <- sr / nano # srn for sum of ranks normalised
  # 20/08/2024: when summing cos, this will simply be average cosine of fingerprints with that annotation
  
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
  
  # 20/08/2024: I think can follow the same logic here when using abscos
  # instead of simply saying "best possible score is n x cos 1.0", we actually look at the data
  # and say "for an annotation with n fingerprints, the best possible sum abscos is this"
  
  ### calculate sum of ranks, including negative ranks
  # this will give a general direction to the enrichment
  # e.g. if all instances are on enriched in negative cosines, we get a negative score here
  # note, a perfect binomial is interesting / a distribution perfectly in the center around cos ~ 0 is not
  # but the two will give score ~ 0
  if(whichRank %in% c('rank', 'rankeq', 'rank0', 'ranks')) {
    srs <- sum(radf[which(radf[,annotationCol]==testAnnotation), 'ranks']) # s for signed
  } else if(whichRank %in% c('cos', 'abscos')) {
    srs <- sum(radf[which(radf[,annotationCol]==testAnnotation), 'cos'])
  }

  # we can also normalise this score
  srsn <- srs / nano # srsn for sum of ranks signed normalised
  # I think there is no point calculating the best possible score here
  # as = 0 or = very negative number or = very positive number are all potentially interesting
  
  
  ### SAMPLING PROCEDURE ###
  # if sufficient number of examples
  # draw number of instances and sum ranks
  
  # if we are given a maximum p-value, after how many more-extreme simulations should we stop?
  # e.g. if we are not interested in p-value above 0.1, there is no point spending time doing 100,000 simulations
  # we can stop as soon as we know that the p-value will be above 0.1
  if(!is.na(maxPval) & maxPval < 1) {
    maxnsim <- ndraws * maxPval
    # e.g. we will do *up to* 100,000 draws and the max p-value is 0.1
    # then as soon as 100,000 * 0.1 = 10,000 draws are more extreme, we know final p-value will be > 0.1
  } else {
    maxnsim <- ndraws
  }
  
  # counter
  # number of simulations which give a more extreme sum of ranks
  # (the more there are, the less evidence for a surprising enrichment there is)
  nsim <- 0
  for(nd in 1:ndraws) {
    
    # draw new `number of instances` ranks
    # and sum them
    dsr <- sum(sample(radf[,whichRank], size=nano)) # draw sum of ranks
    
    # is this sum of ranks more extreme than the real one?
    # if yes, count it
    if(dsr > sr) {
      nsim <- nsim + 1
    }
    
    # should we give up?
    if(nsim >= maxnsim) break
    
  }
  
  ### calculate p-value
  # how many simulations gave a more extreme sum of ranks?
  # approximate p-value is
  pv <- nsim / ndraws
  
  # note, we should not be creating any p-value = 0
  # instead, the minimum p-value depends on the number of draws
  # e.g. if we do N=1000 draws, then the minimum p-value is 1/1000 = 0.001
  if(nsim==0) {
    pv <- 1/ndraws
  }
  
  
  ### test enrichment with KS test
  # ranks where we find the annotation are
  rks <- radf[which(radf[,annotationCol]==testAnnotation), 'rank']
  # note, I think we do not want to use rank0 or rankeq here
  # rank0 is absolute distance from 0 so repeats in the negative
  # i.e. a large rank0 could mean very positive cos OR very negative cos
  # but KS test procedure involves sorting in increasing order, so will mix up positive cos and negative cos
  # between rank and ranks makes no difference (tested)
  ks <- ks.test( rks , radf$rank )
  kspval <- as.numeric(ks$p.value)
  
  
  
  ### return statistics report
  # (one row of the statistics report)
  statres <- data.frame(annotation=testAnnotation, nexamples=nano,
                        sumRanks=sr, sumRanksnorm=srn, bestPos=bsr, sumRanksFracPos=abs(sr/bsr),
                        sumRanksDir=srs, sumRanksDirnorm=srsn,
                        ndraws=ndraws, nhigher=nsim, pval=pv,
                        ksD=as.numeric(ks$statistic), kspval=kspval)
  
  return(statres)
  
}







# function signiBonferonni(...) -------------------------------------------

# given a bunch of p-values, return TRUE or FALSE or NA after checking against threshold Bonferonni-corrected

signiBonferonni <- function(pvs,
                            alphaThr) {
  ### Bonferonni correction of pvals
  # how many p-values did we generate? This is the number of hypotheses we tested
  npvs <- sum(!is.na(pvs))
  # hence, new alpha threshold is
  alpBon <- alphaThr/npvs
  cat('\n \t \t \t \t >>>', npvs, 'hypotheses tested \n')
  # add whether significant or not after Bonferonni correction
  cat('\t \t \t \t \t >>> Bonferonni-corrected alpha threshold:', alpBon,'\n')
  
  return(pvs < alpBon)
}




# function signiBenjaminiHochberg(...) ------------------------------------

signiBenjaminiHochberg <- function(pvs,
                                   alphaThr) {
  
  ### Benjamini-Hochberg procedure for correction
  cat('\t \t \t \t \t >>> Correcting p-values with Benjamini-Hochberg procedure \n')
  
  # how many p-values did we generate? This is the number of hypotheses we tested
  npvs <- sum(!is.na(pvs))
  
  # 1- rank pvals from smallest to largest
  pvord <- order(pvs) # ! record the order so we can return the TRUE/FALSE back in the same order at the end
  pvs <- pvs[pvord]
  
  # 2- for each p-value starting from the smallest
  # threshold is alpha * (rank / total number of tests)
  pvBH <- sapply(1:length(pvs), function(rk) {
    
    # rk because position corresponds to rank, as we sorted the p-values
    # calculate threshold
    alpBH <- alphaThr * (rk / npvs)
    # is the p-val significant or not?
    return(pvs[rk] < alpBH)
    
  })
  
  # 3- as soon as one test is not significant, all the rest is non-significant
  # correct this post-hoc
  
  # in rare cases, all p-values may be significant even after correction, i.e. below BH threshold
  # in which case, we have nothing to correct
  if(all(pvBH[!is.na(pvBH)])) { # if all (except NA) TRUE
    pvBH <- pvBH[order(pvord)] # put the p-values back the original order
    return(pvBH) # return statement closes the function without reading below
  }
  
  # in most cases, some p-values become above threshold (i.e. ns):
  ns1 <- which(!pvBH)[1] # first ns
  # there may be some NAs at the end, which we should leave as NA
  if(length(which(is.na(pvBH))) > 0) { # if there are some NAs at the end
    nsla <- which(is.na(pvBH))[1] - 1
  } else { # if not
    nsla <- length(pvBH)
  }
  # now replace all these by FALSE (i.e. ns)
  pvBH[ns1:nsla] <- FALSE
  
  # put the p-values back the original order
  pvBH <- pvBH[order(pvord)]
  
  return(pvBH)
}
