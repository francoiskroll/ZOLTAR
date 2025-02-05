#####################################################
# ~ predPharma Shiny app: small functions to clean up output tables ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# goal is to make tables clearer for the app, removing not-so-important information
# keep full table for when user downloads


# function cleanDrugsRanked(...) ------------------------------------------

cleanDrugsRanked <- function(vdbr) {
  
  ## remove day Z-scores
  vdbr <- vdbr[ , - which(startsWith(colnames(vdbr), 'day')) ]
  ## remove night Z-scores
  vdbr <- vdbr[ , - which(startsWith(colnames(vdbr), 'night')) ]
  
  ## remove pharma column
  vdbr$pharma <- NULL
  
  ## remove Complexity column
  vdbr$Complexity <- NULL
  
  ## remove rank column
  vdbr$rank <- NULL
  
  ## remove rank0 column
  vdbr$rank0 <- NULL
  
  ## remove cor column
  vdbr$cor <- NULL
  
  ## round cos column
  vdbr$cos <- round(vdbr$cos, digits=2)
  
  ## v6: not using rankeq anymore, so removing column
  # before was: vdbr$rankeq <- round(vdbr$rankeq, digits=2)
  vdbr$rankeq <- NULL
  
  ## re-arrange columns
  vdbr <- vdbr[, c('cos', 'name', 'cleanm', 'cid', 'tid', 'ranks', 'library', 'concentration', 'MolecularWeight', 'inSmall', 'structCluster')]
  
  ## rename columns
  colnames(vdbr) <- c('Cosine', 'Original name', 'Name', 'PubChem CID', 'TTD ID',
                      'Rank from 0', 'Library', 'Concentration', 'Molecular weight', 'Shortlisted', 'Structural cluster')
  
  ## return
  return(vdbr)
}


# function cleanIndications(...) ------------------------------------------

cleanIndications <- function(ind) {
  
  ## remove sumRanksnorm
  ind$sumRanksnorm <- NULL
  
  ## remove sumRanksDir & sumRanksDirnorm
  ind$sumRanksDir <- NULL
  ind$sumRanksDirnorm <- NULL
  
  # make sure ndraws column is displayed as an integer
  ind$ndraws <- as.integer(ind$ndraws)
  ind$nhigher <- as.integer(ind$nhigher)
  
  # round some columns
  ind$sumRanks <- round(ind$sumRanks, digits=2)
  ind$bestPos <- round(ind$bestPos, digits=2)
  ind$sumRanksFracPos <- round(ind$sumRanksFracPos, digits=2)
  ind$ksD <- round(ind$ksD, digits=2)
  ind$kspval <- round(ind$kspval, digits=2)
  
  ## re-name columns
  colnames(ind) <- c('Indication', 'N drugs', 'Sum of cosines', 'Best possible sum of cosines', 'Fraction of best possible', 'N draws', 'N higher',
                     'pval', 'Bon. sign.', 'Ben. sign.', 'KS D', 'KS pval', 'KS Bon. sign.', 'KS Ben. sign.')
  
  ## return
  return(ind)
}


# function cleanTTDtargets(...) -------------------------------------------

cleanTTDtargets <- function(tar) {
  
  ## remove sumRanksnorm
  tar$sumRanksnorm <- NULL
  
  ## remove sumRanksDir & sumRanksDirnorm
  tar$sumRanksDir <- NULL
  tar$sumRanksDirnorm <- NULL
  
  # make sure ndraws & nhigher columns are displayed as an integer
  tar$ndraws <- as.integer(tar$ndraws)
  tar$nhigher <- as.integer(tar$nhigher)
  
  tar$sumRanks <- round(tar$sumRanks, digits=2)
  tar$bestPos <- round(tar$bestPos, digits=2)
  tar$sumRanksFracPos <- round(tar$sumRanksFracPos, digits=2)
  tar$ksD <- round(tar$ksD, digits=2)
  tar$kspval <- round(tar$kspval, digits=2)
  
  # re-arrange columns
  tar <- tar[, c('annotation', 'TARGNAME', 'GENENAME', 'UNIPROID', 'BIOCLASS', 'TARGTYPE',
                 'nexamples', 'sumRanks', 'bestPos', 'sumRanksFracPos', 'ndraws', 'nhigher', 'pval', 'bonferonni', 'benhoch',
                 'ksD', 'kspval', 'KSbonferonni', 'KSbenhoch')]
  
  ## re-name columns
  colnames(tar) <- c('TTD ID', 'Target', 'Gene (human)', 'UniProt name', 'Bioclass', 'Status',
                     'N drugs', 'Sum of cosines', 'Best possible sum of cosines', 'Fraction of best possible', 'N draws', 'N higher',
                     'pval', 'Bon. sign.', 'Ben. sign.', 'KS D', 'KS pval', 'KS Bon. sign.', 'KS Ben. sign.')
  
  ## return
  return(tar)
  
}


# function cleanKEGG(...) -------------------------------------------------

cleanKEGG <- function(keg) {
  
  ## remove sumRanksnorm
  keg$sumRanksnorm <- NULL
  
  ## remove sumRanksDir & sumRanksDirnorm
  keg$sumRanksDir <- NULL
  keg$sumRanksDirnorm <- NULL
  
  # make sure ndraws & nhigher columns are displayed as an integer
  keg$ndraws <- as.integer(keg$ndraws)
  keg$nhigher <- as.integer(keg$nhigher)
  
  keg$sumRanks <- round(keg$sumRanks, digits=2)
  keg$bestPos <- round(keg$bestPos, digits=2)
  keg$sumRanksFracPos <- round(keg$sumRanksFracPos, digits=2)
  keg$ksD <- round(keg$ksD, digits=2)
  keg$kspval <- round(keg$kspval, digits=2)
  
  ## re-name columns
  colnames(keg) <- c('ID', 'KEGG pathway', 'N drug-target', 'Sum of cosines', 'Best possible sum of cosines', 'Fraction of best possible',
                     'N draws', 'N higher', 'pval', 'Bon. sign.', 'Ben. sign.', 'KS D', 'KS pval', 'KS Bon. sign.', 'KS Ben. sign.')
  
  ## return
  return(keg)
}
