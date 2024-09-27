# splitStrTakeNth ---------------------------------------------------------
# split string(s) and take the nth element of each
# e.g. strNthSplit on day1_sunny_23C and day2_rainy_20C
# with split = '_' and nth = 2
# gives sunny, rainy

strNthSplit <- function(stri,
                        split,
                        nth) {
  
  # confirm we are given string(s)
  stri <- as.character(unlist(stri))
  
  as.character(sapply(strsplit(stri, split=split),
                      function(s) {
                        s[nth]
                      }))
}

# function gglegacyFingerprint(...) ---------------------------------------

# expects lFgp (legacy fingerprint) as dataframe created by legacyFingerprint(...)

gglegacyFingerprint <- function(lFgp,
                                onlyGrp=NA,
                                colours=NA,
                                legendOrNo=TRUE,
                                ynameOrNo=TRUE,
                                ytextOrNo=TRUE,
                                xtextOrNo=TRUE,
                                yTitleSize=9,
                                paramNumOrNo=FALSE,
                                paramSize=7,
                                ynumSize=7,
                                nightBgOrNo=TRUE,
                                ymin,
                                ymax,
                                exportOrNo=FALSE,
                                exportPath,
                                width=150,
                                height=100) {
  
  
  ### about settings ###
  
  # if xtextOrNo is FALSE, make sure paramNumOrNo is also set to FALSE to avoid error
  if(!xtextOrNo) {
    paramNumOrNo <- FALSE
  }
  
  
  ### set order of uparam ###
  
  # order which makes most sense to me is:
  paramOrder <- c('averageActivity', 'averageWaking', 'sleep', 'sleepBout', 'sleepLength', 'sleepLatency')
  # set those as levels of parameters
  lFgp$parameter <- factor(lFgp$parameter, levels=paramOrder)
  # re-order the fingerprint following those levels
  lFgp <- lFgp[order(lFgp$parameter),]
  
  # now, order is:
  # parameter1 : night1, night2, day1, day2
  # parameter2 : ..
  # but will look nicer to do "mini-fingerprints",
  # i.e. all day1 parameters, followed by all day2 parameters, etc.
  # set desired order of windows as levels
  lFgp$win <- factor(lFgp$win, levels=c('day1', 'day2', 'night1', 'night2'))
  
  # now we order first by win, then by parameter to get desired order
  lFgp <- lFgp[order(lFgp$win, lFgp$parameter),]
  
  # set uparam column as it is as levels so sure gets plotted this way
  lFgp$uparam <- factor(lFgp$uparam, levels=lFgp$uparam)
  
  
  ### if we have a column (grp) that simply says 'zsco'
  # replace by 'mean fgp' just for aesthetics in the plot
  if('zsco' %in% colnames(lFgp)) {
    colnames(lFgp)[which(colnames(lFgp)=='zsco')] <- 'mean fgp'
  }
  
  
  # below assumes format is e.g. : uparam / win / parameter / fgp1 / fgp2 / fgp3
  # fgp1, fgp2, fgp3 being columns of Z-scores
  lFgpl <- lFgp %>%
    pivot_longer(-c(uparam, win, parameter),
                 names_to='grp',
                 values_to='zsco')
  
  
  ### keep only groups we are plotting ###
  # if not plotting every group:
  if(!is.na(onlyGrp[1])) {
    lFgpl <- lFgpl %>%
      filter(grp %in% onlyGrp)
  }
  

  ### how to connect datapoints in plot ###
  # we want to connect by grp_win:
  lFgpl <- lFgpl %>%
    mutate(grp_win=paste(grp, win, sep='_'), .after='parameter')


  ### colours ###
  # if user did not provide any colours, we use automatic ggplot colours
  if(is.na(colours[1])) {
    colours <- scales::hue_pal()(length(unique(lFgpl$grp)))
  }

  ### set the axes ###
  # set y axis name
  yname <- 'deviation from controls\n(z-score)'

  # below: find last 'day' parameter
  # we want the grey frame to start just after (+0.5)
  # take unique uparam
  xmid <- max(which(startsWith(as.character(unique(lFgpl$uparam)), 'day'))) + 0.5
  

  ### plot ###
  ggFgp <- ggplot(lFgpl, aes(x=uparam, y=zsco, colour=grp, group=grp_win)) +
    geom_hline(yintercept=0, linetype=1, colour='#a7a7a7', linewidth=0.5) +
    geom_point(size=1.5) +
    geom_line() +
    {if(nightBgOrNo) annotate(geom='rect', xmin=xmid, xmax=Inf, ymin=-Inf, ymax=Inf, colour=NA, fill='#1d1d1b', alpha=0.2)} +
    scale_colour_manual(values=colours) +
    theme_minimal() +
    theme(
      panel.grid.minor.y=element_blank(),
      axis.title.x=element_blank(),
      legend.title=element_blank(),
      legend.text=element_text(size=12)) +

    {if(!legendOrNo) theme(legend.position='none')} +

    {if(!ynameOrNo) theme(axis.title.y=element_blank())} +
    {if(ynameOrNo) theme(axis.title.y=element_text(size=yTitleSize, margin=margin(t=0, r=-1.2, b=0, l=0)))} +
    {if(ynameOrNo) ylab(yname)} +

    {if(!ytextOrNo) theme(axis.text.y=element_blank())} +
    {if(ytextOrNo) theme(axis.text.y=element_text(size=ynumSize))} +
    
    # add X axis labels (uparam)
    {if(xtextOrNo & !paramNumOrNo) theme(axis.text.x=element_text(size=paramSize, angle=45, hjust=1))} +
    
    # or X axis labels as parameter numbers
    {if(paramNumOrNo) scale_x_discrete(labels=match(lFgp$parameter, paramOrder))} +
    {if(paramNumOrNo) theme(axis.text.x=element_text(size=paramSize, margin=margin(t=-1.5)))} +
    
    # or turn off X axis labels completely
    {if(!xtextOrNo) theme(axis.text.x=element_blank())} +

    coord_cartesian(ylim=c(ymin, ymax))


  ### export the plot ###
  if(exportOrNo) {
    ggsave(exportPath, ggFgp, width=width, height=height, units='mm', device=cairo_pdf)
  }

  ### return plot ###
  # this way it will display in RStudio
  return(ggFgp)
}


# function ggDrugFgp(...) -------------------------------------------------

# wrapper for ggLegacyFingerprint, but allows to plot drug fingerprints, with or without fingerprint from an experiment (e.g. a knockout)
# note, will lose error bar on fingerprint because drug database does not have the standard deviation

# dnames can be drug name(s) or 'first' or 'last' (assumes you give it ranked drugDb)

ggDrugFgp <- function(drugDb,
                      dnames,
                      legacyFgp=NULL,
                      onlyGrp=NA,
                      colours,
                      legendOrNo,
                      ynameOrNo,
                      ytextOrNo,
                      xtextOrNo,
                      paramNumOrNo,
                      nightBgOrNo,
                      ymin,
                      ymax,
                      exportOrNo,
                      exportPath,
                      width,
                      height) {
  
  
  ### import drugDb
  # if we are given a path
  if(is.character(drugDb)) {
    ddb <- read.csv(drugDb)
  } else { # if not, assume we are given object directly
    ddb <- drugDb
  }
  
  
  ## 05/04/2023, decided when developing predPharam Shiny app to delete parameters daymean_averageWaking & nightmean_averageWaking
  # see rationale in paramsFromMid.R
  ddb$daymean_averageWaking <- NULL
  ddb$nightmean_averageWaking <- NULL
  
  # detect column where Z-scores start
  zcol <- min( which(startsWith(colnames(ddb), c('day', 'night'))) )
  # detect column where Z-scores end
  zcoll <- max( which(startsWith(colnames(ddb), c('day', 'night'))) )
  # (this way if we add columns it does not break the function)
  
  # keep only drugs we want to plot
  # use a loop so we can return some information to user
  rows2take <- as.vector(unlist(sapply(dnames, function(dnm) {
    
    if(dnm=='first') {
      rws <- 1
    } else if(dnm=='last') {
      rws <- nrow(ddb)
    } else {
      rws <- which(ddb$name == dnm)
      if (length(rws)==0) {
        cat('\t \t \t \t >>> No entry for', dnm, '\n')
      } else if (length(rws)>0) {
        cat('\t \t \t \t >>> n =', length(rws), 'entries for', dnm, '\n')
      }
    }
    
    return(rws)
    
  }
  )))
  
  
  # take these rows from drugDb
  dbn <- as.data.frame(t(ddb[rows2take, zcol:zcoll]))
  colnames(dbn) <- ddb[rows2take, 'name']
  
  # will cause an issue if there are duplicate column names
  # which can happen if two fingerprints with the same drug name
  # I think it used to do it automatically below, e.g.   # e.g. Control / Control becomes Control / Control.1
  # but now causes an error
  colnames(dbn) <- make.unique(colnames(dbn))
  
  # match formatting that ggLegacyFingerprint expects
  dbn <- dbn %>%
    add_column(uparam=row.names(dbn), .before=1) %>%
    mutate(parameter=strNthSplit(uparam, '_', 2), .after='uparam') %>%
    mutate(win=strNthSplit(uparam, '_', 1), .after='uparam')
  row.names(dbn) <- NULL
  
  ### keep only Z-scores from legacy fingerprint(s)
  if (!is.null(legacyFgp[1])) {
    
    # note, can work for multiple fingerprints
    # assuming column names of Z-score columns give the group
    
    # first, confirm legacyFgp and drugDb are giving the same parameters
    if(! identical( sort(dbn$uparam) , sort(legacyFgp$uparam )))
      stop('\t \t \t \t >>> Error ggDrugFgp: drugDb and legacyFingerprint do not give the same parameters. \n')
    
    # add fingerprints we find in legacyFgp
    # ! will assume it is every columns not named uparam, win, parameter
    # (we keep uparam for the join below)
    cols2take <- which(! colnames(legacyFgp) %in% c('win', 'parameter'))
    dbn <- right_join(dbn, legacyFgp[cols2take], by='uparam') # join so will work even if parameters not giving in exactly the same order
    # (we checked above all parameters are present)
    
  }
  
  
  ### give cosines of fingerprints comparisons
  fsim <- as.data.frame(lsa::cosine(as.matrix(dbn[4:ncol(dbn)])))
  print(fsim)
  
  
  ### send to gglegacyFingerprint for plotting
  return( gglegacyFingerprint(lFgp=dbn,
                              onlyGrp=onlyGrp,
                              colours=colours,
                              legendOrNo=legendOrNo,
                              ynameOrNo=ynameOrNo,
                              ytextOrNo=ytextOrNo,
                              xtextOrNo=xtextOrNo,
                              paramNumOrNo=paramNumOrNo,
                              nightBgOrNo=nightBgOrNo,
                              ymin=ymin,
                              ymax=ymax,
                              exportOrNo=exportOrNo,
                              exportPath=exportPath,
                              width=width,
                              height=height) )
  
}
