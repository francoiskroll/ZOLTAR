#####################################################
# ~ Shiny app for predictive pharmacology ~
#
# "from behavioural fingerprint to drugs and pathways"
# 
# Francois Kroll 2023
# francois@kroll.be
#####################################################


# packages ----------------------------------------------------------------

library(shiny)
library(tidyr)
library(ggplot2)
library(bslib)
library(DT)
library(tools)

library(callr)
library(shinyjs)


# functions ---------------------------------------------------------------

source('legacyFingerprint.R')
source('gglegacyFingerprint.R')
source('drawEnrich_v5.R')
source('paramsFromMid.R')
source('cleanTables.R')
source('ggEnrich.R')


# settings ----------------------------------------------------------------

Sys.setlocale("LC_ALL","C") # avoids an issue when printing table of ranked drugs, probably because of odd characters in original drug names
# solution StackOverflow question 61656119

ndraws <- 10000
# ndraws <- 2

# set maximum upload to 100 Mb
options(shiny.maxRequestSize = 100*1024^2)
# one middur file is ~ 2 Mb, so 100 Mb should be sufficient for ~ 100 experiments

# seems fairly easy to save uploaded data to Dropbox, if useful at some point
# https://www.r-bloggers.com/2015/07/persistent-data-storage-in-shiny-apps/



# chunks of code for dropdown menus ---------------------------------------

js <- "
$(document).ready(function() {

    var j = 1;

    $('#icon').on('click', function() {
        j = j + 1;
        if (j % 2 === 0) {
            Shiny.setInputValue('show_panel', true);
            $('#icon')[0].innerHTML = $('#icon')[0].innerHTML.replaceAll('-right', '-down');
        } else {
            Shiny.setInputValue('show_panel', false);
            $('#icon')[0].innerHTML = $('#icon')[0].innerHTML.replaceAll('-down', '-right');
        }
    });

    $('#icon').hover(function() {
        $(this).css('cursor', 'pointer');
    });

});
"


# ui ----------------------------------------------------------------------


ui <- fluidPage(
  
  #theme=bslib::bs_theme(bootswatch='sandstone'),
  #theme=bslib::bs_theme(bootswatch='united'), # I quite like it but maybe the font is too childish
  #theme=bslib::bs_theme(bootswatch='simplex'),
  
  # can also set theme manually
  # theme=bslib::bs_theme(
  #   base_font='Rubik'
  # ),
  
  title='ZOLTAR',
  
  tags$head(tags$link(rel="icon", type="image/png", href="favicon.png")),
  
  ### app title
  img(src='zoltar.png', align='top', height=117.1, width=471),
  
  # titlePanel('ZOLTAR'),
  # subtitle
  # p(em('"From a Zebrafish behaviOuraL fingerprinT, I will predict candidate Therapeutics and cAusal pRocesses"'), style='color:grey'),
  
  ### sidebar layout = inputs on left, outputs on right
  sidebarLayout(
    
    ### sidebar
    sidebarPanel(width=4,
      
      ## middur input widget
      fileInput(inputId='mid_drop',
                label='Select or drop your middur.csv file(s)',
                multiple=TRUE,
                accept='.csv') ,
      
      ## genotype input widget
      fileInput(inputId='geno_drop',
                label='Select or drop your genotype.txt file(s)',
                multiple=TRUE,
                accept='.txt') ,
      
      ## treatment group dropdown selection
      selectInput(inputId='treGrp_select',
                  label='Select the treatment group',
                  choices=NULL, # initialise as NULL, we will update based on data given
                  multiple=FALSE) , # can only select one group
      
      ## control group dropdown selection
      selectInput(inputId='conGrp_select',
                  label='Select the control group',
                  choices=NULL, # initialise as NULL, we will update based on data given
                  multiple=FALSE) , # can only select one group
      
      ## GO button
      actionButton(inputId='go_button',
                   label='Go'),
      
      p(''),
      
      ## link to tutorial
      p(a(href="https://github.com/francoiskroll/ZOLTAR#readme",
          HTML("Read instructions"),
          target='_blank')),
      
      ## download sample data
      # downloadLink('mid14_dl', 'Download sample data')
      downloadLink('zip_dl', 'Download sample data'),
      
      ## cite us!
      p(''),
      p('Cite us!'),
      
      ## advanced settings panel
      # from https://stackoverflow.com/questions/77387537/how-can-i-implement-a-conditional-panel-revealed-by-clicking-an-arrowhead/77390635#77390635
      #HTML('<br>'),
      tags$script(HTML(js)),
      
      h5(div(id = "icon", icon("angle-right", style = "color: #7d3040;"),
             "Advanced settings")),
      
      conditionalPanel(condition = "input.show_panel === true",

                       fixedRow(
                         
                         column(5, numericInput('day1_start', 'day1 start', value=24)),
                         column(5, numericInput('day1_stop', 'day1 stop', value=38)),
                         
                         column(5, numericInput('night1_start', 'night1 start', value=38)),
                         column(5, numericInput('night1_stop', 'night1 stop', value=48)),
                         
                         column(5, numericInput('day2_start', 'day2 start', value=48)),
                         column(5, numericInput('day2_stop', 'day2 stop', value=62)),
                         
                         column(5, numericInput('night2_start', 'night2 start', value=62)),
                         column(5, numericInput('night2_stop', 'night2 stop', value=72)),
                       )),
      
      ## need help?
      HTML('<br>'), # empty line
      p('Need help? francois@kroll.be', style='color:#999999')
      
    ),
    
    ### tabs on the right
    mainPanel(
      tabsetPanel(
        
        # tabPanel('Table',
        #          tableOutput('fgp')),
        
        tabPanel('Fingerprint',
                 p(''),
                 downloadButton('ggfgp_dl', 'download pdf'),
                 p(''),
                 plotOutput('ggfgp')),
        
        # tabPanel('Drugs ranked',
        #          p('Drugs are ranked by decreasing cosine.'),
        #          downloadButton('vdbr_dl', 'download'),
        #          h3(paste('Top', showNdrugs)),
        #          tableOutput('topdr'), # topdr is for top X drugs
        #          h3(paste('Bottom', showNdrugs)),
        #          tableOutput('botdr')), # botdr is for bottom X drugs
        
        tabPanel('Drugs ranked',
                 p(''),
                 ## show details panel
                 checkboxInput('showdetails_Drugs', 'Show details', value=FALSE),
                 
                 # conditional panel to show/hide details
                 conditionalPanel(
                   condition = 'input.showdetails_Drugs == true',
                   p('• Click on Cosine column to rank by increasing or decreasing'),
                   p('• Click on a row to reveal fingerprint plot.'),
                   p(''),
                   p(strong(em('• Cosine')), 'cosine similarity score (min. −1, max. +1) between the small molecule fingerprint and the query fingerprint.'),
                   p(strong(em('• Original name')), 'original compound name from Rihel et al. 2010 data.'),
                   p(strong(em('• Name')), 'simplified name.'),
                   p(strong(em('• PubChem CID')), 'PubChem compound ID'),
                   p(strong(em('• TTD ID')), 'Therapeutic Target Database compound ID'),
                   p(strong(em('• Rank from 0')), 'number of random draws which gave a higher sum of ranks than the one observed (Sum of ranks).'),
                   p(strong(em('• Rank eq.')), '= N higher/N draws. Smallest possible p-value is 0.0001, which corresponds to 1 or 0 out of 10,000 random draws giving a more extreme sum of ranks than the observed one. The real p-value is ≤ 0.0001.'),
                   p(strong(em('• Library')), 'source library of the compound.'),
                   p(strong(em('• Concentration')), 'approximate concentration at which the compound was tested.'),
                   p(strong(em('• Molecular weight')), 'compound\'s molecular weight in g/mol.'),
                   p(strong(em('• Shortlisted')), 'whether this compound was labelled as "shortlisted" by Rihel et al., 2010. A compound was shortlisted if it affected at least one behavioural parameter with a large effect size and/or affected the same parameter in the same direction across the two days/nights.'),
                   p(strong(em('• Structural cluster')), 'clusters of compounds with similar structures based on Tanimoto score, see Rihel et al., 2010.'),
                   p('• Data from', a(href='https://www.science.org/doi/abs/10.1126/science.1183090', HTML('Rihel et al., 2010. <em>Science</em>',)), ', reorganised by Kroll et al., 2023.')
                 ), # closes conditional panel
                 downloadButton('vdbr_dl', 'download'),
                 p(''),
                 DTOutput('vdbr_dis')),
        
        tabPanel('Indications',
                 p(''),
                 
                 ## show details panel
                 checkboxInput('showdetails_Indications', 'Show details', value=FALSE),
                 
                 # conditional panel to show/hide details
                 conditionalPanel(
                   condition = 'input.showdetails_Indications == true',
                   p(strong(em('• N examples')), 'number of fingerprints annotated with this Indication. This can include replicate experiments with the same compound.'),
                   p(strong(em('• Sum of ranks')), 'sum of the ranks of the N fingerprints with this Indication.'),
                   p(strong(em('• Best possible sum of ranks')), 'sum of ranks if the N fingerprints (N examples) were at the most extreme positions. For example, if N examples = 4, it would represent the sum of ranks if the fingerprints with this Indication had been at the following positions: top 1 (maximum positive cosine), top 2, before last, last (maximum negative cosine).'),
                   p(strong(em('• Fraction of best possible')), '= Sum of ranks/Best possible sum of ranks, to give a measure of enrichment that is somewhat normalised, i.e. comparable between Indications which have different N examples.'),
                   p(strong(em('• N draws')), 'number of random draws, this is always 10,000.'),
                   p(strong(em('• N higher')), 'number of random draws which gave a higher sum of ranks than the one observed (Sum of ranks).'),
                   p(strong(em('• pval')), '= N higher/N draws. Smallest possible p-value is 0.0001, which corresponds to 1 or 0 out of 10,000 random draws giving a more extreme sum of ranks than the observed one. The real p-value is ≤ 0.0001.'),
                   p(strong(em('• Bon. sign.')), 'whether the p-value remains significant after Bonferroni correction. Alpha threshold is set at 0.05.'),
                   p(strong(em('• Ben. sign.')), 'whether the p-value remains significant after Benjamini-Hochberg correction. Alpha threshold is set at 0.05.'),
                   p(strong(em('• KS D')), 'D statistic of the Kolmogorov–Smirnov (KS) test.'),
                   p(strong(em('• KS pval')), 'p-value statistic of the Kolmogorov–Smirnov (KS) test.'),
                   p(strong(em('• KS Bon. sign.')), 'whether the KS p-value remains significant after Bonferroni correction.'),
                   p(strong(em('• KS Ben. sign.')), 'whether the KS p-value remains significant after Benjamini-Hochberg correction.'),
                   p('• Source of annotations: Therapeutic Target Database'),
                 ), # closes conditional panel
                 downloadButton('ind_dl', 'download'),
                 p(''),
                 uiOutput('ind_waitMessage'),
                 DTOutput('ind_dis')),
        
        tabPanel('Targets',
                 p(''),
                 ## show details panel
                 checkboxInput('showdetails_Targets', 'Show details', value=FALSE),
                 
                 # conditional panel to show/hide details
                 conditionalPanel(
                   condition = 'input.showdetails_Targets == true',
                   p(strong(em('• N examples')), 'number of fingerprints annotated with this Target protein. This can include replicate experiments with the same compound.'),
                   p(strong(em('• Sum of ranks')), 'sum of the ranks of the N fingerprints with this Target protein'),
                   p(strong(em('• Best possible sum of ranks')), 'sum of ranks if the N fingerprints (N examples) were at the most extreme positions. For example, if N examples = 4, it would represent the sum of ranks if the fingerprints with this Target had been at the following positions: top 1 (maximum positive cosine), top 2, before last, last (maximum negative cosine).'),
                   p(strong(em('• Fraction of best possible')), '= Sum of ranks/Best possible sum of ranks, to give a measure of enrichment that is somewhat normalised, i.e. comparable between Targets which have different N examples.'),
                   p(strong(em('• N draws')), 'number of random draws, this is always 10,000.'),
                   p(strong(em('• N higher')), 'number of random draws which gave a higher sum of ranks than the one observed (Sum of ranks).'),
                   p(strong(em('• pval')), '= N higher/N draws. Smallest possible p-value is 0.0001, which corresponds to 1 or 0 out of 10,000 random draws giving a more extreme sum of ranks than the observed one. The real p-value is ≤ 0.0001.'),
                   p(strong(em('• Bon. sign.')), 'whether the p-value remains significant after Bonferroni correction. Alpha threshold is set at 0.05.'),
                   p(strong(em('• Ben. sign.')), 'whether the p-value remains significant after Benjamini-Hochberg correction. Alpha threshold is set at 0.05.'),
                   p(strong(em('• KS D')), 'D statistic of the Kolmogorov–Smirnov (KS) test.'),
                   p(strong(em('• KS pval')), 'p-value statistic of the Kolmogorov–Smirnov (KS) test.'),
                   p(strong(em('• KS Bon. sign.')), 'whether the KS p-value remains significant after Bonferroni correction.'),
                   p(strong(em('• KS Ben. sign.')), 'whether the KS p-value remains significant after Benjamini-Hochberg correction.'),
                   p('• Source of annotations: Therapeutic Target Database'),
                 ), # closes conditional panel
                 downloadButton('ttar_dl', 'download'),
                 p(''),
                 uiOutput('tar_waitMessage'),
                 DTOutput('ttar_dis')),
        
        tabPanel('KEGG pathways',
                 p(''),
                 ## show details panel
                 checkboxInput('showdetails_KEGG', 'Show details', value=FALSE),
                 
                 # conditional panel to show/hide details
                 conditionalPanel(
                   condition = 'input.showdetails_KEGG == true',
                   p(strong(em('• N examples')), 'number of fingerprints annotated with this KEGG pathway. This can include replicate experiments with the same compound.'),
                   p(strong(em('• Sum of ranks')), 'sum of the ranks of the N fingerprints with this KEGG pathway'),
                   p(strong(em('• Best possible sum of ranks')), 'sum of ranks if the N fingerprints (N examples) were at the most extreme positions. For example, if N examples = 4, it would represent the sum of ranks if the fingerprints with this KEGG pathway had been at the following positions: top 1 (maximum positive cosine), top 2, before last, last (maximum negative cosine).'),
                   p(strong(em('• Fraction of best possible')), '= Sum of ranks/Best possible sum of ranks, to give a measure of enrichment that is somewhat normalised, i.e. comparable between KEGG pathways which have different N examples.'),
                   p(strong(em('• N draws')), 'number of random draws, this is always 10,000.'),
                   p(strong(em('• N higher')), 'number of random draws which gave a higher sum of ranks than the one observed (Sum of ranks).'),
                   p(strong(em('• pval')), '= N higher/N draws. Smallest possible p-value is 0.0001, which corresponds to 1 or 0 out of 10,000 random draws giving a more extreme sum of ranks than the observed one. The real p-value is ≤ 0.0001.'),
                   p(strong(em('• Bon. sign.')), 'whether the p-value remains significant after Bonferroni correction. Alpha threshold is set at 0.05.'),
                   p(strong(em('• Ben. sign.')), 'whether the p-value remains significant after Benjamini-Hochberg correction. Alpha threshold is set at 0.05.'),
                   p(strong(em('• KS D')), 'D statistic of the Kolmogorov–Smirnov (KS) test.'),
                   p(strong(em('• KS pval')), 'p-value statistic of the Kolmogorov–Smirnov (KS) test.'),
                   p(strong(em('• KS Bon. sign.')), 'whether the KS p-value remains significant after Bonferroni correction.'),
                   p(strong(em('• KS Ben. sign.')), 'whether the KS p-value remains significant after Benjamini-Hochberg correction.'),
                   p('• Source of annotations: Therapeutic Target Database'),
                 ), # closes conditional panel
                 downloadButton('keg_dl', 'download'),
                 p(''),
                 uiOutput('keg_waitMessage'),
                 DTOutput('keg_dis')),
        
        # 04/04/2023 switching off STITCH for now, needs a bit more work
        # tabPanel('zebrafish STITCH',
        #          p('Source: stitch.embl.de'),
        #          tableOutput('zsti')) ,
        
        # 04/04/2023 switching off STITCH for now, needs a bit more work
        # tabPanel('human STITCH',
        #          p('Source: stitch.embl.de'),
        #          tableOutput('hsti'))
      )
    )
    
  ),
)


# server ------------------------------------------------------------------

# about progress bars:
# currently a bit of a cheap trick: it is started at 0.3 then updated to 1.0 when done
# at least user knows something is going on
# it seems like using a Progress object could make it more precise but requires more work

server <- function(input, output, session) {
  
  ### stop app if tab is closed
  # otherwise, if user stops ongoing run, they cannot open a new one
  # session$onSessionEnded(stopApp)
  
  #### set-up download of sample data
  output$zip_dl <- downloadHandler( # download TTD targets
    filename=function() {
      'sampledata.zip'
    },
    content=function(file) {
      file.copy('sampledata.zip', file)
    }
  )
  
  
  #### update the group selection when user drops the genotype file ####
  # note input$geno_drop is NULL initially
  observeEvent(input$geno_drop, # means: observe the genotype file import widget, if it changes, trigger the code below
               {
                 # remember, user may have given multiple genotype files
                 # so we loop through paths given
                 grpL <- lapply( 1:length(input$geno_drop$datapath) , function(i) {
                   # import the genotype file
                   geno <- importGenotype(input$geno_drop$datapath[i])
                   # get the group names from it
                   grpnms <- colnames(geno)
                   return(grpnms)
                 })
                 # so we have a small list where each slot is a small vector with the groups in that genotype file
                 # as treatment/control groups, user can only choose groups which are present in both genotype files
                 # i.e. the intersection of those character vectors
                 if(length(grpL)>1) {
                   grpnms <- do.call(intersect, grpL)
                 } else {
                   grpnms <- grpL[[1]]
                 }
                 
                 # check that there are more than one group available
                 if(length(grpnms)<2) stop('\t \t \t \t >>> Error: please have at least 2 groups present in all genotype files. \n')
                 
                 # update the selection of treatment / control group so it lists the groups available in the data
                 # selected is the 'pre-fill' option
                 updateSelectInput(session, input='treGrp_select',
                                   choices=grpnms, selected=grpnms[1])
                 updateSelectInput(session, input='conGrp_select',
                                   choices=grpnms, selected=grpnms[2])
               })
  
  
  #### start stuff when user presses Go button ####
  observeEvent(input$go_button,
               {
                 
                 disable('go_button') # supposed to disable go button, but I do not see a difference
                 
                 req(input$mid_drop) # check that middur file is available
                 req(input$geno_drop) # check that genotype file is available
                 
                 
                 ### match middur files & genotype files ###
                 # we will match them by YYMMDD_BX
                 # first, check user gave the same number of middur files and genotype files
                 if(length(input$mid_drop$datapath) != length(input$geno_drop$datapath))
                   stop('\t \t \t \t >>> Error: please give the same number of middur files and genotype files. Each middur file should have its own genotype file. \n')
                 
                 # we will order of middur files are correct order
                 # and sort the genotype files in the same order
                 # get the YYMMDD_BX of all the middur files
                 midymdb <- sapply(1:length(input$mid_drop$datapath), function(i) {
                   filnm <- tools::file_path_sans_ext(input$mid_drop$name[i])
                   return( substr(filnm, 1, 9) )
                 })
                 # check no duplicates
                 if(sum(duplicated(midymdb))>0)
                   stop('\t \t \t \t >>> Error: some of the middur files given have the same YYMMDD_BX. Please edit the file name and start again. \n')
                   
                 # get the YYMMDD_BX of all the genotype files
                 genoymdb <- sapply(1:length(input$mid_drop$datapath), function(i) {
                   filnm <- tools::file_path_sans_ext(input$geno_drop$name[i])
                   return( substr(filnm, 1, 9) )
                 })
                 # check no duplicates
                 if(sum(duplicated(genoymdb))>0)
                   stop('\t \t \t \t >>> Error: some of the genotype files given have the same YYMMDD_BX. Please edit the file name and start again. \n')
                 
                 # now match each middur file with its genotype file
                 mgmatch <- match(midymdb, genoymdb)
                 # this gives the order of the genotype files
                 # any NA is the sign that one middur file is unmatched or one genotype file is unmatched
                 if(sum(is.na(mgmatch))>0)
                   stop('\t \t \t \t Error: please check carefully that each middur file YYMMDD_BX matches one genotype file YYMMDD_BX. \n')
                 
                 # so prepare now the genotype paths in the right order
                 genopaths <- input$geno_drop$datapath[mgmatch]
                 # also record their names in the same order
                 genonms <- input$geno_drop$name[mgmatch]
                 
                 # check
                 sapply(1:length(genopaths), function(i) {
                   cat('\t \t \t \t >>> Middur file', input$mid_drop$name[i], 'matched with genotype file', genonms[i],'\n')
                 })
                 
                 
                 ### calculate fingerprints ###
                 # remember, user may have given multiple middur files
                 # so we loop through paths given
                 fgpL <- lapply( 1:length(input$mid_drop$datapath) , function(i) {
                   
                   # import csv (which assumes it is middur file generated by FramebyFrame)
                   if(file_ext(input$mid_drop$datapath[i])=='csv') {
                     mid <- read.csv(input$mid_drop$datapath[i])
                     
                     # or txt (which assumes it is DATA file generated by Jason's code)
                   } else if(file_ext(input$mid_drop$datapath[i])=='txt') {
                     # if then we call function DATAtoMid which imports and converts format to middur
                     mid <- DATAtoMid(DATApath=input$mid_drop$datapath[i])
                   }
                   
                   # prepare the fingerprint
                   fgp <- legacyFingerprintMid(mid=mid,
                                               genopath=genopaths[i],
                                               treGrp=input$treGrp_select,
                                               conGrp=input$conGrp_select,
                                               nights=c('night1', 'night2'),
                                               days=c('day1', 'day2'),
                                               suns=c(input$day1_start, input$day1_stop,
                                                      input$night1_start, input$night1_stop,
                                                      input$day2_start, input$day2_stop,
                                                      input$night2_start, input$night2_stop))
                   
                   # issue here, function legacyFingerprintMid uses genotype file path to create column as YYMMDD_treGrp
                   # but when running as app, genotype file path is 1.txt, so does not work
                   # correct this here
                  
                   colnames(fgp)[str_detect(colnames(fgp), '.txt')] <- paste(substr(genonms[i], 1, 9), input$treGrp_select, sep='_')
                   
                   # return the fingerprint
                   return(fgp)
                 })
                 
                 
                 # if we had multiple middur files, join fingerprint tables in the list fgpL
                 # (using join_all from plyr)
                 # and calculate mean of all experiments
                 if(length(fgpL)>1) {
                   fgp <- plyr::join_all(fgpL, by=c('uparam', 'parameter', 'win'), type='left')
                   # (avoid loading plyr because it replaces functions from dplyr)
                   
                   fgp <- fgp %>%
                     mutate(zsco = rowMeans(select(., matches('^\\d{6}_\\d{2}')))) # mean of columns that have name YYMMDD_BX
                   # regex is start / 6 digits (YYMMDD) / _ / 2 digits (BX)
                   # column called zsco is now mean of Z-scores
                   
                   # also record fgp as global variable
                   fgp <<- fgp
                   
                 } else { # if only 1 middur file
                   # then we have fingerprint table with column called YYMMDD_BX_treGrp
                   fgp <- fgpL[[1]]
                   # change that column name to simply treGrp
                   colnames(fgp)[which(stringr::str_detect(colnames(fgp), '^\\d{6}_\\d{2}'))] <- 'zsco'
                   # same regex as above
                   
                   # also record as global variable
                   fgp <<- fgp
                 }
                 # make fgp a global variable (<<-) because also used when user clicks on ranked drugs to reveal fingerprint plot
                 
                 
                 ### prepare fingerprint plot ###
                 ggfgp <- gglegacyFingerprint(
                   lFgp=fgp,
                   onlyGrp=NA,
                   colours=NA,
                   legendOrNo=TRUE,
                   ynameOrNo=TRUE,
                   ytextOrNo=TRUE,
                   xtextOrNo=TRUE,
                   yTitleSize=14,
                   paramNumOrNo=FALSE,
                   paramSize=10,
                   ynumSize=12,
                   nightBgOrNo=TRUE,
                   ymin=-3,
                   ymax=3,
                   exportOrNo=FALSE,
                   exportPath=NA,
                   width=NA,
                   height=NA
                 )
                 ## save fingerprint plot
                 # this may appear odd when testing locally because it actually saves the plot as pdf on drive
                 # but when app is deployed online, plot will be saved to server, so user will not see that
                 # then in downloadHandler we copy the file
                 # so really the user is saying "give me a copy of that fingerprint plot saved on the server"
                 ggsave('fingerprint.pdf', plot=ggfgp)
                 
                 
                 ### display fingerprint table ###
                 # will probably delete, it is more useful for debugging
                 output$fgp <- renderTable({
                   return(fgp) # this becomes 'fgp' in ui
                 })
                 
                 
                 ### display fingerprint plot ###
                 output$ggfgp <-  renderPlot({
                   return(ggfgp) # this becomes 'ggfgp' in ui
                 })
                 
                 # set-up the fingerprint plot download
                 output$ggfgp_dl <- downloadHandler(
                   filename=function() {
                     paste('fingerprint.pdf', sep='')
                   },
                   content=function(file) {
                     file.copy('fingerprint.pdf', file)
                   }
                 )
                 
                 
                 ###############################################################
                 ### ranked drugs ###
                 
                 # note about sections: it is better to go: calculate something / display it
                 # then calculate everything and display
                 # this way user sees things appearing while calculations are happening
                 # EDIT: that does not seem to work like that, actually
                 # only displays when done with all calculations, unfortunately
                 
                 ## rank drugs vs fingerprint
                 withProgress(message='ranking drugs', value=0.3, {
                   vdbr <<- rankDrugDb(legacyFgp=fgp, # vdbr is for fingerprint VS drug DB, Ranked
                                       dbPath='drugDb.csv',
                                       metric='cosine')
                   # also used when creating pop-up with ggBarcode, so keep as global variable (<<-)
                   
                   ## prepare the display version of the table
                   vdbr_dis <- cleanDrugsRanked(vdbr)
                   
                   incProgress(1.0)
                 })
                 
                 
                 ## display results
                 
                 # # display top X drugs
                 # output$topdr <- renderTable({ # topdr is for top X drugs
                 #   return(vdbr_dis[1:showNdrugs,]) # this becomes 'topdr' in ui
                 # })
                 # 
                 # # display bottom X drugs
                 # output$botdr <- renderTable({ # botdr is for bottom X drugs
                 #   return(vdbr_dis[(nrow(vdbr_dis)-(showNdrugs-1)):nrow(vdbr_dis),]) # this becomes 'topdr' in ui
                 # })
                 
                 output$vdbr_dis <- renderDT(vdbr_dis,
                                             selection=list(mode='single', target='row'),
                                             rownames=FALSE)
                 
                 # set-up the download
                 output$vdbr_dl <- downloadHandler( # download TTD targets
                   filename=function() {
                     'drugsRanked.csv'
                   },
                   content=function(file) {
                     vroom::vroom_write(vdbr, file, delim=',') # delim = ',' so writes CSV
                     # Note, we download the original vdbr, not the display version
                   }
                 )
                 
                 
                 ###############################################################
                 ### TTD indications ###
                 
                 # 10/12/2023
                 # calculations of enrichment are statistics are done in background processes "asynchronously"
                 # this solves an issue where app could not be reopened when it was processing
                 # a nice side effect is that it displays results as they become available
                 # ideas from
                 # https://stackoverflow.com/questions/30587883/is-it-possible-to-stop-executing-of-r-code-inside-shiny-without-stopping-the-sh/76005099#76005099
                 # https://www.r-bloggers.com/2020/04/asynchronous-background-execution-in-shiny-using-callr/
                 
                 ## calculate enrichment TTD indications
                 # start counter of number of seconds
                 i_ind <- 1
                 
                 # initiate background process
                 indications_bg <- reactive({
                   callr::r_bg(
                     func=drugEnrichment,
                     args=list(vdbr=vdbr,
                               namesPath='compounds.csv',
                               annotationPath='TTDindications.csv',
                               annotation='indications',
                               whichRank='rankeq',
                               minNex=3,
                               ndraws=ndraws,
                               alphaThr=0.05,
                               statsExport=NA),
                     supervise=TRUE
                   )
                 })
                 
                 # collect results, when they are ready
                 observe({
                   if (isolate(indications_bg()$is_alive())==FALSE) {
                     cat('\t \t \t \t >>> Indications done.\n')
                     
                     # empty the wait message
                     output$ind_waitMessage <- renderText({
                       ''
                     })
                     
                     ind <<- indications_bg()$get_result()
                     # also used when creating pop-up with ggBarcode, so keep as global variable (<<-)
                     
                     ## display results
                     # set-up the table
                     
                     ## prepare display version
                     ind_dis <- cleanIndications(ind)
                     
                     output$ind_dis <- renderDT(ind_dis,
                                                selection=list(mode='single', target='row'),
                                                rownames=FALSE)
                     # this becomes 'ind_dis' in ui
                     
                     ## if results not ready
                   } else {
                     # add 1 second to counter
                     i_ind <<- i_ind+1
                     cat('\t \t \t \t >>> Waiting for Indications to finish...\n')
                     
                     # display wait message
                     output$ind_waitMessage <- renderText({
                       paste('Calculations in progress...', i_ind, 'seconds elapsed.')
                     })
                     
                     # try again in 1 sec
                     invalidateLater(1000) # in milliseconds
                   }
                 })
                 
                 ## set-up the download
                 output$ind_dl <- downloadHandler( # download indications
                   filename=function() {
                     'indications.csv'
                   },
                   content=function(file) {
                     vroom::vroom_write(ind, file, delim=',') # delim = ',' so writes CSV
                     # Note, we download the full version, not the display one
                   }
                 )
                 
                 ###############################################################
                 ### TTD targets ###
                 
                 ## calculate enrichment TTD targets
                 # start counter of number of seconds
                 i_tar <- 1
                 
                 # initiate background process
                 targets_bg <- reactive({
                   callr::r_bg(
                     func=drugEnrichment,
                     args=list(vdbr=vdbr,
                               namesPath='compounds.csv',
                               annotationPath='TTDtargets.csv',
                               annotation='TTDtargets',
                               whichRank='rankeq',
                               minNex=3,
                               ndraws=ndraws,
                               alphaThr=0.05,
                               statsExport=NA),
                     supervise=TRUE
                   )
                 })
                 
                 # collect results, when they are ready
                 observe({
                   if (isolate(targets_bg()$is_alive())==FALSE) {
                     cat('\t \t \t \t >>> Targets done.\n')
                     
                     # empty the wait message
                     output$tar_waitMessage <- renderText({
                       ''
                     })
                     
                     ttar <<- targets_bg()$get_result()
                     # also used when creating pop-up with ggBarcode, so keep as global variable (<<-)
                     
                     ## display results
                     # set-up the table
                     
                     ## prepare display version
                     ttar_dis <- cleanTTDtargets(ttar)
                     
                     # TTD targets statistics report
                     output$ttar_dis <- renderDT(ttar_dis,
                                                 selection=list(mode='single', target='row'),
                                                 rownames=FALSE)
                     # this becomes 'tar_dis' in ui
                     
                     ## if results not ready
                   } else {
                     # add 1 second to counter
                     i_tar <<- i_tar+1
                     cat('\t \t \t \t >>> Waiting for Targets to finish...\n')
                     
                     # display wait message
                     output$tar_waitMessage <- renderText({
                       paste('Calculations in progress...', i_tar, 'seconds elapsed.')
                     })
                     
                     # try again in 1 sec
                     invalidateLater(1000) # in milliseconds
                   }
                 })
                 
                 # set-up the download
                 output$ttar_dl <- downloadHandler( # download TTD targets
                   filename=function() {
                     'TTDtargets.csv'
                   },
                   content=function(file) {
                     vroom::vroom_write(ttar, file, delim=',') # delim = ',' so writes CSV
                     # Note, we download the full table, not the display version
                   }
                 )


                 ###############################################################
                 ### KEGG pathways ###

                 ## calculate enrichment KEGG pathways
                 # start counter of number of seconds
                 i_keg <- 1
                 
                 # initiate background process
                 kegg_bg <- reactive({
                   callr::r_bg(
                     func=drugEnrichment,
                     args=list(vdbr=vdbr,
                               namesPath='compounds.csv',
                               annotationPath='TTDkegg.csv',
                               annotation='KEGG',
                               whichRank='rankeq',
                               minNex=3,
                               ndraws=ndraws,
                               alphaThr=0.05,
                               statsExport=NA),
                     supervise=TRUE
                   )
                 })
                 
                 # collect results, when they are ready
                 observe({
                   if (isolate(kegg_bg()$is_alive())==FALSE) {
                     cat('\t \t \t \t >>> KEGG done.\n')
                     
                     # empty the wait message
                     output$keg_waitMessage <- renderText({
                       ''
                     })
                     
                     keg <<- kegg_bg()$get_result()
                     # also used when creating pop-up with ggBarcode, so keep as global variable (<<-)
                     
                     ## display results
                     # set-up the table
                     
                     ## prepare display version
                     keg_dis <- cleanKEGG(keg)
                     
                     # KEGG pathways statistics report
                     output$keg_dis <- renderDT(keg_dis,
                                                selection=list(mode='single', target='row'),
                                                rownames=FALSE)
                     # this becomes 'keg_dis' in ui
                     
                     ## if results not ready
                   } else {
                     # add 1 second to counter
                     i_keg <<- i_keg+1
                     cat('\t \t \t \t >>> Waiting for KEGG to finish...\n')
                     
                     # display wait message
                     output$keg_waitMessage <- renderText({
                       paste('Calculations in progress...', i_keg, 'seconds elapsed.')
                     })
                     
                     # try again in 1 sec
                     invalidateLater(1000) # in milliseconds
                   }
                 })
                 
                 # set-up the download
                 output$keg_dl <- downloadHandler( # download KEGG
                   filename=function() {
                     'KEGGpathways.csv'
                   },
                   content=function(file) {
                     vroom::vroom_write(keg, file, delim=',') # delim = ',' so writes CSV
                     # Note, we download the full table, not the display version
                   }
                 )
                 
                 ###############################################################

  })
  
  
  ###############################################################
  ### clickable drugs ranked table to reveal fingerprint plot ###
  observeEvent(input$vdbr_dis_rows_selected, {
    req(input$vdbr_dis_rows_selected)
    
    # get the index of the selected row
    clickrow <- input$vdbr_dis_rows_selected
    
    # what is the CID of the selected compound?
    cid <- vdbr[clickrow, 'cid']
    
    # search all other rows with that CID and take their compounds names (original names)
    nms2plot <- unique( vdbr[which(vdbr$cid==cid), 'name'] ) # these are the names to plot
    # (it will necessarily include the compound the user clicked on)
    
    # prepare the fingerprint plot
    ggDfp <- ggDrugFgp(drugDb='drugDb.csv',
                       dnames=nms2plot,
                       legacyFgp=fgp,
                       onlyGrp=c('mean fgp', nms2plot), # if user gave multiple experiments, this will only keep 'mean fgp' which is the mean fingerprint
                       colours=NA,
                       legendOrNo=TRUE,
                       ynameOrNo=TRUE,
                       ytextOrNo=TRUE,
                       xtextOrNo=TRUE,
                       paramNumOrNo=TRUE,
                       nightBgOrNo=TRUE,
                       ymin=-3,
                       ymax=3,
                       exportOrNo=FALSE, 
                       exportPath=NA,
                       width=NA,
                       height=NA)
    
    modal <- modalDialog(
      title = 'Barcode plot' ,
      
      # tell user about what we did
      paste(length(which(vdbr$cid==cid)), 'fingerprint(s) for PubChem CID', cid, '\n'),
      
      plotOutput( 'ggDfp' ) ,
      
      easyClose=TRUE,
      footer=modalButton('Close')
    )
    
    showModal(modal)
    
    output$ggDfp <- renderPlot({
      return(ggDfp)
    })
    
  })
  
  
  ###############################################################
  ### clickable Indications table to reveal barcode plot ###
  observeEvent(input$ind_dis_rows_selected, {
    req(input$ind_dis_rows_selected)
    
    # get the index of the selected row
    clickrow <- input$ind_dis_rows_selected
    
    # prepare the barcode plot
    ggBc <- ggBarcode(vdbr=vdbr,
                      namesPath='compounds.csv',
                      annotationPath='TTDindications.csv',
                      annotation='indications',
                      testAnnotation=ind[clickrow , 'annotation'],
                      minScore=NA,
                      barwidth1=2,
                      barwidth2=25,
                      exportPath=NA,
                      width=NA,
                      height=NA)
    
    modal <- modalDialog(
      title = 'Barcode plot' ,
      
      paste(ind[clickrow, 'annotation']) , 
      
      plotOutput( 'ggBc' ) ,
      
      easyClose=TRUE,
      footer=modalButton('Close')
    )
    
    showModal(modal)
    
    output$ggBc <- renderPlot({
      return(ggBc)
    })
    
  })
  
  
  ###############################################################
  ### clickable TTD targets table to reveal barcode plot ###
  observeEvent(input$ttar_dis_rows_selected, {
    req(input$ttar_dis_rows_selected)
    
    # get the index of the selected row
    clickrow <- input$ttar_dis_rows_selected
    
    # prepare the barcode plot
    ggBc <- ggBarcode(vdbr=vdbr,
                      namesPath='compounds.csv',
                      annotationPath='TTDtargets.csv',
                      annotation='TTDtargets',
                      testAnnotation=ttar[clickrow , 'annotation'],
                      minScore=NA,
                      barwidth1=2,
                      barwidth2=25,
                      exportPath=NA,
                      width=NA,
                      height=NA)
    
    modal <- modalDialog(
      title = 'Barcode plot' ,
      
      paste(ttar[clickrow, 'TARGNAME']) , 
      
      plotOutput( 'ggBc' ) ,
      
      easyClose=TRUE,
      footer=modalButton('Close')
    )
    
    showModal(modal)
    
    output$ggBc <- renderPlot({
      return(ggBc)
    })
    
  })
  
  
  ###############################################################
  ### clickable KEGG pathways table to reveal barcode plot ###
  observeEvent(input$keg_dis_rows_selected, {
    req(input$keg_dis_rows_selected)
    
    # get the index of the selected row
    clickrow <- input$keg_dis_rows_selected
    
    # prepare the barcode plot
    ggBc <- ggBarcode(vdbr=vdbr,
                      namesPath='compounds.csv',
                      annotationPath='TTDkegg.csv',
                      annotation='KEGG',
                      testAnnotation=keg[clickrow , 'annotation'],
                      minScore=NA,
                      barwidth1=2,
                      barwidth2=25,
                      exportPath=NA,
                      width=NA,
                      height=NA)
    
    modal <- modalDialog(
      title = 'Barcode plot' ,
      
      paste(keg[clickrow, 'keggname']) , 
      
      plotOutput( 'ggBc' ) ,
      
      easyClose=TRUE,
      footer=modalButton('Close')
    )

    showModal(modal)
    
    output$ggBc <- renderPlot({
      return(ggBc)
    })

  })
  
}


# run the app -------------------------------------------------------------

shinyApp(ui, server)