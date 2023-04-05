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

ndraws <- 2
alphaThr <- 0.2

# set maximum upload to 30 Gb
options(shiny.maxRequestSize = 3000*1024^2)

# ! will face some issues when deployed, likely need Basic plan at 440$/year
# https://support.posit.co/hc/en-us/articles/219449487--How-much-data-can-I-upload-to-shinyapps-io-

# seems fairly easy to save uploaded data to Dropbox, if useful at some point
# https://www.r-bloggers.com/2015/07/persistent-data-storage-in-shiny-apps/


# ui ----------------------------------------------------------------------


ui <- fluidPage(
  
  #theme=bslib::bs_theme(bootswatch='sandstone'),
  #theme=bslib::bs_theme(bootswatch='united'), # I quite like it but maybe the font is too childish
  # theme=bslib::bs_theme(bootswatch='simplex'),
  
  # can also set theme manually
  # theme=bslib::bs_theme(
  #   base_font='Rubik'
  # ),
  
  ### app title
  titlePanel('Predictive pharmacology'),
  # subtitle
  p(em('From behavioural fingerprint to drugs and pathways'), style='color:grey'),
  
  ### sidebar layout = inputs on left, outputs on right
  sidebarLayout(
    
    ### sidebar
    sidebarPanel(
      
      ## middur input widget
      fileInput(inputId='mid_drop',
                label='Select or drop your middur.csv file',
                multiple=TRUE,
                accept='.csv') ,
      
      ## genotype input widget
      fileInput(inputId='geno_drop',
                label='Select or drop your genotype.txt file',
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
                   label='Go')
      
    ),
    
    ### tabs on the right
    mainPanel(
      tabsetPanel(
        
        tabPanel('Table',
                 tableOutput('fgp')),
        
        tabPanel('Fingerprint',
                 p(''),
                 downloadButton('ggfgp_dl', 'download pdf'),
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
                 p('• Click on Cosine column to rank by increasing or decreasing'),
                 p('• Click on a row to reveal fingerprint plot.'),
                 downloadButton('vdbr_dl', 'download'),
                 DTOutput('vdbr_dis')),
        
        tabPanel('Indications',
                 p(''),
                 p('Source: Therapeutic Target Database'),
                 downloadButton('ind_dl', 'download'),
                 DTOutput('ind_dis')) ,
        
        tabPanel('Targets',
                 p(''),
                 p('Source: Therapeutic Target Database'),
                 downloadButton('ttar_dl', 'download'),
                 DTOutput('ttar_dis')) ,
        
        tabPanel('KEGG pathways',
                 p(''),
                 p('Source: Therapeutic Target Database'),
                 downloadButton('keg_dl', 'download'),
                 DTOutput('keg_dis')) ,
        
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
  
  #### update the group selection when user drops the genotype file ####
  # note input$geno_drop is NULL initially
  observeEvent(input$geno_drop, # means: observe the genotype file import widget, if it changes, trigger the code below
               {
                 # import the genotype file
                 geno <- importGenotype(input$geno_drop$datapath)
                 # get the group names from it
                 grpnms <- colnames(geno)
                 
                 # update the selection of treatment / control group so it lists the groups available in the data
                 updateSelectInput(session, input='treGrp_select',
                                   choices=grpnms)
                 updateSelectInput(session, input='conGrp_select',
                                   choices=grpnms)
               })
  
  
  #### start stuff when user presses Go button ####
  observeEvent(input$go_button,
               {
                 req(input$mid_drop) # check that middur file is available
                 req(input$geno_drop) # check that genotype file is available
                 
                 ### calculate fingerprint ###
                 
                 # import middur file
                 mid <- read.csv(input$mid_drop$datapath)
                 
                 fgp <<- legacyFingerprintMid(mid=mid,
                                              genopath=input$geno_drop$datapath,
                                              treGrp=input$treGrp_select,
                                              conGrp=input$conGrp_select,
                                              nights=c('night1', 'night2'),
                                              days=c('day1', 'day2'))
                 
                 ### prepare fingerprint plot ###
                 ggfgp <- gglegacyFingerprint(
                   lFgp=fgp,
                   onlyGrp=NA,
                   colours=NA,
                   legendOrNo=TRUE,
                   ynameOrNo=TRUE,
                   ytextOrNo=TRUE,
                   xtextOrNo=TRUE,
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
                 
                 ## calculate enrichment TTD indications
                 withProgress(message='calculating indications', value=0.3, {
                   ind <<- drugEnrichment(vdbr=vdbr,
                                         namesPath='compounds.csv',
                                         annotationPath='TTDindications.csv',
                                         annotation='indications',
                                         whichRank='rankeq',
                                         minNex=3,
                                         ndraws=ndraws,
                                         alphaThr=alphaThr,
                                         statsExport=NA)
                   # also used when creating pop-up with ggBarcode, so keep as global variable (<<-)
                   
                   ## prepare display version
                   ind_dis <- cleanIndications(ind)
                   
                   incProgress(1.0)
                 })
                 
                 ## display results
                 # set-up the table
                 output$ind_dis <- renderDT(ind_dis,
                                            selection=list(mode='single', target='row'),
                                            rownames=FALSE)
                 # this becomes 'ind_dis' in ui
                 
                 # set-up the download
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
                 withProgress(message='calculating targets', value=0.3, {
                   ttar <<- drugEnrichment(vdbr=vdbr,
                                           namesPath='compounds.csv',
                                           annotationPath='TTDtargets.csv',
                                           annotation='TTDtargets',
                                           whichRank='rankeq',
                                           minNex=3,
                                           ndraws=ndraws,
                                           alphaThr=alphaThr,
                                           statsExport=NA)
                   # also used when creating pop-up with ggBarcode, so keep as global variable (<<-)
                   
                   ## prepare the display version
                   ttar_dis <- cleanTTDtargets(ttar)
                   
                   incProgress(1.0)
                 })
                 
                 # TTD targets statistics report
                 output$ttar_dis <- renderDT(ttar_dis,
                                            selection=list(mode='single', target='row'),
                                            rownames=FALSE)
                 # this becomes 'tar_dis' in ui
                 
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
                 withProgress(message='calculating KEGG pathways', value=0.3, {
                   keg <<- drugEnrichment(vdbr=vdbr,
                                         namesPath='compounds.csv',
                                         annotationPath='TTDkegg.csv',
                                         annotation='KEGG',
                                         whichRank='rankeq',
                                         minNex=3,
                                         ndraws=ndraws,
                                         alphaThr=alphaThr,
                                         statsExport=NA)
                   # also used when creating pop-up with ggBarcode, so keep as global variable (<<-)
                   
                   ## prepare display version
                   keg_dis <- cleanKEGG(keg)
                   
                   incProgress(1.0)
                 })
                 
                 ## display results
                 # set-up the table
                 output$keg_dis <- renderDT(keg_dis, # KEGG pathways statistics report
                                            selection=list(mode='single', target='row'),
                                            rownames=FALSE)
                 # this becomes 'keg_dis' in ui
                 
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
                       onlyGrp=NA,
                       colours=NA,
                       legendOrNo=TRUE,
                       ynameOrNo=TRUE,
                       ytextOrNo=TRUE,
                       xtextOrNo=TRUE,
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