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



# functions ---------------------------------------------------------------

source('legacyFingerprint.R')
source('gglegacyFingerprint.R')
source('drawEnrich_v5.R')



# settings ----------------------------------------------------------------

Sys.setlocale("LC_ALL","C") # avoids an issue when printing table of ranked drugs, probably because of odd characters in original drug names
# solution StackOverflow question 61656119

ndraws <- 10
alphaThr <- 0.05

# options(shiny.maxRequestSize = 30*1024^2) # maximum upload set to 30 Mb


# ui ----------------------------------------------------------------------


ui <- fluidPage(
  
  # theme=bslib::bs_theme(bootswatch='sandstone'),
  #theme=bslib::bs_theme(bootswatch='united'), # maybe the font is too childish
  #theme=bslib::bs_theme(bootswatch='simplex'),
   
  # can also set theme manually
  # theme=bslib::bs_theme(
  #   base_font='Rubik'
  # ),
  
  ### app title
  titlePanel('Predictive pharmacology'),
  # subtitle
  p(em('From behavioural fingerprint to drugs and pathways'), style='color:grey'),
  # h4('Francois Kroll & Jason Rihel'),
  
  ### sidebar layout = inputs on left, outputs on right
  sidebarLayout(
    
    ### sidebar
    sidebarPanel(
      
    ### .mat input widget
      fileInput(inputId='mat_drop',
                label='Select or drop your .mat file',
                multiple=TRUE,
                accept='.mat') ,
      
      selectInput(inputId='treGrp_select',
                  label='Select the treatment group',
                  choices=NULL, # initialise as NULL, we will update based on data given
                  multiple=FALSE) , # can only select one group
      
      selectInput(inputId='conGrp_select',
                  label='Select the control group',
                  choices=NULL, # initialise as NULL, we will update based on data given
                  multiple=FALSE) , # can only select one group
      
      actionButton(inputId='go_button',
                   label='Go')
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel('Table',
                 tableOutput('fgp')),
        
        tabPanel('Fingerprint',
                 plotOutput('ggfgp')),
        
        tabPanel('Drugs ranked',
                 p('Drugs are ranked by decreasing cosine.'),
                 downloadButton('ind_dl', 'download'),
                 h3('Top 20'),
                 tableOutput('topdr'), # topdr is for top X drugs
                 h3('Bottom 20'),
                 tableOutput('botdr')), # botdr is for bottom X drugs
        
        tabPanel('Indications',
                 p('Source: Therapeutic Target Database'),
                 tableOutput('ind')) ,
        
        tabPanel('Targets',
                 p('Source: Therapeutic Target Database'),
                 tableOutput('tar')) ,
        
        tabPanel('KEGG pathways',
                 p('Source: Therapeutic Target Database'),
                 tableOutput('keg')) ,
        
        tabPanel('zebrafish STITCH',
                 p('Source: stitch.embl.de'),
                 tableOutput('zsti')) ,
        
        tabPanel('human STITCH',
                 p('Source: stitch.embl.de'),
                 tableOutput('hsti'))
      )
    )
    
  )
)


# server ------------------------------------------------------------------

# about progress bars:
# currently a bit of a cheap trick: it is started at 0.3 then updated to 1.0 when done
# at least user knows something is going on
# it seems like using a Progress object could make it more precise but requires more work

server <- function(input, output, session) {
  
  #### update the group selection when user drops a .mat file ####
  # note input$mat_drop is NULL initially
  observeEvent(input$mat_drop, # means: observe the .mat import widget, if it changes, trigger the code below
               {
                 # import the .mat file
                 mat <- R.matlab::readMat(input$mat_drop$datapath)$geno
                 
                 # read the group names
                 grpnms <- unlist(mat[,,1]$name)
                 
                 # update the selection of treatment / control group so it lists the groups available in the data
                 updateSelectInput(session, input='treGrp_select',
                                      choices=grpnms)
                 updateSelectInput(session, input='conGrp_select',
                                   choices=grpnms)
               })
  
  #### start stuff when user presses Go button ####
  observeEvent(input$go_button,
               {
                 req(input$mat_drop) # check that mat file is available
                 
                 
                 ### calculate fingerprint ###
                 fgp <- legacyFingerprint(
                   matPath=input$mat_drop$datapath,
                   conGrp=input$conGrp_select,
                   treGrp=input$treGrp_select,
                   nights=c(2,3),
                   days=c(2,3)
                 )
                 
                 
                 ### rank drugs vs fingerprint ###
                 withProgress(message='ranking drugs', value=0.3, {
                   vdbr <- rankDrugDb(legacyFgp=fgp, # vdbr is for fingerprint VS drug DB, Ranked
                                      dbPath='drugDb.csv', 
                                      metric='cosine')
                   incProgress(1.0)
                 })
                 
                 
                 ### calculate enrichment TTD indications ###
                 withProgress(message='calculating indications', value=0.3, {
                   ind <- drugEnrichment(vdbr=vdbr,
                                         namesPath='compounds.csv',
                                         annotationPath='TTDindications.csv',
                                         annotation='indications',
                                         whichRank='rankeq',
                                         minNex=3,
                                         ndraws=ndraws,
                                         alphaThr=alphaThr,
                                         statsExport=NA)
                   incProgress(1.0)
                 })
                 
                 
                 ### calculate enrichment TTD targets ###
                 withProgress(message='calculating targets', value=0.3, {
                   tar <- drugEnrichment(vdbr=vdbr,
                                         namesPath='compounds.csv',
                                         annotationPath='TTDtargets.csv',
                                         annotation='TTDtargets',
                                         whichRank='rankeq',
                                         minNex=3,
                                         ndraws=ndraws,
                                         alphaThr=alphaThr,
                                         statsExport=NA)
                   incProgress(1.0)
                 })
                 
                 
                 ### calculate enrichment KEGG pathways ###
                 withProgress(message='calculating KEGG pathways', value=0.3, {
                   keg <- drugEnrichment(vdbr=vdbr,
                                         namesPath='compounds.csv',
                                         annotationPath='TTDkegg.csv',
                                         annotation='KEGG',
                                         whichRank='rankeq',
                                         minNex=3,
                                         ndraws=ndraws,
                                         alphaThr=alphaThr,
                                         statsExport=NA)
                   incProgress(1.0)
                 })

                 
                 
                 ### fingerprint table ###
                 output$fgp <- renderTable({
                   return(fgp) # this becomes 'fgp' in ui
                 })
                 
                 
                 ### fingerprint plot ###
                 output$ggfgp <-  renderPlot({
                   
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
                   
                   return(ggfgp) # this becomes 'ggfgp' in ui
                   
                 })
                 
                 
                 ### top X drugs ###
                 output$topdr <- renderTable({ # topdr is for top X drugs
                   return(vdbr[1:20,]) # this becomes 'topdr' in ui
                 })
                 
                 ### bottom X drugs ###
                 output$botdr <- renderTable({ # botdr is for bottom X drugs
                   return(vdbr[(nrow(vdbr)-19):nrow(vdbr),]) # this becomes 'topdr' in ui
                 })
                 
                 
                 ### TTD indications ###
                 output$ind <- renderTable({ # indications statistics report
                   return(ind) # this becomes 'ind' in ui
                 })
                 
                 output$ind_dl <- downloadHandler( # download indications
                   filename=function() {
                     'indications.csv'
                   },
                   content=function(file) {
                     vroom::vroom_write(ind, file, delim=',') # delim = ',' so writes CSV
                   }
                 )
                 
                 ### TTD targets ###
                 output$tar <- renderTable({ # TTD targets statistics report
                   return(tar) # this becomes 'tar' in ui
                 })
                 
                 ### KEGG pathways ###
                 output$keg <- renderTable({ # KEGG pathways statistics report
                   return(keg) # this becomes 'keg' in ui
                 })
  })
  
}


# run the app -------------------------------------------------------------

shinyApp(ui, server)