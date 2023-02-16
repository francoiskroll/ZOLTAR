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



# functions ---------------------------------------------------------------

source('legacyFingerprint.R')
source('gglegacyFingerprint.R')



# settings ----------------------------------------------------------------

# options(shiny.maxRequestSize = 30*1024^2) # maximum upload set to 30 Mb


# ui ----------------------------------------------------------------------


ui <- fluidPage(
  
  ### app title
  titlePanel('Predictive pharmacology'),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
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
    
    # # Main panel for displaying outputs ----
    # mainPanel(
    #   
    #   # Output: Data file ----
    #   tableOutput("contents")
    #   
    # ),
    
    mainPanel(
      tabsetPanel(
        tabPanel('Fingerprint',
                 tableOutput('fgp')),
        
        tabPanel('Fingerprint2',
                 plotOutput('ggfgp'))
      )
    )
    
  )
)


# server ------------------------------------------------------------------

server <- function(input, output, session) {
  
  #### update the group selection when user drops a .mat file ####
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
                 output$fgp <- renderTable({
                   
                   # input$mat is NULL initially
                   
                   req(input$mat_drop) # check that mat file is available
                   
                   fgp <- legacyFingerprint(
                     matPath=input$mat_drop$datapath,
                     conGrp=input$conGrp_select,
                     treGrp=input$treGrp_select,
                     nights=c(2,3),
                     days=c(2,3)
                   )
                   
                   return(fgp) # this becomes `fgp` in ui
                   
                 })
                 
                 output$ggfgp <-  renderPlot({
                   
                   fgp <- legacyFingerprint(
                     matPath=input$mat_drop$datapath,
                     conGrp=input$conGrp_select,
                     treGrp=input$treGrp_select,
                     nights=c(2,3),
                     days=c(2,3)
                     )
                   
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
                   
                   return(ggfgp)
                   
                 })
  })
  
  #### make fingerprint plot####
  
}


# run the app -------------------------------------------------------------

shinyApp(ui, server)