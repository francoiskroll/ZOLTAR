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



# functions ---------------------------------------------------------------

source('legacyFingerprint.R')



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
      
      # Input: Select a file ----
      fileInput('mat', 'Select or drop your .mat file',
                multiple=TRUE,
                accept='.mat')
      
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
        tabPanel('Fingerprint', tableOutput('contents')),
        tabPanel('Ranked drugs')
      )
    )
    
  )
)


# server ------------------------------------------------------------------


# Define server logic to read selected file ----
server <- function(input, output) {
  
  output$contents <- renderTable({
    
    # input$mat is NULL initially
    
    req(input$mat) # check that mat file is available
    
    # df <- read.csv(input$file1$datapath,
    #                header = input$header,
    #                sep = input$sep,
    #                quote = input$quote)
    
    df <- R.matlab::readMat(input$mat$datapath)$geno
    
    fgp <- legacyFingerprint(
      matPath=input$mat$datapath,
      conGrp='scr',
      treGrp='sorl1',
      nights=c(2,3),
      days=c(2,3)
    )
    
    return(fgp) # this becomes CONTENTS in ui
    
  })
  
}



# run the app -------------------------------------------------------------

shinyApp(ui, server)