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



# settings ----------------------------------------------------------------

# options(shiny.maxRequestSize = 30*1024^2) # maximum upload set to 30 Mb


ui <- fluidPage(
  
  # App title ----
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
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      tableOutput("contents")
      
    )
    
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    
    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    
  })
  
}
# Run the app ----
shinyApp(ui, server)