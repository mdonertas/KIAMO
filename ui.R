library(shiny)

genelist = readRDS('./data/geneids.rds')$displaygenes

ui <- fluidPage(
  
  theme = bslib::bs_theme(bootswatch = 'journal'),
  
  tags$h1('KIAMO - Killifish Immune Ageing MultiOmics'),
  tags$br(''),
  
  sidebarLayout(
    
    sidebarPanel(width = 3,
      tags$small('Please select (or erase and type) the gene you would like to search for:'),
      tags$br(''),
      selectInput('gene','Gene Symbol', choices = sort(genelist)),
      actionButton('submit', 'GO'),
      
      tags$br(''),
      tags$small('* Click on the tSNE plot to display the cell type annotations.'),
      tags$br(''),
      tags$small('**Normalized expression: Corrected for the library size & log normalized.'),
      tags$br(''),
      tags$small('***All data were generated using GRZ-AD strain Turquoise killifish. Proteomics data include 5 x 7-week-old (labelled young), and 5 x 16-week-old (labelled old) fish, whereas scRNAseq data include 2 x 8-week-old (labelled young) and 3 x 21-week-old (labelled old) fish.')
      
    ),
    
    mainPanel(width = 9,
      
      fluidRow(
        column(5, plotOutput('plasma_box', height = 300)),
        column(5, plotOutput('km_box', height = 300))
      ),
      
      fluidRow(
        column(5, plotOutput('tsne', height = 300, click = 'plot_click')),
        column(7, plotOutput('dotplot', height = 350))
      ),
      
      fluidRow(
        column(5, verbatimTextOutput("info"))
      )
      
    )
    
  )
  
)

shinymanager::secure_app(ui)

