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

      checkboxGroupInput("filtervar", "Filter the gene list based on proteomics expression change with age:",
                         c("Increase with age in kidney marrow" = "km_up",
                           "Decrease with age in kidney marrow" = "km_down",
                           "Significant change with age in kidney marrow" = "km_sig",
                           "Increase with age in plasma" = "pl_up",
                           "Decrease with age in plasma" = "pl_down",
                           "Significant change with age in plasma" = "pl_sig")),
      actionButton('filter', 'Filter'),
      tags$br(''),
      tags$small('* Click on the tSNE plot to display the cell type annotations.'),
      tags$br(''),
      tags$small('* Normalized expression: Corrected for the library size & log normalized.'),
      tags$br(''),
      tags$small('* All data were generated using GRZ-AD strain Turquoise killifish. Proteomics data include 5 x 7-week-old (labelled young), and 5 x 16-week-old (labelled old) fish, whereas scRNAseq data include 2 x 8-week-old (labelled young) and 3 x 21-week-old (labelled old) fish.'),
      tags$br(''),
      tags$small('* For detailed instructions and an example use case, please visit https://github.com/mdonertas/KIAMO/wiki'),
      tags$br(''),
      tags$small('* For detailed information about the data and analysis, please see our paper: Spontaneous onset of cellular markers of inflammation and genome instability during aging in the immune niche of the naturally short-lived turquoise killifish (Nothobranchius furzeri)
Gabriele Morabito, Handan Melike Dönertas, Luca Sperti, Jens Seidel, Aysan Poursadegh, Michael Poeschla, Dario Riccardo Valenzano
bioRxiv 2023.02.06.527346; doi: https://doi.org/10.1101/2023.02.06.527346')
      
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

