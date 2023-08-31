#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for the application
shinyUI(
  navbarPage("apricot",

    # landing page
    tabPanel("Welcome!",
             # enabling styling via the external css
             tags$head(
               tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),

             #content
             tags$div(class = "mainpage-content",
             img(class = "img-logo", src = "logo.png", align = "center", height = 150),
             tags$p(class = "text-highlight", "Welcome to apricot!"),
             tags$ul(
             tags$li("Apricot is a web-based tool for performing differential expression analysis and gene ontology enrichment tests
                    on DIA-NN-generated proteomics datasets"),
             tags$li("To get started, you can read our documentation or head to the 'run analysis' tab to upload and process your dataset"),
             tags$li("To know more about this project and the methodology employed, please check out the 'References' section"),
             tags$li("Please, direct any comment or suggestion to viana.guilherme@lbl.gov")
             ))),

    # main functionalities page
    tabPanel("Run analysis",
      # first screen
      tabPanel("Upload a dataset",
                fileInput(inputId = "upload",
                          label = "Please, upload a DIA-NN output file. You can click the buttom or drag files to the form below",
                          buttonLabel = "Select file..."),
               textOutput("text1"),
               textOutput("text2")
               ),

      # second screen
      tabPanel("",
               sidebarLayout(
                 sidebarPanel(
                   helpText("Please select two samples to perform a comparative analysis", br(),
                   "The first group is the experimental group you're interested in evaluating", br(),
                   "The second group is the reference group to which measure against",br(),
                   "You'll be able to save the current plot by clicking the camera icon on the figure"),
                    selectInput("group1",
                               choices = NULL,
                               label = "select the first group",
                               multiple = FALSE
                              ),
                   selectInput("group2",
                               choices = NULL,
                               label = "select the second group",
                               multiple = FALSE
                   ),
                   actionButton(inputId = "analysis",
                                label = "run analysis!"),
                   br(),
                   br(),
                   sliderInput("pvalue",
                               "Customize the adjusted p-value significant threshold",
                               min = 0,
                               max = 0.05,
                               value = 0.05,
                               step = 0.001),
                   sliderInput("fc",
                               "Customize the Fold Change range",
                               min = -10,
                               max = 10,
                               step = 0.1,
                               value = c(-1, 1)
                               )
                 ),

                 # Show a plot of the generated distribution
                 mainPanel(
                   tabsetPanel(type = "tabs",
                               tabPanel("Volcano plot", plotly::plotlyOutput("volcano")),
                               tabPanel("Gene ontology", plotOutput("ontologyplot")),
                               tabPanel("Results table", DT::dataTableOutput("table"))),
                   tabPanel("Download options",
                            h4("Data download options"),
                            fluidRow(
                              column(width = 4, downloadButton("dl.full", label = "Download full T-test table")),
                              column(width = 4, downloadButton("dl.up", label = "Download T-test significant UP only")),
                              column(width = 4, downloadButton("dl.down", label = "Download T-test significant DOWN only"))
                            )
                   )),
                 )
               )),


    # "about us" screen
    tabPanel("References")

))
