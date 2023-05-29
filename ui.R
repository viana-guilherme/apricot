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
shinyUI(navbarPage("apricot",

    # landing page
    tabPanel("Welcome!",
              h2("Apricot is a tool for automated proteomics differential expression analysis and ontology testing")
             ),

    # main functionalities page
    tabPanel("Run analysis",
      # first screen
      tabPanel("Upload a dataset",
                fileInput(inputId = "upload",
                          label = "Please, upload a file DIA-NN output file",
                          buttonLabel = "Select file..."),
               textOutput("text1"),
               textOutput("text2")
               ),

      # second screen
      tabPanel("Run analysis!",
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
                   conditionalPanel(
                   condition = "input.analysis == true",
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
                 )),

                 # Show a plot of the generated distribution
                 mainPanel(
                   tabsetPanel(type = "tabs",
                               tabPanel("Volcano plot", plotly::plotlyOutput("volcano")),
                               tabPanel("Gene ontology"),
                               tabPanel("Results table", DT::dataTableOutput("table")))
                 )
               )),

      # third screen
      tabPanel("Download options",
               h4("Download options"),
               fluidRow(
                column(width = 4, downloadButton("dl.full", label = "Download full table")),
                column(width = 4, downloadButton("dl.up", label = "Download significant UP only")),
                column(width = 4, downloadButton("dl.down", label = "Download significant DOWN only"))
               )
              )),

    # "about us" screen
    tabPanel("References")

))
