#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
source("utils.R")

# Define server logic

shinyServer(function(input, output, session) {

    # upload the data
    data <- reactive({

      req(input$upload)
      import_diann(file_path = input$upload$datapath)

      })

    # get sample names
    sampleOptions <- reactive({
      unique(data()$Samples)
    })

    # render info about samples on main page
    output$text1 <- renderText({
      paste0("Found ", length(sampleOptions()), " samples in this dataset! they are:")})

    output$text2 <- renderText({
      paste0(sampleOptions())
    })

    ## second tab

    # update selections according to the samples
    observe({
      updateSelectInput(session, "group1",
                        choices = sampleOptions()
      )})

    observe({
      updateSelectInput(session, "group2",
                        choices = setdiff(sampleOptions(), input$group1)
      )})

    # run analysis after button is pressed
    subset_data <-  eventReactive(input$analysis,{
      run_ttest(data(), a = input$group1, b = input$group2) |>
                       dplyr::filter(!is.na(`p.adjusted.BH`))
             })

    # let user know that analysis is running
    observeEvent(input$analysis, {
      showNotification("Running T-test for the selected data! (Please wait a few seconds)", duration = 60, type = "message")
    })

    # enable plot design

    # change the level of significance and expression based on the user's choices
    subset_plotting <- reactive({
      subset_data() |>
        dplyr::mutate(significant = ifelse(`p.adjusted.BH` >= input$pvalue,
                                                          "not significant",
                                                          ifelse(
                                                            `log2FC_A/B` <= input$fc[1], "significant DOWN",
                                                            ifelse(
                                                              `log2FC_A/B` >= input$fc[2], "significant UP",
                                                              "not significant"))),
                      color = case_when(significant == "not significant" ~ "gray",
                                 significant == "significant UP" ~ "orange",
                                 significant == "significant DOWN" ~ "blue"))
    })

    plot_colors <- reactive({subset_plotting() |>
                              dplyr::select(significant, color) |>
                              tibble::deframe()})

    volcanoplot <- reactive({ggplot2::ggplot(subset_plotting(),
                                   aes(angle = Genes,
                                       x = `log2FC_A/B`,
                                       y = -log10(`p.adjusted.BH`),
                                       color = significant,
                                       )) +
      labs(title = glue::glue("Volcano plot - Welch's T-test | {input$group1} over {input$group2}")) +
      geom_point() +
      geom_vline(xintercept = input$fc[1], linetype = "dashed") +
      geom_vline(xintercept = input$fc[2], linetype = "dashed") +
      geom_hline(yintercept = -log10(input$pvalue), linetype = "dashed") +
      scale_color_manual(values = plot_colors())})


    output$volcano <- plotly::renderPlotly({
      plotly::ggplotly(volcanoplot(), tooltip = c("x", "y", "angle"))
    })

    # prepare the current visualized table for download
    toSave <- reactive({subset_plotting() |>
                          dplyr::select(-color)
              })

    output$table <- DT::renderDataTable(toSave())

    output$dl.full <- downloadHandler(
        filename = function() {
          glue::glue("{Sys.Date()}_{input$group1}_over_{input$group2}_fullanalysis.csv")
        },
        content = function(con) {
          write.csv(toSave(), con)
        }
      )

    output$dl.up <- downloadHandler(
      filename = function() {
        glue::glue("{Sys.Date()}_{input$group1}_over_{input$group2}_UP_only.csv")
      },
        content = function(con) {
          toSave_up <- toSave() |>
            dplyr::filter(significant == "significant UP")
            write.csv(toSave_up, con)
      }
    )

    output$dl.down <- downloadHandler(
      filename = function() {
        glue::glue("{Sys.Date()}_{input$group1}_over_{input$group2}_DOWN_only.csv")
      },
      content = function(con) {
        toSave_down <- toSave() |>
          dplyr::filter(significant == "significant DOWN")
          write.csv(toSave_down, con)
      }
    )

})
