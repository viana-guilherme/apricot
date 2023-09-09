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
      base::unique(data()$Samples)
    })

    # render info about samples on main page

    initialMessage <- observeEvent(

      sampleOptions(), {

      optLength <- length(sampleOptions())
      optNames <- stringr::str_flatten_comma(sampleOptions(), last = " and ")

      showModal(modalDialog(
        title = "Opening dataset file...",
        paste0("Found ", optLength, " samples in this dataset! they are: ", optNames),
        easyClose = TRUE))
    })


    ## second tab

    # update selections according to the samples
    observe({
      updateSelectInput(session, "group1",
                        choices = sampleOptions()
      )})

    observe({
      updateSelectInput(session, "group2",
                        choices = base::setdiff(sampleOptions(), input$group1)
      )})

    # run analysis after button is pressed
    subset_data <-  eventReactive(input$analysis,{
      run_ttest(data(), a = input$group1, b = input$group2) |>
                       dplyr::filter(!is.na(`p.adjusted.BH`))
             })

    # let user know that analysis is running
    observeEvent(input$analysis, {
      showModal(modalDialog(
        title = "Running analysis...",
        glue::glue("Running T-test for selected samples! (Please wait a few seconds and close this pop up. The page should refresh soon)."),
        easyClose = TRUE
      ))
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

    # Gene ontology part

    # Running the gene ontology enrichment test based on the selected dataset

    resultsGO <- reactive({
      runEnrichment(dataset = toSave(), significance_threshold = input$pvalue)
      })

    goplot <- reactive({
      resultsGO() |>
        dplyr::filter(Fisher < input$pvalue) |>
        dplyr::mutate(Term = Term |> forcats::fct_relevel(Term)) |>
        dplyr::arrange(GO.Domain, Fisher) |>
        dplyr::mutate(
          plotLabels = glue::glue("{Term} ({GO.ID})"),
          plotLabels = factor(plotLabels, levels = plotLabels)
          ) |>
        ggplot2::ggplot(ggplot2::aes(x = -log10(Fisher), y = rev(plotLabels), fill = GO.Domain)) +
        ggplot2::geom_bar(stat = "identity", color = "black") +
        ggplot2::geom_vline(xintercept = -log10(input$pvalue), linetype = "dashed", color = "gray10") +
        ggplot2::scale_fill_discrete(type = c("#2386A7", "#00A352", "#ED8E28")) +
        ggplot2::labs(x = "-log10(p-value)", y = "", fill = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text = ggplot2::element_text(color = "black", size = 12),
          legend.text = ggplot2::element_text(size = 14)
        )
    })

    # Outputting the plot
    output$ontologyplot <- renderPlot({goplot()}, res = 96)

    #prepare the current visualized table for download
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
