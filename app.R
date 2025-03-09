#
# Data science class
#
# University of Cincinnati/Cincinnati Children's
#
# Demonstrate file upload, creating a signature, submit to iLincs for correlations
#
#

library(shiny)
library(DT)
library(httr)
library(tidyverse)

options(shiny.maxRequestSize = 70 * 1024^2)



ui <- fluidPage(

  # Application title
  titlePanel("File upload"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose TSV File", accept = c("text/tsv", ".tsv")),
      selectInput("variable", "Grouping Variable:", choices = c()),
      selectizeInput("group1", "Group1", choices = NULL, multiple = TRUE),
      selectizeInput("group2", "Group2", choices = NULL, multiple = TRUE),
      selectInput("difffunction", "Differential function:", choices = c("t-test")),
      checkboxInput("usel1000", label = "Use L1000 genes", value = TRUE),
      checkboxInput("top100", label = "select only the top 100 differentially expressed genes", value = TRUE),
      actionButton("compute", "Compute and Submit Signature")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      dataTableOutput("signature_data"),
      hr(),
      dataTableOutput("correlated_data"),
      hr(),
      plotOutput("plot"),
      hr(),
      dataTableOutput("sample_data"),
      hr(),
      # plotOutput("sample_data_heatmap"),
    )
  )
)


server <- function(input, output, session) {
  values <- reactiveValues(data = NULL)

  observe({
    # handle file upload

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.

    inFile <- input$file1

    if (is.null(inFile)) {
      return(NULL)
    }

    isolate({
      file <- (inFile$datapath)

      values$header <- scan(file, nlines = 1, sep = "\t", what = character())
      values$data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names = FALSE)
      names(values$data) <- values$header
      values$header2 <- data.frame(scan(file, skip = 1, nlines = 1, sep = "\t", what = character()))
    })
  })

  # show sample metadata
  output$sample_data <- renderDataTable(
    {
      head(values$data[, 1:10], n = 50)
    },
    caption = "First 50 genes, first 10 samples"
  )
  # heatmap is not performant enough to work, locks up the session
  # output$sample_data_heatmap <- renderPlot({
  #   if (is.null(values$data)) {
  #     return(NULL)
  #   }
  #   d <- values$data
  #   d <- d[!duplicated(d$id), ]
  #   rownames(d) <- d$id
  #   d$id <- NULL
  #   heatmap(as.matrix(d))
  # })

  observe({
    # fill in variable selection from the input file
    updateSelectInput(session, "variable", choices = as.character(values$header2[1, ]))
  })

  # handle UI group values update based on selected variable
  observe({
    if (input$variable != "") {
      updateSelectizeInput(session, "group1", choices = unique(values$header2[-1, values$header2[1, ] == input$variable]), server = TRUE)
      updateSelectizeInput(session, "group2", choices = unique(values$header2[-1, values$header2[1, ] == input$variable]), server = TRUE)
    }
  })

  # Create signature, upload to API, display results
  observeEvent(input$compute, {
    withProgress(message = "Creating signature", value = 0, {
      #
      # Filter into two groups
      #
      group1 <- names(values$data)[values$header2 == input$group1]
      group2 <- names(values$data)[values$header2 == input$group2]

      #
      # ensure the values are numeric
      #
      values$values <- values$data[complete.cases(values$data), -1]
      values$values[] <- lapply(values$values, function(x) {
        as.numeric(as.character(x))
      })

      #
      # select differential function
      #
      incProgress(1 / 3, detail = paste0("Running ", input$difffunction))
      if (input$difffunction == "t-test") {
        diff_result <- as.data.frame(apply(
          values$values, 1,
          function(x) {
            t.test(
              unlist(x[group1], use.names = FALSE),
              unlist(x[group2], use.names = FALSE)
            )$p.value
          }
        ))

        # FIXME: fill in
        # Add t statistic as the measure of differential expression (Value_LogDiffExp)
        diff_result$Value_LogDiffExp <- apply(
          values$values, 1,
          function(x) {
            t.test(
              unlist(x[group1], use.names = FALSE),
              unlist(x[group2], use.names = FALSE)
            )$statistic
          }
        )
      }

      #
      # format signature output
      #
      output_id_column_name <- paste0(values$header[1], "_GeneSymbol")
      diff_result <- data.frame(values$data[values$header[1]], diff_result[, 1:2])
      colnames(diff_result) <- c(output_id_column_name, "Significance_pvalue", "Value_LogDiffExp")



      diff_result <- diff_result[, c(1, 3, 2)]
      # FIXME: fill in
      # choose only L1000 genes
      if (input$usel1000) {
        l1000genes <- read.csv2("L1000.txt", sep = "\t")
        l1000genes <- l1000genes[l1000genes$pr_is_lm == "1", ]$pr_gene_symbol
        diff_result <- diff_result[diff_result[, 1] %in% l1000genes, ]
      }

      # FIXME: fill in
      # choose top 100 most differentially expressed genes
      if (input$top100) {
        diff_result <- head(diff_result[order(diff_result$Value_LogDiffExp, decreasing = TRUE), ], n = 100)
      }
      #
      # show signature in a table
      #
      output$signature_data <- DT::renderDataTable(
        {
          diff_result
        },
        caption = "Signature to submit to iLincs"
      )

      incProgress(1 / 3, detail = paste("Submitting the signature to iLincs"))

      #
      # create temporary csv file to submit into API
      #
      ftemp <- tempfile(pattern = "file", fileext = ".csv", tmpdir = tempdir())
      write.csv(diff_result, ftemp, row.names = F, quote = F)
      cat(ftemp)

      r <- POST("https://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5", body = list(file = upload_file(ftemp)))

      l <- lapply(content(r)$status$concordanceTable, function(x) unlist(x))
      ilincs_result <- data.frame(t(sapply(l, c)))

      #
      # show correlation results
      #
      output$correlated_data <- DT::renderDataTable({
        datatable(ilincs_result, rownames = TRUE, caption = "Correlated signatures")
      })

      ilincs_result$similarity <- as.numeric(ilincs_result$similarity)
      output$plot <- renderPlot({
        ggplot(ilincs_result, aes(x = similarity, y = "")) +
          geom_violin(drop = FALSE) +
          ylab("") +
          title("Distribution of similarity to calculated signature")
      })
    })
  })
}

shinyApp(ui = ui, server = server)
