# This module was created during the St Jude Bio-Hackathon of May 2023 by the team 13.
# author: Christy LaFlamme (christy.laflamme@stjude.org)
# the plot created is adapated from a script originally written by Paras Garg

# Documentation
#' R Shiny module to generate volcano plot
#'
#' @param data A dataframe.
#' @returns A Shiny module.
#' @examples
#' DMR_demo()
#### Library needed #### ----------
library(shiny)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)

#### Function needed to work #### ----------
#' Differentially methylated region (DMR) plot
#'
#' @param data A data frame of B-values from DNA methylation array (eg. Illumina EPIC 850k); must be in ascending order by genomic position
#' @param coordinates Genomic coordinates to plot (e.g. chr16:17562418-17565190)
#' @param chrCol The name of the column containing the chromosome (e.g. chr16)
#' @param startCol The name of the column containing the starting/beginning CpG coordinate (e.g. 17562418)
#' @param endCol The name of the column containing the stopping/ending CpG coordinate (e.g. 17565190)
#' @param intSamples The names of the samples of interest to highlight (e.g. 203866380012_R06C01, 204873630063_R08C01)
#' @param title The title of the plot (if desired)
#' @return plot
#' @export pdf

# Plotting function of the module
plotDMR <- function(data, coordinates, chrCol, startCol, endCol, intSamples, title = "") {
  # get the sample names (column names of data matrix satisfying the following regular expression)
  samples <- str_extract(colnames(data), "............_R0.C0.") # recognizes the sentrix ID format
  samples <- samples[!is.na(samples)]

  # string split the given coordinates to pull the chromosome, start, and stop information
  coordinates <- gsub(",", "", coordinates) # remove commas
  chr <- strsplit(coordinates, split = ":")[[1]][1]
  start <- strsplit(strsplit(coordinates, split = ":")[[1]][2], split = "-")[[1]][1]
  stop <- strsplit(strsplit(coordinates, split = ":")[[1]][2], split = "-")[[1]][2]

  # acquire the row indexes of the data matrix that correspond to the start and stop coordinates

  # cut the data to match with chromosome specificity of the given coordinates
  data <- data[data[, data[[chrCol]] == chr], ]

  # this is with the minimum distance from the given coordiantes to the data +1 / -1 probe
  startPos <- which(abs(data[[startCol]] - as.numeric(start)) == min(abs(data[[startCol]] - as.numeric(start)))) - 1
  endPos <- which(abs(data[[endCol]] - as.numeric(stop)) == min(abs(data[[endCol]] - as.numeric(stop)))) + 1

  # make sure that interest sample list is separated out from commas and remove spaces
  intSamples <- gsub(" ", "", strsplit(intSamples, ",")[[1]])

  # pre-processing of data
  data <- data %>%
    slice(startPos:endPos) %>%
    melt(id.vars = setdiff(colnames(data), samples), variable.name = "sampleID", value.name = "Beta") %>%
    mutate(status = ifelse(sampleID %in% intSamples, "Sample-of-Interest", "normal"))

  # create the line plot
  g <- ggplot(data = data, aes_string(x = startCol, y = "Beta")) +
    geom_line(aes(group = sampleID, size = status, linetype = status, colour = status)) +
    scale_colour_manual(values = c("Sample-of-Interest" = "red", "normal" = "black")) +
    scale_size_manual(values = c("Sample-of-Interest" = 0.9, "normal" = 0.2)) +
    scale_linetype_manual(values = c("Sample-of-Interest" = "solid", "normal" = "dashed")) +
    coord_cartesian(ylim = c(0, 1)) +
    xlab(paste0("Genomic Position at ", chr)) +
    ylab("Methylation Value") +
    labs(title = title) +
    theme_bw()
  return(g)
}

#### UI function of the module #### ----------
DMR_ui <- function(id) {
  ns <- NS(id)
  tagList(
    titlePanel("Interactive DMR plot"),
    textInput(ns("coordinates"), "What genomic coordinates would you like to view? Please use format chrX:start-stop ", value = "chr22:38711890-38714427"), #chr16:17562418-17565190 and chrX:20008363-20009139
    textInput(ns("chrCol"), "What is the name of the column with the chromosome information? (e.g. chr16)", value = "CpG_chrm"),
    textInput(ns("startCol"), "What is the name of the column with the start positions? (e.g. 17562418)", value = "CpG_beg"),
    textInput(ns("endCol"), "What is the name of the column with the end positions? (e.g. 17565190)", value = "CpG_end"),
    textInput(ns("intSamples"), "Are there samples of interest to highlight (e.g. 203866380012_R06C01, 204873630063_R08C01)?", value = ""), #203866380012_R06C01, 204873630063_R08C01
    textInput(ns("title"), "What is the title of the plot?", value = "DMR plot at CSNK1E"), #XYLT1
    verbatimTextOutput(ns("region")),
    plotOutput(ns("plot"))
  )
}

#### Server function of the module #### ----------
DMR_server <- function(id, df) {
  moduleServer(id, function(input, output, session) {
    output$region <- renderText(paste0("The region you have selected is: ", input$coordinates, " and the samples of interest are: ", input$intSamples))

    DMR_plot <- reactive({
      plotDMR(df(), input$coordinates, input$chrCol, input$startCol, input$endCol, input$intSamples, title = input$title)
    })

    output$plot <- renderPlot({
      DMR_plot()
    })

    return(DMR_plot)
  })
}

df <- fread("data-file")

#### Demo function of the module #### ----------
DMR_demo <- function() {
  # df <- fread("/Users/claflamm/OneDrive - St. Jude Children's Research Hospital/Mefford/methylation/data/autosomes.beta.txt.sorted") # example data
  ui <- fluidPage(
    DMR_ui("dmrplot")
  )
  server <- function(input, output, session) {
    DMR_server("dmrplot", reactive({
      df
    }))
  }
  shinyApp(ui, server)
}

DMR_demo()
