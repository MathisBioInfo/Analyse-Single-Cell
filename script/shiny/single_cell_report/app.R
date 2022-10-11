library(shiny)
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

setwd("/home/stagiaire/Documents/Mathis")
load("data/merged_filtered_seurat.RData")
load("data/seurat_filtered.RData")
metadata <- merged_seurat@meta.data
filtered_metadata <- filtered_seurat@meta.data

# Define UI ----
ui <- fluidPage(
    titlePanel("Single Cell Analysis"),
    sidebarLayout(
        sidebarPanel(
            h2("Quality Control"),
            selectInput("select", strong("Choose a variable to explore"), 
                        choices = list("number of cells" = 1, "number of UMI" = 2), selected = 1)
        ),
        mainPanel(
            plotOutput("selected_plot")
        )
    )
    
)

# Define server logic ----
server <- function(input, output) {
    output$quality_plot <- renderPlot({
        metadata %>% 
            ggplot(aes(x=sample, fill=sample)) + 
            geom_bar() +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
            theme(plot.title = element_text(hjust=0.5, face="bold")) +
            ggtitle("NCells")
    })
    
}

# Run the app ----
shinyApp(ui = ui, server = server)