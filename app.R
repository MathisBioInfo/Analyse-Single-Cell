## app.R ##
library(dplyr)
library(ggplot2)
library(shiny)
library(shinydashboard)

options(shiny.maxRequestSize = 1000*1024^2)


ui <- dashboardPage(
    dashboardHeader(title = "Single Cell Analysis"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Quality control and filters", tabName = "QC", icon = icon("filter",lib="glyphicon"))
        )
    ),
    dashboardBody(
        #upload box
        fluidRow(
            box(
                title = "Files",
                fileInput("upload","Upload a Seurat object:",accept = ".RData"))
        ),
        #graph panel
        fluidRow(
                tabBox(
                    tabPanel("nCells",
                        plotOutput("nCells_plot"),
                    ),
                    tabPanel("nUMIs",
                             plotOutput("nUMIs_plot")
                    ),
                    tabPanel("nGenes",
                             plotOutput("nGenes_plot")
                    ),
                    tabPanel("genesExpr",
                             plotOutput("genesExpr_plot")
                    ),
                    tabPanel("nMito",
                             plotOutput("nMito_plot")
                    ),
                    tabPanel("Profile",
                             plotOutput("Profile_plot")
                    )
                ),#tabBox end
        box(
            title = "Filters",
            numericInput('transcript', 'Transcript number', 0, min = 0, max = 10^6, step = 1000), #a ajouter au graph
            numericInput('gene', 'Gene number', 0, min = 0, max = 10^4, step = 100),
            numericInput('ratio', 'Mito ratio', 1, min = 0, max = 1, step = 0.01),
            #verbatimTextOutput("code"),
            actionButton("click", "Apply"))#a compléter
        )#graph panel end
    )#Body end
)#Page end

server <- function(input, output, session) {
    metadata <- reactive({
        req(input$upload)
        file <- input$upload
        ext <- tools::file_ext(file$datapath)
        validate(need(ext == "RData", "Please upload a RData file"))
        load(file$datapath)
        merged_seurat@meta.data
        
        #Test filtrage
        filtered_seurat <- subset(x = merged_seurat, 
                                  subset= (nUMI >= f_transcript()) & 
                                      (nGene >= f_gene()) & 
                                      (mitoRatio < f_ratio())) #test
        filtered_seurat@meta.data
    })
    
    # Variable from the various filters
    t_ratio <- reactive({
        input$ratio
    })
    
    f_ratio <- eventReactive(input$click, {
        input$ratio
    })
    
    t_transcript <- reactive({
        input$transcript
    })
    
    f_transcript <- eventReactive(input$click, {
        input$transcript
    })
    
    t_gene <- reactive({
        input$gene
    })
    
    f_gene <- eventReactive(input$click, {
        input$gene
    })
    
    # Display filtering information
    output$code <- renderPrint({ 
        summary(1:10) #je veux afficher le pourcentage de cellule filtré par les filtres
    })
    
    # Visualize the number of cell counts per sample
    output$nCells_plot <- renderPlot({
        metadata() %>% 
            ggplot(aes(x=sample, fill=sample)) + 
            geom_bar() +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
            theme(plot.title = element_text(hjust=0.5, face="bold")) +
            ggtitle("NCells")
    })
    
    # Visualize the number UMIs/transcripts per cell
    output$nUMIs_plot <- renderPlot({
        metadata() %>% 
            ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
            geom_density(alpha = 0.2) + 
            scale_x_log10() + 
            theme_classic() +
            ylab("Cell density") +
            geom_vline(xintercept = t_transcript())
    })
    
    # Visualize the distribution of genes detected per cell via histogram
    output$nGenes_plot <- renderPlot({
        metadata() %>% 
            ggplot(aes(color=sample, x=nGene, fill= sample)) + 
            geom_density(alpha = 0.2) + 
            theme_classic() +
            scale_x_log10() +
            geom_vline(xintercept = t_gene())
    })
    
    # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
    output$genesExpr_plot <- renderPlot({
        metadata() %>%
            ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
            geom_density(alpha = 0.2) +
            theme_classic()
    })
    
    # Visualize the distribution of mitochondrial gene expression detected per cell
    output$nMito_plot <- renderPlot({
        metadata() %>% 
            ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
            geom_density(alpha = 0.2) + 
            scale_x_log10() + 
            theme_classic() +
            geom_vline(xintercept = t_ratio())
    })
    
    # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
    output$Profile_plot <- renderPlot({
        metadata() %>% 
            ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
            geom_point() + 
            #scale_colour_gradient(low = "gray90", high = "black") +
            scale_colour_gradient2(low = "blue4", high = "red", mid = "white", midpoint = t_ratio()) +
            #stat_smooth(method=lm) +
            scale_x_log10() + 
            scale_y_log10() + 
            theme_classic() +
            facet_wrap(~sample)
    })
}

shinyApp(ui, server)
