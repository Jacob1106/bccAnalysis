#' @title Shiny App for Chip-Seq Analysis
#' @param named_file_list A named list with file paths for narrowPeak files. These files will be imported as GRanges Objects. 
#' @return Returns an Analysis App. Outputs: Overlaps, Unique peak annotations, width analysis and proportion overlap table   
#' @import shiny
#' @import bslib
#' @import bccAnalysis
#' @import GenomicRanges
#' @import rtracklayer
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import org.Hs.eg.db
#' @import ChIPseeker
#' @import tidyverse
#' @import ChIPpeakAnno
#' @import ggvenn
#' @import clusterProfiler
#' @import ggupset
#' @import ReactomePA
#' @import plyranges
#' @import writexl
#' @import readxl
#' @import dplyr
#' @import patchwork
#' @import qpdf
#' @import gridExtra
#' @import gridGraphics
#' @import grid
#' @author Jacob Martin
#' @export



ShinyApp <- function(named_file_list, height = 850) {

    options(shiny.maxRequestSize = 50 * 1024^2)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    ui <- navbarPage("My Genomic App",
        # File Input UI: Contains file input, index assignment and overlap requirements 
        tabPanel("Load Files",
        card(
            card_header("Please Insert a List of GRange Objects in a .rds format"),
            fileInput("GRangeList", label = NULL)
        ),
        card(
            card_header("Select Which objects to use as object1. Objects will be filtered and merged into one."), # nolint
            selectInput("index_object1", "Select integers:",
                choices = 1:10,
                multiple = TRUE)
        ),
        card(
            card_header("Select Which objects to use as object2. Objects will be filtered and merged into one."), # nolint
            selectInput("index_object2", "Select integers:",
                choices = 1:10,
                multiple = TRUE)
        ),
        card(
            card_header("Input the number of Overlaps You require for filtering in Object 1"),
            numericInput("overlaps1", "Input number", value = 0)
        ),
        card(
            card_header("Input the number of Overlaps You require for filtering in Object 2"),
            numericInput("overlaps2", "Input number", value = 0)
        ),
        card(
            card_header("Status"),
            textOutput("FileStatus")
        )
        ),


    # Venn Diagram UI

        tabPanel("Venn Diagram",
            sidebarLayout(
                sidebarPanel(),
        
                mainPanel(
                    card(
                    card_header("Overlap Of Objects"),
                    plotOutput("Venn") 
                    )
                )
            )
        ),

    #Tibble Output: 

        tabPanel("Table of Overlaps:",
            sidebarLayout(
                sidebarPanel(
                    downloadButton("downloadTable", "Download Plot as excel")
                ),

                mainPanel(
                    card(
                        card_header("Overlap Table between Objects"),
                        tableOutput("overlaps_tibble")
                    )
                )
            )
        ),

        #Width Comparison UI 
        tabPanel("Width Comparison Plot",
            sidebarLayout(
                sidebarPanel(
                    downloadButton("downloadPlotsWidth", "Download Plot as PDF")
                ),

                mainPanel(
                    card(
                        card_header("Peak Width Comparison"),
                        plotOutput("width")
                    )
                )
            )
        ),

        #Peak annotation UIs: 3 graphs per Tab. 2 TABS 
        tabPanel("Object 1 Peak Analysis",
            sidebarLayout(
                sidebarPanel(
                    downloadButton("downloadPlotsPeak1", "Download Plots as PDF")

                ),

                mainPanel(
                    card(
                        card_header("Object 1 Peak Analysis"),
                        plotOutput("Object1_Main"), 
                        plotOutput("Object1_Unique"),
                        plotOutput("Object1_Overlap")
                    )
                )
            )
        ),
        tabPanel("Object 2 Peak Analysis",
            sidebarLayout(
                sidebarPanel(
                    downloadButton("downloadPlotsPeak2", "Download Plots as PDF")

                ),

                mainPanel(
                    card(
                        card_header("Object 2 Peak Analysis"),
                        plotOutput("Object2_Main"), 
                        plotOutput("Object2_Unique"),
                        plotOutput("Object2_Overlap")
                    )
                )
            )
        )
    )




    server <- function(input, output) ({
        processedData <- reactive({
            #Checks if a file has been uploaded. 
            req(input$GRangeList)
            # Upstream Filtering and Processing
            datalist <- readRDS(input$GRangeList$datapath)

            # Return All processed Items in a named List: e.g. List(Unique1 = ...., Unique2 = ...., Overlap1 = ...., Overlap2 = ....)

            #Creates a variable with the indexes of the object in the original named list # nolint
            object1_index <- as.integer(input$index_object1)
            object2_index <- as.integer(input$index_object2) #nolint 

            # Creating the two objects: 
            object1 <- datalist[object1_index]
            object2 <- datalist[object2_index] 

            #Merging and Filtering: 
            #Using input$overlaps1 and 2 

            Object1Merged_GRL <- GRangesList(object1)
            object1_M <- overlapGRangeList(Object1Merged_GRL, num_of_overlaps_required = as.integer(input$overlaps1))

            Object2Merged_GRL <- GRangesList(object2)
            object2_M <- overlapGRangeList(Object2Merged_GRL, num_of_overlaps_required = as.integer(input$overlaps2))

            #Filtering By Chromosome: 
            chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
            object1_M <- object1_M[seqnames(object1_M) %in% chromosomes]
            object2_M <- object2_M[seqnames(object2_M) %in% chromosomes]

            #Returning Data as 
            initialData <- list(object1_M = object1_M, object2_M = object2_M)
            return(initialData)  # return processed/reactive data
        })

        #Checking if a file has been loaded 
        output$FileStatus <- renderText({
            if (is.null(input$GRangeList)) "Please upload a file" else "File loaded successfully"
        })

        #Extracting Unique and Overlapped peaks: 
        DifferentialPeaks <- reactive({
            data <- processedData()
            overlaps <- findOverlaps(data$object1_M, data$object2_M)

            #Seperating into unique and overlapped peak sets: 
            #Unique Peaks: 
            uniqueObject1 <- data$object1_M[-queryHits(overlaps)]
            uniqueObject2 <- data$object2_M[-subjectHits(overlaps)]

            #Overlapped: 
            overlapObject1 <- data$object1_M[queryHits(overlaps)]
            overlapObject2 <- data$object2_M[subjectHits(overlaps)]

            #Differential Peak Data output
            DiffPeakList <- list(uniqueObject1 = uniqueObject1, uniqueObject2 = uniqueObject2, overlapObject1 = overlapObject1, overlapObject2 = overlapObject2)
            return(DiffPeakList)
        })

        # Venn diagram Output 
        output$Venn <- renderPlot({
            data1 <- processedData()
            custom_venn_diagram(conditionA = "CHIP", conditionB = "CNR", DataA = data1$object1_M, DataB = data1$object2_M)
        })

    # Overlap Table output 
        output$overlaps_tibble <- renderTable({
            data1 <- processedData()
            data2 <- DifferentialPeaks()

            OverlapTable <- proportionOverlapTibble(object1 = data1$object1_M, object2 = data1$object2_M, overlapobject1 = data2$overlapObject1, overlapobject2 = data2$overlapObject2, object1_unique = data2$uniqueObject1, object2_unique = data2$uniqueObject2)
            OverlapTable
        })
    #Width Diagram Output
        output$width <- renderPlot({
            data1 <- processedData()
            WidthPlot <- compareGRangeWidth(Obj1 = data1$object1_M, Obj2 = data1$object2_M)
            WidthPlot
        })
        PeakAnalysis <- reactive({
            data1 <- processedData()
            data2 <- DifferentialPeaks()
            peakAnno_Object1 <- annotatePeak(data1$object1_M, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
            peakAnno_Object2 <- annotatePeak(data1$object2_M, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")

            #Unique Peak Annotations: 

            unique_peak_anno_Object1 <- annotatePeak(data2$uniqueObject1, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
            unique_peak_anno_Object2 <- annotatePeak(data2$uniqueObject2, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")

            overlapRegions1 <- annotatePeak(data2$overlapObject1, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
            overlapRegions2 <- annotatePeak(data2$overlapObject2, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")

            #Output List 
            outputList <- list(peakAnno_Object1 = peakAnno_Object1, peakAnno_Object2 = peakAnno_Object2, unique_peak_anno_Object1 = unique_peak_anno_Object1, unique_peak_anno_Object2 = unique_peak_anno_Object2, overlapRegions1 = overlapRegions1, overlapRegions2 = overlapRegions2)
            return(outputList)
        })

        output$Object1_Main <- renderPlot({
            data <- PeakAnalysis()
            plotAnnoBar(data$peakAnno_Object1)
        })

        output$Object1_Unique <- renderPlot({
            data <- PeakAnalysis()
            plotAnnoBar(data$unique_peak_anno_Object1)
        })

        output$Object1_Overlap <- renderPlot({
            data <- PeakAnalysis()
            plotAnnoBar(data$overlapRegions1)
        })

        output$Object2_Main <- renderPlot({
            data <- PeakAnalysis()
            plotAnnoBar(data$peakAnno_Object2)
        })

        output$Object2_Unique <- renderPlot({
            data <- PeakAnalysis()
            plotAnnoBar(data$unique_peak_anno_Object2)
        })

        output$Object2_Overlap <- renderPlot({
            data <- PeakAnalysis()
            plotAnnoBar(data$overlapRegions2)
        })

        output$downloadPlotsPeak1 <- downloadHandler(
            filename = function() {
                paste("plots-", Sys.Date(), ".pdf", sep = "")
            },
            content = function(file) {
                pdf(file, width = 12, height = 4)  
            
                data <- PeakAnalysis() 
                
                plot1 <- plotAnnoBar(data$peakAnno_Object1)
                plot2 <- plotAnnoBar(data$unique_peak_anno_Object1)
                plot3 <- plotAnnoBar(data$overlapRegions1)
                
                gridExtra::grid.arrange(plot1, plot2, plot3, ncol = 3)
                
                dev.off()
            }
        )
        output$downloadPlotsPeak2 <- downloadHandler(
            filename = function() {
                paste("plots-", Sys.Date(), ".pdf", sep = "")
            },
            content = function(file) {
                pdf(file, width = 12, height = 4)  
            
                data <- PeakAnalysis() 
                
                plot1 <- plotAnnoBar(data$peakAnno_Object1)
                plot2 <- plotAnnoBar(data$unique_peak_anno_Object1)
                plot3 <- plotAnnoBar(data$overlapRegions1)
                
                gridExtra::grid.arrange(plot1, plot2, plot3, ncol = 3)
                
                dev.off()
            }
        )

        output$downloadTable <- downloadHandler(
            filename = function() {
                paste("overlap_table-", Sys.Date(), ".xlsx", sep = "")
            },
            content = function(file) {
                data1 <- processedData()
                data2 <- DifferentialPeaks()

        
                overlap_table <- proportionOverlapTibble(object1 = data1$object1_M, object2 = data1$object2_M, overlapobject1 = data2$overlapObject1, overlapobject2 = data2$overlapObject2, object1_unique = data2$uniqueObject1, object2_unique = data2$uniqueObject2)

                writexl::write_xlsx(overlap_table, path = file)
            }
        )

        output$downloadPlotsWidth <- downloadHandler(
            filename = function() {
                paste("widthPlot-", Sys.Date(), ".pdf", sep = "")
            },
            content = function(file) {
                pdf(file, width = 12, height = 4)
                data1 <- processedData()
                compareGRangeWidth(Obj1 = data1$object1_M, Obj2 = data1$object2_M)
                dev.off()
            }
        )
    })

shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
}

#Example: 
# data <- processeData()
# data1 <- data$object1_M









