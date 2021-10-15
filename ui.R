library(shiny)
library(plotly)
library(ggplot2)
library(shinyjs)
library(DT)
library(shinyBS)
library(Seurat)
library(SeuratData)
library(cowplot)
library(SeuratWrappers)
library(dplyr)
library(pbmcapply)
library(harmony)
library(rliger)
library(reshape2)
library(shinydashboard)
library(shinyalert)
library(shinyFiles)
library(shinyWidgets)

shiny_one_panel = fluidPage(
    titlePanel("Seurat analysis of scRNAseq data"),
    hr(),

    fluidRow(
        ##------Sidebar---------
        column(3,
               h4('Load Data:'),
               wellPanel(
                   fileInput(inputId = 'tpmFiles',
                             label = "Gene expression file",
                             multiple = FALSE,
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".Robj")),
                   fileInput(inputId = 'cellAnnoFiles',
                             label = "Metadata",
                             multiple = FALSE,
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain")),
                   checkboxInput(inputId = 'norm',
                                 label = "Normalise?",
                                 value = TRUE),
                   fluidRow(
                       column(6,
                              numericInput(inputId = "min.genes",
                                           label = "Min. genes",
                                           value = 200,
                                           min = 1)
                       ),
                       column(6,
                              numericInput(inputId = "min.cells",
                                           label = "Min. cells",
                                           value = 3,
                                           min = 1)
                       )
                   ),
                   textInput(inputId = "projName",
                             label = "Project Name",
                             value = "Seurat_analysis"),
                   fluidRow(
                       column(6,
                              actionButton("loadButton", "Create  Seurat Object", icon = icon("hand-o-right"))
                       ),
                       column(6,
                              actionButton("reset", "Reset Data", icon = icon("repeat"))
                       )
                   )
               ),

               ##------Plot download---------
               h4("Export to PDF:"),
               wellPanel(
                   ## Conditional panel for different plots
                   conditionalPanel(" input.tabs == 'Data Integration using Seurat' ",
                                    actionButton("PDFl", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'Data Integration using Harmony' ",
                                    actionButton("PDFm", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'Data Integration using RLiger' ",
                                    actionButton("PDFn", "Download", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'scMultiomics using Seurat' ",
                                    actionButton("PDFo", "Download", icon = icon("download"))
                   ),
                   ## ensure no spill over in button text
                   tags$head(
                       tags$style(HTML('
                                   .btn {
                                   white-space: normal;
                                   }'
                       )
                       )
                   ),
                   ## Conditional is separate from pdf options
                   hr(),
                   fluidRow(
                       column(6,
                              sliderInput(inputId="pdf_w", label = "PDF width(in):",
                                          min=3, max=20, value=8, width=100, ticks=F)
                       ),
                       column(6,
                              sliderInput(inputId="pdf_h", label = "PDF height(in):",
                                          min=3, max=20, value=8, width=100, ticks=F)
                       )),

                   #actionButton("OpenDir", "Open download folder", icon = icon("folder"))
               ),

               ##------Save Data---------
               hr(),
               actionButton("saveButton", "Save Data", icon = icon("hand-o-right")),

               hr(),
               h4(tags$a(href="mailto:Chen_Jinmiao@immunol.a-star.edu.sg?subject=[cytof-question]",
                         "Contact Us")),
               imageOutput("logo", height = "60px")
        ),
        ##------Main area---------
        column(9,
               tabsetPanel(type = "pills", id = "tabs",

                           ##------Data Integration using Seurat---------

                           tabPanel("Data Integration using Seurat", fluidPage(
                               hr(),
                               fluidRow(
                                   #br(),
                                   fluidRow(
                                       fluidRow(
                                       column(3,
                                              br(),
                                              actionButton("doIntg", "Running Data Integration", icon = icon("hand-pointer-o")),
                                              br(),
                                       ),
                                       fluidRow(
                                           uiOutput("ident.swap1")
                                       ),
                                       #selectInput("dim.used",
                                       #            label = "Dimensions",
                                       #            choices = c(1:50)
                                       #),
                                       #br(),
                                       #fluidRow(
                                       column(9,
                                              br(),
                                              actionButton("runPCA", "Run PCA", icon = icon("hand-pointer-o")),
                                              br(),
                                              #column(6,
                                              br(),
                                              plotlyOutput("PCAplot_a", width = "100%"),
                                              br(),
                                              plotlyOutput("PCAplot_b", width = "100%"),
                                              br(),
                                              plotlyOutput("PCAplot_c", width = "100%")),
                                       br(),
                                       #),
                                       #),
                                       #selectInput("assays2",
                                       #            label = "Default Assay",
                                       #            choices = c("RNA", "Integrated")
                                       #),
                                       #column(6,
                                       #       br(),
                                       #       actionButton("doCluster", "Find Clusters", icon = icon("hand-pointer-o")),
                                       #       br(),
                                       #column(6,
                                       #       br(),
                                       #       plotlyOutput("Clusterplot", width = "100%"),
                                       #br(),
                                       #),
                                       column(9,
                                              br(),
                                              actionButton("runUMAP", "Run UMAP", icon = icon("hand-pointer-o")),
                                              #textOutput("Intg.done"),
                                              br(),
                                              #),
                                              #column(6,
                                              br(),
                                              plotlyOutput("UMAPplot_a", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot_b", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot_c", width = "100%")),
                                       br(),
                                       #),
                                       column(9,
                                              br(),
                                              actionButton("runTSNE", "Run TSNE", icon = icon("hand-pointer-o")),
                                              #textOutput("Intg.done"),
                                              br(),
                                              br(),
                                              plotlyOutput("TSNEplot_a", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot_b", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot_c", width = "100%")),
                                       br(),
                                       ),
                                   ))
                           )),


                           ##------Data Integration using Harmony---------

                           tabPanel("Data Integration using Harmony", fluidPage(
                               hr(),
                               fluidRow(
                                   #br(),
                                   fluidRow(
                                       fluidRow(
                                           column(3,
                                                  br(),
                                                  actionButton("doIntg1", "Running Data Integration", icon = icon("hand-pointer-o")),
                                                  br(),
                                           ),
                                           fluidRow(
                                               uiOutput("ident.swap2")
                                           ),
                                           column(9,
                                                  br(),
                                                  actionButton("runPCA1", "Run PCA", icon = icon("hand-pointer-o")),
                                                  br(),
                                                  #column(6,
                                                  br(),
                                                  plotlyOutput("PCAplot1_a", width = "100%"),
                                                  br(),
                                                  plotlyOutput("PCAplot1_b", width = "100%"),
                                                  br(),
                                                  plotlyOutput("PCAplot1_c", width = "100%")),
                                           br(),
                                           column(9,
                                                  br(),
                                                  actionButton("runUMAP1", "Run UMAP", icon = icon("hand-pointer-o")),
                                                  br(),
                                                  br(),
                                                  plotlyOutput("UMAPplot1_a", width = "100%"),
                                                  br(),
                                                  plotlyOutput("UMAPplot1_b", width = "100%"),
                                                  br(),
                                                  plotlyOutput("UMAPplot1_c", width = "100%")),
                                           br(),
                                           column(9,
                                                  br(),
                                                  actionButton("runTSNE1", "Run TSNE", icon = icon("hand-pointer-o")),
                                                  br(),
                                                  br(),
                                                  plotlyOutput("TSNEplot1_a", width = "100%"),
                                                  br(),
                                                  plotlyOutput("TSNEplot1_b", width = "100%"),
                                                  br(),
                                                  plotlyOutput("TSNEplot1_c", width = "100%")),
                                           br(),
                                       ),
                                   ))
                           )),

                           ##------Data Integration using RLiger---------

                           tabPanel("Data Integration using RLiger", fluidPage(
                               hr(),
                               fluidRow(
                                   #br(),
                                   fluidRow(
                                       fluidRow(
                                           column(3,
                                                  br(),
                                                  actionButton("doIntg2", "Running Data Integration", icon = icon("hand-pointer-o")),
                                                  br(),
                                           ),
                                           fluidRow(
                                               uiOutput("ident.swap3")
                                           ),
                                           column(9,
                                                  br(),
                                                  actionButton("runUMAP2", "Run UMAP", icon = icon("hand-pointer-o")),
                                                  br(),
                                                  br(),
                                                  plotlyOutput("UMAPplot2_a", width = "100%"),
                                                  br(),
                                                  plotlyOutput("UMAPplot2_b", width = "100%")),
                                            br(),
                                           column(9,
                                                  br(),
                                                  actionButton("runTSNE2", "Run TSNE", icon = icon("hand-pointer-o")),
                                                  br(),
                                                  br(),
                                                  plotlyOutput("TSNEplot2_a", width = "100%"),
                                                  br(),
                                                  plotlyOutput("TSNEplot2_b", width = "100%")),
                                            br(),
                                       ),
                                   ))
                           )),


                           ##------END---------
               )
        )
    )
)

dbHeader <- dashboardHeader(title = "cytofkit2")
dbHeader$children[[2]]$children <-  tags$a(href='https://github.com/JinmiaoChenLab',
                                           tags$img(src='https://avatars1.githubusercontent.com/u/8896007?s=400&u=b0029c2e64f405ea0a46d311239b674a430ec77c&v=4'
                                                    ,height='60',width='60', align='left')
                                           , tags$div('cytofkit2', style='color:white;font-family:arial rounded MT bold'))

dashboardPage(skin = "yellow",
              dbHeader,
              dashboardSidebar(
                  sidebarMenu(id = "sbm",
                              menuItem(tags$p(style = "display:inline;font-size: 20px;", "Seurat"), tabName = "seurat", icon = icon('cog'))


                  )# end of sidebarMenu
              ),#end of dashboardSidebar
              dashboardBody(
                  #includeCSS("www/custom.css")
                  useShinyalert()
                  , shinyjs::useShinyjs()
                  , tabItem(
                      tabName = "seurat"
                      , shiny_one_panel
                  ) # End of tabItem

              )# end of dashboard body
)# end of dashboard page
