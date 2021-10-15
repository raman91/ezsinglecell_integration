## max data size
options(shiny.maxRequestSize = 1024^10)
options(shiny.launch.browser = T)

shinyServer(function(input, output, session) {
  v <- reactiveValues(scData = NULL,
                      scData1 = NULL,
                      scData2 = NULL,
                      scData3 = NULL,
                      scData4 = NULL,
                      scData5 = NULL,
                      scData6 = NULL,
                      scData7 = NULL,
                      scData8 = NULL,
                      scData9 = NULL,
                      scData10 = NULL,
                      scData11 = NULL,
                      scData12 = NULL,
                      scDatat = NULL,
                      idents = NULL,
                      isPCAdone = NULL,
                      isUMAPdone = NULL,
                      isTSNEdone = NULL,
                      isPCAdone1 = NULL,
                      isUMAPdone1 = NULL,
                      isTSNEdone1 = NULL,
                      isPCAdone2 = NULL,
                      isUMAPdone2 = NULL,
                      isTSNEdone2 = NULL,
                      isPCAdone3 = NULL,
                      isUMAPdone3 = NULL,
                      isTSNEdone3 = NULL,
                      isUMAPdone4 = NULL,
                      isTSNEdone4 = NULL,
                      isVisdone = NULL,
                      isClusterdone = NULL,
                      isDataIntegration = NULL,
                      pcGenes = NULL,
                      plotlySelection = NULL,
                      ips.markers = NULL)
  #celltypes <- NULL
  prePlot <- function(){
    while(names(dev.cur()) != "null device"){
      dev.off()
    }
  }
  observe({
    #s <- event_data("plotly_selected")
    #cells <- s[["key"]]
    v$plotlySelection <- event_data("plotly_selected")[["key"]]
  })
  ##-------------------Side Panel-------------------

  normMethod <- NULL

  output$name.field <- renderUI({
    if(is.null(input$cellAnnoFiles)){
      numericInput(inputId = "field",
                   label = "Field",
                   value = 1,
                   min = 1)
    }else{
      annoFile <- input$cellAnnoFiles
      anno.data <- read.table(annoFile$datapath[1], header = T,
                              sep = "\t", stringsAsFactors = FALSE)
      groupings <- colnames(anno.data)
      selectInput("groupby",
                  label = "Group by:",
                  choices = groupings)
    }
  })

  observeEvent(input$loadButton, {
    tpmFiles <- input$tpmFiles
    annoFile <- input$cellAnnoFiles
    names.field <- input$field
    if (is.null(tpmFiles)){
      v$scData <- NULL
    }else{
      withProgress(message="Loading and Processing Data...", value=0, {
        print(tpmFiles$datapath)
        print(tpmFiles$name)
        print(file.exists(paste(tpmFiles$datapath[1], "/", tpmFiles$name[1], sep="")))
        exp.data <- read.table(tpmFiles$datapath,
                               sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)
        #additional.ident1 <- NULL
        if(!is.null(annoFile)){
          anno.data <- read.table(annoFile$datapath[1], header = T,
                                  sep = "\t", stringsAsFactors = FALSE, row.names=1)
          #to.append <- apply(anno.data1, 1, paste, collapse = "_")
          #colnames(exp.data1) <- to.append
          #names.field1 <- match(input$groupby, colnames(anno.data1))
          #additional.ident1 <- data.frame(data.frame(anno.data1[,-1], row.names = to.append))
          #additional.ident1[] <- lapply(additional.ident1, factor)
        }
        incProgress(0.5, "Creating Seurat Object")
        sObj <- CreateSeuratObject(exp.data,
                                   meta.data = anno.data,
                                   project = input$projName,
                                   names.field = names.field,
                                   names.delim = input$delim,
                                   is.expr = input$expThres,
                                   normalization.method = normMethod,
                                   min.genes = input$min.genes,
                                   min.cells = input$min.cells)
        #mito.genes <- grep("^MT-", rownames(sObj@assays$RNA@data), ignore.case = TRUE, value = TRUE)
        sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
        #incProgress(0.5, "Adding metadata")
        #sObj <- AddMetaData(sObj, percent.mt, "percent.mt")
        #if(!is.null(additional.ident1)){
        #  sObj1 <- AddMetaData(sObj1, additional.ident1)
        #}
        v$scData <- sObj
      })
    }
  })
  dir.create("Seurat_results")
  #})

  observeEvent(input$reset, {
    session$reload()
    print("Reset done")
  })

  observeEvent(input$saveButton, {
    if(!is.null(input$tpmFiles)){
      withProgress(message="Saving Results...", value=0, {
        print(getwd())
        dir.create("Seurat_results")
        resultDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results")
        filename <- paste0(resultDir, .Platform$file.sep, v$scData@project.name, "_", Sys.Date())
        sObj <- v$scData
        save(sObj, file= paste0(resultDir, .Platform$file.sep, sObj@project.name, "_", Sys.Date(), ".Robj"))
      })
      ## open the results directory
      opendir(resultDir)
    }
  })

  output$ident.swap <- renderUI({
    if(is.null(v$scData)){
      return(NULL)
    }else{
      groupings1 <- names(v$scData@meta.data[,!names(v$scData@meta.data) %in% c("nFeature_RNA", "nCount_RNA", "percent.mt")])
      tagList(
        h4("Set current identity:"),
        fluidRow(
          column(6,
                 selectInput("active.ident", label = NULL,
                             choices = groupings1)
          ),
          column(6,
                 actionButton("swap.ident",label = NULL, icon = icon("arrow-right"))
          )
        )

      )
    }
  })

  observeEvent(input$swap.ident, {
    v$scData <- SetIdent(v$scData, value = as.character(v$scData@meta.data[,input$active.ident]))
  })

  output$logo <- renderImage({
    return(list(
      src = "inst/extdata/logo.png",
      contentType = "image/png",
      alt = "Singapore Immunology Network"
    ))
  }, deleteFile = FALSE)

  opendir <- function(dir = getwd()){
    if (.Platform['OS.type'] == "windows"){
      shell.exec(dir)
    } else {
      system(paste(Sys.getenv("R_BROWSER"), dir))
    }
  }

  observeEvent(input$OpenDir, {
    resultDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results")
    if(!dir.exists(resultDir)){
      dir.create("Seurat_results")
    }
    pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
    if(dir.exists(pdfDir)){
      opendir(pdfDir)
    }else{
      warning("No reports created yet!")
      dir.create(pdfDir)
    }
  })

  ##---------------Data Integration using Seurat-------------------

  observeEvent(input$doIntg, {
    withProgress(message = "Running Data Integration...", value = 0.3, {
      v$scData.list <- SplitObject(v$scData, split.by = "batchlb")
      print(v$scData)
      print(v$scData.list)
      v$scData.list <- pbmclapply(mc.cores = 20, X = v$scData.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      })

      features <- SelectIntegrationFeatures(object.list = v$scData.list)
      print(features)
      v$scData.anchors <- FindIntegrationAnchors(object.list = v$scData.list, anchor.features = features)
      v$scData.combined <- IntegrateData(anchorset = v$scData.anchors)
      print(v$scData.anchors)
      print(v$scData.combined)

      DefaultAssay(v$scData.combined) <- "integrated"
      v$scData.combined <- ScaleData(v$scData.combined, verbose = FALSE)
      print(v$scData.combined)
    })
  })

  observeEvent(input$runPCA, {
    withProgress(message = "Running PCA...", value = 0,{
      incProgress(0.5, message = "Running PCA...")
      v$scData.combined <- RunPCA(v$scData.combined, verbose = FALSE)
      print(v$scData.combined[["pca"]], dims = 1:5, nfeatures = 5)
      v$isPCAdone1 <- TRUE
      PCA_plot1a <- DimPlot(v$scData.combined, reduction = "pca", label = T)
      PCA_plot1b <- DimPlot(v$scData.combined, reduction = "pca", label = T, group.by = 'batch')
      PCA_plot1c <- DimPlot(v$scData.combined, reduction = "pca", label = T, group.by = 'celltype')
      print(PCA_plot1a)
      print(PCA_plot1b)
      print(PCA_plot1c)
    })
  })

  output$PCAplot_a <- renderPlotly({
    if(is.null(v$isPCAdone1)){
      plotly_empty()
    }else{
      withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "pca", label = T)
      })
    }
  })

  output$PCAplot_b <- renderPlotly({
    if(is.null(v$isPCAdone1)){
      plotly_empty()
    }else{
      withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "pca", label = T, group.by = 'batch')
      })
    }
  })

  output$PCAplot_c <- renderPlotly({
    if(is.null(v$isPCAdone1)){
      plotly_empty()
    }else{
      withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "pca", label = T, group.by = 'celltype')
      })
    }
  })

  observeEvent(input$runUMAP, {
    withProgress(message = "Running UMAP...", value = 0,{
      incProgress(0.5, message = "Running UMAP...")
      v$scData.combined <- FindNeighbors(v$scData.combined, reduction = "pca", dims = 1:30)
      v$scData.combined <- FindClusters(v$scData.combined, resolution = 0.5)
      v$scData.combined <- RunUMAP(v$scData.combined, reduction = "pca", dims = 1:30)
      v$isUMAPdone1 <- TRUE
      UMAP_plot1a <- DimPlot(v$scData.combined, reduction = "umap", label = T)
      UMAP_plot1b <- DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'batch')
      UMAP_plot1c <- DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'celltype')
      print(UMAP_plot1a)
      print(UMAP_plot1b)
      print(UMAP_plot1c)
    })
  })

  output$UMAPplot_a <- renderPlotly({
    if(is.null(v$isUMAPdone1)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "umap", label = T)
      })
    }
  })

  output$UMAPplot_b <- renderPlotly({
    if(is.null(v$isUMAPdone1)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'batch')
      })
    }
  })

  output$UMAPplot_c <- renderPlotly({
    if(is.null(v$isUMAPdone1)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'celltype')
      })
    }
  })

  observeEvent(input$runTSNE, {
    withProgress(message = "Running TSNE...", value = 0,{
      incProgress(0.5, message = "Running TSNE...")
      v$scData.combined <- RunTSNE(v$scData.combined, reduction = "pca", dims = 1:30)
      v$isTSNEdone1 <- TRUE
      TSNE_plot1a <- DimPlot(v$scData.combined, reduction = "tsne", label = T)
      TSNE_plot1b <- DimPlot(v$scData.combined, reduction = "tsne", label = T, group.by = 'batch')
      TSNE_plot1c <- DimPlot(v$scData.combined, reduction = "tsne", label = T, group.by = 'celltype')
      print(TSNE_plot1a)
      print(TSNE_plot1b)
      print(TSNE_plot1c)
    })
  })

  output$TSNEplot_a <- renderPlotly({
    if(is.null(v$isTSNEdone1)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "tsne", label = T)
      })
    }
  })

  output$TSNEplot_b <- renderPlotly({
    if(is.null(v$isTSNEdone1)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "tsne", label = T, group.by = 'batch')
      })
    }
  })

  output$TSNEplot_c <- renderPlotly({
    if(is.null(v$isTSNEdone1)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "tsne", label = T, group.by = 'celltype')
      })
    }
  })

  output$ident.swap1 <- renderUI({
    if(is.null(v$scData.combined)){
      return(NULL)
    }else{
      groupings1 <- names(v$scData.combined@meta.data[,!names(v$scData.combined@meta.data) %in% c("nGene", "nUMI", "percent.mito")])
      tagList(
        h4("Set current identity:"),
        fluidRow(
          column(3,
                 selectInput("active.ident1", label = NULL,
                             choices = groupings1)
          ),
          column(3,
                 actionButton("swap.ident1",label = NULL, icon = icon("arrow-right"))
          )
        )

      )
    }
  })

  observeEvent(input$PDFl, {
    if(!is.null(v$scData.combined) ){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"Data_integration_Seurat_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"Dimension_reduction_plots_Seurat_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        PCA_plot1a <- DimPlot(v$scData.combined, reduction = "pca", label = T)
        PCA_plot1b <- DimPlot(v$scData.combined, reduction = "pca", label = T, group.by = 'type')
        PCA_plot1c <- DimPlot(v$scData.combined, reduction = "pca", label = T, group.by = 'celltype')
        UMAP_plot1a <- DimPlot(v$scData.combined, reduction = "umap", label = T)
        UMAP_plot1b <- DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'type')
        UMAP_plot1c <- DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'celltype')
        TSNE_plot1a <- DimPlot(v$scData.combined, reduction = "tsne", label = T)
        TSNE_plot1b <- DimPlot(v$scData.combined, reduction = "tsne", label = T, group.by = 'type')
        TSNE_plot1c <- DimPlot(v$scData.combined, reduction = "tsne", label = T, group.by = 'celltype')
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        print(PCA_plot1a)
        print(PCA_plot1b)
        print(PCA_plot1c)
        print(UMAP_plot1a)
        print(UMAP_plot1b)
        print(UMAP_plot1c)
        print(TSNE_plot1a)
        print(TSNE_plot1b)
        print(TSNE_plot1c)
        dev.off()
      })
      withProgress(message="Downloading Dimension_reduction coordinates...", value=0.6, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2a <- paste0(pdfDir, .Platform$file.sep,"seurat_pca_", Sys.Date(), ".txt")
        filename2b <- paste0(pdfDir, .Platform$file.sep,"seurat_umap_", Sys.Date(), ".txt")
        filename2c <- paste0(pdfDir, .Platform$file.sep,"seurat_tsne_", Sys.Date(), ".txt")
        i = 0
        while(file.exists(filename2a)){
          filename2a <- paste0(pdfDir, .Platform$file.sep,"seurat_pca_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        while(file.exists(filename2b)){
          filename2b <- paste0(pdfDir, .Platform$file.sep,"seurat_umap_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        while(file.exists(filename2c)){
          filename2c <- paste0(pdfDir, .Platform$file.sep,"seurat_tsne_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        write.csv(v$scData.combined@reductions$pca@cell.embeddings, file = filename2a)
        write.csv(v$scData.combined@reductions$umap@cell.embeddings, file = filename2b)
        write.csv(v$scData.combined@reductions$tsne@cell.embeddings, file = filename2c)
      })
    }
  })

  ##---------------Data Integration using Harmony-------------------

  observeEvent(input$doIntg1, {
    withProgress(message = "Running Data Integration...", value = 0.3, {
      #v$scData <- merge(v$scData4, y = v$scData5, add.cell.ids = c("1", "2"), project = input$projName)
      #v$scData[['type']] <- ifelse(startsWith(v$scData6@assays$RNA@data@Dimnames[[2]], '1'), '1', '2')
      v$scData.list <- SplitObject(v$scData, split.by = "batch")
      print(v$scData)
      print(v$scData.list)
      v$scData.list <- pbmclapply(mc.cores = 20, X = v$scData.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      })

      features <- SelectIntegrationFeatures(object.list = v$scData.list)
      print(features)
      v$scData.anchors <- FindIntegrationAnchors(object.list = v$scData.list, anchor.features = features)
      v$scData.combined <- IntegrateData(anchorset = v$scData.anchors)
      print(v$scData.anchors)
      print(v$scData.combined)

      DefaultAssay(v$scData.combined) <- "integrated"
      v$scData.combined <- ScaleData(v$scData.combined, verbose = FALSE)
      print(v$scData.combined)
    })
  })

  observeEvent(input$runPCA1, {
    withProgress(message = "Running PCA...", value = 0,{
      incProgress(0.5, message = "Running PCA...")
      v$scData.combined <- RunPCA(v$scData.combined, verbose = FALSE)
      print(v$scData.combined[["pca"]], dims = 1:5, nfeatures = 5)
      v$isPCAdone2 <- TRUE
      v$scData.combined <- RunHarmony(v$scData.combined, "batch")
      PCA_plot2a <- DimPlot(v$scData.combined, reduction = "harmony", label = T)
      PCA_plot2b <- DimPlot(v$scData.combined, reduction = "harmony", group.by = 'batch')
      PCA_plot2c <- DimPlot(v$scData.combined, reduction = "harmony", group.by = 'celltype')
      print(PCA_plot2a)
      print(PCA_plot2b)
      print(PCA_plot2c)
    })
  })

  output$PCAplot1_a <- renderPlotly({
    if(is.null(v$isPCAdone2)){
      plotly_empty()
    }else{
      withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "harmony", label = T)
      })
    }
  })

  output$PCAplot1_b <- renderPlotly({
    if(is.null(v$isPCAdone2)){
      plotly_empty()
    }else{
      withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "harmony", label = T, group.by = 'batch')
      })
    }
  })

  output$PCAplot1_c <- renderPlotly({
    if(is.null(v$isPCAdone2)){
      plotly_empty()
    }else{
      withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "harmony", label = T, group.by = 'celltype')
      })
    }
  })

  observeEvent(input$runUMAP1, {
    withProgress(message = "Running UMAP...", value = 0,{
      incProgress(0.5, message = "Running UMAP...")
      v$scData.combined <- FindNeighbors(v$scData.combined, reduction = "harmony", dims = 1:30)
      v$scData.combined <- FindClusters(v$scData.combined, resolution = 0.5)
      v$scData.combined <- RunUMAP(v$scData.combined, reduction = "harmony", dims = 1:30)
      v$isUMAPdone2 <- TRUE
      UMAP_plot2a <- DimPlot(v$scData.combined, reduction = "umap", label = T)
      UMAP_plot2b <- DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'batch')
      UMAP_plot2c <- DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'celltype')
      print(UMAP_plot2a)
      print(UMAP_plot2b)
      print(UMAP_plot2c)
    })
  })

  output$UMAPplot1_a <- renderPlotly({
    if(is.null(v$isUMAPdone2)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "umap", label = T)
      })
    }
  })

  output$UMAPplot1_b <- renderPlotly({
    if(is.null(v$isUMAPdone2)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'batch')
      })
    }
  })

  output$UMAPplot1_c <- renderPlotly({
    if(is.null(v$isUMAPdone2)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'celltype')
      })
    }
  })

  observeEvent(input$runTSNE1, {
    withProgress(message = "Running TSNE...", value = 0,{
      incProgress(0.5, message = "Running TSNE...")
      v$scData.combined <- RunTSNE(v$scData.combined, reduction = "harmony", dims = 1:30)
      v$isTSNEdone2 <- TRUE
      TSNE_plot2a <- DimPlot(v$scData.combined, reduction = "tsne", label = T)
      TSNE_plot2b <- DimPlot(v$scData.combined, reduction = "tsne", label = T, group.by = 'batch')
      TSNE_plot2c <- DimPlot(v$scData.combined, reduction = "tsne", label = T, group.by = 'celltype')
      print(TSNE_plot2a)
      print(TSNE_plot2b)
      print(TSNE_plot2c)
    })
  })

  output$TSNEplot1_a <- renderPlotly({
    if(is.null(v$isTSNEdone2)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "tsne", label = T)
      })
    }
  })

  output$TSNEplot1_b <- renderPlotly({
    if(is.null(v$isTSNEdone2)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "tsne", label = T, group.by = 'batch')
      })
    }
  })

  output$TSNEplot1_c <- renderPlotly({
    if(is.null(v$isTSNEdone2)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scData.combined, reduction = "tsne", label = T, group.by = 'celltype')
      })
    }
  })

  output$ident.swap2 <- renderUI({
    if(is.null(v$scData6.combined)){
      return(NULL)
    }else{
      groupings2 <- names(v$scData.combined@meta.data[,!names(v$scData.combined@meta.data) %in% c("nGene", "nUMI", "percent.mito")])
      tagList(
        h4("Set current identity:"),
        fluidRow(
          column(3,
                 selectInput("active.ident2", label = NULL,
                             choices = groupings2)
          ),
          column(3,
                 actionButton("swap.ident2",label = NULL, icon = icon("arrow-right"))
          )
        )

      )
    }
  })

  observeEvent(input$swap.ident2, {
    v$scData.combined <- SetIdent(v$scData.combined, ident.use = as.character(v$scData.combined@meta.data[,input$active.ident2]))
  })

  observeEvent(input$PDFm, {
    if(!is.null(v$scData.combined) ){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"Data_integration_Harmony_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"Dimension_reduction_plots_Harmony_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        PCA_plot2a <- DimPlot(v$scData.combined, reduction = "pca", label = T)
        PCA_plot2b <- DimPlot(v$scData.combined, reduction = "pca", label = T, group.by = 'batchlb')
        PCA_plot2c <- DimPlot(v$scData.combined, reduction = "pca", label = T, group.by = 'celltype')
        UMAP_plot2a <- DimPlot(v$scData.combined, reduction = "umap", label = T)
        UMAP_plot2b <- DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'batchlb')
        UMAP_plot2c <- DimPlot(v$scData.combined, reduction = "umap", label = T, group.by = 'celltype')
        TSNE_plot2a <- DimPlot(v$scData6.combined, reduction = "tsne", label = T)
        TSNE_plot2b <- DimPlot(v$scData6.combined, reduction = "tsne", label = T, group.by = 'batchlb')
        TSNE_plot2c <- DimPlot(v$scData6.combined, reduction = "tsne", label = T, group.by = 'celltype')
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        print(PCA_plot2a)
        print(PCA_plot2b)
        print(PCA_plot2c)
        print(UMAP_plot2a)
        print(UMAP_plot2b)
        print(UMAP_plot2c)
        print(TSNE_plot2a)
        print(TSNE_plot2b)
        print(TSNE_plot2c)
        dev.off()
      })
      withProgress(message="Downloading Dimension_reduction coordinates...", value=0.6, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2a <- paste0(pdfDir, .Platform$file.sep,"harmony_pca_", Sys.Date(), ".txt")
        filename2b <- paste0(pdfDir, .Platform$file.sep,"harmony_umap_", Sys.Date(), ".txt")
        filename2c <- paste0(pdfDir, .Platform$file.sep,"harmony_tsne_", Sys.Date(), ".txt")
        i = 0
        while(file.exists(filename2a)){
          filename2a <- paste0(pdfDir, .Platform$file.sep,"harmony_pca_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        while(file.exists(filename2b)){
          filename2b <- paste0(pdfDir, .Platform$file.sep,"harmony_umap_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        while(file.exists(filename2c)){
          filename2c <- paste0(pdfDir, .Platform$file.sep,"harmony_tsne_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        write.csv(v$scData.combined@reductions$pca@cell.embeddings, file = filename2a)
        write.csv(v$scData.combined@reductions$umap@cell.embeddings, file = filename2b)
        write.csv(v$scData.combined@reductions$tsne@cell.embeddings, file = filename2c)
      })
    }
  })

  ##---------------Data Integration using RLiger-------------------

  observeEvent(input$doIntg2, {
    withProgress(message = "Running Data Integration...", value = 0.3, {
      #v$scData16 <- rbind(v$scData30, v$scData31)
      #v$scData9 <- merge(v$scData7, y = v$scData8, add.cell.ids = c("1", "2"), project = input$projName)
      #v$scData9[['type']] <- ifelse(startsWith(v$scData9@assays$RNA@data@Dimnames[[2]], '1'), '1', '2')
      v$scData.list <- SplitObject(v$scData, split.by = "batchlb")
      print(v$scData)
      print(v$scData.list)
      #v$scData9.list <- pbmclapply(mc.cores = 20, X = v$scData9.list, FUN = function(x) {
      #  x <-  SCTransform(x, method = "glmGamPoi", vars.to.regress = "percent.mt")
      #})

      v$scData.list <- lapply(X = v$scData.list, FUN = SCTransform)

      features <- SelectIntegrationFeatures(object.list = v$scData.list, nfeatures = 3000)
      print(features)
      v$scData.list <- PrepSCTIntegration(object.list = v$scData.list, anchor.features = features)
      v$scData.anchors <- FindIntegrationAnchors(object.list = v$scData.list, normalization.method = "SCT", anchor.features = features)
      v$scData.combined <- IntegrateData(anchorset = v$scData.anchors, normalization.method = "SCT")
      print(v$scData.anchors)
      print(v$scData.combined)

      DefaultAssay(v$scData.combined) <- "integrated"
      #v$scData9.combined <- ScaleData(v$scData9.combined, verbose = FALSE)
      print(v$scData.combined)
      v$scData.combined@meta.data -> v$scData16
      #v$scData1 <- list(b1 = GetAssayData(subset(v$scData.combined, subset = batch == '1'), slot="counts", assay='SCT'), b2 = GetAssayData(subset(v$scData.combined, subset = batch == '2'), slot="counts", assay='SCT'))
      v$scData1 <- list(b1 = GetAssayData(subset(v$scData.combined, subset = batch == '1'), slot="counts", assay='SCT'), b2 = GetAssayData(subset(v$scData.combined, subset = batch == '2'), slot="counts", assay='SCT'), b3 = GetAssayData(subset(v$scData.combined, subset = batch == '3'), slot="counts", assay='SCT'), b4 = GetAssayData(subset(v$scData.combined, subset = batch == '4'), slot="counts", assay='SCT'), b5 = GetAssayData(subset(v$scData.combined, subset = batch == '5'), slot="counts", assay='SCT'))
      #v$scData13 <- list(Tcell = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Tcell'), slot="counts", assay='SCT'), NK = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'NK'), slot="counts", assay='SCT'), Bcell = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Bcell'), slot="counts", assay='SCT'), Monocyte = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Monocyte'), slot="counts", assay='SCT'), Macrophage = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Macrophage'), slot="counts", assay='SCT'), Dendritic = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Dendritic'), slot="counts", assay='SCT'), Endothelial = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Endothelial'), slot="counts", assay='SCT'), SmoothMuscle = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'SmoothMuscle'), slot="counts", assay='SCT'), Neutrophil = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Neutrophil'), slot="counts", assay='SCT'), Stromal = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Stromal'), slot="counts", assay='SCT'), Epithelial = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Epithelial'), slot="counts", assay='SCT'))
    })
  })

  #observeEvent(input$runPCA2, {
  # withProgress(message = "Running PCA...", value = 0,{
  #    incProgress(0.5, message = "Running PCA...")
  #    v$scData9.combined <- RunOptimizeALS(v$scData9.combined, k = 30, lambda = 5, split.by = "batch")
  #    v$scData9.combined <- RunQuantileNorm(v$scData9.combined, split.by = "batch")
  #    v$scData9.combined <- RunPCA(v$scData9.combined, verbose = FALSE)
  #    print(v$scData9.combined[["pca"]], dims = 1:5, nfeatures = 5)
  #    v$isPCAdone3 <- TRUE
  #    PCA_plot3a <- DimPlot(v$scData9.combined, reduction = "pca", label = T)
  #    PCA_plot3b <- DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'batch')
  #    PCA_plot3c <- DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'orig.ident')
  #    print(PCA_plot3a)
  #    print(PCA_plot3b)
  #    print(PCA_plot3c)
  #  })
  # })

  #output$PCAplot2_a <- renderPlotly({
  #  if(is.null(v$isPCAdone3)){
  #    plotly_empty()
  #  }else{
  #    withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
  #      DimPlot(v$scData9.combined, reduction = "pca", label = T)
  #    })
  #  }
  # })

  #output$PCAplot2_b <- renderPlotly({
  #  if(is.null(v$isPCAdone3)){
  #    plotly_empty()
  #  }else{
  #    withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
  #      DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'batch')
  #    })
  #  }
  #})

  #output$PCAplot2_c <- renderPlotly({
  #  if(is.null(v$isPCAdone3)){
  #    plotly_empty()
  #  }else{
  #    withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
  #      DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'orig.ident')
  #    })
  #  }
  # })


  observeEvent(input$runUMAP2, {
    withProgress(message = "Running UMAP...", value = 0,{
      incProgress(0.5, message = "Running UMAP...")
      v$scData11 <- rliger::createLiger(v$scData1)
      v$scData11 <- rliger::normalize(v$scData11)
      v$scData11 <- rliger::selectGenes(v$scData11, var.thresh = 0.2, do.plot = F)
      v$scData11 <- rliger::scaleNotCenter(v$scData11)
      v$scData11 <- rliger::online_iNMF(v$scData11, k = 20, miniBatch_size = 5000, max.epochs = 5)
      v$scData11 <- rliger::quantile_norm(v$scData11)
      v$scData11 <- rliger::runUMAP(v$scData11)
      v$scData17 <- data.frame(v$scData11@tsne.coords)
      v$scData18 <- rownames(v$scData17)
      #gsub('1_','', v$scData18) -> v$scData19
      #gsub('2_','', v$scData19) -> v$scData20
      v$scData16 <- v$scData16[v$scData18,]
      v$scData11@clusters <- as.factor(v$scData16[,"celltype"])
      names(v$scData11@clusters) <- rownames(v$scData16)
      #plotByDatasetAndCluster(v$scData21, pt.size = 1, text.size = 0)

      #v$scData14 <- rliger::createLiger(v$scData13)
      #v$scData14 <- rliger::normalize(v$scData14)
      #v$scData14 <- rliger::selectGenes(v$scData14, var.thresh = 0.2, do.plot = F)
      #v$scData14 <- rliger::scaleNotCenter(v$scData14)
      #v$scData14 <- rliger::online_iNMF(v$scData14, k = 20, miniBatch_size = 5000, max.epochs = 5)
      #v$scData14 <- rliger::quantile_norm(v$scData14)
      #v$scData14 <- rliger::runUMAP(v$scData14)
      #v$scData11 <- ligerToSeurat(v$scData11, nms = names(v$scData11@H), renormalize = TRUE, use.liger.genes = TRUE, by.dataset = FALSE)
      #v$scData9.combined <- FindNeighbors(v$scData9.combined, reduction = "iNMF", dims = 1:30)
      #v$scData9.combined <- FindClusters(v$scData9.combined, resolution = 0.5)
      #v$scData9.combined <- RunUMAP(v$scData9.combined, dims = 1:30, reduction = "iNMF")
      v$isUMAPdone3 <- TRUE
      UMAP_plot3a <- rliger::plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"))[[1]]
      UMAP_plot3b <- rliger::plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"))[[2]]
      #UMAP_plot3b <- DimPlot(v$scData9.combined, reduction = "umap", label = T, group.by = 'batch')
      #UMAP_plot3c <- DimPlot(v$scData9.combined, reduction = "umap", label = T, group.by = 'orig.ident')
      print(UMAP_plot3a)
      print(UMAP_plot3b)
      #print(UMAP_plot3c)
    })
  })

  output$UMAPplot2_a <- renderPlotly({
    if(is.null(v$isUMAPdone3)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"), return.plots = T)[[1]]
      })
    }
  })

  output$UMAPplot2_b <- renderPlotly({
    if(is.null(v$isUMAPdone3)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"), return.plots = T)[[2]]
      })
    }
  })

  #output$UMAPplot2_c <- renderPlotly({
  #  if(is.null(v$isUMAPdone3)){
  #    plotly_empty()
  #  }else{
  #    withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
  #      DimPlot(v$scData9.combined, reduction = "umap", label = T, group.by = 'orig.ident')
  #    })
  #  }
  # })

  observeEvent(input$runTSNE2, {
    withProgress(message = "Running TSNE...", value = 0,{
      incProgress(0.5, message = "Running TSNE...")
      #v$scData9.combined <- RunTSNE(v$scData9.combined, dims = 1:ncol(v$scData9.combined[["iNMF"]]), reduction = "iNMF", tsne.method = "Rtsne", check_duplicates = FALSE)
      v$scData12 <- rliger::runTSNE(v$scData11)
      v$scData22 <- data.frame(v$scData12@tsne.coords)
      v$scData23 <- rownames(v$scData22)
      #gsub('1_','', v$scData23) -> v$scData24
      #gsub('2_','', v$scData24) -> v$scData25
      v$scData16 <- v$scData16[v$scData23,]
      v$scData12@clusters <- as.factor(v$scData16[,"celltype"])
      names(v$scData12@clusters) <- rownames(v$scData16)
      #plotByDatasetAndCluster(v$scData25, pt.size = 1, text.size = 0)
      #v$scData15 <- rliger::runTSNE(v$scData14)
      v$isTSNEdone3 <- TRUE
      TSNE_plot3a <- rliger::plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"))[[1]]
      TSNE_plot3b <- rliger::plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"))[[2]]
      #TSNE_plot3c <- DimPlot(v$scData9.combined, reduction = "tsne", label = T, group.by = 'orig.ident')
      print(TSNE_plot3a)
      print(TSNE_plot3b)
      #print(TSNE_plot3c)
    })
  })

  output$TSNEplot2_a <- renderPlotly({
    if(is.null(v$isTSNEdone3)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"), return.plots = T)[[1]]
      })
    }
  })

  output$TSNEplot2_b <- renderPlotly({
    if(is.null(v$isTSNEdone3)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"), return.plots = T)[[2]]
      })
    }
  })

  #output$TSNEplot2_c <- renderPlotly({
  #  if(is.null(v$isTSNEdone3)){
  #    plotly_empty()
  #  }else{
  #    withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
  #      DimPlot(v$scData9.combined, reduction = "tsne", label = T, group.by = 'orig.ident')
  #    })
  #  }
  # })

  output$ident.swap3 <- renderUI({
    if(is.null(v$scData.combined)){
      return(NULL)
    }else{
      groupings3 <- names(v$scData.combined@meta.data[,!names(v$scData.combined@meta.data) %in% c("nGene", "nUMI", "percent.mito")])
      tagList(
        h4("Set current identity:"),
        fluidRow(
          column(6,
                 selectInput("active.ident3", label = NULL,
                             choices = groupings3)
          ),
          column(6,
                 actionButton("swap.ident3",label = NULL, icon = icon("arrow-right"))
          )
        )

      )
    }
  })

  observeEvent(input$swap.ident3, {
    v$scData.combined <- SetIdent(v$scData.combined, ident.use = as.character(v$scData.combined@meta.data[,input$active.ident3]))
  })

  observeEvent(input$PDFn, {
    if(!is.null(v$scData.combined) ){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"Data_integration_RLiger_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"Dimension_reduction_plots_RLiger_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        #PCA_plot3a <- DimPlot(v$scData9.combined, reduction = "pca", label = T)
        #PCA_plot3b <- DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'batch')
        #PCA_plot3c <- DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'orig.ident')
        UMAP_plot3a <- plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"))[[1]]
        UMAP_plot3b <- plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"))[[2]]
        #UMAP_plot3c <- DimPlot(v$scData9.combined, reduction = "umap", label = T, group.by = 'orig.ident')
        TSNE_plot3a <- plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"))[[1]]
        TSNE_plot3b <- plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"))[[2]]
        #TSNE_plot3c <- DimPlot(v$scData9.combined, reduction = "tsne", label = T, group.by = 'orig.ident')
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        #print(PCA_plot3a)
        #print(PCA_plot3b)
        #print(PCA_plot3c)
        print(UMAP_plot3a)
        print(UMAP_plot3b)
        #print(UMAP_plot3c)
        print(TSNE_plot3a)
        print(TSNE_plot3b)
        #print(TSNE_plot3c)
        dev.off()
      })
      withProgress(message="Downloading Dimension_reduction coordinates...", value=0.6, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2a <- paste0(pdfDir, .Platform$file.sep,"rliger_pca_", Sys.Date(), ".txt")
        filename2b <- paste0(pdfDir, .Platform$file.sep,"rliger_umap_", Sys.Date(), ".txt")
        filename2c <- paste0(pdfDir, .Platform$file.sep,"rliger_tsne_", Sys.Date(), ".txt")
        i = 0
        while(file.exists(filename2a)){
          filename2a <- paste0(pdfDir, .Platform$file.sep,"rliger_pca_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        while(file.exists(filename2b)){
          filename2b <- paste0(pdfDir, .Platform$file.sep,"rliger_umap_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        while(file.exists(filename2c)){
          filename2c <- paste0(pdfDir, .Platform$file.sep,"rliger_tsne_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        write.csv(v$scData.combined@reductions$pca@cell.embeddings, file = filename2a)
        write.csv(v$scData.combined@reductions$umap@cell.embeddings, file = filename2b)
        write.csv(v$scData.combined@reductions$tsne@cell.embeddings, file = filename2c)
      })
    }
  })

  ##---------------scMultiomics using Seurat-------------------

  observeEvent(input$loadButton7, {
    tpmFiles7 <- input$tpmFiles7
    #annoFile3 <- input$cellAnnoFiles3
    #names.field3 <- input$field3
    if (is.null(tpmFiles7)){
      v$scDatat <- NULL
    }else{
      withProgress(message="Loading and Processing Data...", value=0, {
        print(tpmFiles7$datapath)
        print(tpmFiles7$name)
        print(file.exists(paste(tpmFiles7$datapath[1], "/", tpmFiles7$name[1], sep="")))
        exp.data7 <- Read10X_h5(tpmFiles7$datapath)
        exp.data7.rna <- exp.data7$`Gene Expression`
        exp.data7.ab <- exp.data7$`Antibody Capture`
        #additional.ident <- NULL
        incProgress(0.5, "Creating Seurat Object")
        sObj7 <- CreateSeuratObject(exp.data7.rna)
        #mito.genes <- grep("^MT-", rownames(sObj@assays$RNA@data), ignore.case = TRUE, value = TRUE)
        sObj7[["percent.mt"]] <- PercentageFeatureSet(sObj7, pattern = "^MT-")
        #incProgress(0.5, "Adding metadata")
        #sObj <- AddMetaData(sObj, percent.mt, "percent.mt")
        #if(!is.null(additional.ident)){
        #  sObj3 <- AddMetaData(sObj3, additional.ident)
        #}
        v$scDatat <- sObj7
        v$scDatab <- CreateAssayObject(counts = exp.data7.ab)
        v$scDatat[["ADT"]] <- v$scDatab
        print(Assays(v$scDatat))
      })
    }
  })

  observeEvent(input$runUMAP3, {
    withProgress(message = "Running UMAP...", value = 0,{
      incProgress(0.5, message = "Running UMAP...")
      DefaultAssay(v$scDatat) <- "RNA"
      v$scDatat <- NormalizeData(v$scDatat)
      v$scDatat <- FindVariableFeatures(v$scDatat)
      v$scDatat <- ScaleData(v$scDatat)
      v$scDatat <- RunPCA(v$scDatat, verbose = FALSE)
      v$scDatat <- FindNeighbors(v$scDatat, dims = 1:input$dim.used1)
      v$scDatat <- FindClusters(v$scDatat, resolution = input$clus.res1, verbose = FALSE)
      v$scDatat <- RunUMAP(v$scDatat, dims = 1:input$dim.used1)
      v$isUMAPdone4 <- TRUE
      UMAP_plot4a <- DimPlot(v$scDatat, label = TRUE, reduction = "umap")
      print(UMAP_plot4a)
    })
  })

  output$UMAPplot4_a <- renderPlotly({
    if(is.null(v$isUMAPdone4)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scDatat, reduction = "umap", label = T)
      })
    }
  })

  observeEvent(input$runTSNE3, {
    withProgress(message = "Running TSNE...", value = 0,{
      incProgress(0.5, message = "Running TSNE...")
      v$scDatat <- RunTSNE(v$scDatat, dims = 1:input$dim.used1)
      v$isTSNEdone4 <- TRUE
      TSNE_plot4a <- DimPlot(v$scDatat,  label = TRUE, reduction = "tsne")
      print(TSNE_plot4a)
    })
  })

  output$TSNEplot4_a <- renderPlotly({
    if(is.null(v$isTSNEdone4)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
        DimPlot(v$scDatat,  label = TRUE, reduction = "tsne")
      })
    }
  })

  observeEvent(input$Vis3, {
    withProgress(message = "Running TSNE...", value = 0,{
      incProgress(0.5, message = "Running TSNE...")
      DefaultAssay(v$scDatat) <- "ADT"
      #v$scDatat <- NormalizeData(v$scDatat, normalization.method = "CLR", margin = 2)
      #DefaultAssay(v$scDatat) <- "RNA"
      #v$scDatat <- NormalizeData(v$scDatat, normalization.method = "CLR", margin = 2, assay = "ADT")
      #DefaultAssay(v$scDatat) <- "ADT"
      #p1 <- FeaturePlot(v$scDatat, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
      #DefaultAssay(v$scDatat) <- "RNA"
      #p2 <- FeaturePlot(v$scDatat, "CD19") + ggtitle("CD19 RNA")
      #v$isVisdone <- TRUE
      #print(p1)
      #print(p2)
    })
  })

  output$vis.gene.select <- renderUI({
    if(is.null(v$scDatat)){
      return(NULL)
    }else{
      selectInput("vis.gene", label = "Gene to visualise",
                  choices = rownames(v$scDatat))
    }
  })

  output$vis.plot <- renderPlotly({
    if(is.null(v$scDatat)){
      return(NULL)
    }else{
      withProgress(message="Generating DEG Plot...", value=0, {
        FeaturePlot(v$scDatat, input$vis.gene, cols = c("lightgrey", "darkgreen"))
      })
    }
  })

  output$vis.gene.select1 <- renderUI({
    if(is.null(v$scDatat)){
      return(NULL)
    }else{
      selectInput("vis.gene1", label = "Gene to visualise",
                  choices = rownames(v$scDatat[["ADT"]]))
    }
  })

  output$vis1.plot <- renderPlotly({
    if(is.null(v$scDatat)){
      return(NULL)
    }else{
      withProgress(message="Generating DEG Plot...", value=0, {
        FeaturePlot(v$scDatat, input$vis.gene1, cols = c("lightgrey", "darkgreen"))
      })
    }
  })

  #output$Vis_a <- renderPlotly({
  #  if(is.null(v$isVisdone)){
  #    plotly_empty()
  #  }else{
  #    withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
  #      FeaturePlot(v$scDatat, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
  #    })
  #  }
  #})

  #output$Vis_b <- renderPlotly({
  #  if(is.null(v$isVisdone)){
  #    plotly_empty()
  #  }else{
  #    withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
  #      FeaturePlot(v$scDatat, "CD19") + ggtitle("CD19 RNA")
  #    })
  #  }
  #})

  observeEvent(input$PDFm, {
    if(!is.null(v$scData6.combined) ){
      withProgress(message="Downloading plot PDF files...", value=0, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2 <- paste0(pdfDir, .Platform$file.sep,"Data_integration_Harmony_", Sys.Date(), ".pdf")
        i = 0
        while(file.exists(filename2)){
          filename2 <- paste0(pdfDir, .Platform$file.sep,"Dimension_reduction_plots_Harmony_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
          i = i + 1;
        }
        PCA_plot2a <- DimPlot(v$scData6.combined, reduction = "pca", label = T)
        PCA_plot2b <- DimPlot(v$scData6.combined, reduction = "pca", label = T, group.by = 'type')
        PCA_plot2c <- DimPlot(v$scData6.combined, reduction = "pca", label = T, group.by = 'orig.ident')
        UMAP_plot2a <- DimPlot(v$scData6.combined, reduction = "umap", label = T)
        UMAP_plot2b <- DimPlot(v$scData6.combined, reduction = "umap", label = T, group.by = 'type')
        UMAP_plot2c <- DimPlot(v$scData6.combined, reduction = "umap", label = T, group.by = 'orig.ident')
        TSNE_plot2a <- DimPlot(v$scData6.combined, reduction = "tsne", label = T)
        TSNE_plot2b <- DimPlot(v$scData6.combined, reduction = "tsne", label = T, group.by = 'type')
        TSNE_plot2c <- DimPlot(v$scData6.combined, reduction = "tsne", label = T, group.by = 'orig.ident')
        prePlot()
        pdf(filename2,
            width=as.numeric(input$pdf_w),
            height=as.numeric(input$pdf_h))
        print(PCA_plot2a)
        print(PCA_plot2b)
        print(PCA_plot2c)
        print(UMAP_plot2a)
        print(UMAP_plot2b)
        print(UMAP_plot2c)
        print(TSNE_plot2a)
        print(TSNE_plot2b)
        print(TSNE_plot2c)
        dev.off()
      })
      withProgress(message="Downloading Dimension_reduction coordinates...", value=0.6, {
        print(getwd())
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(!dir.exists(pdfDir)){
          dir.create(pdfDir)
        }
        filename2a <- paste0(pdfDir, .Platform$file.sep,"harmony_pca_", Sys.Date(), ".txt")
        filename2b <- paste0(pdfDir, .Platform$file.sep,"harmony_umap_", Sys.Date(), ".txt")
        filename2c <- paste0(pdfDir, .Platform$file.sep,"harmony_tsne_", Sys.Date(), ".txt")
        i = 0
        while(file.exists(filename2a)){
          filename2a <- paste0(pdfDir, .Platform$file.sep,"harmony_pca_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        while(file.exists(filename2b)){
          filename2b <- paste0(pdfDir, .Platform$file.sep,"harmony_umap_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        while(file.exists(filename2c)){
          filename2c <- paste0(pdfDir, .Platform$file.sep,"harmony_tsne_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
          i = i + 1;
        }
        write.csv(v$scData6.combined@reductions$pca@cell.embeddings, file = filename2a)
        write.csv(v$scData6.combined@reductions$umap@cell.embeddings, file = filename2b)
        write.csv(v$scData6.combined@reductions$tsne@cell.embeddings, file = filename2c)
      })
    }
  })



  ##---------------Summary tab

  ##------Clean up when ending session----
  session$onSessionEnded(function(){
    prePlot()
  })
})


