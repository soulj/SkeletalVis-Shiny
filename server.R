library(shiny)
library(DT)
library(visNetwork)
library(squash)
library(feather)
library(readr)
library(dplyr)
library(gplots)
library(shinyWidgets)
library(RankProd)
library(enrichR)
library(shinyjs)

cat(file=stderr(),"\n")

addResourcePath("microarrayQC_html", "www/microarrayQC.html/")
addResourcePath("plots", "www/plots/")

cat(file=stderr(),  "\n")

function(input, output, session) {
  
  
  getJaccardHist <-function(ID,jaccardZScores) {
    zscores <- getSimilarityTable()[,"jaccard",F]
    colnames(zscores)<-"V1"
    ggplot(zscores,aes(V1)) +  geom_histogram(bins = 50,fill="lightgrey") + cowplot::theme_cowplot(font_size = 24) +scale_x_continuous(name = "Signed Jaccard coefficient") + scale_y_continuous(name = "Count")
  }
  
  
  getCosineHist <-function(ID,cosineZScores) {
    zscores <- getSimilarityTable()[,"cosine",F]
    colnames(zscores)<-"V1"
    ggplot(zscores,aes(V1)) +  geom_histogram(bins = 50,fill="lightgrey") + cowplot::theme_cowplot(font_size = 24) +scale_x_continuous(name = "Cosine coefficient") + scale_y_continuous(name = "Count")
  }
  
  
  getJaccardZscoreHist <-function(ID,jaccardZScores) {
    zscores <- getJaccardSim()[,3,F]
    zscores[is.nan(zscores[,1]),] <- NA
    colnames(zscores)<-"V1"
    ggplot(zscores,aes(V1)) +  geom_histogram(bins = 50,fill="lightgrey") + cowplot::theme_cowplot(font_size = 24) +scale_x_continuous(name = "z-score") + scale_y_continuous(name = "Count")
  }
  
  
  getChrDirZscoreHist <-function(ID,cosineZScores) {
    zscores <- getJaccardSim()[,4,F]
    zscores[is.nan(zscores[,1]),] <- NA
    colnames(zscores)<-"V1"
    ggplot(zscores,aes(V1)) +  geom_histogram(bins = 50,fill="lightgrey") + cowplot::theme_cowplot(font_size = 24) +scale_x_continuous(name = "z-score") + scale_y_continuous(name = "Count")
  }
  
  
  
  #generic null table for when no experiment/comparison is selected
  makeNullTable <- function(){
    pathways<-as.data.frame("Select an experiment and a comparison from the Choose Data tab first")
    colnames(pathways)=""
    pathways <-
      datatable(
        pathways,
        rownames = F,
        options = list(bInfo = FALSE,dom = 'rti')) %>% formatStyle(columns = 1,color = 'red')
    
    pathways
  }
  

  
  
  #gene search button
  getGeneFoldChangeTable <- eventReactive(input$geneSearch, {
    
    shinyjs::disable("geneSearch")
    
    
    geneName <- input$GeneName
    
    if (input$speciesSelectGene != "Human"){
      geneName <- convertGenes(geneName, tolower(input$speciesSelectGene))
    }
    
    print(head(geneName))

    
    if (is.na(geneName) | nrow(foldChangeTable[ foldChangeTable$ID == geneName,])==0){
      dt<-as.data.frame(paste("No entry for",input$GeneName,"for", input$speciesSelectGene,sep=" "))
      colnames(dt)=""
      dt <-
        datatable(
          dt,
          rownames = F,
          options = list(bInfo = FALSE,dom = 'rti')
        )
    } else{ 
      
      dt <- foldChangeTable[ foldChangeTable$ID == geneName,]
      dt<-stack(dt[,!colnames(dt)=="ID" ])
      colnames(dt) <- c("log2 FoldChange","Comparison")
      dt$Accession<-accessions[,1]
      dt$adj.Pval <- t(pvalTable[ pvalTable$GeneName==geneName,as.character(dt$Comparison)])[,1]
      print(head(dt$adj.Pval))
      dt$Comparison <- accessions[,"comparisonsText"]
      dt<-merge(dt,expTable,by.x="Accession",by.y="ID")
      dt<-dt[,c(2,4,5,3,1,6:13)]
      dt$Accession <- sapply( dt$Accession,createIDLink)
      dt$PMID<- sapply(dt$PMID,createPMIDLink)
      #dt<-na.omit(dt)
      
      dt[,1] <- signif(dt[,1], digits = 3)
      dt[,2] <- signif(dt[,2], digits = 3)
      
      dt<-dt[order(dt[,1],decreasing = T),]
      brks<-seq(from = min(dt[,1],na.rm = T),to=max(dt[,1],na.rm = T),length.out=100)
      clrs <- round(seq(255, 40, length.out = length(brks) + 1), 0) %>% {paste0("rgb(255,", ., ",", ., ")")}
      
      maxVal<-max(abs(dt[,1]),na.rm = T) + 0.001
      minVal<- -maxVal
      
      breaksList = seq(minVal, maxVal, by = 0.001)
      colfunc <- colorRampPalette(c("green", "white","red"))
      map <- makecmap(unique(dt[,1]),n = 3,breaks = breaksList,symm = T,colFn = colfunc)
      # clrs <- round(seq(255, 40, length.out = length(brks) + 1), 0) %>% {paste0("rgb(255,", ., ",", ., ")")}
      
      colours<-  cmap(unique(dt[,1]), map = map)
      style <- styleEqual(unique(dt[,1]), colours)

      dt<-datatable(
        dt,
        rownames = F,
        extensions = 'Scroller',
        options = list(
          scrollY = 300,
          scrollX = 500,
          scroller = TRUE,
          pageLength = 50,
          lengthChange = TRUE,
          bInfo = FALSE,
          order = list(0, 'desc'),buttons = c('copy', 'csv')
        ),
        selection = list(mode = "single"),
        filter = list(position = 'top', clear = FALSE, plain=TRUE)
       ,escape=F)
       dt <- formatStyle(dt,1, backgroundColor = style)

    }
    shinyjs::enable("geneSearch")
    dt
  })
  
  
  #gene fold change output
  output$geneFoldChangeTable <- DT::renderDataTable({ 
    dt <- getGeneFoldChangeTable()
    }
  )
  
  output$value <- renderText({ input$GeneName })
  
  
  cat(file=stderr(),"\n")
  
  #Load the comparisonTable
  getComparisonTable<-reactive({
    if (!is.null(input$expTable_row_last_clicked)){
      comparisons<-read_delim(paste0("data/", getAccession(),"/", getAccession(),"_comparisons.txt"),"\t", escape_double = FALSE, trim_ws = TRUE)
      comparisons
    }}
  )
  
  cat(file=stderr(),"\n")
  getSubnetworkTable <- reactive({
    subnetworkTable.string <-read.delim(paste0("data/", getAccession(), "/", getAccession(), "_summaryTable_",v$comparisonRow ,".txt"))
    subnetworkTable.biogrid <-read.delim(paste0("data/", getAccession(), "/", getAccession(), "_biogridsummaryTable_",v$comparisonRow ,".txt"))
    subnetworkTable <- rbind(subnetworkTable.string,subnetworkTable.biogrid)
    subnetworkTable$Source <- c(rep("StringDB",nrow(subnetworkTable.string)),rep("BioGrid",nrow(subnetworkTable.biogrid)))
    subnetworkTable <- subnetworkTable[,c(1,2,5,3:4)]
    subnetworkTable
    
  })
  
  #get the accession from the selected row in the experiment table
  getAccession<-renderText({
    # expTable<-getExpData()
    id = input$expTable_row_last_clicked
    if (!is.null(id)){
      accession<-expTable[id,"ID"]
      accession<-as.character(accession)
      accession
    }
    
  })
  cat(file=stderr(),"\n")
  
  #links for experiment table
  createIDLink <- function(ID) {
    accession = gsub("(.*)_.*", "\\1",ID)
    if ( grepl("GSE",accession)){
    as.character(tags$a(href=paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", accession), target="_blank",ID))
    } else {
      as.character(tags$a(href=paste0("https://www.ebi.ac.uk/arrayexpress/experiments/",accession),target="_blank", ID))
    }
  }
  
  createPMIDLink <- function(val) {
      if(!is.na(val)){
        as.character(tags$a(href=paste0("https://www.ncbi.nlm.nih.gov/pubmed/", val), target="_blank",val))
      } else{
        NA
      }
  }
  
  #make experiment table
  output$expTable = DT::renderDataTable({
    # expTable <- getExpData()
      expTable$ID <- sapply(expTable$ID,createIDLink)
      expTable$PMID <- sapply(expTable$PMID,createPMIDLink)
      datatable(
        expTable,
        rownames = F,
        extensions = 'Scroller',
        options = list(
          scrollY = 300,
          scrollX = 1000,
          scroller = TRUE,
          #scrollX=T,
          columnDefs=list(list(width="150px",targets=c(0,1)),list(width="250px",targets=c(2:4))),
          autoWidth = TRUE,
          pageLength = 10,
          lengthChange = TRUE,
          bInfo = FALSE,deferRender=T
        ),
        selection = list(mode = "single"),
        filter = list(position = 'top', clear = FALSE),
        escape = F)
  })
  
  comparisonTable <- reactive(getComparisonTable())
  
  
  # make comparisonTable
  output$comparisons = DT::renderDataTable({
    if ({
      id = input$expTable_row_last_clicked
      !is.null(id)
    }){
      comparisons<-comparisonTable()
      comparisons<-datatable(comparisons,rownames = F,
                             options = list(scrollY = "200px",
                                            lengthChange = FALSE,
                                            bInfo = FALSE,searching = FALSE, paging =FALSE                             ),
                             selection = list(mode = "single"),callback = JS("table.on('click.dt', 'td', function() {
            var row_=table.cell(this).index().row;
            var col=table.cell(this).index().column;
            var rnd= Math.random();
            var data = [row_, col, rnd];
           Shiny.onInputChange('rows',data );
    });")
      )
      comparisons
    }
  })
  
  #UIs
  output$AccessionText = renderUI({
    id = input$expTable_row_last_clicked
    if (is.null(id)){
      div(h4("Click a row in the table to select an Experiment"),style="color:red")
    }else {
      div(h5(getAccession()),downloadButton("downloadData", label = "Download Experiment"))
    }
  })
  
  output$ComparisonText = renderUI({
    id = v$comparisonRow
    if (is.null(id)){
      div(h4("Click a row in the table to select a Comparison"),style="color:red")
    }else {
      id<-v$comparisonRow
      comparison<-comparisonTable()[id,]
      comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep=" ")
      div(h5(comparison))
    }
  })
  
  output$reminderText = renderUI({
    id = v$comparisonRow
    if (!is.null(id) & input$tabset =="Choose Data"){
      div(h4("Step 3: Explore the analysed data through the tabs at the top"),style="color:green")
    }
  })
  
  output$subNetReminderText = renderUI({
    id = v$comparisonRow
    if (!is.null(id) & !is.null(v$subnetworkRow) & input$tabset =="subnet" & input$activesubnet =="Summary Table"){
      div(h4("View the selected subnetwork in the Subnetwork tab"),style="color:green")
    }
  })
  
  output$sharedResponseReminderText = renderUI({
    id = v$comparisonRow
    if (!is.null(id) & !is.null(v$simSummaryRow) & input$tabset =="response" & input$sharedResponse =="Summary"){
      if (input$responseType=="Pairwise"){
      div(h4("View the experimental overlap in the Gene Overlap tab"),style="color:green")
      }
    }
  })
  
  output$geneSearchUI = renderUI({
    id<-input$search
    if (!is.null(id) & id == "gene"){
      div(textInput("GeneName", "Enter a Gene Name"),selectInput("speciesSelectGene", "Select species", c("Human","Mouse","Rat","Zebrafish","Cow","Horse","Pig")),
          actionButton("geneSearchExample", "Example"),actionButton("geneSearch", "Search"))
    }
  })
  
  output$sigSearchUI = renderUI({
    id<-input$search
    if (!is.null(id) & id == "sig"){
      div(textAreaInput("UpRegulated", label = "Enter up-regulated gene names",resize="vertical"),
          textAreaInput("DownRegulated", "Enter down-regulated gene names",resize="vertical"),
          textAreaInput("Background", "Enter the background/all expressed genes (optional)",resize="vertical"),
          selectInput("speciesSelectSig", "Select species", c("Human","Mouse","Rat","Zebrafish","Cow","Horse","Pig")),
          actionButton("loadSigExample", "Example"),actionButton("sigSearch", "Search"))
    }
    
  })
  
  
  output$expHelp = renderUI({
    div(icon("fas fa-question-circle"))
      })
  output$fcHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$pathwayHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$GOHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$GOslimHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$GOMDSHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$TFHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$mimicDrugHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$reverseDrugHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$networkHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$subnetworkHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  output$chrDirHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  
  output$geneSearchHelp = renderUI({
    id<-input$search
    if (!is.null(id) & id == "gene"){
    div(icon("fas fa-question-circle"))
    }
  })
  
  output$sigSearchHelp = renderUI({
    id<-input$search
    if (!is.null(id) & id == "sig"){
      div(icon("fas fa-question-circle"))
    }
  })

  output$responseHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$responseHelpMulti = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$geneOverlapHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$chrDirOverlapHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  output$sigOverlapHelp = renderUI({
    div(icon("fas fa-question-circle"))
  })
  
  

  #sliders for fold change thresholds
  output$thresholdPicker <- renderUI({
    id = v$comparisonRow
    if (!is.null(id) & input$tabset %in% c("Pathways","TF Enrichment","go","drug")){
      fcOnly =expTable[input$expTable_row_last_clicked,"foldChangeOnly"]
      print(fcOnly)
      if (fcOnly=="FALSE"){
        pickerInput("threshold", h3("Differential Expression Threshold"),
                choices = c("FoldChange 1 padj 0.05"="_1_0.05","FoldChange 1.5 padj 0.05"="_1.5_0.05","FoldChange 2 padj 0.05"="_2_0.05","FoldChange 1.5 padj 1"="_1.5_1","FoldChange 2 padj 1"="_2_1")
                ,selected=v$thresholdPickerPval)
      } else{
        pickerInput("threshold", h3("Differential Expression Threshold"),
                    choices = c("FoldChange 1.5"="_1.5_1","FoldChange 2"="_2_1")
                    ,selected=v$thresholdPickerFC
        )
                    
      }
    }
  })
  
  
  
  #sliders for p-value threshold
  output$thresholdSlider <- renderUI({

    id = v$comparisonRow
    if (!is.null(id) & input$sharedResponseMulti %in% c("Gene Overlap","multiGeneOverlap") & !is.null(v$similaritySummaryMultiRow)){
      output = tagList()
      output[[1]] = sliderInput("thresholdSlider","Choose p-value threshold",0,1,value = 0.05)
      output[[2]] = selectInput("databaseSelect","Choose Database",choices = databases)
      output[[3]] = actionButton(inputId = "pathwaySubmit","Pathway Analysis")
      output
    }

  })

  output$thresholdSliderHelp <- renderUI({
    id = v$comparisonRow
    if (!is.null(id) & input$sharedResponseMulti %in% c("Gene Overlap","multiGeneOverlap") & !is.null(v$similaritySummaryMultiRow)){
    div(icon("fas fa-question-circle"))
    }
  })
  
  #sliders for vote counting method
  output$voteSlider <- renderUI({
    
    id = v$comparisonRow
    if (!is.null(id) & input$sharedResponseMulti %in% c("ChrDir Overlap") & !is.null(v$similaritySummaryMultiRow)){
      output = tagList()
      output[[1]] = sliderInput("voteSlider","Choose vote threshold",1,plotReady$chrDirCols,value = 2,step=1)
      output[[2]] = selectInput("databaseSelectChrDir","Choose Database",choices = databases)
      output[[3]] = actionButton(inputId = "pathwaySubmitChrDir","Pathway Analysis")
      output
      
      
    }
  })
  
  output$voteSliderHelp <- renderUI({
    id = v$comparisonRow
    if (!is.null(id) & input$sharedResponseMulti %in% c("ChrDir Overlap") & !is.null(v$similaritySummaryMultiRow)){
      div(icon("fas fa-question-circle"))
    }
  })
  
  #sliders for vote counting method
  output$voteSliderSig <- renderUI({
    
    id = v$comparisonRow
    if (!is.null(id) & input$sharedResponseMulti %in% c("Sig Overlap") & !is.null(v$similaritySummaryMultiRow)){
      output = tagList()
      output[[1]] = sliderInput("voteSliderSig","Choose vote threshold",1,plotReady$sigCols,value = 2,step=1)
      output[[2]] = selectInput("databaseSelectSig","Choose Database",choices = databases)
      output[[3]] = actionButton(inputId = "pathwaySubmitSig","Pathway Analysis")
	    output
      
    }
  })
  
  output$voteSliderSigHelp <- renderUI({
    id = v$comparisonRow
    if (!is.null(id) & input$sharedResponseMulti %in% c("Sig Overlap") & !is.null(v$similaritySummaryMultiRow)){
      div(icon("fas fa-question-circle"))
    }
  })
  
  
  
  #data download handler
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(getAccession(),'.zip', sep='')
    },
    content = function(con) {
      dataDir <- paste0("data/",getAccession())
      wwwDir <- paste0("www/microarrayQC.html/",getAccession())
      dataFiles<-list.files(path = dataDir, recursive=TRUE,full.names = T)
      wwwFiles <- list.files(path = wwwDir, recursive=TRUE,full.names = T)
      plotDir <- paste0("www/plots/",getAccession())
      plotFiles <- list.files(path = plotDir, pattern="GO.MDS", recursive=TRUE,full.names = T)
      
      downloadDir <- makeDownloadDirectory(dataFiles,wwwFiles,plotFiles)
      
      files <- c(dataFiles,wwwFiles,plotFiles)
      print(files)
      zip(con,files)
    },contentType = "application/zip"
  )

  makeDownloadDirectory <- function(dataFiles,wwwFiles,plotFiles){
    downloadDir <- basename(tempfile())
    dir.create(downloadDir)
    
    #for every comparison make a directory
    comparisons <- getComparisonTable()
    comparisons <- paste(comparisons$Numerator,comparisons$Denominator,sep="vs")
    sapply(comparisons,function(x) dir.create(paste(downloadDir,x,sep="/")))
    
    files <- c(dataFiles,wwwFiles,plotFiles)
    files<-files[ -grep("zip",files)]
    #collect the comparison files
    comparisonFiles <- lapply(seq_along(comparisons),function(x){
      temp <- files[grep(paste0("_",x,"|plot",x),files)]
    })
    
    #get the files into the correct directory
    mapply(function(comp,files,downloadDir){
      
      for(file in files){
        file.copy(file,paste(downloadDir,comp,basename(file),sep="/"))
      }
    }, comparisons,comparisonFiles,MoreArgs=list(downloadDir=downloadDir))
    
    #add the experiment level files
    tempFiles <- files[grep(paste("PCA","comparisons","sampletable","microarrayQC",sep="|"),files,ignore.case = T)]
    for(tempFile in tempFiles){
      file.copy(tempFile,paste(downloadDir,basename(tempFile),sep="/"))
    }
    
    return(downloadDir)
    
  }
  
  
  #create the QC panel
  output$QC = renderUI({
    id = input$expTable_row_last_clicked
    if (!is.null(id)){
      fileName <- paste0(
        "www/microarrayQC.html/",
        getAccession(),
        "/microarrayQC_html.html"
      )
      if(file.exists(fileName)){
        QC <-
          tags$iframe(
            # seamless = "seamless",
            src =paste0(
              "microarrayQC.html/",
              getAccession(),
              "/microarrayQC_html.html"
            ),
            height = 800,
            width = "100%",
            scrolling = "yes",
            frameBorder = 0
          )
        QC
      }
      else{
        h5("QC not available",style="color:red")
        
      }
    }
  })
  
  #create the PCA panel
  output$PCA = renderImage(
    list(
      src = paste0("data/", getAccession(), "/", getAccession(), "_PCA.png"),
      contentType = 'image/png',
      width = 1000,
      height = 800,
      alt = "PCA"
    ),
    deleteFile = F
  )
  
  #load fold change table
  getFoldChangesTable = reactive({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID)){
      fc <- read_delim(paste0("data/", getAccession(), "/", getAccession(), "_FC_",comparisonID,".txt"),"\t", escape_double = FALSE, trim_ws = TRUE)
      fc <- dplyr::mutate_if(fc,is.numeric, signif, 3)
      fc <- arrange_at(fc,2,desc)
      fc <- fc[fc[,1] != " ",]
      
      if (ncol(fc)==2){
        colnames(fc)[1:2] <- c("Gene Symbol","log2 Fold Change")
      } else {
        colnames(fc)[1:3] <- c("Gene Symbol","log2 Fold Change","adj. p-value")
      }
      fc
    }
  })
  
  
  #create the foldchange table
  output$foldchanges = DT::renderDataTable({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID)){
      comparison<-comparisonTable()[v$comparisonRow,]
      comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep="_")
      fc <- getFoldChangesTable()
      fc <-
        datatable(
          fc,
          rownames = F,
          extensions = c('Buttons','Scroller'), options = list(
            deferRender = TRUE,
            scrollY = 500,
            scroller = TRUE,
            pageLength = 50,
            lengthChange = FALSE,
            bInfo = FALSE,
            dom = 'frtipB',
            buttons = list('copy',list(extend='csv',filename=paste(getAccession(),comparison,"FC",sep="_")))
            ,deferRender=TRUE
          ),
          selection = list(mode = "single"),
          filter = "top"
        )
      fc
    } else{
      makeNullTable()
    }
    
  })
  
  
  
  #create the pathway table
  output$pathways = DT::renderDataTable({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID) & !is.null(input$threshold)){
      pathways<-paste0("data/", getAccession(), "/", getAccession(), "_pathways_",comparisonID,input$threshold,".txt")
      print(pathways)
      if(file.size(pathways)<5){
        pathways<-as.data.frame("No significant pathways")
        colnames(pathways)=""
        pathways <-
          datatable(
            pathways,
            rownames = F,
            options = list(bInfo = FALSE,dom = 'rti')
          )
      } else {
        comparison<-comparisonTable()[v$comparisonRow,]
        comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep="_")
        pathways <- read.delim(pathways)[,-5]
        pathways[, 3:ncol(pathways)] <-
          signif(pathways[, 3:ncol(pathways)], digits = 3)
        colnames(pathways)[4] <- "Percentage Cover"
        pathways <-
          datatable(
            pathways,
            rownames = F,
            options = list(
              pageLength = 50,
              scrollY = "500px",
              lengthChange = FALSE,
              bInfo = FALSE,
              dom = 'frtipB',
              buttons = list('copy',list(extend='csv',filename=paste(getAccession(),comparison,"pathways",sep="_"))),

              bPaginate = FALSE,
              deferRender=TRUE
            ),
            selection = list(mode = "single"),
            filter = list(position = 'top', clear = FALSE)
          )
      }
      pathways
    } else {
      makeNullTable()
    }
  })
  
  #load characteristic direction table
  getChrDirTable <- reactive({
    chrdirs<-paste0("data/", getAccession(), "/", getAccession(), "_chrDirTable.txt")
    if(file.size(chrdirs)<5){return(NA)}
    chrdirs <- read.delim(chrdirs)
    if ("gene_name" %in% colnames(chrdirs)){
      chrdirs <- chrdirs[ chrdirs$comparisonNumber == v$comparisonRow,c("gene_name","chrDir"),]
    } else {
      chrdirs <- chrdirs[ chrdirs$comparisonNumber == v$comparisonRow,c("ID","chrDir"),]
    }
    colnames(chrdirs)<-c("Gene Symbol","Characterstic Direction")
    chrdirs[, 2] <- chrdirs[, 2]*-1
    chrdirs[, 2] <- signif(chrdirs[, 2], digits = 3)
    chrdirs <- chrdirs[order(chrdirs[,2],decreasing = T),]
    chrdirs
  })
  
  
  
  #create the chrDir table
  output$chrDir = DT::renderDataTable({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID)){
      chrdirs <- getChrDirTable()
      if(is.null(dim(chrdirs))){
        chrdirs<-as.data.frame("Characteristic direction analysis only possible with experimental replicates")
        colnames(chrdirs)=""
        chrdirs <-
          datatable(
            chrdirs,
            rownames = F,
            options = list(bInfo = FALSE,dom = 'rti')
          )
      } else {
        
        comparison<-comparisonTable()[v$comparisonRow,]
        comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep="_")
        
        chrdirs <-
          datatable(
            chrdirs,
            rownames = F,
            options = list(
              pageLength = 50,
              scrollY = "500px",
              lengthChange = FALSE,
              bInfo = FALSE,
              dom = 'frtipB',
              buttons = list('copy',list(extend='csv',filename=paste(getAccession(),comparison,"chrDirTable",sep="_"))),
              
              bPaginate = FALSE,referRender=TRUE
            ),
            selection = list(mode = "single"),
            filter = "top"
          )
      }
      chrdirs
    } else {
      makeNullTable()
    }
  })
  
  #create the GO enrichment table
  output$goTable = DT::renderDataTable({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID) & !is.null(input$threshold)){
      terms <- paste0("data/", getAccession(), "/", getAccession(), "_goterms_",comparisonID,input$threshold,".txt")
      if(file.size(terms)<5){
        terms<-as.data.frame("No significant GO Terms")
        colnames(terms)=""
        terms <-
          datatable(
            terms,
            rownames = F,
            options = list(bInfo = FALSE,dom = 'rti')
          )
      } else {
        comparison<-comparisonTable()[v$comparisonRow,]
        comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep="_")
        terms <-read.delim(terms)
        terms[, 4:5] <-
        signif(terms[,4:5], digits = 3)
        terms[,5]<-terms[,5]*100
        colnames(terms)[5] <- "Percentage Cover"
         terms <-
            datatable(
            cbind(' ' = '&oplus;', terms),escape = FALSE,
            rownames = F,
            options = list(lengthMenu = c(10,50,100,nrow(terms)), order = list(4, 'asc'),
                         columnDefs = list(
                           list(visible = FALSE, targets = 3),
                           list(orderable = FALSE, className = 'details-control', targets = 0)
                         ),
                         pageLength = 50,
                         scrollY = "500px",
                         lengthChange = FALSE,
                         bInfo = FALSE,
                         dom = 'frtipB',
                         buttons = list('copy',list(extend='csv',filename=paste(getAccession(),comparison,"GOTerms",sep="_"))),
                         bPaginate = FALSE,referRender=TRUE
            ),
            callback = JS("
                table.column(1).nodes().to$().css({cursor: 'pointer'});
                var format = function(d) {
                return '<div style=\"background-color:#eee; padding: .5em;\">' +
                '<b>'+\"Genes: \" + '</b>'+ d[3] +  '</div>';
                };
                table.on('click', 'td.details-control', function() {
                var td = $(this), row = table.row(td.closest('tr'));
                if (row.child.isShown()) {
                row.child.hide();
                td.html('&oplus;');
                } else {
                row.child(format(row.data())).show();
                td.html('&CircleMinus;');
                }
                });"
           ),
           selection = list(mode = "single")
            #          filter = "top"
          )
      }
      terms
    } else {
      makeNullTable()
    }
  })
  
  #create the GO Reduced enrichment table
  output$goReducedTable = DT::renderDataTable({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID) & !is.null(input$threshold)){
      terms <- paste0("data/", getAccession(), "/", getAccession(), "_goterms_reduced_",comparisonID,input$threshold,".txt")
      if(file.size(terms)<5){
        terms<-as.data.frame("No significant GO Terms")
        colnames(terms)=""
        terms <-
          datatable(
            terms,
            rownames = F,
            options = list(bInfo = FALSE,dom = 'rti',buttons = c('copy', 'csv'))
          )
      } else {
        comparison<-comparisonTable()[v$comparisonRow,]
        comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep="_")
        terms <-read.delim(terms)
        terms[, 4:5] <-
          signif(terms[,4:5], digits = 3)
        terms[,5]<-terms[,5]*100
        terms <-
          datatable(
            cbind(' ' = '&oplus;', terms),escape = FALSE,
            rownames = F,
            options = list(lengthMenu = c(10,50,100,nrow(terms)), order = list(4, 'asc'),
                           columnDefs = list(
                             list(visible = FALSE, targets = 3),
                             list(orderable = FALSE, className = 'details-control', targets = 0)
                           ),
                           pageLength = 50,
                           scrollY = "500px",
                           lengthChange = FALSE,
                           bInfo = FALSE,
                           dom = 'frtipB',
                           buttons = list('copy',list(extend='csv',filename=paste(getAccession(),comparison,"ReducedGOTerms",sep="_"))),
                           bPaginate = FALSE,referRender=TRUE
            ),
            callback = JS("
                table.column(1).nodes().to$().css({cursor: 'pointer'});
                var format = function(d) {
                return '<div style=\"background-color:#eee; padding: .5em;\">' +
                '<b>'+\"Genes: \" + '</b>'+ d[3] +  '</div>';
                };
                table.on('click', 'td.details-control', function() {
                var td = $(this), row = table.row(td.closest('tr'));
                if (row.child.isShown()) {
                row.child.hide();
                td.html('&oplus;');
                } else {
                row.child(format(row.data())).show();
                td.html('&CircleMinus;');
                }
                });"
            ),
            selection = list(mode = "single")
            #          filter = "top"
          )
      }
      terms
    } else {
      makeNullTable()
    }
  })
  
  
  
  #create the GO MDS plot
  output$goMDS = renderUI({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID) & !is.null(input$threshold)){
      go.mds <- paste0("plots/", getAccession(), "/plot",comparisonID,input$threshold,"_html/","GO.MDS_html.html")
   		if(file.info(paste0("www/",go.mds))$size != 0){
        go.mds <-
        tags$iframe(
          seamless = "seamless",
          src = go.mds,
          height = 800,
          width = "100%",
          scrolling = "yes",
          frameBorder = 0
        )
      go.mds
      } else{
        h5("GO MDS not available",style="color:red")
        
      }
      
    }
  })
  
      
  #create the drug mimic table
  output$mimicTable = DT::renderDataTable({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID) & !is.null(input$threshold)){
      mimic <-paste0("data/", getAccession(), "/", getAccession(), "_drugsMimic_",comparisonID,input$threshold,".txt")
      if(file.size(mimic)<5){
        mimic<-as.data.frame("No significant drugs")
        colnames(mimic)=""
        mimic <-
          datatable(
            mimic,
            rownames = F,
            options = list(bInfo = FALSE,dom = 'rti')
          )
      } else {
        comparison<-comparisonTable()[v$comparisonRow,]
        comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep="_")
        mimic <-read.delim(mimic)
        mimic <-read.delim(paste0("data/", getAccession(), "/", getAccession(), "_drugsMimic_",comparisonID,input$threshold,".txt"))
        mimic <-
          datatable(
            cbind(' ' = '&oplus;', mimic),escape = FALSE,
            rownames = F,
            options = list(lengthMenu = c(10,50,100,nrow(mimic)), order = list(6, 'desc'),
                           columnDefs = list(
                             list(visible = FALSE, targets = c(8:10)),
                             list(orderable = FALSE, className = 'details-control', targets = 0)
                           ),
                           pageLength = 50,
                           scrollY = "500px",
                           lengthChange = FALSE,
                           bInfo = FALSE,
                           dom = 'frtipB',
                           buttons = list('copy',list(extend='csv',filename=paste(getAccession(),comparison,"MimicDrugs",sep="_"))),
                           bPaginate = FALSE,referRender=TRUE
            ),
            callback = JS("
                  table.column(1).nodes().to$().css({cursor: 'pointer'});
                  var format = function(d) {
                  return '<div style=\"background-color:#eee; padding: .5em;\">' +
                  '<b>'+\"Up Overlap: \" + '</b>'+ d[8] + '<BR/>' +'<b>'+ \"Down Overlap: \" +'</b>'+ d[9] + '<BR/>' + '<b>' + \"Targets: \" +'</b>'+ d[10] +  '</div>';
                  };
                  table.on('click', 'td.details-control', function() {
                  var td = $(this), row = table.row(td.closest('tr'));
                  if (row.child.isShown()) {
                  row.child.hide();
                  td.html('&oplus;');
                  } else {
                  row.child(format(row.data())).show();
                  td.html('&CircleMinus;');
                  }
                  });"
            ),
            selection = list(mode = "single")
            #          filter = "top"
          )
      }
        mimic
    }else {
      makeNullTable()
    }
  })
  
  #create the drug reverse table
  output$reverseTable = DT::renderDataTable({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID) & !is.null(input$threshold)){
      reverse <-paste0("data/", getAccession(), "/", getAccession(), "_drugsReverse_",comparisonID,input$threshold,".txt")
      if(file.size(reverse)<5){
        reverse<-as.data.frame("No significant drugs")
        colnames(reverse)=""
        reverse <-
          datatable(
            reverse,
            rownames = F,
            options = list(bInfo = FALSE,dom = 'rti')
          )
      } else {
        comparison<-comparisonTable()[v$comparisonRow,]
        comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep="_")
        reverse <-read.delim(reverse)
        reverse <-read.delim(paste0("data/", getAccession(), "/", getAccession(), "_drugsReverse_",comparisonID,input$threshold,".txt"))
        reverse <-
          datatable(
            cbind(' ' = '&oplus;', reverse),escape = FALSE,
            rownames = F,
            options = list(lengthMenu = c(10,50,100,nrow(reverse)), order = list(6, 'desc'),
                           columnDefs = list(
                             list(visible = FALSE, targets = c(8:10)),
                             list(orderable = FALSE, className = 'details-control', targets = 0)
                           ),
                           pageLength = 50,
                           scrollY = "500px",
                           lengthChange = FALSE,
                           bInfo = FALSE,
                           dom = 'frtipB',
                           buttons = list('copy',list(extend='csv',filename=paste(getAccession(),comparison,"ReverseDrugs",sep="_"))),
                           bPaginate = FALSE,referRender=TRUE
            ),
            callback = JS("
                  table.column(1).nodes().to$().css({cursor: 'pointer'});
                  var format = function(d) {
                  return '<div style=\"background-color:#eee; padding: .5em;\">' +
                  '<b>'+\"Up Overlap: \" + '</b>'+ d[8] + '<BR/>' +'<b>'+ \"Down Overlap: \" +'</b>'+ d[9] + '<BR/>' + '<b>' + \"Targets: \" +'</b>'+ d[10] +  '</div>';
                  };
                  table.on('click', 'td.details-control', function() {
                  var td = $(this), row = table.row(td.closest('tr'));
                  if (row.child.isShown()) {
                  row.child.hide();
                  td.html('&oplus;');
                  } else {
                  row.child(format(row.data())).show();
                  td.html('&CircleMinus;');
                  }
                  });"
            ),
            selection = list(mode = "single")
            #          filter = "top"
          )
      }
      reverse
    }else {
      makeNullTable()
    }
  })
  
  #create the TF enrichment table
  output$TFs = DT::renderDataTable({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    
    species <-expTable[expTable$ID==getAccession(),"Species"]
    
    if (!species %in% c("Human","Mouse")){
      TFs<-as.data.frame("TF enrichment only possible for Human and Mouse experiments")
      colnames(TFs)=""
      TFs <-
        datatable(
          TFs,
          rownames = F,
          options = list(bInfo = FALSE,dom = 'rti')
        )
    } else if (!is.null(id) & !is.null(comparisonID) & !is.null(input$threshold)){
      TFs <- paste0("data/", getAccession(), "/", getAccession(), "_TFs_",comparisonID,input$threshold,".txt")
      
      if(file.size(TFs)<5 | !file.exists(TFs)){
        TFs<-as.data.frame("No significant TFs")
        colnames(TFs)=""
        TFs <-
          datatable(
            TFs,
            rownames = F,
            options = list(bInfo = FALSE,dom = 'rti')
          )
      } else {
        comparison<-comparisonTable()[v$comparisonRow,]
        comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep="_")
        TFs <-
          read.delim(paste0("data/", getAccession(), "/", getAccession(), "_TFs_",comparisonID,input$threshold,".txt"))[, c(-1,-8)]
        TFs[,7] <- gsub(";"," ",TFs[,7])
        TFs<-TFs[ order(TFs[,2],decreasing=T),]
        TFs <-
          datatable(
            cbind(' ' = '&oplus;', TFs),escape = FALSE,
            rownames = F,
            options = list(lengthMenu = c(10,50,100,nrow(TFs)), order = list(2, 'desc'),
                           columnDefs = list(
                             list(visible = FALSE, targets = c(7)),
                             list(orderable = FALSE, className = 'details-control', targets = 0)
                           ),
                           pageLength = 50,
                           scrollY = "500px",
                           lengthChange = FALSE,
                           bInfo = FALSE,
                           dom = 'frtipB',
                           buttons = list('copy',list(extend='csv',filename=paste(getAccession(),comparison,input$threshold,"_TFs",sep="_"))),
                           bPaginate = FALSE,referRender=TRUE
            ),
            callback = JS("
                  table.column(1).nodes().to$().css({cursor: 'pointer'});
                  var format = function(d) {
                  return '<div style=\"background-color:#eee; padding: .5em;\">' +
                  d[7] +  '</div>';
                  };
                  table.on('click', 'td.details-control', function() {
                  var td = $(this), row = table.row(td.closest('tr'));
                  if (row.child.isShown()) {
                  row.child.hide();
                  td.html('&oplus;');
                  } else {
                  row.child(format(row.data())).show();
                  td.html('&CircleMinus;');
                  }
                  });"
            ),
            selection = list(mode = "single")
            #          filter = "top"
          )
      }
      TFs
    }else {
      makeNullTable()
    }
  })
  
  #create the subnetworkTable
  output$subnetworkTable =  DT::renderDataTable({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID)){
      comparison<-comparisonTable()[v$comparisonRow,]
      comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep="_")
      subnetworkTable <-getSubnetworkTable()
	  colnames(subnetworkTable)[4:5] <- c("P-Value","Top GO Term")
      
      subnetworkTable <-
        datatable(
          subnetworkTable,
          rownames = F,
          options = list(
            pageLength = 50,
            scrollY = "500px",
            lengthChange = FALSE,
            order = list(3, 'asc'),
            bInfo = FALSE,
            dom = 'frtipB',
            buttons = list('copy',list(extend='csv',filename=paste(getAccession(),comparison,"subnetworks",sep="_"))),
            bPaginate = FALSE,referRender=TRUE
          ),
          selection = list(mode = "single"),
          filter = list(position = 'top', clear = FALSE),callback = JS("table.on('click.dt', 'td', function() {
            var row_=table.cell(this).index().row;
            var col=table.cell(this).index().column;
            var rnd= Math.random();
            var data = [row_, col, rnd];
           Shiny.onInputChange('subRows',data );
    });")
        )
      subnetworkTable <-formatSignif(subnetworkTable, columns = 4, digits = 3)
      subnetworkTable
    }else {
      makeNullTable()
    }
  }
  )
  
  #create the subnetwork info div
  output$subnetworkInfo<-renderUI({
    if(input$tabset=="subnet" & !is.null(v$comparisonRow)){
      if(is.null(v$subnetworkRow)){
        div(h4("Subnetwork:",h4("Click on a row in the table to select a subnetwork to visualise",style="color:red")))
      } else{
        data<-getSubnetworkTable()
        div(h4("Subnetwork:"),h5("Subnetwork No. :",input$subnetworkTable_row_last_clicked,HTML("<br>"),"Function:",data[input$subnetworkTable_row_last_clicked,5]))
      }
      
    }
    
    
  })
  
  #create the comparison info div
  output$compareInfo<-renderUI({
    if(input$tabset=="response" & input$toptabset=="comparisons" & !is.null(v$comparisonRow)){
      print(input$tabset)
      if (input$responseType=="Pairwise"){
        if(is.null(v$simSummaryRow)){
          div(h4("Compare:",h4("Select another experiment to see the gene overlap",style="color:red")))
        } else{
          rowID = input$similaritySummary_row_last_clicked
          ID = similaritySummary()[rowID,"ID"]
          accession = gsub("(.*)_.*", "\\1",ID)
          comparison = similaritySummary()[rowID,"Comparison"]
          div(h4("Compare:"),h5("Experiment. :",accession,HTML("<br>"),"Comparison:",comparison))
        }
      } else {
        if(is.null(v$similaritySummaryMultiRow)){
          div(h4("Select other experiments to calculate consensus signatures",style="color:red"))
        }
        
      }
    }
  })

  
  
  #interactive subnetwork
  output$subnetwork = renderVisNetwork({
    comparisonID = v$comparisonRow
    if (!is.null(v$subnetworkRow) & !is.null(comparisonID)) {
      load(paste0(
        "data/",
        getAccession(),
        "/",
        getAccession(),
        "_networks_",comparisonID,".rdata"
      ))
      networks.string<-networks
      load(paste0(
        "data/",
        getAccession(),
        "/",
        getAccession(),
        "_biogridnetworks_",comparisonID,".rdata"
      ))
      networks <- c(networks.string,networks)
      maxValue <- signif(max(abs(unlist(sapply(networks,function(x) as.numeric(levels(x$x$nodes$title)))))),2)
      if(!file.exists(paste0("www/legend/tempNode",maxValue,".png"))){
      colfunc <- colorRampPalette(c("green", "white","red"))
      color.bar(colfunc(100),-maxValue)
      }
      
      networks[[input$subnetworkTable_row_last_clicked]] %>% visEdges(smooth = FALSE) %>%   visPhysics(enabled = F) %>% 
        
      visLegend(main = list(text = "log2 FoldChange"),addNodes = list(shape = "image",
                                image = paste0("legend/tempNode",maxValue,".png"),shapeProperties=list(useImageSize=T)),  useGroups = T,position="left",width=0.1) %>%
        
        visExport(
          type = "png",
          name = "subnetwork",
          float = "middle",
          label = "Save network",
          background = "white",
          style = "color: #ffffff;background-color: #2c3e50;border-color: #2c3e50;    display: inline-block;margin-bottom: 0;
          font-weight: normal;
          text-align: center;
          vertical-align: middle;
          -ms-touch-action: manipulation;
          touch-action: manipulation;
          cursor: pointer;
          background-image: none;
          border: 1px solid transparent;
          white-space: nowrap;
          padding: 10px 15px;
          font-size: 15px;
          line-height: 1.42857143;
          border-radius: 4px;"
        )
    }
  })
  
  color.bar <- function(lut, minVal, maxVal=-minVal, nticks=11, ticks=pretty(c(minVal, maxVal)), title='') {
    
    png(filename = paste0("www/legend/tempNode",maxVal,".png"),width = 300,height=500)
    scale = (length(lut)-1)/(maxVal-minVal)
    plot(c(0,10), c(minVal,maxVal), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1,lwd = 2,lwd.ticks = 2,cex.axis = 1.5)
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + minVal
      rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
    dev.off()
  }
  

  
  #create the cosineSim plot
  output$cosineSim = renderPlot({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID)){
      ID <- paste0(getAccession(),"_",comparisonID)
      cosinePlot <- getCosineHist(ID,cosineZScores)
      cosinePlot
    }
      
    
  })
  
  #create the jaccardSim plot
  output$jaccardSim = renderPlot({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID)){
      
      ID <- paste0(getAccession(),"_",comparisonID)
      
      jaccardPlot <- getJaccardHist(ID,jaccardZScores)
      jaccardPlot
    }
    
    
  })
  
  #create the jaccard plot for the sig search
  output$sigJaccZscore = renderPlot({
      jaccardPlot <- getJaccardZscoreHist()
      jaccardPlot
    })
  
  #create the chrdir plot for the sig search
  output$chrDirZscore = renderPlot({
    chrDirPlot <- getChrDirZscoreHist()
    chrDirPlot
  })
  
  #load the similarity calculations
  getSimilarityTable <- reactive({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID)){
      ID <- paste0(getAccession(),"_",comparisonID)
     
      dt <- read.delim(paste0("data/similarity/", ID,".txt"))
    }
  })
  
  
  #get the similarity summary
  similaritySummary <- reactive({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID)){
      dt <- getSimilarityTable()
      remove <- which(accessions$combined == paste0(getAccession(),"_",comparisonID))
      dt$ID<-accessions[-remove,1]
      dt$comparison <- accessions[match(dt$comparisonID,accessions$combined),"comparisonsText"]
      dt<-merge(dt,expTable,by.x="ID",by.y="ID")
      dt$ID = NULL
      
      dt<-dt[,c(2,3,4,5,7,6,1,8:12)]
      colnames(dt)[6:7] <- c("Comparison","ID")
      dt<-dt[order(dt[,1],decreasing = T),]
      dt[,1] <- signif(dt[,1], digits = 3)
      dt[,2] <- signif(dt[,2], digits = 3)
      dt[,3] <- signif(dt[,3], digits = 3)
      dt[,4] <- signif(dt[,4], digits = 3)
      colnames(dt)[1:4] <- c("Signed Jaccard (FC)","Cosine","Signed Jaccard (Sig)","Signed Jaccard (ChrDir)")
      dt
    }
  })
  
  #out the similarity summary
  output$similaritySummary = renderDataTable({
      id = input$expTable_row_last_clicked
      comparisonID = v$comparisonRow
      if (!is.null(id) & !is.null(comparisonID)){
        comparison<-comparisonTable()[v$comparisonRow,]
        comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep="_")
        dt <- similaritySummary()
        
        maxVal<-max(abs(dt[,1:4]),na.rm = T) + 0.001
        minVal<- -maxVal
        breaksList = seq(minVal, maxVal, by = 0.001)
        colfunc <- colorRampPalette(c("green", "white","red"))
        map <- makecmap(dt[,1],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
        colours<-  cmap(dt[,1], map = map)
        style.jacc <- styleEqual(dt[,1], colours)
        map <- makecmap(dt[,2],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
        colours<-  cmap(dt[,2], map = map)
        style.cosine <- styleEqual(dt[,2], colours)
        map <- makecmap(dt[,3],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
        colours<-  cmap(dt[,3], map = map)
        style.sigjacc <- styleEqual(dt[,3], colours)
        map <- makecmap(dt[,4],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
        colours<-  cmap(dt[,4], map = map)
        style.chrDir <- styleEqual(dt[,4], colours)
        
        dt<-datatable(
          dt,
          rownames = F,
          extensions = c('Scroller',"Buttons"),
          options = list(
            scrollY = 300,
            scrollX = 500,
            scroller = TRUE,
            pageLength = 50,
            lengthChange = TRUE,
            bInfo = FALSE,
            order = list(0, 'desc'),
            buttons = list('copy',list(extend='csv',filename=paste(getAccession(),comparison,"similarity",sep="_"))),
            dom = 'frtipB',referRender=TRUE
          ),
          selection = list(mode = "single"),
          filter = list(position = 'top', clear = FALSE, plain=TRUE),callback = JS("table.on('click.dt', 'td', function() {
              var row_=table.cell(this).index().row;
              var col=table.cell(this).index().column;
              var rnd= Math.random();
              var data = [row_, col, rnd];
             Shiny.onInputChange('simSummaryRows',data );
      });")
        )
        dt <- formatStyle(dt,1, backgroundColor = style.jacc)
        dt <- formatStyle(dt,2, backgroundColor = style.cosine)
        dt <- formatStyle(dt,3, backgroundColor = style.sigjacc)
        dt <- formatStyle(dt,4, backgroundColor = style.chrDir)
        dt
      }else {
        makeNullTable()
      }
  })
  
  #out the similarity summary
  output$similaritySummaryMulti = renderDataTable({
    id = input$expTable_row_last_clicked
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID)){
      comparison<-comparisonTable()[v$comparisonRow,]
      comparison<-paste(comparison["Numerator"],"vs",comparison["Denominator"],sep="_")
      dt <- similaritySummary()
      
      maxVal<-max(abs(dt[,1:4]),na.rm = T) + 0.001
      minVal<- -maxVal
      breaksList = seq(minVal, maxVal, by = 0.001)
      colfunc <- colorRampPalette(c("green", "white","red"))
      map <- makecmap(dt[,1],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
      colours<-  cmap(dt[,1], map = map)
      style.jacc <- styleEqual(dt[,1], colours)
      map <- makecmap(dt[,2],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
      colours<-  cmap(dt[,2], map = map)
      style.cosine <- styleEqual(dt[,2], colours)
      map <- makecmap(dt[,3],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
      colours<-  cmap(dt[,3], map = map)
      style.sigjacc <- styleEqual(dt[,3], colours)
      map <- makecmap(dt[,4],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
      colours<-  cmap(dt[,4], map = map)
      style.chrDir <- styleEqual(dt[,4], colours)
      
      dt<-datatable(
        dt,
        rownames = F,
        extensions = c('Scroller',"Buttons","Select"),
        options = list(
          scrollY = 300,
          scrollX = 500,
          scroller = TRUE,
          pageLength = 50,
          lengthChange = TRUE,
          bInfo = FALSE,
          order = list(0, 'desc'),
          buttons = list('copy',list(extend='csv',filename=paste(getAccession(),comparison,"similarity",sep="_"))),
          dom = 'frtipB',referRender=TRUE
        ),
        selection = list(mode = "multiple"),
        filter = list(position = 'top', clear = FALSE, plain=TRUE),callback = JS("table.on('click.dt', 'td', function() {
              var row_=table.cell(this).index().row;
              var col=table.cell(this).index().column;
              var rnd= Math.random();
              var data = [row_, col, rnd];
             Shiny.onInputChange('similaritySummaryMultiRow',data );
      });
      table.on('click.dt', 'tr', function() {
        if ($(this).hasClass('selected'))  {
        $(this).toggleClass('row_selected');  
      } else if(table.$('tr.selected').length >5){
        $(this).toggleClass('selected');
        }
    });")
      )
      dt <- formatStyle(dt,1, backgroundColor = style.jacc)
      dt <- formatStyle(dt,2, backgroundColor = style.cosine)
      dt <- formatStyle(dt,3, backgroundColor = style.sigjacc)
      dt <- formatStyle(dt,4, backgroundColor = style.chrDir)
      dt
    }else {
      makeNullTable()
    }
  })
  
 
  rankProductTable <- reactive({ 
  id = getAccession()
  comparisonID = v$comparisonRow
  if (!is.null(id) & !is.null(comparisonID) & !is.null(v$similaritySummaryMultiRow)){
    print("hello")
    IDs = as.character(similaritySummary()[input$similaritySummaryMulti_rows_selected,"ID"])
    IDs <- c(paste(id,comparisonID,sep="_"),IDs)
    fc <- as.data.frame(getFoldChangeSignatures(IDs))
    
    fc[,2:ncol(fc)] <- signif(fc[,2:ncol(fc)], digits = 3)
    
    fc <- fc[apply(fc[,-1],1,function(x) length(which(is.na(x)))<(length(x)*0.5)),]
    
    rankProd <- RankProducts(data = fc[,-1],cl = rep(as.factor(1),ncol(fc)-1))$pfp[,2:1]
    print(head(rankProd))
    colnames(rankProd) <- c("UpRegulated PFP","DownRegulated PFP")
    rankProd[, 1:2] <- signif(rankProd[, 1:2], digits = 3)
    fc <- cbind(fc,rankProd)
    fc
  }
    })
  
  
  output$sharedFC = renderDataTable({
    id = getAccession()
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID) & !is.null(v$similaritySummaryMultiRow)){
      
      shinyjs::hide("sharedFC")
      shinyjs::hide("consensusPathway")
      shinyjs::show("loadingTable")

      
      fc <- rankProductTable()
      
      maxVal<-max(abs(fc[,c(-1,-ncol(fc),-ncol(fc)-1)]),na.rm = T) + 0.001
      minVal<- -maxVal
      breaksList = seq(minVal, maxVal,length.out = 100)
      colfunc <- colorRampPalette(c("green", "white","red"))
      map <- makecmap(fc[,-1],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
      values <- unique(unlist(c(fc[,-1])))
      colours <-  cmap(values, map = map)
      style <- styleEqual(values, colours)
      
      


      dt<-datatable(
        fc,
        rownames = F,
        extensions = 'Scroller',
        caption = tags$caption(
          style = 'caption-side: centre; color: black;text-align: centre;', tags$b("Fold Changes")),
        options = list(
          scrollY = 300,
          scrollX = 500,
          scroller = FALSE,
          pageLength = 50,
          lengthChange = FALSE,
          bInfo = FALSE,
          buttons = c('copy', 'csv'),dom = 'frtipB',
          order = list(ncol(fc)-1, 'asc'),
          deferRender=TRUE
        ),
        selection = list(mode = "single"),
        filter = list(position = 'top', clear = FALSE)
        ,escape=F)
      #dt <- formatStyle(dt,2:ncol(fc)-2, backgroundColor = style)
      shinyjs::hide("loadingTable")
      shinyjs::show("sharedFC")
      shinyjs::show("consensusPathway")
      
    
    dt
    }
    
  })
  
  output$sharedChrDir = renderDataTable({
    id = getAccession()
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID) & !is.null(v$similaritySummaryMultiRow)){
      
      res <- getChrDirSigs()
      
      if(is.na(res)){
        dt<-as.data.frame("No characteristic direction signatures available")
        colnames(dt)=""
        dt <-
          datatable(
            dt,
            rownames = F,
            options = list(bInfo = FALSE,dom = 'rti')) %>% formatStyle(columns = 1,color = 'red')
        
      } else {
  
        dt<-datatable(
          res,
          rownames = F,
          caption = tags$caption(
            style = 'caption-side: centre; color: black;text-align: centre;', tags$b("Characteristic Direction Signatures")),
          extensions = 'Scroller',
          options = list(
            scrollY = 300,
            scrollX = 500,
            scroller = TRUE,
            pageLength = 50,
            lengthChange = TRUE,
            bInfo = FALSE,
            buttons = c('copy', 'csv'),referRender=TRUE
          ),
          selection = list(mode = "single"),
          filter = list(position = 'top', clear = FALSE)
          ,escape=F)
        
        backgroundColor = styleEqual(c(0, 1,-1), c('white', 'red', 'green'))
        dt <- formatStyle(dt,1:ncol(res), backgroundColor = backgroundColor)
      }
        
      dt
    }
    
  })
  
  output$sharedSigs = renderDataTable({
    id = getAccession()
    comparisonID = v$comparisonRow
    if (!is.null(id) & !is.null(comparisonID) & !is.null(v$similaritySummaryMultiRow)){
      
      res <- getSigs()
      
      if(is.na(res)){
        dt<-as.data.frame("No differential expression signatures available")
        colnames(dt)=""
        dt <-
          datatable(
            dt,
            rownames = F,
            options = list(bInfo = FALSE,dom = 'rti')) %>% formatStyle(columns = 1,color = 'red')
        
      } else{
        dt<-datatable(
          res,
          rownames = F,
          caption = tags$caption(
            style = 'caption-side: centre; color: black;text-align: centre;', tags$b("Differentially Expressed Genes")),
          extensions = 'Scroller',
          options = list(
            scrollY = 300,
            scrollX = 500,
            scroller = TRUE,
            pageLength = 50,
            lengthChange = TRUE,
            bInfo = FALSE,
            buttons = c('copy', 'csv'),referRender=TRUE
          ),
          selection = list(mode = "single"),
          filter = list(position = 'top', clear = FALSE)
          ,escape=F)
        
        backgroundColor = styleEqual(c(0, 1,-1), c('white', 'red', 'green'))
        dt <- formatStyle(dt,1:ncol(res), backgroundColor = backgroundColor)
      }
        
        dt
    }
    
  })
  
  
  getconsensusPathwayTable <- eventReactive(input$pathwaySubmit, {
    
    id = getAccession()
    comparisonID = v$comparisonRow
    
    if (!is.null(id) & !is.null(comparisonID) & !is.null(v$similaritySummaryMultiRow)){
    
      print("disabling")
      shinyjs::disable("pathwaySubmit")
      shinyjs::show("hiddenLoad")
      shinyjs::hide("consensusPathway")
      plotReady$ok <- FALSE
      
    
      print(input$thresholdSlider)
      threshold <- input$thresholdSlider
      database <- input$databaseSelect
      rp <- rankProductTable()
      genes <- na.omit(c(rp[rp[,ncol(rp)]<=threshold,1],rp[rp[,ncol(rp)-1]<=threshold,1]))
      
      if(length(genes)==0){
        plotReady$ok <- TRUE
        shinyjs::hide("hiddenLoad")
        shinyjs::show("consensusPathway")
        return(NULL)
      }
      
      if(length(genes)>5000){
        plotReady$ok <- TRUE
        shinyjs::hide("hiddenLoad")
        shinyjs::show("consensusPathway")
        return(tooLargeInputError())
      }
      
      pathways <- try(enrichr(genes,database)[[1]])
      if(is(pathways, "try-error") || nrow(pathways) == 0) {
        plotReady$ok <- TRUE
        shinyjs::hide("hiddenLoad")
        shinyjs::show("consensusPathway")
        return(NULL)
      }
      
      pathways <-pathways[,c(-3,-5,-6)]
      pathways[,3:5] <- signif(pathways[,3:5], digits = 3)
      pathways[,6] <- gsub(";"," ",pathways[,6])
      
      plotReady$ok <- TRUE
      shinyjs::hide("hiddenLoad")
      shinyjs::show("consensusPathway")
      pathways
    }
  })
    
    getconsensusPathwayTableChrDir <- eventReactive(input$pathwaySubmitChrDir, {
      
      id = getAccession()
      comparisonID = v$comparisonRow
      
      if (!is.null(id) & !is.null(comparisonID) & !is.null(v$similaritySummaryMultiRow)){
      
      print("disablingChrDirButton")
      shinyjs::disable("pathwaySubmitChrDir")
      shinyjs::show("chrDirPathwayLoad")
      shinyjs::hide("consensusPathwayChrDir")
      plotReady$chrDir <- FALSE
      
      print("hello")
      print(input$voteSlider)
      threshold <- input$voteSlider
      database <- input$databaseSelectChrDir
      chrDirTable <- getChrDirSigs()
      genes <- as.character(na.omit(chrDirTable[ apply(chrDirTable[,-ncol(chrDirTable)],1,function(x) length(which(x==1))>=threshold | length(which(x==-1))>=threshold),ncol(chrDirTable)]))
      print(head(genes))

      if(length(genes)==0){
        plotReady$chrDir <- TRUE
        # shinyjs::hide("hiddenLoadChrDir")
        # shinyjs::show("consensusPathwayChrDir")
        shinyjs::show("pathwaySubmitChrDir")
        return(NULL)
      }
      
      if(length(genes)>5000){
        plotReady$chrDir <- TRUE
        # shinyjs::hide("hiddenLoadChrDir")
        # shinyjs::show("consensusPathwayChrDir")
        shinyjs::show("pathwaySubmitChrDir")
        return(tooLargeInputError())
      }
      
      pathways <- try(enrichr(genes,database)[[1]])
      if(is(pathways, "try-error")  || nrow(pathways) == 0 ) {
        plotReady$chrDir <- TRUE
        shinyjs::hide("hiddenLoadChrDir")
        shinyjs::show("consensusPathwayChrDir")
        return(NULL)
      }
      
      pathways <-pathways[,c(-3,-5,-6)]
      pathways[,3:5] <- signif(pathways[,3:5], digits = 3)
      pathways[,6] <- gsub(";"," ",pathways[,6])

      plotReady$chrDir <- TRUE
      #shinyjs::hide("chrDirPathwayLoad")
      shinyjs::show("pathwaySubmitChrDir")
      pathways
      }
  })
    
  getconsensusPathwayTableSig <- eventReactive(input$pathwaySubmitSig, {
    
    id = getAccession()
    comparisonID = v$comparisonRow
    
    if (!is.null(id) & !is.null(comparisonID) & !is.null(v$similaritySummaryMultiRow)){
      
      print("disablingSigButton")
      shinyjs::disable("pathwaySubmitSig")
      shinyjs::show("sigPathwayLoad")
      shinyjs::hide("consensusPathwaySig")
      plotReady$sig <- FALSE
      
      print("hello")
      print(input$voteSliderSig)
      threshold <- input$voteSliderSig
      database <- input$databaseSelectSig
      sigTable <- getSigs()
      genes <- as.character(sigTable[ apply(sigTable[,-ncol(sigTable)],1,function(x) length(which(x==1))>=threshold | length(which(x==-1))>=threshold),ncol(sigTable)])
      print(genes)
      
      if(length(genes)==0){
        plotReady$sig <- TRUE
        shinyjs::hide("hiddenLoad")
        shinyjs::show("consensusPathway")
        return(NULL)
      }
      
      if(length(genes)>5000){
        plotReady$sig <- TRUE
        shinyjs::hide("hiddenLoad")
        shinyjs::show("consensusPathway")
        return(tooLargeInputError())
      }
      
      pathways <- try(enrichr(genes,database)[[1]])
      if(is(pathways, "try-error")  || nrow(pathways) == 0) {
        plotReady$sig <- TRUE
        shinyjs::hide("hiddenLoad")
        shinyjs::show("consensusPathway")
        return(NULL)
      }
      
      pathways <-pathways[,c(-3,-5,-6)]
      pathways[,3:5] <- signif(pathways[,3:5], digits = 3)
      pathways[,6] <- gsub(";"," ",pathways[,6])
      
      plotReady$sig <- TRUE
      #shinyjs::hide("chrDirPathwayLoad")
      shinyjs::show("pathwaySubmitSig")
      pathways
    }
    })
  
  
    output$consensusPathway = renderDataTable({
      
      #getSlider Value
      if(plotReady$ok){
        shinyjs::enable("pathwaySubmit")
        terms <- getconsensusPathwayTable()
        shinyjs::hide("hiddenLoad")
        
        if (inherits(terms,"tooLargeInputError")){
          
          terms<-as.data.frame("Gene list length too large. Increase the threshold")
          colnames(terms)=""
          terms <-
            datatable(
              terms,
              rownames = F,
              options = list(bInfo = FALSE,dom = 'rti')
            ) %>% formatStyle(columns = 1,color = 'red')
        } else if (!is.null(terms)){
          
          terms <-
            datatable(
              cbind(' ' = '&oplus;', terms),escape = FALSE,
              rownames = F,
              caption = tags$caption(
                style = 'caption-side: centre; color: black;text-align: centre;', tags$b("Consensus Pathways")),
              options = list(lengthMenu = c(10,50,100,nrow(terms)), order = list(3, 'asc'),
                             columnDefs = list(
                               list(visible = FALSE, targets = 6),
                               list(orderable = FALSE, className = 'details-control', targets = 0)
                             ),
                             pageLength = 50,
                             scrollY = "500px",
                             lengthChange = FALSE,
                             bInfo = FALSE,
                             dom = 'frtipB',
                             buttons = list('copy',list(extend='csv',filename="enrichR_Results.txt")),
                             bPaginate = FALSE, referRender=TRUE
              ),
              callback = JS("
                  table.column(1).nodes().to$().css({cursor: 'pointer'});
                  var format = function(d) {
                  return '<div style=\"background-color:#eee; padding: .5em;\">' +
                  '<b>'+\"Genes: \" + '</b>'+ d[6] +  '</div>';
                  };
                  table.on('click', 'td.details-control', function() {
                  var td = $(this), row = table.row(td.closest('tr'));
                  if (row.child.isShown()) {
                  row.child.hide();
                  td.html('&oplus;');
                  } else {
                  row.child(format(row.data())).show();
                  td.html('&CircleMinus;');
                  }
                  });"
              ),
              selection = list(mode = "single")
              #          filter = "top"
            )
        } else{
          terms<-as.data.frame("No significant pathways")
          colnames(terms)=""
          terms <-
            datatable(
              terms,
              rownames = F,
              options = list(bInfo = FALSE,dom = 'rti')
            ) %>% formatStyle(columns = 1,color = 'red')
        }
      shinyjs::show("consensusPathway")
      terms
        
        
      }
    })
    
    output$consensusPathwayChrDir = renderDataTable({
      
      #getSlider Value
      if(plotReady$chrDir){
        terms <- getconsensusPathwayTableChrDir()
        shinyjs::enable("pathwaySubmitChrDir")
        shinyjs::hide("chrDirPathwayLoad")

        if (inherits(terms,"tooLargeInputError")){
          
          terms<-as.data.frame("Gene list length too large. Increase the threshold")
          colnames(terms)=""
          terms <-
            datatable(
              terms,
              rownames = F,
              options = list(bInfo = FALSE,dom = 'rti')
            ) %>% formatStyle(columns = 1,color = 'red')
        } else if (!is.null(terms)){
          
          terms <-
            datatable(
              cbind(' ' = '&oplus;', terms),escape = FALSE,
              rownames = F,
              caption = tags$caption(
                style = 'caption-side: centre; color: black;text-align: centre;', tags$b("Consensus Pathways")),
              options = list(lengthMenu = c(10,50,100,nrow(terms)), order = list(3, 'asc'),
                             columnDefs = list(
                               list(visible = FALSE, targets = 6),
                               list(orderable = FALSE, className = 'details-control', targets = 0)
                             ),
                             pageLength = 50,
                             scrollY = "500px",
                             lengthChange = FALSE,
                             bInfo = FALSE,
                             dom = 'frtipB',
                             buttons = list('copy',list(extend='csv',filename="enrichR_Results.txt")),
                             bPaginate = FALSE, referRender=TRUE
              ),
              callback = JS("
                            table.column(1).nodes().to$().css({cursor: 'pointer'});
                            var format = function(d) {
                            return '<div style=\"background-color:#eee; padding: .5em;\">' +
                            '<b>'+\"Genes: \" + '</b>'+ d[6] +  '</div>';
                  };
                  table.on('click', 'td.details-control', function() {
                  var td = $(this), row = table.row(td.closest('tr'));
                  if (row.child.isShown()) {
                  row.child.hide();
                  td.html('&oplus;');
                  } else {
                  row.child(format(row.data())).show();
                  td.html('&CircleMinus;');
                  }
                  });"
              ),
              selection = list(mode = "single")
              #          filter = "top"
            )
        } else{
          terms<-as.data.frame("No significant pathways")
          colnames(terms)=""
          terms <-
            datatable(
              terms,
              rownames = F,
              options = list(bInfo = FALSE,dom = 'rti')
            ) %>% formatStyle(columns = 1,color = 'red')
        }
        shinyjs::show("consensusPathwayChrDir")
        terms
        
        
      }
    })
    
    output$consensusPathwaySig= renderDataTable({
      
      #getSlider Value
      if(plotReady$sig){
        terms <- getconsensusPathwayTableSig()
        shinyjs::enable("pathwaySubmitSig")
        shinyjs::hide("sigPathwayLoad")
        if (inherits(terms,"tooLargeInputError")){
          
          terms<-as.data.frame("Gene list length too large. Increase the threshold")
          colnames(terms)=""
          terms <-
            datatable(
              terms,
              rownames = F,
              options = list(bInfo = FALSE,dom = 'rti')
            ) %>% formatStyle(columns = 1,color = 'red')
        } else if (!is.null(terms)){
          
          terms <-
            datatable(
              cbind(' ' = '&oplus;', terms),escape = FALSE,
              rownames = F,
              caption = tags$caption(
                style = 'caption-side: centre; color: black;text-align: centre;', tags$b("Consensus Pathways")),
              options = list(lengthMenu = c(10,50,100,nrow(terms)), order = list(3, 'asc'),
                             columnDefs = list(
                               list(visible = FALSE, targets = 6),
                               list(orderable = FALSE, className = 'details-control', targets = 0)
                             ),
                             pageLength = 50,
                             scrollY = "500px",
                             lengthChange = FALSE,
                             bInfo = FALSE,
                             dom = 'frtipB',
                             buttons = list('copy',list(extend='csv',filename="enrichR_Results.txt")),
                             bPaginate = FALSE, referRender=TRUE
              ),
              callback = JS("
                            table.column(1).nodes().to$().css({cursor: 'pointer'});
                            var format = function(d) {
                            return '<div style=\"background-color:#eee; padding: .5em;\">' +
                            '<b>'+\"Genes: \" + '</b>'+ d[6] +  '</div>';
                  };
                  table.on('click', 'td.details-control', function() {
                  var td = $(this), row = table.row(td.closest('tr'));
                  if (row.child.isShown()) {
                  row.child.hide();
                  td.html('&oplus;');
                  } else {
                  row.child(format(row.data())).show();
                  td.html('&CircleMinus;');
                  }
                  });"
              ),
              selection = list(mode = "single")
              #          filter = "top"
            )
        } else{
          terms<-as.data.frame("No significant pathways")
          colnames(terms)=""
          terms <-
            datatable(
              terms,
              rownames = F,
              options = list(bInfo = FALSE,dom = 'rti')
            ) %>% formatStyle(columns = 1,color = 'red')
        }
        shinyjs::show("consensusPathwaySig")
        terms
        
        
      }
    })
  
  
  
  #get GOTerms given an ID
  getGOTerms <- function(ID) {
    accession = gsub("(.*)_.*", "\\1",ID)
    comparisonID = gsub(".*_(.*)", "\\1",ID) 
    terms <- paste0("data/", accession, "/", accession, "_goterms_",comparisonID,".txt")
  if(file.size(terms)<5){
    return(as.data.frame(NA))
  } else {
    terms <- read.delim(terms)
    return(terms)
  }
  }
  
  #get TFTerms given an ID
  getTFs <- function(ID) {
    accession = gsub("(.*)_.*", "\\1",ID)
    comparisonID = gsub(".*_(.*)", "\\1",ID) 
    tfs <- paste0("data/", accession, "/", accession, "_TFs_",comparisonID,".txt")
    print(tfs)
    if(file.size(tfs)<5){
      return(as.data.frame(NA))
    } else {
      tfs <- read.delim(tfs)[, c(-1,-8)]
      return(tfs)
    }
  }
  
  #get foldchanges given IDs
  getFoldChangeSignatures <- function(IDs) {
    accession = gsub("(.*)_.*", "\\1",IDs)
    comparisonID = gsub(".*_(.*)", "\\1",IDs) 
    fc <- foldChangeTable[,c("ID",IDs)]
    fc <- fc[rowSums(fc[,-1])!=0,]
  }
  
  #get chrdir sigs given an ID
  getChrDirSigs <- reactive({
    id = getAccession()
    comparisonID = v$comparisonRow
    IDs = as.character(similaritySummary()[input$similaritySummaryMulti_rows_selected,"ID"])
    IDs <- c(paste(id,comparisonID,sep="_"),IDs)
    accession = gsub("(.*)_.*", "\\1",IDs)
    comparisonID = gsub(".*_(.*)", "\\1",IDs)
    chrDirs <- chrDirsList[IDs]
    columns <- as.data.frame(rbind_list(chrDirs))
    if(ncol(columns)==0){
      return(NA)
    }
    columns<-unique(unlist(columns[,1]))
    res<-as.data.frame(t(do.call(rbind,lapply(chrDirs,function(x) ifelse(columns %in% x[x[,2]>0,1],1,ifelse(columns %in% x[x[,2]<0,1],-1,0))))))
    res$Gene <- columns
    res <- res[order(rowSums(abs(res[,-ncol(res)])),decreasing = T),]
    
    plotReady$chrDirCols <- ncol(res)-1
    
    return(res)
    
  })
  
  #get fc sigs given an ID
  getSigs <- reactive({
    id = getAccession()
    comparisonID = v$comparisonRow
    IDs = as.character(similaritySummary()[input$similaritySummaryMulti_rows_selected,"ID"])
    IDs <- c(paste(id,comparisonID,sep="_"),IDs)
    accession = gsub("(.*)_.*", "\\1",IDs)
    comparisonID = gsub(".*_(.*)", "\\1",IDs)
    up <- upSigs[IDs]
    down <- downSigs[IDs]
    columns<-unique(c(unlist(up),unlist(down)))
    print(head(columns))
    res<-as.data.frame(t(do.call(rbind,mapply(function(Up,Down) ifelse(columns %in% Up,1,ifelse(columns %in% Down,-1,0)),up,down,SIMPLIFY = F))))
    print(head(res))
    res$Gene <- columns
    print(head(res))
    if(nrow(res)==0){
      return(NA)
    }
    
  res <- res[order(rowSums(abs(res[,-ncol(res)])),decreasing = T),]

	plotReady$sigCols <- ncol(res)-1
    
    return(res)
    
  })
  
    
  #get the selected comparison fold change table
  getComparisonFoldChangeTable = reactive({
    rowID = input$similaritySummary_row_last_clicked
    if (!is.null(rowID)){
      ID = similaritySummary()[rowID,"ID"]
      accession = gsub("(.*)_.*", "\\1",ID)
      comparisonID = gsub(".*_(.*)", "\\1",ID) 
      fc <- read.delim(paste0("data/", accession, "/", accession, "_FC_",comparisonID,".txt"))
    }
  })
  
  #get the selected comparison chrDir table
  getComparisonChrDirTable = reactive({
    rowID = input$similaritySummary_row_last_clicked
    if (!is.null(rowID)){
      ID = similaritySummary()[rowID,"ID"]
      accession = gsub("(.*)_.*", "\\1",ID)
      comparisonID = gsub(".*_(.*)", "\\1",ID) 
      chrdirs<-paste0("data/", accession, "/", accession, "_chrDirTable.txt")
      if(file.size(chrdirs)>5){
        chrdirs <- read.delim(chrdirs)
        if ("gene_name" %in% colnames(chrdirs)){
          chrdirs <- chrdirs[ chrdirs$comparisonNumber == comparisonID,c("gene_name","chrDir"),]
        } else {
          chrdirs <- chrdirs[ chrdirs$comparisonNumber == comparisonID,c("ID","chrDir"),]
        }
        colnames(chrdirs)<-c("Gene Symbol","Characterstic Direction")
        chrdirs[, 2] <- chrdirs[, 2]*-1
        chrdirs
      }
    }
  })
  
  #get the selected sigResponse comparison fold change table
  getResponseComparisonFoldChangeTable = reactive({
    rowID = input$sigSummary_row_last_clicked
    if (!is.null(rowID)){
      ID = getJaccardSim()[rowID,"ID"]
      accession = gsub("(.*)_.*", "\\1",ID)
      comparisonID = gsub(".*_(.*)", "\\1",ID) 
      fc <- read.delim(paste0("data/", accession, "/", accession, "_FC_",comparisonID,".txt"))
    }
  })
  
  #get the selected sigResponse comparison fold change table
  getResponseComparisonChrDirTable = reactive({
    rowID = input$sigSummary_row_last_clicked
    if (!is.null(rowID)){
      ID = getJaccardSim()[rowID,"ID"]
      accession = gsub("(.*)_.*", "\\1",ID)
      comparisonID = gsub(".*_(.*)", "\\1",ID) 
      chrdirs<-paste0("data/", accession, "/", accession, "_chrDirTable.txt")
      if(file.size(chrdirs)<5){return(NA)}
      chrdirs <- read.delim(chrdirs)
      if ("gene_name" %in% colnames(chrdirs)){
        chrdirs <- chrdirs[ chrdirs$comparisonNumber == comparisonID,c("gene_name","chrDir"),]
      } else {
        chrdirs <- chrdirs[ chrdirs$comparisonNumber == comparisonID,c("ID","chrDir"),]
      }
      colnames(chrdirs)<-c("Gene Symbol","Characterstic Direction")
      chrdirs[, 2] <- chrdirs[, 2]*-1
      chrdirs[, 2] <- signif(chrdirs[, 2], digits = 3)
      chrdirs <- chrdirs[order(chrdirs[,2],decreasing = T),]
      chrdirs
    }
  })
  
  #map gene symbols to human
  convertGenes <- function(genes,species){
    
    map <- human2otherspecies[ human2otherspecies$species==species,]
   
    genes <- unique(map[ match(genes,map[,"queryGene"]),"humanGene"])
    
    print(head(genes))
    
    return(genes)
  }
  
  #find the overlap of the fold-change gene signatures
  makeVenn = reactive ({
    
    rowID = input$similaritySummary_row_last_clicked
    if (!is.null(rowID)){
        
      queryFC <- as.data.frame(getFoldChangesTable())
      comparisonFC <- as.data.frame(getComparisonFoldChangeTable())
      
      queryFCUp <- na.omit(queryFC[ queryFC[,2] >= log2(1.5),])
      queryFCDown <- na.omit(queryFC[ queryFC[,2] <= log2(1/1.5),])
      
      comparisonFCUp <- na.omit(comparisonFC[ comparisonFC[,2] >= log2(1.5),])
      comparisonFCDown <- na.omit(comparisonFC[ comparisonFC[,2] <= log2(1/1.5),])
      
      querySpecies <- expTable[expTable$ID==getAccession(),"Species"]
      comparisonID =  gsub("(.*)_.*", "\\1",similaritySummary()[input$similaritySummary_row_last_clicked,"ID"])
      comparisonSpecies <- expTable[expTable$ID==comparisonID,"Species"]
      querySpecies <- tolower(pull(querySpecies, Species))
      comparisonSpecies <- tolower(pull(comparisonSpecies, Species))
      
      cat(file=stderr(), dim(queryFCUp),"\n")
      cat(file=stderr(), querySpecies,"\n")
      cat(file=stderr(), dim(comparisonFCUp),"\n")
      cat(file=stderr(), comparisonSpecies,"\n")
      if (querySpecies!="human"){
        queryFCUp.genes <-convertGenes(queryFCUp[,1],querySpecies)
        queryFCDown.genes <-convertGenes(queryFCDown[,1],querySpecies)
      } else{
        queryFCUp.genes <-queryFCUp[,1]
        queryFCDown.genes <-queryFCDown[,1]
      }
      
      
      if (comparisonSpecies!="human"){
        comparisonFCUp.genes <-convertGenes(comparisonFCUp[,1],comparisonSpecies)
        comparisonFCDown.genes <-convertGenes(comparisonFCDown[,1],comparisonSpecies)
      } else{
        comparisonFCUp.genes <-comparisonFCUp[,1]
        comparisonFCDown.genes <-comparisonFCDown[,1]
      }
      cat(file=stderr(), dim(queryFCUp.genes),"\n")
      cat(file=stderr(), dim(comparisonFCUp.genes),"\n")
      cat(file=stderr(), head(queryFCUp.genes),"\n")
      cat(file=stderr(), head(comparisonFCUp.genes),"\n")
      
      
      upup <- paste(sort(intersect( queryFCUp.genes, comparisonFCUp.genes)),collapse=" ")
      downdown <- paste(sort(intersect( queryFCDown.genes, comparisonFCDown.genes)),collapse=" ")
      updown <- paste(sort(intersect( queryFCUp.genes, comparisonFCDown.genes)),collapse=" ")
      downup <- paste(sort(intersect( queryFCDown.genes ,comparisonFCUp.genes)),collapse=" ")
      
      venn <- data.frame(Overlap=c("Query:Up Comparison:Up","Query:Down Comparison:Down","Query:Up Comparison:Down","Query:Down Comparison:Up"), Genes=c(upup,downdown,updown,downup))
      
      if(ncol(queryFC)==3 & ncol(comparisonFC)==3 ){
        queryFCUp <- na.omit(queryFCUp[ queryFCUp[,3] <=0.05,])
        queryFCDown <- na.omit(queryFCDown[ queryFCDown[,3] <=0.05,])
        
        comparisonFCUp <- na.omit(comparisonFCUp[ comparisonFCUp[,3] <=0.05,])
        comparisonFCDown <- na.omit(comparisonFCDown[ comparisonFCDown[,3] <=0.05,])
        
        if (querySpecies!="human"){
          queryFCUp.genes <-convertGenes(queryFCUp[,1],querySpecies)
          queryFCDown.genes <-convertGenes(queryFCDown[,1],querySpecies)
        } else{
          queryFCUp.genes <-queryFCUp[,1]
          queryFCDown.genes <-queryFCDown[,1]
        }
        
        if (comparisonSpecies!="human"){
          comparisonFCUp.genes <-convertGenes(comparisonFCUp[,1],comparisonSpecies)
          comparisonFCDown.genes <-convertGenes(comparisonFCDown[,1],comparisonSpecies)
        } else{
          comparisonFCUp.genes <-comparisonFCUp[,1]
          comparisonFCDown.genes <-comparisonFCDown[,1]
        }
        
        upup <- paste(sort(intersect( queryFCUp.genes, comparisonFCUp.genes)),collapse=" ")
        downdown <- paste(sort(intersect( queryFCDown.genes, comparisonFCDown.genes)),collapse=" ")
        updown <- paste(sort(intersect( queryFCUp.genes, comparisonFCDown.genes)),collapse=" ")
        downup <- paste(sort(intersect( queryFCDown.genes ,comparisonFCUp.genes)),collapse=" ")
        venn.pval <- data.frame(Overlap=c("Query:Up Comparison:Up","Query:Down Comparison:Down","Query:Up Comparison:Down","Query:Down Comparison:Up"), Genes=c(upup,downdown,updown,downup))
        
      } else {
        venn.pval = as.data.frame(NULL)
      }
      return(list(venn=venn,venn.pval=venn.pval))
    }
    
        
    
  })
  
  #find the overlap of the chrDir gene signatures
  makeVennChrDir = reactive ({
    
    rowID = input$similaritySummary_row_last_clicked
    if (!is.null(rowID)){
      
      queryFC <- na.omit(as.data.frame(getChrDirTable()))
      comparisonFC <- na.omit(as.data.frame(getComparisonChrDirTable()))
      
      if(nrow(queryFC)==0 | nrow(comparisonFC)==0){return(data.frame(NULL))}
      
      querySpecies <- expTable[expTable$ID==getAccession(),"Species"]
      comparisonID =  gsub("(.*)_.*", "\\1",similaritySummary()[input$similaritySummary_row_last_clicked,"ID"])
      comparisonSpecies <- expTable[expTable$ID==comparisonID,"Species"]
      querySpecies <- tolower(pull(querySpecies, Species))
      comparisonSpecies <- tolower(pull(comparisonSpecies, Species))
      
      queryFCUp <- na.omit(queryFC[ queryFC[,2] > 0,])
      queryFCDown <- na.omit(queryFC[ queryFC[,2] < 0,])
      
      comparisonFCUp <- na.omit(comparisonFC[ comparisonFC[,2] > 0,])
      comparisonFCDown <- na.omit(comparisonFC[ comparisonFC[,2] < 0,])
      
      if (querySpecies!="human"){
        queryFCUp.genes <-convertGenes(queryFCUp[,1],querySpecies)
        queryFCDown.genes <-convertGenes(queryFCDown[,1],querySpecies)
      } else{
        queryFCUp.genes <-queryFCUp[,1]
        queryFCDown.genes <-queryFCDown[,1]
      }
      
      if (comparisonSpecies!="human"){
        comparisonFCUp.genes <-convertGenes(comparisonFCUp[,1],comparisonSpecies)
        comparisonFCDown.genes <-convertGenes(comparisonFCDown[,1],comparisonSpecies)
      } else{
        comparisonFCUp.genes <-comparisonFCUp[,1]
        comparisonFCDown.genes <-comparisonFCDown[,1]
      }
      
      upup <- paste(sort(intersect( queryFCUp.genes, comparisonFCUp.genes)),collapse=" ")
      downdown <- paste(sort(intersect( queryFCDown.genes, comparisonFCDown.genes)),collapse=" ")
      updown <- paste(sort(intersect( queryFCUp.genes, comparisonFCDown.genes)),collapse=" ")
      downup <- paste(sort(intersect( queryFCDown.genes ,comparisonFCUp.genes)),collapse=" ")
      
      venn <- data.frame(Overlap=c("Query:Postive Comparison:Positive","Query:Negative Comparison:Negative","Query:Positive Comparison:Negative","Query:Negative Comparison:Postive"), Genes=c(upup,downdown,updown,downup))
      
      return(venn)
    }
    
    
    
  })
  
  #output the gene overlap
  output$geneOverlap<- renderDataTable({
      if(!is.null(v$simSummaryRow)){
        venn <- makeVenn()
        if(!is.null(venn)){
          overlapTable <- datatable(venn[[1]],rownames = F,
          options = list(scrollY = "400px",
                         lengthChange = FALSE,
                         bInfo = FALSE,searching = FALSE, paging =FALSE,buttons = c('copy', 'csv') ),
          selection = list(mode = "single"))
        }
      }
    })
  
  #output the gene overlap
  output$geneOverlapSig<- renderDataTable({ 
    
    venn <- makeVenn()
    if(!is.null(venn) & !is.null(venn[[2]])){
    overlapTable <- datatable(venn[[2]],rownames = F,
                              options = list(scrollY = "400px",
                                             lengthChange = FALSE,
                                             bInfo = FALSE,searching = FALSE, paging =FALSE,buttons = c('copy', 'csv')),
                              selection = list(mode = "single"))
    }
  })
  
  #output the gene overlap
  output$geneOverlapChrDir<- renderDataTable({ 
    
    venn <- makeVennChrDir()
    if(!is.null(venn)){
      overlapTable <- datatable(venn,rownames = F,
                                options = list(scrollY = "400px",
                                               lengthChange = FALSE,
                                               bInfo = FALSE,searching = FALSE, paging =FALSE,buttons = c('copy', 'csv')),
                                selection = list(mode = "single"))
    }
  })
  
  #calculate jaccard similarity
  jaccard<-function(set1,set2){
    I <- length(intersect(set1,set2))
    S <- I/(length(set1)+length(set2)-I)
    return(S)
  }
  
  #calculate signed jaccard similarity
  signedJaccard<-function(set1Up,set1Down,set2Up,set2Down){
    signedJaccard<-(jaccard(set1Up,set2Up)+jaccard(set1Down,set2Down)-jaccard(set1Up,set2Down)-jaccard(set1Down,set2Up))/2
    return(signedJaccard)
  }
  
  #calculate signed jaccard similarity for differentially expressed genes
  signedJaccardSigCreate<-function(i,queryUp,queryDown,foldChangeListUp,foldChangeListDown,background){
    
    set1Up<-foldChangeListUp[[i]]
    set1Down<-foldChangeListDown[[i]]
    
    expressed <- foldChangeTable[foldChangeTable[,i] !=0,"ID"]
    
    queryUp <- queryUp[ queryUp %in% expressed$ID]
    queryDown <- queryDown[ queryDown %in% expressed$ID]
    
    if(length(background)>0){
      set1Up <- set1Up[set1Up %in% background]
      set1Down <- set1Down[set1Down %in% background]
    }
    
    return(signedJaccard(set1Up,set1Down,queryUp,queryDown))
  }
  
  #calculate signed jaccard similarity for chrDir
  signedJaccardCreateChrDir<-function(i,queryUp,queryDown,chrDirs,background){

    if (is.na( chrDirs[i])){
      return(NA)
    } 
    
    chrDir1 <- chrDirs[[i]]
    
    set1Up <- chrDir1[chrDir1$chrDir<0,1]
    set1Down<-chrDir1[chrDir1$chrDir>0,1]
    
    expressed <- foldChangeTable[foldChangeTable[,i] !=0,"ID"]
    
    queryUp <- queryUp[ queryUp %in% expressed$ID]
    queryDown <- queryDown[ queryDown %in% expressed$ID]
    
    
    return(signedJaccard(set1Up,set1Down,queryUp,queryDown))
  }
  
  #similarity search button
  getJaccardSim <- eventReactive(input$sigSearch, {
    
    shinyjs::disable("sigSearch")
    
    upRegulated <- unlist(strsplit(input$UpRegulated,"[ \t\r\n]"))
    downRegulated <- unlist(strsplit(input$DownRegulated,"[ \t\r\n]"))
    background <- unlist(strsplit(input$Background,"[ \t\r\n]"))
    
    
    if (input$speciesSelectSig != "Human"){
      upRegulated <- convertGenes(upRegulated,tolower(input$speciesSelectSig))
      downRegulated <- convertGenes(downRegulated,tolower(input$speciesSelectSig))
      background <- convertGenes(background,tolower(input$speciesSelectSig))
    }
          
    sigJaccards <- sapply(seq_along(upSigs),signedJaccardSigCreate,upRegulated,downRegulated,upSigs,downSigs,background)
    
    chrDirs <- sapply(seq_along(chrDirsList),signedJaccardCreateChrDir,upRegulated,downRegulated,chrDirsList,background)
    
    
    
    sigJaccards <- data.frame(ID=names(upSigs),sigJaccard=sigJaccards,chrDir=chrDirs)
    sigJaccards$accessions <- gsub("(.*)_.*", "\\1",sigJaccards$ID)
    sigJaccards <- merge(sigJaccards,expTable,by.x="accessions",by.y="ID")
    sigJaccards <- sigJaccards[,-1]
    sigJaccards$Comparison <- accessions[ match(sigJaccards$ID,accessions$combined),"comparisonsText"]
   
    sigJaccards <- sigJaccards[ ,c(2,3,4,14,1,5:13)]
    sigJaccards <- sigJaccards[ order(sigJaccards[,1]),]
    
    sigJaccards <- tibble::add_column(sigJaccards, as.vector(scale(sigJaccards[,1])),.after = 2)
    sigJaccards <- tibble::add_column(sigJaccards, as.vector(scale(sigJaccards[,2])),.after = 3)
    
    shinyjs::enable("sigSearch")

    sigJaccards
    
    
  })
  
  #sim output
  output$sigSummary <- DT::renderDataTable({ 
    dt <- getJaccardSim()
 
    if (!all(na.omit(dt[,1])==0)){
      
      dt[,1] <- signif(dt[,1], digits = 3)
      dt[,2] <- signif(dt[,2], digits = 3)
      dt[,3] <- signif(dt[,3], digits = 3)
      dt[,4] <- signif(dt[,4], digits = 3)
      
      maxVal<-max(abs(c(dt[,1],dt[,2])),na.rm = T) + 0.001
      minVal<- -maxVal
      breaksList = seq(minVal, maxVal, by = 0.001)
      colfunc <- colorRampPalette(c("green", "white","red"))
      map <- makecmap(dt[,1],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
      colours<-  cmap(dt[,1], map = map)
      style.sigjacc <- styleEqual(dt[,1], colours)
      map <- makecmap(dt[,2],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
      colours<-  cmap(dt[,2], map = map)
      style.chrDir <- styleEqual(dt[,2], colours)
      
      maxVal<-max(abs(c(dt[,3],dt[,4])),na.rm = T) + 0.001
      minVal<- -maxVal
      breaksList = seq(minVal, maxVal, by = 0.001)
      colfunc <- colorRampPalette(c("green", "white","red"))
      map <- makecmap(dt[,3],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
      colours<-  cmap(dt[,3], map = map)
      style.sigjaccZscore <- styleEqual(dt[,3], colours)
      map <- makecmap(dt[,4],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
      colours<-  cmap(dt[,4], map = map)
      style.chrDirZscore <- styleEqual(dt[,4], colours)
      
      colnames(dt)[1:4] <- c("Signed Jaccard (sig)","Signed Jaccard (ChrDir)","Signed Jaccard (sig) zscore", "Signed Jaccard (ChrDir) zscore")
      
      dt<-datatable(
        dt,
        rownames = F,
        caption = tags$caption(
          style = 'caption-side: centre; color: black;text-align: centre;', tags$b("Search Summary")),
        extensions = c("Buttons",'Scroller'),
        options = list(
          scrollY = 300,
          scrollX = 500,
          scroller = F,
          pageLength = -1,
          lengthChange = T,
          bInfo = F,
          dom = 'frtB',
          lengthMenu = list( c(10, 20, -1), c(10, 20, "All")),
          order = list(0, 'desc'),
          buttons = c("copy", "csv"),referRender=TRUE
        ),
        selection = list(mode = "single"),
        filter = list(position = 'top', clear = FALSE, plain=TRUE))
        dt <- formatStyle(dt,1, backgroundColor = style.sigjacc)
        dt <- formatStyle(dt,2, backgroundColor = style.chrDir)
        dt <- formatStyle(dt,3, backgroundColor = style.sigjaccZscore)
        dt <- formatStyle(dt,4, backgroundColor = style.chrDirZscore)
      dt
    } else {
      dt<-as.data.frame("No matches found")
      colnames(dt)=""
      dt <-
        datatable(
          dt,
          rownames = F,
          options = list(bInfo = FALSE,dom = 'rti')) %>% formatStyle(columns = 1,color = 'red')
      
      dt
    }
    
  },server = T
  )
  
  #sig search gene overlap
  makeQueryResponseVenn = reactive ({
    
    rowID = input$sigSummary_row_last_clicked
    if (!is.null(rowID)){
      
      comparisonFC <- as.data.frame(getResponseComparisonFoldChangeTable())
     
      
      if(ncol(comparisonFC)==3 ){
        
      comparisonFCUp <- na.omit(comparisonFC[ comparisonFC[,2] >= log2(1.5),])
      comparisonFCDown <- na.omit(comparisonFC[ comparisonFC[,2] <= log2(1/1.5),])
      comparisonFCUp <- na.omit(comparisonFCUp[ comparisonFCUp[,3] <=0.05,])
      comparisonFCDown <- na.omit(comparisonFCDown[ comparisonFCDown[,3] <=0.05,])

      
      queryFCUp.genes <- unlist(strsplit(input$UpRegulated,"[ \t\r\n]"))
      queryFCDown.genes <- unlist(strsplit(input$DownRegulated,"[ \t\r\n]"))
      
      if (input$speciesSelectSig!="Human"){
        queryFCUp.genes <- convertGenes(queryFCUp.genes,tolower(input$speciesSelectSig))
        queryFCDown.genes <- convertGenes(queryFCDown.genes,tolower(input$speciesSelectSig))
      }
      comparisonID =  gsub("(.*)_.*", "\\1",getJaccardSim()[input$sigSummary_row_last_clicked,"ID"])
      comparisonSpecies <- expTable[expTable$ID==comparisonID,"Species"]
      comparisonSpecies <- tolower(pull(comparisonSpecies, Species))
      
      if (comparisonSpecies!="human"){
        comparisonFCUp.genes <-convertGenes(comparisonFCUp[,1],comparisonSpecies)
        comparisonFCDown.genes <-convertGenes(comparisonFCDown[,1],comparisonSpecies)
      } else{
        comparisonFCUp.genes <-comparisonFCUp[,1]
        comparisonFCDown.genes <-comparisonFCDown[,1]
      }

      
      upup <- paste(sort(intersect( queryFCUp.genes, comparisonFCUp.genes)),collapse=" ")
      downdown <- paste(sort(intersect( queryFCDown.genes, comparisonFCDown.genes)),collapse=" ")
      updown <- paste(sort(intersect( queryFCUp.genes, comparisonFCDown.genes)),collapse=" ")
      downup <- paste(sort(intersect( queryFCDown.genes ,comparisonFCUp.genes)),collapse=" ")

      venn <- data.frame(Overlap=c("Query:Up Comparison:Up","Query:Down Comparison:Down","Query:Up Comparison:Down","Query:Down Comparison:Up"), Genes=c(upup,downdown,updown,downup))
    
      return(venn=venn)
      }
    }
    
  })
  
  #sig search gene overlap
  makeQueryResponseChrDirVenn = reactive ({
    
    rowID = input$sigSummary_row_last_clicked
    if (!is.null(rowID)){
      
      comparisonFC <- as.data.frame(getResponseComparisonChrDirTable())
      
      comparisonFCUp <- na.omit(comparisonFC[ comparisonFC[,2] > 0,])
      comparisonFCDown <- na.omit(comparisonFC[ comparisonFC[,2] < 0,])
      
      queryFCUp.genes <- unlist(strsplit(input$UpRegulated,"[ \t\r\n]"))
      queryFCDown.genes <- unlist(strsplit(input$DownRegulated,"[ \t\r\n]"))
      
      if (input$speciesSelectSig!="Human"){
        queryFCUp.genes <- convertGenes(queryFCUp.genes,tolower(input$speciesSelectSig))
        queryFCDown.genes <- convertGenes(queryFCDown.genes,tolower(input$speciesSelectSig))
      }
      comparisonID =  gsub("(.*)_.*", "\\1",getJaccardSim()[input$sigSummary_row_last_clicked,"ID"])
      comparisonSpecies <- expTable[expTable$ID==comparisonID,"Species"]
      comparisonSpecies <- tolower(pull(comparisonSpecies, Species))
      
      if (comparisonSpecies!="human"){
        comparisonFCUp.genes <-convertGenes(comparisonFCUp[,1],comparisonSpecies)
        comparisonFCDown.genes <-convertGenes(comparisonFCDown[,1],comparisonSpecies)
      } else{
        comparisonFCUp.genes <-comparisonFCUp[,1]
        comparisonFCDown.genes <-comparisonFCDown[,1]
      }
      
      upup <- paste(sort(intersect( queryFCUp.genes, comparisonFCUp.genes)),collapse=" ")
      downdown <- paste(sort(intersect( queryFCDown.genes, comparisonFCDown.genes)),collapse=" ")
      updown <- paste(sort(intersect( queryFCUp.genes, comparisonFCDown.genes)),collapse=" ")
      downup <- paste(sort(intersect( queryFCDown.genes ,comparisonFCUp.genes)),collapse=" ")
      
      venn <- data.frame(Overlap=c("Query:Postive Comparison:Positive","Query:Negative Comparison:Negative","Query:Positive Comparison:Negative","Query:Negative Comparison:Postive"), Genes=c(upup,downdown,updown,downup))
      
     return(venn=venn)
      }
    })
  
  #output sig search gene overlap
  output$sigOverlap<- renderDataTable({ 
    
    venn <- makeQueryResponseVenn()
    if(!is.null(venn)){
      overlapTable <- datatable(venn,rownames = F,caption = tags$caption(
        style = 'caption-side: centre; color: black;text-align: centre;', tags$b("Differentially Expressed Genes")),
                                options = list(scrollY = "400px",
                                               lengthChange = FALSE,
                                               dom = 'rtiB',
                                               extensions = c("Buttons",'Scroller'),
                                               bInfo = FALSE,searching = FALSE, paging =FALSE ,buttons = c('copy', 'csv')),
                                selection = list(mode = "single"))
    }
  })
  
  #output sig search gene overlap
  output$sigChrDirOverlap<- renderDataTable({ 
    
    venn <- makeQueryResponseChrDirVenn()
    if(!is.null(venn)){
      overlapTable <- datatable(venn,rownames = F,caption=tags$caption(
        style = 'caption-side: centre; color: black;text-align: centre;', tags$b("Characteristic Direcion")),
                                options = list(scrollY = "400px",
                                               extensions = c("Buttons",'Scroller'),
                                               lengthChange = FALSE,
                                               dom = 'rtiB',
                                               bInfo = FALSE,searching = FALSE, paging =FALSE ,buttons = c('copy', 'csv')),
                                selection = list(mode = "single"))
    }
  })
  
  
  #when the comparisons table changes reset the selected rows
  v <- reactiveValues(comparisonRow = NULL, subnetworkRow=NULL, simSummaryRow=NULL,similaritySummaryMultiRow=NULL,thresholdPickerPval="_1.5_0.05",thesholdFC="_1.5")
  
  #reset the values when experiment changes
  observeEvent(input$expTable_row_last_clicked, {
    v$comparisonRow <- NULL
    v$subnetworkRow <- NULL
    v$simSummaryRow <- NULL
    v$similaritySummaryMultiRow <- NULL
    updateTabsetPanel(session, "activesubnet",selected = "Summary Table")
    updateTabsetPanel(session, "sharedResponse",selected = "Summary")
    updateTabsetPanel(session, "sharedResponseMulti",selected = "SummaryMulti")
    
    
  })
  
  #set values on click
  observeEvent( input$rows,{
    v$comparisonRow <- input$comparisons_row_last_clicked
  })
  
  observeEvent(input$subRows, {
    v$subnetworkRow <- input$subnetworkTable_row_last_clicked
  })
  
  observeEvent(input$simSummaryRows, {
    v$simSummaryRow <- input$similaritySummary_rows_selected
  })
  
  observeEvent(input$similaritySummaryMultiRow, {
    v$similaritySummaryMultiRow <- input$similaritySummaryMulti_rows_selected
  })
  
  observeEvent(input$threshold, {
    if(expTable[input$expTable_row_last_clicked,"foldChangeOnly"]==FALSE){
      v$thresholdPickerPval <- input$threshold
    } else {
      v$thresholdPickerFC <- input$threshold
    }
    
  })
  
  #geneSearch example button
  observeEvent(input$geneSearchExample, {
    updateTextInput(session,inputId = "GeneName",value = "DDIT3")
  })
  
  #sig example button
  observeEvent(input$loadSigExample, {
    updateTextAreaInput(session,inputId = "UpRegulated",value = "STMN2\nABCB5\nTHBS4\nMMP13\nC21orf37\nEGFL6\nCOL16A1\nGPR158\nAK5\nAK5\nATP6V0D2\nALU2\nDIRAS2\nLOC101929504\nDNASE2B\nSTXBP5L\nLHX2\nPSIP1\nNYAP2\nCTNND2\nFNDC1\nEN1\nMIR31HG\nXLOC_006820\nFRMPD4\nLOC101929450\nFGFR2\nLPAR3\nCOL12A1\nLOC613266\nGLDC\nESYT3\nKIAA1549L\nNSG1\nLCTL\nTRIM67\nPPEF1\nDPP4\nLRRC8E\nIRX1\nCTHRC1\nTM4SF19\nESYT3\nKPNA7\nFRAS1\nCALCR\nCALHM3\nTOX3\nFAP\nDGKI\nPTPRD\nIGDCC4\nGLT8D2\nIGSF3\nC11orf87\nESPNL\nHEY1\nDCSTAMP\nGALNT5\nLOC101927619\nSATB2\nGJA1\nSLC9B2\nPLEKHG4B\nLINC00673\nHMSD\nGPR68\nCALCR\nPRSS12\nTHBS2\nSP6\nPLEKHA5\nUST\nSGMS2\nMMP9\nMYO1B\nFKBP7\nIGSF3\nSHC3\nXLOC_005452\nSPNS2\nCASC14\nLINC00605\nGNPTAB\nARMC4\nMMP11\nLAMP3\nLOC100132705\nPRMT8\nSRPX\nST8SIA1\nFERMT1\nCUBN\nTDRD3\nPPIC\nGFRA1\nCLEC4A\nRUNX2\nMAGED4B\nSRD5A1\nVAV2\nNKD2\nSCD5\nPRSS8\nFRMD6\nIL21R\nLINC01057\nTP53I3\nCD276\nARRDC4\nTPBG\nRBFOX2\nTRAF3IP2-AS1\nZNF815P\nTPMT\nATP6AP2\nAARS")
    updateTextAreaInput(session,inputId = "DownRegulated",value = "RAPGEF1\nWDR74\nPRMT5\nPRKD2\nWASF2\nLOC102725378\nTSPAN14\nPEA15\nTNIP1\nMFNG\nRCSD1\nMGAT1\nCNPPD1\nPPP5C\nPECAM1\nABLIM3\nHSPB1\nPEA15\nSAMD1\nTM4SF1\nSYNPO\nPXN\nDSP\nFLT4\nAFAP1L1\nABCG1\nC10orf54\nDACH1\nNR1H2\nTSPAN14\nIQCA1\nFAM107A\nSLPI\nXLOC_010730\nSLC16A13\nLOC100128343\nIQCA1\nFBXO31\nGRAP\nTP53TG3C\nLYVE1\nPRR26\nKCNQ3\nXLOC_001496\nTNFRSF1B\nINHBB\nALU1\nC6orf25\nAQP10\nFFAR2\nTRBV28\nHBA2\nXLOC_014512\nHBA2\nCLEC1B\nPF4V1\nMMP25\nLOC102724484\nSMIM24\nTUBB1\nHSJ1167H4\nTTTY16\nLINC00570\nCMTM2\nCXCR2\nHBD\nGYPB\nSEC14L3\nFLJ46249\nCLEC1B\nCLEC1B\nGATA1\nPDZK1IP1\nNFE2\nXLOC_013489\nPF4\nXLOC_000346\nALAS2\nS100P\nHBQ1\nPPBP\nHBA2\nHBM\nHEMGN\nS100A12\nS100A12\nHBG1")
  })
  
  plotReady <- reactiveValues(ok = TRUE,chrDir=TRUE,sig=TRUE,chrDirCols=1,sigCols=1)
  
  #call condition class
  condition <- function(subclass, message, call = sys.call(-1), ...) {
    structure(
      class = c(subclass, "condition"),
      list(message = message, call = call, ...)
    )
  }
  
  tooLargeInputError <- function() {
    msg <- paste0("Too large input")
    condition(c("tooLargeInputError", "error"),
              message = msg, 
              text = text
    )
  }
  
  
  
}
