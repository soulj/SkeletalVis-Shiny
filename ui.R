library(shiny)
library(visNetwork)
library(plotly)
library(shinyBS)

navbarPage(
  title = "SkeletalVis",
  theme = "bootstrap.min.css",
  
  tabPanel(
    title = "Home",
    icon = icon("home"),
    div(
      style = "width:800px; margin:0 auto;",
      div(
        class = "jumbotron",
        h2(
          "Welcome to ",
          tags$b("SkeletalVis"),
          ": a user friendly web application for exploration of skeletal disease related expression datasets."
        )
      ),
      div(
        class = "jumbotron",
        h2("How does it work?"),
        h3(strong("SkeletalVis"), " is simple to use."),
        h3(
          icon("search"),
          strong("Explore"),
          " and look for similarities between existing gene expression datasets"
        ),
        h3(
          icon("cogs"),
          strong("Compare"),
          " by gene across all the expression datasets or enter your own differentially expressed genes to search for similar experiments"
        )
      ),
      div(
        class = "jumbotron",
        h2("What data does SkeletalVis contain?"),
        h3(
          "SkeletalVis uses public microarray and RNA-Seq expression data from ArrayExpress and GEO"
        ),
        h3(
          "These repositories are searched for skeletal disease experiments which are then consistently analysed"
        )
      ),
      div(
        class = "jumbotron",
        h3("Created by the Pathways and biological systems modelling group"),
        h3(
          "Funded by the European Community’s Seventh Framework Programme grant (602300)"
        ),
        img(
          src = "images/sybil.png",
          max.width = "20",
          width = "20%",
          height = "auto",
          style = "display: block;
          margin-left: auto;
          margin-right: auto"
        )
      )
    )
    
  ),
  
  #Exploration Panel
  tabPanel(
    title = "Explore",
    tags$style(
      type = 'text/css',
      ".navbar-nav {font-size: 20px};.main-header .logo {font-size: 22px}; "
    ),
    icon = icon("search"),
    
    
    
    sidebarLayout(
      position = "right",
      sidebarPanel(
        h4("Experiment"),
        uiOutput("AccessionText"),
        h4("Comparison"),
        uiOutput("ComparisonText"),
        uiOutput("reminderText"),
        uiOutput("subnetworkInfo"),
        uiOutput("subNetReminderText"),
        uiOutput("compareInfo"),
        uiOutput("sharedResponseReminderText"),
        width = 2
        
      ),
      
      mainPanel(
        width = 10,
        tags$style(".popover{max-width: 100%;}"),
        tabsetPanel(
          id = "toptabset",
          tabPanel(
            "Comparisons",
            value = "comparisons",
            icon = icon("table"),
            tabsetPanel(
              id = "tabset",
              tabPanel(
                "Choose Data",
                div(h3("Step 1: Select an Experiment"), style = "display:inline-block"),
                uiOutput("expHelp", style = "display:inline-block"),
                bsPopover(
                  "expHelp",
                  title = "Select an Experiment",
                  content = "The table can be explored with the search box and sorted by clicking the columns. <br> Click a row to select the experiment and view the comparisons <br> Then select a comparison to load the data",
                  placement = "bottom",
                  trigger = "hover",
                  options = NULL
                ),
                conditionalPanel("$('#expTable').hasClass('recalculating')",
                                 tags$div(h3(
                                   'Loading the data ... '
                                 ))),
                DT::dataTableOutput("expTable"),
                shiny::HTML("<h3>Step 2: Choose a Comparison</h3>"),
                DT::dataTableOutput("comparisons"),
                icon = icon("table")
              ),
              tabPanel(
                "Fold Changes",
                div(h3("Differential Expression"), style = "display:inline-block"),
                uiOutput("fcHelp", style = "display:inline-block"),
                bsPopover(
                  "fcHelp",
                  title = "Differential Expression Analysis",
                  content = "Table shows the differential expression analysis for the selected experiment and comparison. <br> Table can be searched and sorted",
                  placement = "bottom",
                  trigger = "hover",
                  options = NULL
                ),
                DT::dataTableOutput("foldchanges"),
                icon = icon("table")
              ),
              tabPanel(
                "Chr Dir",
                div(h3("Characteristic Direction"), style = "display:inline-block"),
                uiOutput("chrDirHelp", style = "display:inline-block"),
                bsPopover(
                  "chrDirHelp",
                  title = "Characteristic Direction",
                  content = "Table shows the Characteristic Direction analysis for the selected experiment and comparison. <br> Table can be searched and sorted",
                  placement = "bottom",
                  trigger = "hover",
                  options = NULL
                ),
                DT::dataTableOutput("chrDir"),
                icon = icon("table")
              ),
              
              tabPanel(
                "Pathways",
                div(h3("Differentially Regulated Pathways"), style = "display:inline-block"),
                uiOutput("pathwayHelp", style = "display:inline-block"),
                bsPopover(
                  "pathwayHelp",
                  title = "Differentially regulated pathways",
                  content = "Table shows the significant (adjusted p-value ≤ 0.05) pathways. <br> The genes in each pathway, the significance and the percentage enrichment are shown. <br> Table can be searched and sorted",
                  placement = "bottom",
                  trigger = "hover",
                  options = NULL
                ),
                DT::dataTableOutput("pathways"),
                icon =
                  icon("table")
              ),
              tabPanel(
                "GO Enrichment",
                value = "go",
                icon = icon("table"),
                tabsetPanel(
                  tabPanel(
                    "GO Enrichment",
                    div(h3("Differentially Enriched GO Terms"), style = "display:inline-block"),
                    uiOutput("GOHelp", style = "display:inline-block"),
                    bsPopover(
                      "GOHelp",
                      title = "Differentially enriched GO terms",
                      content = "Table shows the significant (adjusted p-value ≤ 0.05) GO Terms. <br> The differentially expressed genes corresponding to each enriched GOTerm can be viewed by selecting the + symbol on the left of the table <br> The statistical significance and the percentage enrichment are shown <br> Table can be searched and sorted",
                      placement = "bottom",
                      trigger = "hover",
                      options = NULL
                    ),
                    DT::dataTableOutput("goTable")
                  ),
                  tabPanel(
                    "Reduced GO Enrichment",
                    div(h3("Reduced Differentially Enriched GO Terms"), style = "display:inline-block"),
                    uiOutput("GOslimHelp", style = "display:inline-block"),
                    bsPopover(
                      "GOslimHelp",
                      title = "Reduced differentially enriched GO terms",
                      content = "Table shows the significant (adjusted p-value ≤ 0.05) GO Terms with redundancy reduction <br> The differentially expressed genes corresponding to each enriched GOTerm can be viewed by selecting the + symbol on the left of the table <br> The significance and the percentage enrichment are shown. <br> Table can be searched and sorted",
                      placement = "bottom",
                      trigger = "hover",
                      options = NULL
                    ),
                    DT::dataTableOutput("goReducedTable")
                  ),
                  tabPanel(
                    "GO MDS",
                    div(h3("Reduced GO Term MDS"), style = "display:inline-block"),
                    uiOutput("GOMDSHelp", style = "display:inline-block"),
                    bsPopover(
                      "GOMDSHelp",
                      title = "Differentially regulated pathways",
                      content = "Plot shows the MDS of the differentially enriched GO Terms using semantic similarity to group similar terms <br>The colour corresponds to the GO Term enrichment p-value",
                      placement = "bottom",
                      trigger = "hover",
                      options = NULL
                    ),
                    htmlOutput("goMDS")
                  )
                )
              ),
              tabPanel(
                "TF Enrichment",
                div(h3("Transcription Factor Enrichment"), style = "display:inline-block"),
                uiOutput("TFHelp", style = "display:inline-block"),
                bsPopover(
                  "TFHelp",
                  title = "Differentially regulated pathways",
                  content = "Table shows the significantly enriched (adjusted p-value ≤ 0.05) TFs <br> The differentially expressed genes regulated by each TF can be viewed by selecting the + symbol on the left of the table. <br> Table can be searched and sorted for a target gene/TF of interest",
                  placement = "bottom",
                  trigger = "hover",
                  options = NULL
                ),
                DT::dataTableOutput("TFs"),
                icon = icon("table")
              ),
              tabPanel(
                "Drug Enrichment",
                value = "drug",
                icon = icon("table"),
                tabsetPanel(
                  tabPanel(
                    "Mimic",
                    div(h3("Mimic Drug Enrichment"), style = "display:inline-block"),
                    uiOutput("mimicDrugHelp", style = "display:inline-block"),
                    bsPopover(
                      "mimicDrugHelp",
                      title = "Mimic Drugs",
                      content = "Table shows the drugs that give the most similar transcriptomic response <br> The drug name and predicted drug targets, differentially expressed genes that overlap with the drug pertubation, the percentage overlap are shown <br> Table can be searched and sorted for a target gene/TF of interest",
                      placement = "bottom",
                      trigger = "hover",
                      options = NULL
                    ),
                    DT::dataTableOutput("mimicTable")
                  ),
                  tabPanel(
                    "Reverse",
                    div(h3("Reverse Drug Enrichment"), style = "display:inline-block"),
                    uiOutput("reverseDrugHelp", style = "display:inline-block"),
                    bsPopover(
                      "reverseDrugHelp",
                      title = "Reverse Drugs",
                      content = "Table shows the drugs that give the opposite transcriptomic response i.e could be used to reverse the experimental pertubation <br> The drug name and predicted drug targets, differentially expressed genes that overlap with the drug pertubation, the percentage overlap are shown <br> Table can be searched and sorted for a target gene/TF of interest",
                      placement = "bottom",
                      trigger = "hover",
                      options = NULL
                    ),
                    DT::dataTableOutput("reverseTable")
                  )
                )
              ),
              tabPanel(
                "Active Sub-networks",
                value = "subnet",
                icon = icon("dot-circle-o"),
                tabsetPanel(
                  id = "activesubnet",
                  tabPanel(
                    "Summary Table",
                    div(h3("Active Subnetwork Summary"), style = "display:inline-block"),
                    uiOutput("networkHelp", style = "display:inline-block"),
                    bsPopover(
                      "networkHelp",
                      title = "Active Subnetwork",
                      content = "Table shows the significant active subnetworks and the top GO enrichment for that network <br> Click on a row to select that network and view it in the subnetwork tab",
                      placement = "bottom",
                      trigger = "hover",
                      options = NULL
                    ),
                    DT::dataTableOutput("subnetworkTable")
                  ),
                  tabPanel(
                    "Subnetwork",
                    div(h3("Subnetwork Visualisation"), style = "display:inline-block"),
                    uiOutput("subnetworkHelp", style = "display:inline-block"),
                    bsPopover(
                      "subnetworkHelp",
                      title = "Differentially regulated subnetwork",
                      content = "View the subnetwork. Nodes can be selected and dragged. Colour corresponds to the fold-change. Hover over a node to view the fold change",
                      placement = "bottom",
                      trigger = "hover",
                      options = NULL
                    ),
                    visNetworkOutput("subnetwork", width = "100%", height = "650px")
                  )
                )
              ),
              tabPanel(
                "Shared Response",
                icon = icon("table"),
                value = "response",
                tabsetPanel(
                  id = "sharedResponse",
                  tabPanel(
                    "Summary",
                    div(h3("Shared Response Summary"), style = "display:inline-block"),
                    uiOutput("responseHelp", style = "display:inline-block"),
                    bsPopover(
                      "responseHelp",
                      title = "Response summary",
                      content = "Table shows the similarity of the loaded experiment to all the other experiments <br> Click on a row to select that experiment and view the overlapping genes in the comparison tab",
                      placement = "bottom",
                      trigger = "hover",
                      options = NULL
                    ),
                    DT::dataTableOutput("similaritySummary")
                  ),
                  tabPanel(
                    "Gene Overlap",
                    h4("FoldChange Overlap"),
                    DT::dataTableOutput("geneOverlap"),
                    h4("Significant Overlap"),
                    DT::dataTableOutput("geneOverlapSig"),
                    h4("Charecteristic Direction Overlap"),
                    DT::dataTableOutput("geneOverlapChrDir")
                  ),
                  tabPanel("Cosine Plot", plotOutput("cosineSim")),
                  tabPanel("Signed Jaccard Plot",  plotOutput("jaccardSim"))
                )
              )
              
              
              
            )
          ),
          tabPanel("Quality Control", htmlOutput("QC"), icon =
                     icon("check-square-o")),
          tabPanel("PCA",  imageOutput("PCA"), icon = icon("line-chart"))
          
        )
        
        
      )
    )
  ),
  #Compare panel
  tabPanel(
    title = "Compare",
    icon = icon("cogs"),
    sidebarLayout(
      position = "left",
      sidebarPanel(
        uiOutput("geneSearchUI"),
        popify(
          uiOutput("geneSearchHelp"),
          title = "Compare by Gene",
          content = "Enter a gene and select the species to search the fold change tables"
        ),
        uiOutput("sigSearchUI"),
        popify(
          uiOutput("sigSearchHelp"),
          title = "Compare by Signature",
          content = "Enter up and down regulated genes to compare your transcriptomic signature to the database. Optionally enter the background of all expressed genes in your experiment to improve the accuracy of the similarity calculations"
        )
      ),
      mainPanel(tabsetPanel(
        id = "search",
        tabPanel(
          title = "Compare by gene",
          value = "gene",
          icon = icon("table"),
          conditionalPanel(
            "$('#geneFoldChangeTable').hasClass('recalculating')",
            tags$div(h2('Searching ... '))
          ),
          DT::dataTableOutput("geneFoldChangeTable")
        ),
        tabPanel(
          title = "Compare by signature",
          value = "sig",
          icon = icon("table"),
          tabsetPanel(tabPanel(
            title="Search Results", icon = icon("table"),
          conditionalPanel("$('#sigSummary').hasClass('recalculating')",
                           tags$div(h2('Searching ... '))),
          DT::dataTableOutput("sigSummary"),
          DT::dataTableOutput("sigOverlap"),
          DT::dataTableOutput("sigChrDirOverlap")
        ),
        tabPanel("Signed Jaccard Plot", plotOutput("sigJaccZscore")),
        tabPanel("Characteristic Direction Plot",  plotOutput("chrDirZscore"))
        ))
      ))
    )
  ),
  tabPanel(
    tags$style(type = 'text/css', 'body { overflow-y: scroll; }'),
    title = "Help",
    icon = icon("question"),
    sidebarLayout(
      sidebarPanel(
        style = "position: fixed; width=200px; overflow: visible;white-space: nowrap;",
        HTML("<a href='#ToolTips'><h3>Display tooltips</h3></a>"),
        HTML("<a href='#Explore'><h3>Explore</h3></a>"),
        HTML("<a href='#Step1'><h4>Search and choose data</h4></a>"),
        HTML(
          "<a href='#Step2'><h4>View the experiment metadata and QC</h4></a>"
        ),
        HTML(
          "<a href='#Step3'><h4>View the expression response and downstream analysis</h4></a>"
        ),
        HTML("<a href='#DiffExp'><h4>Differential Expression</h4></a>"),
        HTML(
          "<a href='#charDirect'><h4>Characteristic Direction</h4></a>"
        ),
        HTML("<a href='#DiffPath'><h4>Differential Pathways</h4></a>"),
        HTML(
          "<a href='#DiffTF'><h4>Differential Transcription Factors</h4></a>"
        ),
        HTML("<a href='#DiffGO'><h4>Enriched GO Terms</h4></a>"),
        HTML("<a href='#DiffDrug'><h4>Enriched Drugs</h4></a>"),
        HTML("<a href='#Subnet'><h4>Active Subnetworks</h4></a>"),
        HTML("<a href='#Shared'><h4>Shared Responses</h4></a>"),
        HTML("<a href='#Compare'><h3>Compare</h3></a>"),
        HTML("<a href='#CompareGene'><h4>Compare by Gene</h4></a>"),
        HTML("<a href='#CompareSig'><h4>Compare by Signature</h4></a>")
      ),
      mainPanel(
        h1("Display tooltips", align = "center", id = "ToolTips"),
        h4(
          "Hover over the",
          icon("fas fa-question-circle"),
          "icons to display tool tips in the app",
          align = "left"
        ),
        img(src = "images/tooltips.png", style =
              'border:1px solid #000000'),
        h1("Explore", align = "center", id =
             "Explore"),
        h4(
          "The explore section of the app allows browsing of pre-analysed transcriptomics data and comparison between these datasets.\n"
        ),
        img(src = "images/explore.png"),
        
        h2("Step 1: Search and choose data", align = "left", id =
             "Step1") ,
        h4(
          "The experiment table shows the analysed datasets available to explore.\n"
        ),
        img(
          src = "images/expTable.png",
          style = 'border:1px solid #000000',
          max.width = "100%",
          width = "100%",
          height = "auto"
        ),
        h4(
          "To help find an experiment of interest the whole table can be searched or a column can be sorted/searched."
        ),
        img(
          src = "images/expTableSearched.png",
          style = 'border:1px solid #000000',
          max.width = "100%",
          width = "100%",
          height = "auto"
        ),
        h4("Click on a row in the experiment table to select that experiment"),
        h4(
          "The comparisons in the experiment will be displayed in the comparisons table"
        ),
        h4(
          "A comparison to view can be chosen by selecting a row in the comparison table."
        ),
        img(
          src = "images/comparisonsTable",
          style = 'border:1px solid #000000',
          max.width = "100%",
          width = "100%",
          height = "auto"
        ),
        h4(
          "The selected experiment and comparison will be displayed in the left hand information panel."
        ),
        img(
          src = "images/infoPanel",
          style = 'border:1px solid #000000',
          max.width = "30%",
          width = "50%",
          height = "30%",
          align = "centre"
        ),
        h4("The selected experiment and comparison data will be now be loaded"),
        
        h2(
          "Step 2: View the experiment metadata, quality control and PCA plots",
          align = "left",
          id = "Step2"
        ),
        h4(
          "The QC data and the PCA plot after normalisation and batch effect correction can be viewed by clicking the corresponding tabs."
        ),
        h4("The QC allows assessment of the data quality"),
        img(
          src = "images/qc.png",
          style = 'border:1px solid #000000',
          max.width = "50%",
          width = "50%",
          height = "auto"
        ),
        h4(
          "The PCA shows the seperation of the samples by the top two components explaining the variation after any batch effect correction"
        ),
        img(
          src = "images/pca.png",
          style = 'border:1px solid #000000',
          max.width = "50%",
          width = "50%",
          height = "auto"
        ),
        
        h2(
          "Step 3: View the expression response and downstream analysis",
          align = "left",
          id = "Step3"
        ),
        h4(
          "Select selected expression response can be explored by examining the fold-change, pathway, gene ontology enrichments, transcription factor predictions and enriched drugs."
        ),
        h3("Differential Expression", align = "left", id =
             "DiffExp"),
        h4(
          "The differential expression table shows the gene names, log2 fold changes and if there are experimental replicates the Benjamini-Hochberg adjusted p-values."
        ),
        h4(
          "Differential expression analysis is performed using DESeq2 for RNA-seq and limma for microarray."
        ),
        h4(
          "The table can be searched using the search box, sorted by clicking the columns and filtered using the column filters."
        ),
        h4(
          "The copy and csv buttons at the bottom of the table allow export of visible filtered table"
        ),
        img(
          src = "images/fc.png",
          style = 'border:1px solid #000000',
          max.width = "100%",
          width = "100%",
          height = "auto"
        ),
        h3("Characteristic Direction", align = "left", id =
             "charDirect"),
        h4("The characteristic direction method gives the genes that best seperate the samples in the comparison and has been shown to be more sensitive than limma/DESeq2 in identifying differentially expressed genes. The top 500 genes are shown which can be considered a differential expression signature of the expression response. Please see"),
        tags$a(href="http://doi.org/10.1186/1471-2105-15-79","http://doi.org/10.1186/1471-2105-15-79",target="_blank"),
        h4("for details of the method"),
        h4(
          "The table can be searched using the search box, sorted by clicking the columns and filtered using the column filters."
        ),
        img(
          src = "images/chrDir.png",
          style = 'border:1px solid #000000',
          max.width = "50%",
          width = "50%",
          height = "auto"
        ),
        h3("Differential Pathways", align = "left", id =
             "DiffPath"),
        h4(
          "The pathway table shows the differentially regulated pathways, the differentially expressed in those pathways, the adjusted p-value and the percentage pathway coverage"
        ),
        h4(
          "The pathway enrichment is performed differentially expressed genes with a threshold of 1.5 fold change and an adjusted pvalue ≤0.05 (if applicable)"
        ),
        h4(
          "The table can be searched using the search box, sorted by clicking the columns and filtered using the column filters."
        ),
        h4(
          "The copy and csv buttons at the bottom of the table allow export of visible filtered table"
        ),
        img(
          src = "images/pathways.png",
          style = 'border:1px solid #000000',
          max.width = "100%",
          width = "100%",
          height = "auto"
        ),
        h3("Enriched GO Terms", align = "left", id =
             "DiffGO"),
        h4(
          "The GO Term enrichment table shows the differentially enriched Gene Ontology Biological Process terms, the adjusted p-value and the percentage of differentially expressed genes out of all genes annotated to that GO Term"
        ),
        h4(
          "Click on the left hand + symbol to view the differentially expressed genes annotated to that GO Term"
        ),
        h4(
          "The GO enrichment is performed differentially expressed genes with a threshold of 1.5 fold change and an adjusted pvalue ≤0.05 (if applicable)"
        ),
        h4(
          "The table can be searched using the search box, sorted by clicking the columns and filtered using the column filters."
        ),
        h4(
          "The copy and csv buttons at the bottom of the table allow export of visible filtered table"
        ),
        h4(
          "The Reduced GO Term enrichment table similarly shows GO Enrichment results after redundancy reduction based on the semantic similarity of the GO Terms"
        ),
        img(
          src = "images/go.png",
          style = 'border:1px solid #000000',
          max.width = "100%",
          width = "100%",
          height = "auto"
        ),
        h4(
          "The GO MDS plots shows an overview of the reduced GO Terms, seperated based on their semantic simiarity so similar terms are grouped together"
        ),
        h4("Hover over a point to show the name of the GO term"),
        h4("The points are coloured by adjusted p-value"),
        img(
          src = "images/gomds.png",
          style = 'border:1px solid #000000',
          max.width = "50%",
          width = "50%",
          height = "auto"
        ),
        h3(
          "Enriched Transcription Factors",
          align = "left",
          id = "DiffTF"
        ),
        h4(
          "The transcription factor (TF) enrichment table shows the transcription factor motifs predicted to significantly regulate the differentially expressed genes using RcisTarget"
        ),
        h4("The NES gives the significance of the enrichment for each motif"),
        h4(
          "The known TF corresponding to motifs are given in the TF_direct column"
        ),
        h4(
          "Inferred TF annotations through motif simialrity are given in the TF_indirect column"
        ),
        h4(
          "Click on the left hand + symbol to view the differentially expressed genes predicted to be regulated by that TF"
        ),
        h4(
          "The TF enrichment is performed using differentially expressed genes with a threshold of 1.5 fold change and an adjusted pvalue ≤0.05 (if applicable)"
        ),
        img(
          src = "images/TF.png",
          style = 'border:1px solid #000000',
          max.width = "100%",
          width = "100%",
          height = "auto"
        ),
        h3("Enriched Drugs", align = "left", id =
             "DiffDrug"),
        h4(
          "The drug enrichment table shows the drugs with overlapping transcriptomic signatures from the L1000CDS2 tool"
        ),
        h4(
          "The drug enrichment is performed using differentially expressed genes with a threshold of 1.5 fold change and an adjusted pvalue ≤0.05 (if applicable)"
        ),
        h4(
          "The search score is the overlap between the input and drug signature differentially expressed genes divided by the total number of genes"
        ),
        
        h4(
          "Drugs can either mimic the transcriptional response or have the opposite response (reverse)"
        ),
        h4(
          "Click on the left hand + symbol to view the overlapping differentially expressed genes and the predicted drug targets"
        ),
        img(
          src = "images/drug.png",
          style = 'border:1px solid #000000',
          max.width = "100%",
          width = "100%",
          height = "auto"
        ),
        h3("Active Subnetworks", align = "left", id =
             "Subnet"),
        h4(
          "The subnetwork table shows a summary of the significant subnetworks (de novo pathways) identified using the GIGA algorithm"
        ),
        h4(
          "This algorithm uses all the differential expression data without a threshold"
        ),
        h4(
          "The signifcance of the pathway and the top enriched GO term are shown"
        ),
        img(
          src = "images/subnetwork.png",
          style = 'border:1px solid #000000',
          max.width = "100%",
          width = "100%",
          height = "auto"
        ),
        h4(
          "Click on a row to load that subnetwork, the information will be shown and it can now be viewed in the subnetwork tab"
        ),
        img(
          src = "images/subnetworkinfo.png",
          style = 'border:1px solid #000000',
          max.width = "30%",
          width = "30%",
          height = "auto"
        ),
        h4(
          "The subnetwork can be zoomed, moved around the canvas and each node moved to improve the layout"
        ),
        h4("Hover over a node to view the fold change of that node"),
        img(
          src = "images/subnetworkvis.png",
          style = 'border:1px solid #000000',
          max.width = "50%",
          width = "50%",
          height = "auto"
        ),
        h3("Shared Response", align = "left", id =
             "Shared"),
        h4(
          "The shared response summary allows comparison of the selected dataset against the other datasets in the database"
        ),
        h4(
          "The cosine similarity is a measure of correlation of the fold-changes"
        ),
        h4(
          "Signed jaccard coefficient is the overlap of up and down differentially expressed genes in the query and comparison datasets"
        ),
        h4(
          "Differentially expressed genes are defined using either just a fold-change threshold (1.5 fold) or using both fold-change (1.5 fold) and statistical significance (FDR ≤ 0.05"
        ),
        h4(
          "Select a row to choose a comparison dataset for further detail in the comparison tab"
        ),
        img(
          src = "images/sharedResponse.png",
          style = 'border:1px solid #000000',
          max.width = "50%",
          width = "50%",
          height = "auto"
        ),
        h4(
          "This will show the up and down differentially expressed genes in the query dataset that overlap the selected comparison"
        ),
        img(
          src = "images/sharedGenes.png",
          style = 'border:1px solid #000000',
          max.width = "40%",
          width = "40%",
          height = "auto"
        ),
        h4(
          "The cosine and jaccard plots show the distribution of scores over the database"
        ),
        
        h1("Compare", align = "center", id =
             "Compare"),
        h4(
          "The compare section of the app allows searching of the database by a gene or by a transcriptomics signature"
        ),
        img(src = "images/compare.png"),
        h2("Compare by Gene", align = "left", id =
             "CompareGene"),
        h4(
          "To search all the fold change tables enter a gene symbol and select the appropriate species"
        ),
        img(
          src = "images/searchGene.png",
          style = 'border:1px solid #000000',
          max.width = "25",
          width = "25%",
          height = "auto"
        ),
        h4(
          "The fold change of that gene in every comparison in all the experiments will be shown"
        ),
        h4("Datasets of interest can then be fully explored in the Explore tab"),
        img(
          src = "images/compareGene.png",
          style = 'border:1px solid #000000',
          max.width = "100%",
          width = "100%",
          height = "auto"
        ),
        h2("Compare by signature", align = "left", id =
             "CompareSig"),
        h4(
          "To search all the fold change tables enter the differentially expressed up and down regulated genes and select the appropriate species"
        ),
        img(
          src = "images/inputSig.png",
          style = 'border:1px solid #000000',
          max.width = "30%",
          width = "30%",
          height = "auto"
        ),
        h4(
          "The similarity to every other comparison in all the experiments will be shown using the signed jaccard measure for significant differentially expressed genes and characterstic direction signatures"
        ),
        h4("Select a row to see the gene overlap for that comparison"),
        img(
          src = "images/sigTable.png",
          style = 'border:1px solid #000000',
          max.width = "100%",
          width = "100%",
          height = "auto"
        )
      )
    )
  )
  
  
)
