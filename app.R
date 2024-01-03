rm(list = ls())

list.of.packages <- c('shiny', 'rmarkdown', 'seqinr', 'shinydashboard', 'dashboardthemes', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable', 'rgl',
                      'DT', 'htmltools', 'biomformat', 'phangorn', 'bios2mds', 'zip', 'dplyr', 'forestplot', 'fossil', 'picante',
                      'entropart', 'dirmult', 'robustbase', 'erer', 'BiasedUrn', 'CompQuadForm', 'robCompositions', 
                      'GUniFrac', 'ecodist', 'MiRKAT', 'gridExtra', 'ggplot2', 'patchwork', 'ggthemes', 'DiagrammeR', 'stringr', 
                      'devtools', 'reticulate', 'remotes', 'gridGraphics', 'compositions', 'ICSNP', 'xtable', 'BiocManager', 'PMCMRplus', 'vegan', 'bda', 'mediation', 'survival', 'survminer', 'coin', 'dplyr')

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(!require('fBasics')) install.packages("fBasics", type = "binary")
if(!require('BiocParallel')) BiocManager::install("BiocParallel") 
if(!require('phyloseq')) remotes::install_github('joey711/phyloseq')
if(!require('biomformat')) remotes::install_github('joey711/biomformat') 
if(!require('DACT')) remotes::install_github("https://github.com/zhonghualiu/DACT")

library(seqinr)
library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(tidyverse)
library(phyloseq)
library(plotly)
library(shinyWidgets)
library(shinyjs)
library(googleVis)
library(xtable)
library(DT)
library(DiagrammeR)
library(htmltools)
library(biomformat)
library(phangorn)
library(bios2mds)
library(zip)
library(fresh)
library(bda)
library(devtools)
library(mediation)
library(DACT)
library(survival)
library(survminer)
library(coin)
library(dplyr)

source("Source/MiMedSurv.Non.Med.R")
source("Source/MiMedSurv.Data.Upload.R")
source("Source/MiMedSurv.Alpha.Mediation.R")
source("Source/MiMedSurv.Beta.Mediation.R")
source("Source/MiMedSurv.Taxa.Mediation.R")

# Comments ---------------------------

{
  TITLE = p(strong("MiMedSurv: A Unified Cloud Platform for Microbiome Causal Mediation Analysis with Survival Responses", style = "font-size:18pt"))
  HOME_COMMENT = p(strong("MiMedSurv", style = "font-size:13pt"), "(microbiome mediation analysis with survival responses) is a unified web cloud computing platform for comprehensive microbiome mediation analysis with survival (i.e., time-to-event) responses. MiMedSurv surveys the roles of the microbiome as a mediator between a treatment (e.g., medical intervention, environmental exposure) and a survival response on the host’s health or disease (e.g., time-to-disease, time-to-cure). The main features of MiMedSurv are as follows. First, MiMedSurv conducts basic exploratory non-mediational survival analysis, not involving microbiome, to survey the disparity in survival time between medical treatments (e.g., treatment vs. placebo, new treatment vs. old treatment) / environmental exposures (e.g., rural vs. urban, smoking vs. non-smoking). (see Non-Mediational Analysis). Second, MiMedSurv identifies the mediating roles of the microbiome in various aspects (see Mediational Analysis): (i) as a microbial ecosystem using ecological measures (e.g., alpha- and beta-diversity indices) (see Community-level Analysis) and (ii) as individual microbial taxa in various hierarchies (e.g., phyla, classes, orders, families, genera, species) (see Taxonomy-level Analysis). We also stress that MiMedSurv can conduct covariate-adjusted analysis to control for potential confounding factors (e.g., age, sex) to enhance the causality of the results especially for observational studies. MiMedSurv also provides user-friendly data preprocessing and analytic modules and makes nice visualizations.", style = "font-size:13pt")
  
  HOME_COMMENT2 = p(strong("URLs:"), " Web server: ", tags$a(href = "http://mimedsurv.micloud.kr", "http://mimedsurv.micloud.kr"), 
                    "; GitHub repository (local implementation):", 
                    tags$a(href = "https://github.com/yj7599/mimedsurvgit", "https://github.com/yj7599/mimedsurvgit"), style = "font-size:13pt")
  HOME_COMMENT3 = p(strong("Maintainer:"), " Hyojung Jang (", tags$a(href = "hyojung.jang@northwestern.edu", "hyojung.jang@northwestern.edu"), 
                    ")", style = "font-size:13pt")
  HOME_COMMENT4 = p(strong("Reference:"), " Jang, H., Koh, H. MiMedSurv: A Unified Web Cloud Computing Platform for Microbiome Causal Mediation Analysis with Survival Responses (In Review)", style = "font-size:13pt")
  
  INPUT_PHYLOSEQ_COMMENT1 = p(strong("Notice:"), "This should be an '.Rdata' or '.rds' file, and the data should be in the 'phyloseq' format", br(), 
                              "(see ", a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"), ")",
                              br(), br(),
                              "(1) To fully perform all community-level analyses for all non-phylogenetic and phylogenetic alpha- and beta-diversity indices as well as all taxonomy-level analyses, users should upload a feature table, a taxonomic table, a phylogenetic tree, and metadata. ", 
                              br(), br(),
                              "(2) To perform community-level analyses for only non-phylogenetic alpha- and beta-diversity indices as well as all taxonomy-level analyses, users can upload only a feature table, a taxonomic table, and metadata. ", br(), br(),
                              
                              strong("Details:"), br(), br(),
                              strong("Feature table: "), "It should contain counts, where rows are features (OTUs or ASVs) and columns are subjects
                              (row names are feature IDs and column names are subject IDs).", br(), br(),
                              strong("Taxonomic table: "), "It should contain taxonomic names, where rows are features and columns are seven
                              taxonomic ranks (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family',
                              'Genus', 'Species' or 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(), br(),
                              strong("Metadata/sample information:"), "It should contain variables for the subjects about host phenotypes, medical
                              interventions, disease status or environmental/behavioral factors, where rows are subjects and columns are variables
                              (row names are subject IDs, and column names are variable names).", br(), br(),
                              strong("Phylogenetic tree: "), "It should be a rooted tree. Otherwise, MiMedSurv automatically roots the tree through
                              midpoint rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs.", br(), br(),
                              "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                              The subjects should be matched and identical between feature table and metadata/sample information. 
                              MiMedSurv will analyze only the matched features and subjects.", style = "font-size:11pt")
  
  INPUT_PHYLOSEQ_COMMENT2 = p("You can download example microbiome data 'FTMP.Rdata', 'FTM.Rdata', in the 'phyloseq'
                              format. For more details about 'phyloseq', see", 
                              a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"), br(), br(), 
                              "> setwd('/yourdatadirectory/')", br(), br(), 
                              "> load(file = 'FTMP.Rdata')", br(), br(), 
                              "> library(phyloseq)", br(), br(), 
                              " > otu.tab <- otu_table(biom)", br(), 
                              " > tax.tab <- tax_table(biom)", br(), 
                              " > tree <- phy_tree(biom)", br(), 
                              " > sam.dat <- sample_data(biom)", br(), br(), 
                              "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree,
                              and the subjects are matched and identical between feature table and metadata/sample information using following 
                              code.", br(), br(), 
                              " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                              " > identical(colnames(otu.tab), rownames(sam.dat))", 
                              br(), br(),
                              strong("Reference:"), "Zhang XS, Li J, Krautkramer KA, Badri M, Battaglia T, Borbet TC et al. Antibiotic-induced acceleration of type 1 diabetes alters maturation of innate intestinal immunity. eLife 2018:7:e37816.", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT = p(strong("Notice:"), br(), br(),
                                    "(1) To fully perform all the community-level analyses for all the abundance-based and phylogenetic alpha- and beta-diversity indices as well as the taxonomic analyses, you have to upload all the four data components,
                                    feature (OTU or ASV) table, taxonomic table, metadata/sample information, and phylogenetic tree.", br(), br(),
                                    "(2) To perform the community-level analyses for only abundance-based alpha-diversity indices as well
                                    as the taxonomic analyses, you can upload only the three data components, feature (OTU or ASV) table, 
                                    taxonomic table, and metadata/sample information.",
                                    br(), br(), br(),
                                    strong("Details:"), br(), br(),
                                    strong("Feature table: "), "It should contain counts, where rows are features (OTUs or ASVs) and columns are
                                    subjects (row names are feature IDs and column names are subject IDs).", br(), br(),
                                    strong("Taxonomic table: "), "It should contain taxonomic names, where rows are features and columns are seven
                                    taxonomic ranks (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family',
                                    'Genus', 'Species' or 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(), br(),
                                    strong("Metadata/sample information: "), "It should contain variables for the subjects about host phenotypes,
                                    medical interventions, disease status or environmental/behavioral factors, where rows are subjects and columns
                                    are variables (row names are subject IDs, and column names are variable names).", br(), br(),
                                    strong("Phylogenetic tree: "), "It should be a rooted tree. Otherwise, MiMedSurv automatically roots the tree
                                    through midpoint rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs.",
                                    br(), br(),
                                    "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                                    The subjects should be matched and identical between feature table and metadata/sample information. MiMedSurv will
                                    analyze only the matched features and subjects.", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT2 = p("You can download example microbiome data 'oral.zip'. This zip file contains four data components, feature
                                     table (otu.tab.txt), taxonomic table (tax.tab.txt), metadata/sample information (sam.dat.txt), and 
                                     phylogenetic tree (tree.tre).", br(), br(),
                                     "> setwd('/yourdatadirectory/')", br(), br(), 
                                     "> otu.tab <- read.table(file = 'otu.tab.txt', check.names = FALSE) ", br(), 
                                     "> tax.tab <- read.table(file = 'tax.tab.txt', check.names = FALSE)", br(), 
                                     " > sam.dat <- read.table(file = 'sam.dat.txt', check.names = FALSE) ", br(),
                                     "> tree <- read.tree(file = 'tree.tre')", br(), br(), 
                                     "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic
                                     tree, and the subjects are matched and identical between feature table and metadata/sample information using
                                     following code.", br(), br(), 
                                     " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                                     " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                                     " > identical(colnames(otu.tab), rownames(sam.dat))", br(), br(),
                                     strong("Reference: "), "Park B, Koh H, Patatanian M, Reyes-Caballero H, Zhao N, Meinert J, et al. 
                                     The mediating roles of the oral microbiome in saliva and subgingival sites between e-cigarette smoking and
                                     gingival inflammation. BMC Microbiology. 2023;23(35):1-18.", style = "font-size:11pt")
  
  QC_KINGDOM_COMMENT = p("A microbial kingdom to be analyzed. Default is 'Bacteria' for 16S data. Alternatively, you can type 'Fungi' for ITS data 
                         or any other kingdom of interest for shotgun metagenomic data.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT1 = p("Remove subjects that have low library sizes (total read counts). Default is 3,000.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT2 = p("Library size: The total read count per subject.", style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT1 = p("Remove features (OTUs or ASVs) that have low mean relative abundances (Unit: %). Default is 0.02%.",style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT2 = p("Mean proportion: The average of relative abundances (i.e., proportions) per feature.", style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT1 = p('Remove taxonomic names in the taxonomic table that are completely matched with the specified character strings. 
                            Multiple character strings should be separated by a comma. Default is "", "metagenome", "gut metagenome", "mouse gut metagenome".',
                            style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT2 = p('Remove taxonomic names in the taxonomic table that are partially matched with the specified character strings (i.e., taxonomic names that contain 
                            the specified character strings). Multiple character strings should be separated by a comma. Default is "uncultured", "incertae", "Incertae",
                            "unidentified", "unclassified", "unknown".',style = "font-size:11pt")
  
  ALPHA_COMMENT = p("Calculate alpha-diversity indices: Richness (Observed), Shannon (Shannon, 1948), Simpson (Simpson, 1949), Inverse Simpson (Simpson, 1949), 
                    Fisher (Fisher et al., 1943), Chao1 (Chao, 1984), ACE (Chao and Lee, 1992), ICE (Lee and Chao, 1994), PD (Faith, 1992).")
  ALPHA_REFERENCES = p("1. Chao A, Lee S. Estimating the number of classes via sample coverage. J Am Stat Assoc. 1992:87:210-217.", br(),
                       "2. Chao A. Non-parametric estimation of the number of classes in a population. Scand J Stat. 1984:11:265-270.", br(),
                       "3. Faith DP. Conservation evaluation and phylogenetic diversity. Biol Conserv. 1992:61:1-10.", br(),
                       "4. Fisher RA, Corbet AS, Williams CB. The relation between the number of species and the number of individuals 
                       in a random sample of an animal population. J Anim Ecol. 1943:12:42-58.", br(),
                       "5. Lee S, Chao A. Estimating population size via sample coverage for closed capture-recapture models. Biometrics. 1994:50:1:88-97.", br(),
                       "6. Shannon CE. A mathematical theory of communication. Bell Syst Tech J. 1948:27:379-423 & 623-656.", br(),
                       "7. Simpson EH. Measurement of diversity. Nature 1949:163:688.", br())
  BETA_COMMENT = p("Calculate beta-diversity indices: Jaccard dissimilarity (Jaccard, 1912), Bray-Curtis dissimilarity (Bray and Curtis, 1957), Unweighted UniFrac distance 
                   (Lozupone and Knight, 2005), Generalized UniFrac distance (Chen et al., 2012), Weighted UniFrac distance (Lozupone et al., 2007).")
  BETA_REFERENCES = p("1. Bray JR, Curtis JT. An ordination of the upland forest communities of Southern Wisconsin. Ecol Monogr. 1957;27(32549).", br(),
                      "2. Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD., et al. Associating microbiome composition with environmental 
                      covariates using generalized UniFrac distances. Bioinformatics. 2012;28(16):2106-13.", br(),
                      "3. Jaccard P. The distribution of the flora in the alpine zone. New Phytol. 1912;11(2):37-50.", br(),
                      "4. Lozupone CA, Hamady M, Kelley ST, Knight R. Quantitative and qualitative β-diversity measures lead to 
                      different insights into factors that structure microbial communities. Appl Environ Microbiol. 2007;73(5):1576-85.", br(),
                      "5. Lozupone CA, Knight R. UniFrac: A new phylogenetic method for comparing microbial communities. Appl Environ Microbiol. 2005;71(12):8228-35.")
  DATA_TRANSFORM_COMMENT = p("Normalize the data into the CLR (centered log ratio) (Aitchison, 1982) format for each taxonomic rank (phylum, class, order, family, genus, species).")
  DATA_TRANSFORM_REFERENCE = p("1. Aitchison J. The statistical analysis of compositional data. J R Statist Soc B. 1982;44:2:139-77")
}


# UI ---------------------------------
ui <- dashboardPage(
  
  title = "MiMedSurv",
  
  dashboardHeader(title = span(TITLE, style = "float:left;font-size: 20px"), titleWidth = "100%"),
  dashboardSidebar(
    width = 240,

    chooseSliderSkin("Flat"),
    tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
    sidebarMenu(id = "side_menu",
                menuItem(strong("Home"), tabName = "home"),
                menuItem(strong("Data Processing"),
                         menuSubItem(strong("Data Input"), tabName = "step1", icon = icon("mouse")),
                         menuSubItem(strong("Quality Control"), tabName = "step2", icon = icon("chart-bar")),
                         menuSubItem(strong("Data Transformation"), tabName = "step3", icon = icon("calculator"))
                ),
                
                menuItem(strong("Non-Mediational Analysis"), tabName = "nonMed"),
                
                menuItem(strong("Mediational Analysis"), tabName = "surv-med",
                         menuItem(strong("Community-level Analysis"),  
                                  menuSubItem(strong("Alpha Diversity"), tabName = "alphaMed", icon = icon("font")),
                                  menuSubItem(strong("Beta Diversity"), tabName = "betaMed", icon = icon("bold"))),
                         menuItem(strong("Taxonomy-level Analysis"), tabName = "taxaMed")
                ))), 
  
  dashboardBody(
    customTheme, 
    tags$head(tags$style(HTML(".content { padding-top: 2px;}"))),
    
    tags$head(tags$style(HTML('.progress-bar {background-color: #ff7f50;}'))),
    
    tags$style(".nav-tabs-custom .nav-tabs li.active {border-top-color: #ff7f50;}
                 .shiny-notification {height: 80px; width: 600px; position:absolute; right: 0px; bottom: 0px; font-size: 11pt;}"),
    
    
    tags$script(src = "fileInput_text.js"),
    useShinyjs(),
    tabItems(
      
      ###### HOME ######
      
      tabItem(tabName = "home",
              div(id = "homepage", br(), HOME_COMMENT,
                  div(tags$img(src="mimedsurv_workflow.png", height="50%", width="50%",), style =  "text-align: center;"),
                  br(),
                  HOME_COMMENT2, HOME_COMMENT3, HOME_COMMENT4
              )),
      
      ##### DATA INPUT #####
      tabItem(tabName = "step1", br(), 
              fluidRow(column(width = 6, style='padding-left:+20px',
                              box(
                                width = NULL, status = "primary", solidHeader = TRUE,
                                title = strong("Data Input", style = "color:white"),
                                selectInput("inputOption", h4(strong("Data Type")), c("Choose one" = "", "Phyloseq", "Individual Data"), width = '30%'),
                                div(id = "optionsInfo", tags$p("You can choose phyloseq or individual data.", style = "font-size:11pt"), style = "margin-top: -15px"),
                                uiOutput("moreOptions"))),
                       column(width = 6, style='padding-left:0px', uiOutput("addDownloadinfo")))),
      
      
      ##### QC #####
      tabItem(tabName = "step2", br(), 
              sidebarLayout(
                position = "left",
                sidebarPanel(width = 3,
                             textInput("kingdom", h4(strong("Kingdom")), value = "Bacteria"),
                             
                             QC_KINGDOM_COMMENT,
                             
                             # Control the size of the numbers on the slider
                             tags$style(type = 'text/css', '#slider1 .irs-grid-text {font-size: 1px}'),
                             tags$style(type = 'text/css', '#slider2 .irs-grid-text {font-size: 1px}'), 
                             
                             sliderInput("slider1", h4(strong("Library size")), min=0, max=10000, value = 3000, step = 1000),
                             
                             QC_LIBRARY_SIZE_COMMENT1,
                             QC_LIBRARY_SIZE_COMMENT2,
                             
                             sliderInput("slider2", h4(strong("Mean proportion")), min = 0, max = 0.1, value = 0.02, step = 0.001,  post  = " %"),
                             
                             QC_MEAN_PROP_COMMENT1,
                             QC_MEAN_PROP_COMMENT2,
                             
                             br(),
                             p(" ", style = "margin-bottom: -20px;"),
                             
                             shinyjs::hidden(
                               shiny::div(id = "error_taxa", 
                                          h4(strong("Erroneous taxonomic names")),
                                          textInput("rem.str", label = "Complete match", value = ""),
                                          QC_TAXA_NAME_COMMENT1,
                                          
                                          textInput("part.rem.str", label = "Partial match", value = ""),
                                          QC_TAXA_NAME_COMMENT2,
                               )),
                             
                             actionButton("run", (strong("Run!")), class = "btn-info",
                                          style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50"), br(), br(),
                             uiOutput("moreControls")),
               
                mainPanel(width = 9,
                          fluidRow(width = 12,
                                   color = "primary", solidHeader = TRUE, 
                                   valueBoxOutput("sample_Size", width = 3),
                                   valueBoxOutput("OTUs_Size", width = 3),
                                   valueBoxOutput("phyla", width = 3),
                                   valueBoxOutput("classes", width = 3)),
                          fluidRow(width = 12, 
                                   status = "primary", solidHeader = TRUE,
                                   valueBoxOutput("orders", width = 3),
                                   valueBoxOutput("families", width = 3),
                                   valueBoxOutput("genera", width = 3),
                                   valueBoxOutput("species", width = 3)),
                          fluidRow(style = "position:relative",
                                   tabBox(width = 6, title = strong("Library Size", style = "color:black"), 
                                          tabPanel("Histogram",
                                                   plotlyOutput("hist"),
                                                   sliderInput("binwidth", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                   chooseSliderSkin("Round", color = "#ff7f50")),
                                          tabPanel("Box Plot", 
                                                   plotlyOutput("boxplot"))),
                                   tabBox(width = 6, title = strong("Mean Proportion", style = "color:black"), 
                                          tabPanel("Histogram",
                                                   plotlyOutput("hist2"),
                                                   sliderInput("binwidth2", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                   chooseSliderSkin("Round", color = "#ff7f50")),
                                          tabPanel("Box Plot",
                                                   plotlyOutput("boxplot2"))))))),
      
      ##### DATA TRANSFORMATION ##### 
      tabItem(tabName = "step3", br(),  
              fluidRow(column(width = 6, style = 'padding-left:+20px',
                              box(title = strong("Data Transformation", style = "color:white"), width = NULL, status = "primary", solidHeader = TRUE,
                                  ALPHA_COMMENT, 
                                  BETA_COMMENT,
                                  DATA_TRANSFORM_COMMENT,
                                  actionButton("divCalcRun", strong("Run!"), class = "btn-info"), 
                                  p(" ", style = "margin-bottom: +10px;"),
                                  p(strong("Attention:"), "You have to click this Run button to perform following community-level (alpha-and beta-diversity) analyses and taxonomy-level analyses.")
                              ),
                              
                              uiOutput("divCalcDownload"),
                              uiOutput("datTransDownload")),
                       column(width = 6, style='padding-left:0px',
                              box(title = strong("References", style = "color:white"), width = NULL, status = "primary", solidHeader = TRUE,
                                  p("Alpha Diversity", style = "font-size:12pt"),
                                  ALPHA_REFERENCES,
                                  p("Beta Diversity", style = "font-size:12pt"),
                                  BETA_REFERENCES, 
                                  p("Taxonomic Abundance", style = "font-size:12pt"),
                                  DATA_TRANSFORM_REFERENCE)))),
      
      ##### Non-Mediational Analysis #### 
      tabItem(tabName = "nonMed", br(), 
              tabPanel(
                title = "Non-Mediational Analysis",
                sidebarLayout(
                  position = "left",
                  sidebarPanel(width = 3, 
                               uiOutput("non.survtime"),
                               uiOutput("non.censor"),
                               uiOutput("non.treat"),
                               uiOutput("non.treat.rename"), 
                               uiOutput("non.cov"),
                               uiOutput("non.method"),
                               uiOutput("non.downloadTable"),
                               br(), 
                               
                               uiOutput("non.references")),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("non.display_results"))))
              )),
      
      ##### Mediational Analysis #####
      tabItem(tabName = "alphaMed", br(), 
              tabPanel(
                title = "Alpha Diversity",
                sidebarLayout(
                  position = "left",
                  sidebarPanel(width = 3,
                               uiOutput("outvars"),
                               uiOutput("censorvar"),
                               uiOutput("treatvars"),
                               uiOutput("outcome_types"),
                               br(), 
                               uiOutput("alpha.downloadTable"),
                               uiOutput("alpha.references")),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("alpha.display_results")))))),
      
      tabItem(tabName = "betaMed", br(), 
              tabPanel(
                title = "Beta Diversity",
                sidebarLayout(
                  position = "left",
                  sidebarPanel(width = 3,
                               uiOutput("beta.outvars"),
                               uiOutput("beta.censorvar"),
                               uiOutput("beta.treatvar"),
                               uiOutput("beta.rename.treatvar"), 
                               uiOutput("beta.outcome_types"),br(), 
                               uiOutput("beta.downloadTable"),
                               uiOutput("beta.references")),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("beta.display_results")))))),
      
      tabItem(tabName = "taxaMed", br(), 
              tabPanel(
                title = "Taxonomic Analysis",
                sidebarLayout( 
                  sidebarPanel(width = 3,
                               uiOutput("taxa_outcome"),
                               uiOutput("taxa_censor"),
                               uiOutput("taxa_treat"), 
                               p(" ", style = "margin-bottom: -10px;"),
                               uiOutput("taxa_int"),
                               uiOutput("taxa_cov"),  
                               p(" ", style = "margin-bottom: -10px;"),
                               uiOutput("taxa_method"),
                               uiOutput("taxa_reg_med"),
                               uiOutput("downloadTable_taxa"),
                               uiOutput("taxa_references")),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     div(style='height:2000px;overflow-y: scroll;', uiOutput("taxa_display")),
                                     br(),
                                     uiOutput("taxa_display_dend"),
                                     uiOutput("chat_explanation")
                            ),
                            br(), br())))
      )
    )
  )
)

# Server -----------------------------

server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  
  ###### Load Example Data ######
  
  sub.biom <- readRDS("Data/val_physeq.rds")
  
  biom <- sub.biom
  ori.biom <- biom
  otu.tab <- otu_table(ori.biom)
  tax.tab <- tax_table(ori.biom)
  tree <- phy_tree(ori.biom)
  sam.dat <- sample_data(ori.biom)
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("biom",".Rdata", sep = "")
    },
    content = function(file1) {
      save(biom, file = file1)
    })
  
  output$downloadZip <- downloadHandler(
    filename = function() {
      paste("biom",".zip", sep = "")
    },
    content <- function(fname) {
      temp <- setwd(tempdir())
      on.exit(setwd(temp))
      dataFiles = c("otu.tab.txt", "tax.tab.txt", "sam.dat.txt" ,"tree.tre")
      write.table(otu.tab, "otu.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(tax.tab, "tax.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(sam.dat, "sam.dat.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.tree(tree, "tree.tre")
      zip(zipfile=fname, files=dataFiles)
    })
  
  env <- new.env()
  nm <- load(file = "Data/FTMP.Rdata", env)[1]
  FTMP <- env[[nm]]
  
  env_2 <- new.env()
  nm_2 <- load(file = "Data/FTM.Rdata", env_2)[1]
  FTM <- env_2[[nm_2]]
  
  B.otu.tab <- otu_table(FTMP)
  B.tax.tab <- tax_table(FTMP)
  B.tree <- phy_tree(FTMP)
  B.sam.dat <- sample_data(FTMP)
  
  output$downloadData.1 <- downloadHandler(
    filename = function() {
      paste("FTMP.Rdata", sep = "")
    },
    content = function(file1) {
      save(FTMP, file = file1)
    })
  
  output$downloadData.2 <- downloadHandler(
    filename = function() {
      paste("FTM.Rdata", sep = "")
    },
    content = function(file1) {
      save(FTM, file = file1)
    })
  
  output$downloadData.1.Ind <- downloadHandler(
    filename = function() {
      paste("oral",".zip", sep = "")
    },
    content <- function(fname) {
      temp <- setwd(tempdir())
      on.exit(setwd(temp))
      dataFiles = c("Oral.otu.tab.txt", "Oral.tax.tab.txt", "Oral.sam.dat.txt" ,"Oral.tree.tre")
      write.table(B.otu.tab, "Oral.otu.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(B.tax.tab, "Oral.tax.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(B.sam.dat, "Oral.sam.dat.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.tree(B.tree, "Oral.tree.tre")
      zip(zipfile=fname, files=dataFiles)
    })
  
  
  ##### Define Variables #####
  infile = reactiveValues(biom = NULL, qc_biom = NULL, rare_biom = NULL)
  ds.Ks <- reactiveValues(res = NULL)
  chooseData = reactiveValues(sam.dat = NULL, mon.sin.rev.bin.con = NULL, treat_vars = NULL, out_vars = NULL, alpha.div = NULL,
                              alpha.div.rare = NULL, alpha.div.qc = NULL, taxa.out = NULL, tax.tab = NULL)
  is.results = reactiveValues(treatResult = NULL, outResult = NULL, beta.treatResult = NULL, beta.outResult = NULL)
  is.results.long = reactiveValues(treatResult = NULL, outResult = NULL)
  multi.test = reactiveValues(boolval = FALSE)
  multi.test.long = reactiveValues(boolval = FALSE)
  
  alpha.categos <- reactiveValues(cat1 = NULL, cat2 = NULL, cat3 = NULL, cat4 = NULL)
  alpha.categos.long <- reactiveValues(cat1 = NULL, cat2 = NULL)
  alpha.data.results = reactiveValues(table.out = NULL, data.q.out = NULL, table.p.out = NULL)
  alpha.results = reactiveValues(bin.var = NULL, alpha_div = NULL, alpha.bin.sum.out = NULL)
  alpha.results.cont = reactiveValues(alpha.con.out = NULL, alpha.table.out = NULL)
  alpha.reg.results = reactiveValues(bin.var = NULL, cov.var = NULL, alpha.div = NULL)
  alpha.resultslong = reactiveValues(bin.var = NULL, id.var = NULL, alpha_div = NULL, alpha.bin.sum.out = NULL)
  alpha.noncovs_res = reactiveValues(con.var = NULL, id.var = NULL, alpha.div = NULL)
  data.alphaBin_res = reactiveValues(table.output = NULL, data.output = NULL, alpha.bin.sum.out = NULL, table.p_outbin = NULL)
  data.results.cont_long = reactiveValues(table.out = NULL, data.q.out = NULL, table_p.out = NULL, alpha.table.out = NULL)
  
  beta.data.results = reactiveValues(data.q.out = NULL)
  beta.results = reactiveValues(result = NULL)
  beta.resultscont = reactiveValues(beta.cont.out = NULL)
  beta.data.results_long = reactiveValues(beta.bin.out = NULL)
  beta.resultscon_long = reactiveValues(beta.con.out = NULL)
  beta.categos <- reactiveValues(cat1 = NULL, cat2 = NULL, cat3 = NULL, cat4 = NULL)
  beta.categos.con <- reactiveValues(cat1 = NULL, cat2 = NULL)
  beta.down.results <- reactiveValues(CS = NULL, LONG = NULL)
  
  taxa.categos <- reactiveValues(cat1 = NULL, cat2 = NULL)
  taxa.data.results = reactiveValues(data.q.out = NULL)
  taxa.results = reactiveValues(bin.var = NULL, cov.var = NULL, id.var = NULL, taxa = NULL, taxa.bin.sum.out = NULL, con.var = NULL, taxa.con.sum.out = NULL, lib.size = NULL)
  taxa.types = reactiveValues(dataType = NULL, regression = NULL)
  taxa.outputs = reactiveValues(DAoutput = NULL, DAoutput_or = NULL, DAoutputlong = NULL)
  
  rcol = reactiveValues(selected = "lightblue")
  
  ##### DATA INPUT #####
  ## Input Options ##
  observeEvent(input$inputOption, {
    observe({
      if (input$inputOption == "Phyloseq") {
        
        shinyjs::hide(id = "optionsInfo")     
        
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("phyloseqData", strong("Please upload your 'phyloseq' data (.Rdata, .rds)", style = "color:black; border-color: #ff7f50"), 
                      accept = c(".Rdata", ".rds"), width = '80%'), div(style = "margin-top: -15px"),
            actionButton('Load_Phyloseq_Data', 'Upload', class = "btn-info", 
                         style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50"), br(),br(),
            shinyjs::hidden(
              shiny::div(id = "phyloseqUpload_error",
                         shiny::tags$p("Please upload a Rdata file!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_PHYLOSEQ_COMMENT1
          )
        })
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data", style = "color:white"), width = NULL, status = "primary", solidHeader = TRUE,
                p("4 Components (feature table, taxonomic table, metadata/sample information, and phylogenetic tree)"),
                downloadButton("downloadData.1", "FTMP.Rdata", width = '30%', 
                               style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50"),
                br(), br(), 
                
                p("3 Components (feature table, taxonomic table, and metadata/sample information)"),
                downloadButton("downloadData.2", "FTM.Rdata", width = '30%', 
                               style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50"),
                br(), br(),
                
                INPUT_PHYLOSEQ_COMMENT2
            )
          )
        })
      } 
      else if (input$inputOption == "Individual Data") {
        shinyjs::hide(id = "optionsInfo")
        
        
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("otuTable", strong("Please upload your feature (OTU or ASV) table (.txt, .csv, .biom)", style = "color:black"), 
                      accept = c(".txt", ".csv", ".biom"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("taxTable", strong("Please upload your taxonomic table (.txt, .tsv)", style = "color:black"), 
                      accept = c(".txt", ".tsv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("samData", strong("Please upload your metadata/sample information (.txt, .csv)", style = "color:black"), 
                      accept = c(".txt", ".csv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("tree", strong("Please upload your phylogenetic tree (.tre, .nwk)", style = "color:black"), 
                      accept = c(".tre", ".nwk"), width = '80%'), div(style = "margin-top: -15px"),
            actionButton('Load_Individual_Data', 'Upload', class = "btn-info",
                         style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50"), br(),br(),
            shinyjs::hidden(
              shiny::div(id = "textfilesUpload_error",
                         shiny::tags$p("Please upload txt and tre files!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_INDIVIDUAL_DATA_COMMENT
          )
        })
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data", style = "color:white"), width = NULL, status = "primary", solidHeader = TRUE,
                p("Oral Microbiome"), 
                downloadButton("downloadData.1.Ind", "oral.zip", width = '30%', 
                               style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50"),
                
                br(),br(),
                INPUT_INDIVIDUAL_DATA_COMMENT2
            )
          )
        })
      }
    })
  }, ignoreInit = TRUE, once = TRUE, ignoreNULL = TRUE)
  
  
  ## Enable or Disable an Input Element ##
  observe({
    
    toggleState("Load_Phyloseq_Data", !is.null(input$phyloseqData))
    toggleState("Load_Individual_Data", 
                !(is.null(input$otuTable) | is.null(input$samData) ))
    toggleState("run", !is.null(infile$biom))
    toggleState("skip", !is.null(infile$biom))
    toggleState("slider1", !is.null(infile$biom))
    toggleState("slider2", !is.null(infile$biom))
    toggleState("kingdom", !is.null(infile$biom))
    toggleState("binwidth", !is.null(infile$biom))
    toggleState("binwidth2", !is.null(infile$biom))
    toggleState("divCalcRun", !is.null(infile$rare_biom))
    
    if(!is.null(input$phyloseqData)){
      toggleState("divCalcRun", !is.null(access(infile$rare_biom, "sample_data")))
      toggleState("divCalcRun", !is.null(access(infile$rare_biom, "otu_table")))
      toggleState("divCalcRun", !is.null(access(infile$rare_biom, "tax_table")))
    }else{
      toggleState("divCalcRun", !is.null(input$otuTable))
      toggleState("divCalcRun", !is.null(input$samData))
      toggleState("divCalcRun", !is.null(input$taxTable))
    }
  })
  
  observeEvent(input$Load_Phyloseq_Data, {
    
    if (!is.null(input$phyloseqData)) {
      
      
      dataInfile  = reactive({
        phyloseq.data = input$phyloseqData
        ext <- tools::file_ext(phyloseq.data$datapath)
        
        req(phyloseq.data)
        if (ext == "Rdata") {
          
          phyloseq.dataPath = phyloseq.data$datapath
          e = new.env()
          name <- load(phyloseq.dataPath, envir = e)
          data <- e[[name]]
          
          is.taxa <- NULL 
          try(is.taxa <- tax_table(data), silent = TRUE)
          
          if (!is.null(is.taxa)){
            shinyjs::show("error_taxa")
          }else{
            shinyjs::hide("error_taxa")
          }
          shinyjs::disable("Load_Individual_Data")
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          
          try(colnames(tax_table(data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), silent = TRUE)
          
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          return(data)
        } 
        else if (ext == "rds") {
          phyloseq.dataPath = phyloseq.data$datapath
          data <- readRDS(phyloseq.dataPath)
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          
          try(colnames(tax_table(data)) <-  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), silent = TRUE) 
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          
          return(data)
        } 
        else {
          shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade")
          shinyjs::delay(5000, shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade"))
          return(NULL)
        }
      })
    } 
    else {
      return(NULL)
    }
    
    if (is.null(dataInfile)) {
      infile$biom <- NULL
      infile$qc_biom <- NULL
      infile$rare_biom = NULL
    } 
    else {
      infile$biom <- dataInfile()
      infile$qc_biom <- dataInfile()
      infile$rare_biom = NULL
    }
    
    updateTabsetPanel(session, "side_menu",
                      selected = "step2")
    rcol$selected = "rgba(86, 180, 233, 0.6)"
    
    if (!is.null(infile$biom)) QC$resume()
  })
  
  observeEvent(input$Load_Individual_Data, {
    shinyjs::disable("Load_Individual_Data")
    
    if (!is.null(input$taxTable)){
      shinyjs::show("error_taxa")
    }else{
      shinyjs::hide("error_taxa")
    }
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        incProgress(3/10, message = "File Check")
        
        
        if (!is.null(input$otuTable) & !is.null(input$taxTable) & !is.null(input$samData) & !is.null(input$tree)) {

          dataInfile  = reactive({
            otu.table = input$otuTable
            ext1 <- tools::file_ext(otu.table$datapath)
            
            tax.table = input$taxTable
            ext2 <- tools::file_ext(tax.table$datapath)
            
            sam.data = input$samData
            ext3 <- tools::file_ext(sam.data$datapath)
            
            tree.data = input$tree
            ext4 <- tools::file_ext(tree.data$datapath)
            
            req(otu.table, tax.table, sam.data, tree.data)
            if ((ext1 == "txt"| ext1 == "csv" | ext1 == "biom") & (ext2 == "txt" | ext2 == "tsv") &
                (ext3 == "txt" | ext3 == "csv") & (ext4 == "tre" | ext4 == "nwk")) {
              otu.table.path = otu.table$datapath
              tax.table.path = tax.table$datapath
              sam.data.path = sam.data$datapath
              tree.data.path = tree.data$datapath
              
              if (ext1 == "txt") {
                otu.tab <- read.table(otu.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext1 == "csv") {
                otu.tab <- read.csv(otu.table.path, check.names = FALSE)
                rownames(otu.tab) = otu.tab[,1];otu.tab = otu.tab[,-1]
              } else if (ext1 == "biom") {
                biom <- read_biom(otu.table.path)
                otu.tab <- as.matrix(biom_data(biom))
              }
              
              if (ext2 == "txt") {
                tax.tab <- read.table(tax.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext2 == "tsv") {
                tax.tab <- read.table(tax.table.path, header=TRUE, sep="\t")
                tax.tab = preprocess.tax.tab(tax.tab)
              }
              
              if (ext3 == "txt") {
                sam.dat <- read.table(sam.data.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext3 == "csv") {
                sam.dat <- read.csv(sam.data.path, check.names = FALSE)
                rownames(sam.dat) = sam.dat[,1]
                sam.dat = sam.dat[,-1]
              }
              
              if (ext4 == "tre") {
                tree <- read.tree(file = tree.data.path)
              } else if (ext4 == "nwk") {
                tree <- read.tree(file = tree.data.path) 
              }
              
              otu.tab <- otu_table(otu.tab, taxa_are_rows = TRUE)
              tax.tab <- tax_table(as.matrix(tax.tab))
              sam.dat <- sample_data(sam.dat)
              tree <- phy_tree(tree)
              
              if (sum(colnames(otu.tab) %in% rownames(sam.dat)) < sum(rownames(otu.tab) %in% rownames(sam.dat))) {
                otu.tab = t(otu.tab)
              }
              
              incProgress(3/10, message = "Validating")
              validate(
                if (biom.check.samples(otu.tab, sam.dat)) {
                  if (biom.check.otu(otu.tab, tax.tab, tree)) {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data. And
                                        there is no common OTUs among OTU/feature table, taxonomic table and tree tip labels"),
                                     type = "error")
                  } else {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data"),
                                     type = "error")
                  }
                } else if (biom.check.otu(otu.tab, tax.tab, tree)) {
                  showNotification(h4("Error: There is no common OTUs among OTU/feature table, taxonomic table and tree tip labels"),
                                   type = "error")
                } else {
                  NULL
                }
              )
              
              incProgress(1/10, message = "Merging")
              biomData <- merge_phyloseq(otu.tab, tax.tab, sam.dat, tree)
              return(biomData)
            } else {
              shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(5000, shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade"))
              return(NULL)
            }
          })
        } else if(!is.null(input$otuTable) & is.null(input$taxTable) & !is.null(input$samData) & is.null(input$tree))   #only sam.dat / otu.tab
        {
          
          dataInfile  = reactive({
            otu.table = input$otuTable
            ext1 <- tools::file_ext(otu.table$datapath)
            
            sam.data = input$samData
            ext3 <- tools::file_ext(sam.data$datapath)
            
            if ((ext1 == "txt"| ext1 == "csv" | ext1 == "biom") & (ext3 == "txt" | ext3 == "csv") ) {
              otu.table.path = otu.table$datapath
              sam.data.path = sam.data$datapath
              
              if (ext1 == "txt") {
                otu.tab <- read.table(otu.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext1 == "csv") {
                otu.tab <- read.csv(otu.table.path, check.names = FALSE)
                rownames(otu.tab) = otu.tab[,1];otu.tab = otu.tab[,-1]
              } else if (ext1 == "biom") {
                biom <- read_biom(otu.table.path)
                otu.tab <- as.matrix(biom_data(biom))
              }
              
              
              if (ext3 == "txt") {
                sam.dat <- read.table(sam.data.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext3 == "csv") {
                sam.dat <- read.csv(sam.data.path, check.names = FALSE)
                rownames(sam.dat) = sam.dat[,1]
                sam.dat = sam.dat[,-1]
              }
              
              otu.tab <- otu_table(otu.tab, taxa_are_rows = TRUE)
              sam.dat <- sample_data(sam.dat)
              
              if (sum(colnames(otu.tab) %in% rownames(sam.dat)) < sum(rownames(otu.tab) %in% rownames(sam.dat))) {
                otu.tab = t(otu.tab)
              }
              
              incProgress(3/10, message = "Validating")
              
              validate(
                if (biom.check.samples(otu.tab, sam.dat)) {
                  showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data"),
                                   type = "error")
                } else {
                  NULL
                }
              )
              
              incProgress(1/10, message = "Merging")
              biomData <- merge_phyloseq(otu.tab, sam.dat)
              return(biomData)
            } else {
              shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(5000, shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade"))
              return(NULL)
            }
          })
          
        }
        else if(!is.null(input$otuTable) & !is.null(input$taxTable) & !is.null(input$samData) & is.null(input$tree))   #only sam.dat / otu.tab
        {
          
          dataInfile  = reactive({
            otu.table = input$otuTable
            ext1 <- tools::file_ext(otu.table$datapath)
            
            tax.table = input$taxTable
            ext2 <- tools::file_ext(tax.table$datapath)
            
            sam.data = input$samData
            ext3 <- tools::file_ext(sam.data$datapath)
            
            req(otu.table, tax.table, sam.data)
            
            if ((ext1 == "txt"| ext1 == "csv" | ext1 == "biom") & (ext2 == "txt" | ext2 == "tsv") &(ext3 == "txt" | ext3 == "csv") ) {
              otu.table.path = otu.table$datapath
              tax.table.path = tax.table$datapath
              sam.data.path = sam.data$datapath
              
              if (ext1 == "txt") {
                otu.tab <- read.table(otu.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext1 == "csv") {
                otu.tab <- read.csv(otu.table.path, check.names = FALSE)
                rownames(otu.tab) = otu.tab[,1];otu.tab = otu.tab[,-1]
              } else if (ext1 == "biom") {
                biom <- read_biom(otu.table.path)
                otu.tab <- as.matrix(biom_data(biom))
              }
              
              if (ext2 == "txt") {
                tax.tab <- read.table(tax.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext2 == "tsv") {
                tax.tab <- read.table(tax.table.path, header=TRUE, sep="\t")
                tax.tab = preprocess.tax.tab(tax.tab)
              }
              
              if (ext3 == "txt") {
                sam.dat <- read.table(sam.data.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext3 == "csv") {
                sam.dat <- read.csv(sam.data.path, check.names = FALSE)
                rownames(sam.dat) = sam.dat[,1]
                sam.dat = sam.dat[,-1]
              }
              
              otu.tab <- otu_table(otu.tab, taxa_are_rows = TRUE)
              tax.tab <- tax_table(as.matrix(tax.tab))
              sam.dat <- sample_data(sam.dat)
              
              if (sum(colnames(otu.tab) %in% rownames(sam.dat)) < sum(rownames(otu.tab) %in% rownames(sam.dat))) {
                otu.tab = t(otu.tab)
              }
              
              incProgress(3/10, message = "Validating")
              validate(
                if (biom.check.samples(otu.tab, sam.dat)) {
                  if (biom.check.otu.without.tree(otu.tab, tax.tab)) {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data. And
                                        there is no common OTUs among OTU/feature table, taxonomic table and tree tip labels"),
                                     type = "error")
                  } else {
                    showNotification(h4("Error: There is no common OTUs among OTU/feature table and Sample Data"),
                                     type = "error")
                  }
                } else if (biom.check.otu.without.tree(otu.tab, tax.tab)) {
                  showNotification(h4("Error: There is no common OTUs among OTU/feature table and taxonomic table"),
                                   type = "error")
                } else {
                  NULL
                }
              )
              
              incProgress(1/10, message = "Merging")
              biomData <- merge_phyloseq(otu.tab, tax.tab, sam.dat)
              return(biomData)
            } else {
              shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(5000, shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade"))
              return(NULL)
            }
          })
          
        }
        else if(!is.null(input$otuTable) & is.null(input$taxTable) & !is.null(input$samData) & !is.null(input$tree))
        {
          dataInfile  = reactive({
            otu.table = input$otuTable
            ext1 <- tools::file_ext(otu.table$datapath)
            
            sam.data = input$samData
            ext3 <- tools::file_ext(sam.data$datapath)
            
            tree.data = input$tree
            ext4 <- tools::file_ext(tree.data$datapath)
            
            req(otu.table, sam.data, tree.data)
            
            if ((ext1 == "txt"| ext1 == "csv" | ext1 == "biom")  &
                (ext3 == "txt" | ext3 == "csv") & (ext4 == "tre" | ext4 == "nwk")) {
              otu.table.path = otu.table$datapath
              sam.data.path = sam.data$datapath
              tree.data.path = tree.data$datapath
              
              if (ext1 == "txt") {
                otu.tab <- read.table(otu.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext1 == "csv") {
                otu.tab <- read.csv(otu.table.path, check.names = FALSE)
                rownames(otu.tab) = otu.tab[,1];otu.tab = otu.tab[,-1]
              } else if (ext1 == "biom") {
                biom <- read_biom(otu.table.path)
                otu.tab <- as.matrix(biom_data(biom))
              }
              
              if (ext3 == "txt") {
                sam.dat <- read.table(sam.data.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext3 == "csv") {
                sam.dat <- read.csv(sam.data.path, check.names = FALSE)
                rownames(sam.dat) = sam.dat[,1]
                sam.dat = sam.dat[,-1]
              }
              
              if (ext4 == "tre") {
                tree <- read.tree(file = tree.data.path)
              } else if (ext4 == "nwk") {
                tree <- read.tree(file = tree.data.path) 
              }
              
              otu.tab <- otu_table(otu.tab, taxa_are_rows = TRUE)
              sam.dat <- sample_data(sam.dat)
              tree <- phy_tree(tree)
              
              if (sum(colnames(otu.tab) %in% rownames(sam.dat)) < sum(rownames(otu.tab) %in% rownames(sam.dat))) {
                otu.tab = t(otu.tab)
              }
              
              incProgress(3/10, message = "Validating")
              validate(
                if (biom.check.samples(otu.tab, sam.dat)) {
                  if (biom.check.otu.without.tax.tab(otu.tab, tree)) {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data. And
                                        there is no common OTUs among OTU/feature tableand tree tip labels"),
                                     type = "error")
                  } else {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data"),
                                     type = "error")
                  }
                } else if (biom.check.otu.without.tax.tab(otu.tab, tree)) {
                  showNotification(h4("Error: There is no common OTUs among OTU/feature tableand tree tip labels"),
                                   type = "error")
                } else {
                  NULL
                }
              )
              
              incProgress(1/10, message = "Merging")
              biomData <- merge_phyloseq(otu.tab, tax.tab, sam.dat, tree)
              return(biomData)
            } else {
              shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(5000, shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade"))
              return(NULL)
            }
          })
          
        }else {
          showNotification(h4("Error: OTU/feature table and Sample Data are required to proceed with data processing"),
                           type = "error")
        }
        
        if (is.null(dataInfile)) {
          infile$biom <- NULL
          infile$qc_biom <- NULL
          infile$rare_biom = NULL
        } else {
          infile$biom <- dataInfile()
          infile$qc_biom <- dataInfile()
          infile$rare_biom = NULL
        }
        
        updateTabsetPanel(session, "side_menu",
                          selected = "step2")
        rcol$selected = "rgba(86, 180, 233, 0.6)"
        
        if (!is.null(infile$biom)) QC$resume()
      })
    shinyjs::enable("Load_Individual_Data")
  })
  
  ##### QC #####
  # This reactive expression stores the input data from either the individual data or phyloseq data
  QC = observe(suspended = T,{
    taxa.results$lib.size <- lib.size.func(infile$biom)$lib.size
    
    dat_to_check <<- infile$biom
    is.sam.dat <- NULL 
    try(is.sam.dat <- sample_data(dat_to_check), silent = TRUE)
    
    if(is.null(access(infile$biom, "tax_table"))){
      showNotification(h4("There is no taxonomic table. Therefore, no outcome is available for taxonomic analysis."), duration = 4, type = "message")  
    }else if(is.null(is.sam.dat)){
      showNotification(h4("There is no sample data. As the minimum requirements are sample data and feature table, no outcome is available!"), duration = 4, type = "message")  
    }else if (is.null(access(infile$biom, "otu_table"))){
      showNotification(h4("There is no feature table. As the minimum requirements are sample data and feature table, no outcome is available!"), duration = 4, type = "message")  
    }
    
    # Plots graphs using example infile data
    output$hist <- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      plot_ly(x = ~lib_size, nbinsx = input$binwidth,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Library Size", zeroline = FALSE))
    })
    
    output$hist2 <- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      plot_ly(x = ~mean_prop, nbinsx = input$binwidth2,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Mean Proportion", zeroline = FALSE))
    })
    
    output$boxplot<- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      
      plot_ly(x = ~lib_size, type = "box", notched=TRUE, name = "Library Size",
              color = ~"lib_size", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    output$boxplot2<- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      
      plot_ly(x = ~mean_prop, type = "box", notched=TRUE, name = "Mean Proportion",
              color = ~"mean_prop", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    ## Number of Taxonomic Rank for biom before QC
    
    num_tax.rank <- NULL 
    tax.tab <- NULL 
    try(tax.tab <- tax_table(infile$qc_biom), silent = TRUE)
    
    if(!is.null(tax.tab)){
      num_tax.rank <- reactive({
        num.tax.rank(tax.tab)
      })
    }
    
    ## Fills value boxes using example biom data
    output$sample_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.sams), style = "font-size: 75%;"),
        "Sample Size", icon = icon("user-circle"), color = "yellow")
    })
    
    output$OTUs_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.otus), style = "font-size: 75%;"),
        "Number of Features", icon = icon("dna"), color = "aqua")
    })
    
    if(!is.null(num_tax.rank)){
      
      shinyjs::show("phyla")
      shinyjs::show("classes")
      shinyjs::show("orders")
      shinyjs::show("families")
      shinyjs::show("genera")
      shinyjs::show("species")
      
      output$phyla <- renderValueBox({
        num.phyla = num_tax.rank()[1]
        valueBox(
          value = tags$p(paste0(num.phyla), style = "font-size: 75%;"),
          "Number of Phyla", icon = icon("sitemap"), color = "orange")
      })
      
      output$classes <- renderValueBox({
        num.classes = num_tax.rank()[2]
        valueBox(
          value = tags$p(paste0(num.classes), style = "font-size: 75%;"),
          "Number of Classes", icon = icon("sitemap"), color = "purple")
      })
      
      output$orders <- renderValueBox({
        num.orders = num_tax.rank()[3]
        valueBox(
          value = tags$p(paste0(num.orders), style = "font-size: 75%;"),
          "Number of Orders", icon = icon("sitemap"), color = "blue")
      })
      
      output$families <- renderValueBox({
        num.families = num_tax.rank()[4]
        valueBox(
          value = tags$p(paste0(num.families), style = "font-size: 75%;"),
          "Number of Families", icon = icon("sitemap"), color = "red")
      })
      
      output$genera <- renderValueBox({
        num.genera = num_tax.rank()[5]
        valueBox(
          value = tags$p(paste0(num.genera), style = "font-size: 75%;"),
          "Number of Genera", icon = icon("sitemap"), color = "light-blue")
      })
      
      output$species <- renderValueBox({
        num.species = num_tax.rank()[6]
        valueBox(
          value = tags$p(paste0(num.species), style = "font-size: 75%;"),
          "Number of Species", icon = icon("sitemap"), color = "teal" )
      })
    }else{
      shinyjs::hide("phyla")
      shinyjs::hide("classes")
      shinyjs::hide("orders")
      shinyjs::hide("families")
      shinyjs::hide("genera")
      shinyjs::hide("species")
    }
    
    ## This event handler checks whether there is an input file and updates the slider options accordingly
    maxi.slider1 = as.numeric(lib.size.func(infile$qc_biom)$lib.size.sum["3rd quartile"])
    max.mean.prop = as.numeric(mean.prop.func(infile$qc_biom)$mean.prop.sum["3rd quartile"])
    maxi.slider2 = round(max.mean.prop, digits = 6)
    
    if (maxi.slider2 < 2e-05) {
      maxi.slider2 = 2e-05
    }
    
    updateSliderInput(session, "slider1", min = 0, max = round(maxi.slider1,-3))
    updateSliderInput(session, "slider2", min = 0, max = maxi.slider2*100)
  })
  
  ##### Non-Mediational Analysis #####

  observeEvent(chooseData$sam.dat, {
    out_taxa_options <- surv.pri.func(chooseData$sam.dat, is.mon.sin.rev.bin.con(chooseData$sam.dat))[[2]]
    censor_taxa_options <- surv.pri.func(chooseData$sam.dat, is.mon.sin.rev.bin.con(chooseData$sam.dat))[[1]]

    output$non.survtime <- renderUI({
      tagList(
        h4(strong("Survival time", style = "color:black")),
        p(" ", style = "margin-bottom: +15px;"),
        selectInput("non.outvar", label = NULL,
                    c("Choose one" = "", out_taxa_options), "T1Dweek", width = '80%'),
        p(" ", style = "margin-bottom: -10px;"),
        p("A variable for follow-up time (required).",
          style = "font-size:10pt")
      )
    })

    output$non.censor <- renderUI({
      tagList(
        p(" ", style = "margin-top: 25px;"),
        h4(strong("Censored/event", style = "color:black")),
        p(" ", style = "margin-bottom: +15px;"),
        selectInput("non.censorvar", label = NULL,
                    c("Choose one" = "", censor_taxa_options), "T1D", width = '80%'),
        p(" ", style = "margin-bottom: -10px;"),
        p("A binary indicator for censored (0) or event (1) (required).",
          style = "font-size:10pt")
      )
    })

    observeEvent(input$non.outvar, {
      observeEvent(input$non.censorvar, {
        output$non.treat <- renderUI({
          
          tagList(
            p(" ", style = "margin-top: 25px;"),
            h4(strong("Treatment variable", style = "color:black")),
            p(" ", style = "margin-bottom: +15px;"),

            selectInput("non.treatvar", label = NULL,
                        c("Choose one" = "", select.outcome.func(chooseData$sam.dat, colnames(chooseData$sam.dat), c(input$non.censorvar, input$non.outvar))), selected = "Antibiotics", width = '80%'),
            p(" ", style = "margin-bottom: -10px;"),
            p("Environmental, behavioral or medical exposures (e.g., diet, residence, smoking, preterm birth, delivery mode, antibiotic/probiotic use).",
              style = "font-size:10pt"))
        })
        
        #is.results$treatResult = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$treatvar)
        ntselected.treat_vars = select.outcome.func(chooseData$sam.dat, chooseData$treat_vars, input$non.treatvar)
      })
    })

    observeEvent(input$non.outvar, {
      observeEvent(input$non.censorvar, {
        observeEvent(input$non.treatvar, {
          
          non.med.treat_type <- is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$non.treatvar) 
          

          if(non.med.treat_type == "Binary"){
            non.cat.1 <<- beta.bin.cat.func(chooseData$sam.dat, input$non.treatvar)[1]
            non.cat.2 <<- beta.bin.cat.func(chooseData$sam.dat, input$non.treatvar)[2]
            
            output$non.treat.rename <- renderUI({
              tagList(
                p(" ", style = "margin-top: 25px;"),
                h4(strong("Rename treatment variable?", style = "color:black")),
                p(" ", style = "margin-bottom: +15px;"),
                textInput("non.med.rename_bin_var1", label = (paste0("Reference: ", non.cat.1)), value = non.cat.1, width = '80%'),
                textInput("non.med.rename_bin_var2", label = (paste0("Comparison: ", non.cat.2)), value = non.cat.2, width = '80%'),
                p(" ", style = "margin-bottom: -10px;"),
                p("You can rename the categories of treatment variable. MiMedSurv keeps up to 8 characters on graphs.", style = "font-size:10pt"))
            })
          }else {
            shinyjs::hide("non.treat.rename")
          }
          
          if (input$non.outvar %in% chooseData$out_vars) {
            
            ntselected.vars <<- select.covariates.func.surv(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$non.treatvar, input$non.outvar, input$non.censorvar)
            
            is.results$outResult = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$non.outvar)
            
            output$non.cov <- renderUI({
              tagList(
                p(" ", style = "margin-top: 25px;"),
                h4(strong("Covariate(s)", style = "color:black")),
                p(" ", style = "margin-bottom: +15px;"),
                prettyRadioButtons("non.covariate",label = NULL,
                                   icon = icon("check"),
                                   animation = "jelly", c("None", "Covariate(s)"), selected = "None", width = '70%'),
                p(" ", style = "margin-bottom: -10px;"),
                p("Potential confounders (e.g., age, gender) to be adjusted for.", style = "font-size:10pt"),
                
                shinyjs::hidden(
                  shiny::div(id = "non.covariates_variables", style = "margin-left: 2%",
                             prettyCheckboxGroup("non.covariateOptions"," Please select covariate(s)",
                                                 choices = ntselected.vars, width = '70%'))))
            })
            
            observeEvent(input$non.covariate,{
              if (input$non.covariate == "Covariate(s)") {
                output$non.method <- renderUI({
                  tagList(
                    p(" ", style = "margin-top: 25px;"),
                    h4(strong("Analytic method", style = "color:black")),
                    p(" ", style = "margin-bottom: +15px;"),
                    prettyRadioButtons("non.chooseMethod",label = NULL, icon = icon("check"),
                                       animation = "jelly", "Cox model",
                                       selected = "Cox model", width = '70%'),
                    actionButton("runbtn_nonmed", (strong("Run!")), class = "btn-info",
                                 style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50")
                  )
                })
                
                shinyjs::hide("non.choosetest")
                shinyjs::show("non.covariates_variables")
              }
              else if (input$non.covariate == "None") {
                
                output$non.method <- renderUI({
                  tagList(
                    p(" ", style = "margin-top: 25px;"),
                    h4(strong("Analytic method", style = "color:black")),
                    p(" ", style = "margin-bottom: +15px;"),
                    prettyRadioButtons("non.chooseMethod",label = NULL, icon = icon("check"),
                                       animation = "jelly", c("Kaplan-Meier analysis (default)", "Cox model"),
                                       selected = "Kaplan-Meier analysis (default)", width = '70%'),
                    uiOutput("non.choosetest"),
                    
                    actionButton("runbtn_nonmed", (strong("Run!")), class = "btn-info",
                                 style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50")
                  )
                })
                shinyjs::hide("non.covariates_variables")
              }
            })
            
            output$non.choosetest <- renderUI({
              tagList(
                p(" ", style = "margin-top: -20px;"),
                selectInput("non.test", label = h4(strong("", style = "color:black")),
                            choices = c("Log-rank test", "Wilcoxon test"),
                            selected = "Log-rank test", width = '80%'),
              )
            })
            
            observeEvent(input$non.chooseMethod,{
              if(input$non.chooseMethod == "Kaplan-Meier analysis (default)") {
                shinyjs::show("non.choosetest")
              }
              else if(input$non.chooseMethod == "Cox model"){
                shinyjs::hide("non.choosetest")
              }
            })
          }
        })
      })
    })
  })
  ##### DIVERSITY CALCULATION #####
  
  observeEvent(chooseData$alpha.div, {
    
    ##### ALPHA DIVERSITY #####
    
    out_taxa_options <- surv.pri.func(chooseData$sam.dat, is.mon.sin.rev.bin.con(chooseData$sam.dat))[[2]]
    censor_taxa_options <- surv.pri.func(chooseData$sam.dat, is.mon.sin.rev.bin.con(chooseData$sam.dat))[[1]]
    
    
    ## Outcome Variable UI
    output$outvars <- renderUI({
      tagList(
        h4(strong("Survival time", style = "color:black")),
        p(" ", style = "margin-bottom: +15px;"),
        selectInput("outvar", label = NULL,
                    c("Choose one" = "", out_taxa_options), "T1Dweek", width = '80%'),
        p(" ", style = "margin-bottom: -10px;"),
        p("A variable for follow-up time (required).",
          style = "font-size:10pt"),
        
      )
    })
    
    output$censorvar <- renderUI({
      tagList(
        p(" ", style = "margin-top: 25px;"),
        h4(strong("Censored/event", style = "color:black")),
        p(" ", style = "margin-bottom: +15px;"),
        selectInput("censorvar", label = NULL,
                    c("Choose one" = "", censor_taxa_options), "T1D", width = '80%'),
        p(" ", style = "margin-bottom: -10px;"),
        p("A binary indicator for censored (0) or event (1) (required).",
          style = "font-size:10pt")
      )
    })
    
    observeEvent(input$outvar, {
      observeEvent(input$censorvar, {
        output$treatvars <- renderUI({
          tagList(
            p(" ", style = "margin-top: 25px;"),
            h4(strong("Treatment variable", style = "color:black")),
            p(" ", style = "margin-bottom: +15px;"),
            
            selectInput("treatvar", label = NULL,
                        c("Choose one" = "", select.outcome.func(chooseData$sam.dat, colnames(chooseData$sam.dat), c(input$censorvar, input$outvar))), selected = "Antibiotics", width = '80%'),
            p(" ", style = "margin-bottom: -10px;"),
            p("Environmental, behavioral or medical exposures (e.g., diet, residence, smoking, preterm birth, delivery mode, antibiotic/probiotic use).",
              style = "font-size:10pt"))
        })
        
        ntselected.treat_vars = select.outcome.func(chooseData$sam.dat, chooseData$treat_vars, input$treatvar)
        
      }) 
    }) 
    
    output$outcome_types <- renderUI({
      tagList(
        uiOutput("covariates"),
        uiOutput("chooseMethodType"),
        actionButton("runbtn_alpha", (strong("Run!")), class = "btn-info",
                     style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50")
      )
    })
    
    observeEvent(input$outvar,{
      observeEvent(input$censorvar, {
        observeEvent(input$treatvar,{
          # User selects whether to use rarefied or non rarefied biom data for outcome variable
          if (input$outvar %in% chooseData$out_vars) {
            
            ntselected.vars <<- select.covariates.func.surv (chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$treatvar, input$outvar, input$censorvar)
            
            is.results$outResult = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$outvar)
            
            if (is.results$outResult == "Binary") {
              output$covariates <- renderUI({
                tagList(
                  p(" ", style = "margin-top: 25px;"),
                  h4(strong("Covariate(s)", style = "color:black")),
                  p(" ", style = "margin-bottom: +15px;"),
                  prettyRadioButtons("covariate",label = NULL, 
                                     icon = icon("check"),
                                     animation = "jelly", c("None", "Covariate(s)"), selected = "None", width = '70%'),
                  p(" ", style = "margin-bottom: -10px;"),
                  p("Potential confounders (e.g., age, gender) to be adjusted for.", style = "font-size:10pt"),
                  
                  shinyjs::hidden(
                    shiny::div(id = "covariates_variables", style = "margin-left: 2%",
                               prettyCheckboxGroup("covariateOptions"," Please select covariate(s)", 
                                                   choices = ntselected.vars, width = '70%'))))
              })
              
              ## Covariates Effect
              observeEvent(input$covariate,{
                if (input$covariate == "Covariate(s)") {
                  shinyjs::show("covariates_variables")
                } 
                else if (input$covariate == "None") {
                  shinyjs::hide("covariates_variables")
                }
              })
              
              output$chooseMethodType <- renderUI({
                tagList(
                  p(" ", style = "margin-top: 25px;"),
                  h4(strong("Analytic method", style = "color:black")),
                  p(" ", style = "margin-bottom: +15px;"),
                  prettyRadioButtons("chooseMethod",label = NULL, icon = icon("check"),
                                     animation = "jelly", c("Imai Method"), 
                                     selected = "Imai Method", width = '70%'))
              })
              
            }else if (is.results$outResult == "Continuous") {
              
              output$covariates <- renderUI({
                tagList(
                  p(" ", style = "margin-top: 25px;"),
                  h4(strong("Covariate(s)", style = "color:black")),
                  p(" ", style = "margin-bottom: +15px;"),
                  prettyRadioButtons("covariate",label = NULL, 
                                     icon = icon("check"),
                                     animation = "jelly", c("None", "Covariate(s)"), selected = "None", width = '70%'),
                  p(" ", style = "margin-bottom: -10px;"),
                  p("Potential confounders (e.g., age, gender) to be adjusted for.", style = "font-size:10pt"),
                  
                  shinyjs::hidden(
                    shiny::div(id = "covariates_variables", style = "margin-left: 2%",
                               prettyCheckboxGroup("covariateOptions"," Please select covariate(s)", 
                                                   choices = ntselected.vars, width = '70%'))))
              })
              
              ## Covariates Effect
              observeEvent(input$covariate,{
                if (input$covariate == "Covariate(s)") {
                  
                  shinyjs::show("covariates_variables")
                  
                  output$chooseMethodType <- renderUI({
                    tagList(
                      p(" ", style = "margin-top: 25px;"),
                      h4(strong("Analytic method", style = "color:black")),
                      p(" ", style = "margin-bottom: +15px;"),
                      prettyRadioButtons("chooseMethod",label = NULL, 
                                         icon = icon("check"),
                                         animation = "jelly", c("Imai Method"), 
                                         selected = "Imai Method", width = '70%'))
                  })
                  
                }else{
                  shinyjs::hide("covariates_variables")
                  
                  output$chooseMethodType <- renderUI({
                    tagList(
                      p(" ", style = "margin-top: 25px;"),
                      h4(strong("Analytic method", style = "color:black")),
                      p(" ", style = "margin-bottom: +15px;"),
                      prettyRadioButtons("chooseMethod",label = NULL, 
                                         icon = icon("check"),
                                         animation = "jelly", c("Imai Method"), 
                                         selected = "Imai Method", width = '70%'))
                  })
                }
              })
            }
          }
        })
        
      })
    })
    
    ##### BETA DIVERSITY #####
    
    out_taxa_options <- surv.pri.func(chooseData$sam.dat, is.mon.sin.rev.bin.con(chooseData$sam.dat))[[2]]
    censor_taxa_options <- surv.pri.func(chooseData$sam.dat, is.mon.sin.rev.bin.con(chooseData$sam.dat))[[1]]
    
    ## Outcome Variable UI
    output$beta.outvars <- renderUI({
      tagList(
        h4(strong("Survival time", style = "color:black")),
        p(" ", style = "margin-bottom: +15px;"),
        selectInput("beta.outvar", label = NULL,
                    c("Choose one" = "", out_taxa_options), "T1Dweek", width = '80%'),
        p(" ", style = "margin-bottom: -10px;"),
        p("A variable for follow-up time (required).",
          style = "font-size:10pt")
      )
    })
    
    output$beta.censorvar <- renderUI({
      tagList(
        p(" ", style = "margin-top: 25px;"),
        h4(strong("Censored/event", style = "color:black")),
        p(" ", style = "margin-bottom: +15px;"),
        selectInput("beta.censorvar", label = NULL,
                    c("Choose one" = "", censor_taxa_options), "T1D", width = '80%'),
        p(" ", style = "margin-bottom: -10px;"),
        p("A binary indicator for censored (0) or event (1) (required).",
          style = "font-size:10pt")
      )
    })
    
    observeEvent(input$beta.outvar, {
      observeEvent(input$beta.censorvar, {
        
        output$beta.treatvar <- renderUI({
          tagList(
            p(" ", style = "margin-top: 25px;"),
            h4(strong("Treatment variable", style = "color:black")),
            p(" ", style = "margin-bottom: +15px;"),
            
            selectInput("beta.treatvar", label = NULL,
                        c("Choose one" = "", select.outcome.func(chooseData$sam.dat, colnames(chooseData$sam.dat), c(input$beta.censorvar, input$beta.outvar))), selected = "Antibiotics", width = '80%'),
            p(" ", style = "margin-bottom: -10px;"),
            p("Environmental, behavioral or medical exposures (e.g., diet, residence, smoking, preterm birth, delivery mode, antibiotic/probiotic use).",
              style = "font-size:10pt"))
        })
        
        observeEvent(input$beta.treatvar, {
          beta.treat_type <- is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.treatvar) 
          
          if(beta.treat_type == "Binary"){
            beta.categos$cat1 = beta.bin.cat.func(chooseData$sam.dat, input$beta.treatvar)[1]
            beta.categos$cat2 = beta.bin.cat.func(chooseData$sam.dat, input$beta.treatvar)[2]
            
            output$beta.rename.treatvar <- renderUI({
              tagList(
                p(" ", style = "margin-top: 25px;"),
                h4(strong("Rename treatment variable?", style = "color:black")),
                p(" ", style = "margin-bottom: +15px;"),
                textInput("beta.rename_bin_var1", label = (paste0("Reference: ",beta.categos$cat1)), value = beta.categos$cat1, width = '80%'),
                textInput("beta.rename_bin_var2", label = (paste0("Comparison: ",beta.categos$cat2)), value = beta.categos$cat2, width = '80%'),
                p(" ", style = "margin-bottom: -10px;"),
                p("You can rename the categories of treatment variable. MiMedSurv keeps up to 8 characters on graphs.", style = "font-size:10pt"))
            })
          }else if (beta.treat_type == "Continuous"){
            output$beta.rename.treatvar <- renderUI({
              tagList(
                p(" ", style = "margin-top: 25px;"),
                h4(strong("Rename treatment variable?", style = "color:black")),
                p(" ", style = "margin-bottom: +15px;"),
                textInput("beta.rename_con_var1", label = NULL, value = input$beta.treatvar_cross, width = '80%'),
                p(" ", style = "margin-bottom: -10px;"),
                p("You can rename the outcome variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:10pt"))
            })
          }
        
        ntselected.treat_vars = select.outcome.func(chooseData$sam.dat, chooseData$treat_vars, input$beta.treatvar)
      }) 
    }) 
    
    output$beta.outcome_types <- renderUI({
      tagList(
        uiOutput("beta.covariates"),
        uiOutput("beta.chooseMethodType"),
        actionButton("runbtn_beta", (strong("Run!")), class = "btn-info",
                     style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50")
      )
    })
    
    observeEvent(input$beta.outvar,{
      observeEvent(input$beta.treatvar, {
        
        observeEvent(input$beta.censorvar, {
        
        # User selects whether to use rarefied or non rarefied biom data for outcome variable
        if (input$beta.outvar %in% chooseData$out_vars) {
          
          ntselected.vars.beta <<- select.covariates.func.surv(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.treatvar, input$beta.outvar, input$beta.censorvar)
          
          beta.out.result = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.outvar)
          
          if (beta.out.result == "Binary") {
            
            ## Covariates UI
            output$beta.covariates <- renderUI({
              tagList(
                p(" ", style = "margin-top: 25px;"),
                h4(strong("Covariate(s)", style = "color:black")),
                p(" ", style = "margin-bottom: +15px;"),
                prettyRadioButtons("beta.covariate",label = NULL, 
                                   icon = icon("check"),
                                   animation = "jelly", c("None", "Covariate(s)"), selected = "None", width = '70%'),
                p(" ", style = "margin-bottom: -10px;"),
                p("Potential confounders (e.g., age, gender) to be adjusted for.", style = "font-size:10pt"),
                
                shinyjs::hidden(
                  shiny::div(id = "beta.covariates_variables", style = "margin-left: 2%",
                             prettyCheckboxGroup("beta.covariateOptions"," Please select covariate(s)", 
                                                 choices = ntselected.vars.beta, width = '70%'))))
            })
            
            ## Covariates Effect
            observeEvent(input$beta.covariate,{
              if (input$beta.covariate == "Covariate(s)") {
                shinyjs::show("covariates_variables")
              } 
              else if (input$beta.covariate == "None") {
                shinyjs::hide("covariates_variables")
              }
            })
            
            ## Mediation Analysis Type UI
            output$beta.chooseMethodType <- renderUI({
              tagList(
                p(" ", style = "margin-top: 25px;"),
                h4(strong("Analytic method", style = "color:black")),
                p(" ", style = "margin-bottom: +15px;"),
                prettyRadioButtons("beta.chooseMethod",label = NULL, icon = icon("check"),
                                   animation = "jelly", c("DACT"), 
                                   selected = "DACT", width = '70%'))
            })
          }else if (beta.out.result == "Continuous") {
            
            ## Covariates UI
            output$beta.covariates <- renderUI({
              tagList(
                p(" ", style = "margin-top: 25px;"),
                h4(strong("Covariate(s)", style = "color:black")),
                p(" ", style = "margin-bottom: +15px;"),
                prettyRadioButtons("beta.covariate",label = NULL, 
                                   icon = icon("check"),
                                   animation = "jelly", c("None", "Covariate(s)"), selected = "None", width = '70%'),
                p(" ", style = "margin-bottom: -10px;"),
                p("Potential confounders (e.g., age, gender) to be adjusted for.", style = "font-size:10pt"),
                
                shinyjs::hidden(
                  shiny::div(id = "beta.covariates_variables", style = "margin-left: 2%",
                             prettyCheckboxGroup("beta.covariateOptions"," Please select covariate(s)", 
                                                 choices = ntselected.vars.beta, width = '70%'))))
            })
            
            ## Covariates Effect
            observeEvent(input$beta.covariate,{
              if (input$beta.covariate == "Covariate(s)") {
                
                shinyjs::show("beta.covariates_variables")
                
                ## Mediation Analysis Type UI
                output$beta.chooseMethodType <- renderUI({
                  tagList(
                    p(" ", style = "margin-top: 25px;"),
                    h4(strong("Analytic method", style = "color:black")),
                    p(" ", style = "margin-bottom: +15px;"),
                    prettyRadioButtons("beta.chooseMethod",label = NULL, 
                                       icon = icon("check"),
                                       animation = "jelly", c("DACT"), 
                                       selected = "DACT", width = '70%'))
                })
              }else{
                shinyjs::hide("beta.covariates_variables")
                output$beta.chooseMethodType <- renderUI({
                  tagList(
                    p(" ", style = "margin-top: 25px;"),
                    h4(strong("Analytic method", style = "color:black")),
                    p(" ", style = "margin-bottom: +15px;"),
                    prettyRadioButtons("beta.chooseMethod",label = NULL, 
                                       icon = icon("check"),
                                       animation = "jelly", c("DACT"), 
                                       selected = "DACT", width = '70%'))
                })
              }
            })
          }
        }
      })
      })
    })
    })
  })
  
  observeEvent(chooseData$taxa.out, {
    
    ######################################
    ######################################
    ######## Taxonomic Analysis ##########
    ######################################
    ###################################### 
    
    out_taxa_options <- surv.pri.func(chooseData$sam.dat, is.mon.sin.rev.bin.con(chooseData$sam.dat))[[2]]
    censor_taxa_options <- surv.pri.func(chooseData$sam.dat, is.mon.sin.rev.bin.con(chooseData$sam.dat))[[1]]
    
    output$taxa_outcome <- renderUI({
      tagList(
        
        h4(strong("Survival time", style = "color:black; font-size:12.5pt")), 
        selectInput("outcome_taxa", label = NULL, 
                    choices = sort(out_taxa_options), selected = "T1Dweek", width = '80%'),
        p(" ", style = "margin-bottom: -10px;"),
        p("A variable for follow-up time (required).", style = "font-size:10pt"),
        p(" ", style = "margin-bottom: +22px;"),
      )
    })
    
    output$taxa_censor <- renderUI({ 
      tagList(
        h4(strong("Censored/event", style = "color:black; font-size:12.5pt")), 
        selectInput("censor_taxa", label = NULL, 
                    choices = sort(censor_taxa_options), selected = "T1D", width = '80%'),
        p(" ", style = "margin-bottom: -10px;"),
        p("A binary indicator for censored (0) or event (1) (required).", style = "font-size:10pt"),
        p(" ", style = "margin-bottom: +22px;"),
      )
    })
    
    taxa.types$dataType = "clr"
    
    observeEvent(input$outcome_taxa, {
      observeEvent(input$outcome_taxa, {
        output$taxa_treat <- renderUI({
          
          treat.option <- select.outcome.func(chooseData$sam.dat, colnames(chooseData$sam.dat), c(input$censor_taxa, input$outcome_taxa))
          
          tagList(
            h4(strong("Treatment variable", style = "color:black; font-size:12.5pt")),
            selectInput("treat_taxa", label = NULL,
                        choices = sort(treat.option), selected = "Antibiotics", width = '80%'),
            p(" ", style = "margin-bottom: -10px;"),
            p("Environmental, behavioral or medical exposures (e.g., diet, residence, smoking, preterm birth, delivery mode, antibiotic/probiotic use).", style = "font-size:10pt"), 
            p(" ", style = "margin-bottom: +30px;"),
          ) 
        })
        
        observeEvent(input$treat_taxa, {
          ntselected.treat_vars = cov.func_2(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, c(input$treat_taxa, input$outcome_taxa, input$censor_taxa))
          out.ind <- which(names(chooseData$sam.dat) == input$outcome_taxa)
          
          
          output$taxa_cov <- renderUI({
            tagList(
              h4(strong("Covariate(s)", style = "color:black; font-size:12.5pt")),
              p(" ", style = "margin-bottom: -10px;"),
              prettyRadioButtons("cov_taxa",label = NULL, status = "primary", icon = icon("check"),
                                 animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
              p(" ", style = "margin-bottom: -10px;"),
              p("Potential confounders (e.g., age, gender) to be adjusted for.", style = "font-size:10pt"),
              
              shinyjs::hidden(
                shiny::div(id = "covariates_t", style = "margin-left: 2%",
                           prettyCheckboxGroup("covariates_taxa"," Please select covariate(s)", status = "primary",
                                               choices = sort(ntselected.treat_vars), width = '70%'))),
              
              p(" ", style = "margin-bottom: +25px;"),
              
              h4(strong("Taxonomic ranks", style = "color:black; font-size:12.5pt")),
              p(" ", style = "margin-bottom: -10px;"),
              prettyRadioButtons("include_species_taxa", label = NULL, status = "primary", animation = "jelly",
                                 c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                 icon = icon("check"), width = '80%'),
              p(" ", style = "margin-bottom: -10px;"),
              p("Taxonomic ranks to be analyzed. ‘Phylum - Genus’ (default) can be used for 16S data; ‘Phylum - Species’ can be used for shotgun metagenomic data. ", style = "font-size:10pt"),
              p(" ", style = "margin-bottom: +30px;")
              
            )
          })
        })
      })
    }) 
    
    observeEvent(input$cov_taxa,{
      if (input$cov_taxa == "Covariate(s)") {
        shinyjs::show("covariates_t")
      } 
      else if (input$cov_taxa == "None") {
        shinyjs::hide("covariates_t")
      }
    })
    
    observeEvent(input$cov_taxa, {
      observeEvent(input$outcome_taxa, {
        output$taxa_method <- renderUI({
          tagList(
            h4(strong("Analytic method", style = "color:black; font-size:12.5pt")),
            prettyRadioButtons("choose_taxa", label = NULL,
                               choices = c("Imai Method"), selected = "Imai Method", icon = icon("check"), width = '80%')
          )
        })
      })
    })
    
    observeEvent(input$choose_taxa, {
      observeEvent(input$outcome_taxa, {
        output$taxa_reg_med <- renderUI({
          tagList(
            actionButton("runbtn_taxa", (strong("Run!")), class = "btn-info",
                         style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50")
          )
        })
      })
    })
  }) 
  
  ######################################
  ##########  RUN BUTTONS   ############
  ######################################
  
  ##### QC #####
  observeEvent(input$run, {
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        
        incProgress(1/10, message = "Data Trimming in progress")
        
        if (nchar(input$part.rem.str) == 0) {
          rem.tax.complete <- rem.tax.d
          rem.tax.partial <- rem.tax.str.d
        } else {
          rem.tax.complete <- unique(c(unlist(strsplit(input$rem.str, split = ",")), rem.tax.d))
          rem.tax.partial <- unique(c(unlist(strsplit(input$part.rem.str, split = ",")), rem.tax.str.d))
        }
        
        tax.tab <- NULL 
        try(tax.tab <- tax_table(infile$biom), silent = TRUE) 
        
        tree <- NULL
        try(tree <- phy_tree(infile$biom), silent = TRUE) 
        
        if(!is.null(tax.tab)){
          if (input$kingdom != "all") {
            ind <- is.element(tax.tab[,1], input$kingdom)
            validate(
              if (sum(ind) == 0) {
                showNotification(h4(paste("Error: Please select valid Kingdom. Available kingdoms are:", 
                                          paste(c(na.omit(unique(tax.tab[,1])) ,"and all"), collapse = ", "))),
                                 type = "error")
              } else {
                NULL
              }
            )
          }
        }
        
        shinyjs::disable("run")
        shinyjs::disable("slider1")
        shinyjs::disable("slider2")
        shinyjs::disable("kingdom")
        shinyjs::disable("skip")
        shinyjs::disable("binwidth")
        shinyjs::disable("binwidth2")
        
        rcol$selected = "rgba(255, 0, 0, 0.6)"
        
        if (is.null(tax.tab) & !is.null(tree)){
          
          infile$qc_biom = biom.clean.no.tax.tab(infile$biom, 
                                                 input$kingdom, 
                                                 lib.size.cut.off = input$slider1, 
                                                 mean.prop.cut.off = input$slider2/100,
                                                 rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial)
          
          infile$qc_biomNA = biom.cleanSNA.no.tax.tab(infile$biom,
                                                      input$kingdom,
                                                      lib.size.cut.off = input$slider1,
                                                      mean.prop.cut.off = input$slider2/100,
                                                      rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
                                                      surv.time = chooseData$surv.Time.select,
                                                      censor = chooseData$censor.select,
                                                      follow = chooseData$follow.time)
          
        }
        
        else if (!is.null(tax.tab) & is.null(tree)){
    
          infile$qc_biom = biom.clean.no.tree(infile$biom, 
                                              input$kingdom, 
                                              lib.size.cut.off = input$slider1, 
                                              mean.prop.cut.off = input$slider2/100,
                                              rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial)
          
          
          infile$qc_biomNA = biom.cleanSNA.no.tree(infile$biom,
                                                   input$kingdom,
                                                   lib.size.cut.off = input$slider1,
                                                   mean.prop.cut.off = input$slider2/100,
                                                   rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
                                                   surv.time = chooseData$surv.Time.select,
                                                   censor = chooseData$censor.select,
                                                   follow = chooseData$follow.time)
          
          
        }else if (is.null(tax.tab) & is.null(tree)){
          
          infile$qc_biom = biom.clean.no.both(infile$biom, 
                                              input$kingdom, 
                                              lib.size.cut.off = input$slider1, 
                                              mean.prop.cut.off = input$slider2/100,
                                              rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial)
          
          infile$qc_biomNA = biom.cleanSNA.no.both(infile$biom,
                                                   input$kingdom,
                                                   lib.size.cut.off = input$slider1,
                                                   mean.prop.cut.off = input$slider2/100,
                                                   rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
                                                   surv.time = chooseData$surv.Time.select,
                                                   censor = chooseData$censor.select,
                                                   follow = chooseData$follow.time)
        }else{
          infile$qc_biom = biom.clean(infile$biom, 
                                      input$kingdom, 
                                      lib.size.cut.off = input$slider1, 
                                      mean.prop.cut.off = input$slider2/100,
                                      rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial)
          
          infile$qc_biomNA = biom.cleanSNA(infile$biom,
                                           input$kingdom,
                                           lib.size.cut.off = input$slider1,
                                           mean.prop.cut.off = input$slider2/100,
                                           rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
                                           surv.time = chooseData$surv.Time.select,
                                           censor = chooseData$censor.select,
                                           follow = chooseData$follow.time)
        }
        
        
        incProgress(3/10, message = "Rarefying in progress")
        lib_size.sum = lib.size.func(infile$qc_biom)$lib.size.sum
        infile$rare_biom = rarefy.func(infile$qc_biom, 
                                       cut.off = lib_size.sum["Minimum"],
                                       multi.rarefy = 1)
        
        infile$rare_biomNA = rarefy.func(infile$qc_biomNA, 
                                         cut.off = lib_size.sum["Minimum"],
                                         multi.rarefy = 1)
        
        incProgress(2/10, message = "Saving File in progress")
        
        chooseData$sam.dat = sample_data(infile$qc_biom)
        chooseData$mon.sin.rev.bin.con = is.mon.sin.rev.bin.con(chooseData$sam.dat)
        chooseData$treat_vars = extract.bin.con.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con)
        chooseData$out_vars = extract.bin.con.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con)
        
        if (!is.null(input$taxTable)){
          chooseData$tax.tab = tax_table(infile$rare_biom)
          chooseData$tax.tabNA = tax_table(infile$rare_biomNA)
        }
        
        output$moreControls <- renderUI({
          tagList(
            box(title = strong("Download Data", style = "color:white"), width = NULL, status = "primary", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:13pt"),
                h5("Data after Quality Control"),
                downloadButton("downloadData2", "Download", width = '50%', 
                               style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50"),br(),
                h5("Data after Quality Control and Rarefaction"),
                downloadButton("downloadData3", "Download", width = '50%', 
                               style="color: #000000; background-color: #FFFFFF; border-color: #ff7f50"),br(),
                p("For your reference, you can download the data files above for the phyloseq object (biom.after.qc) after QC and
                      (rare.biom.after.qc) after QC and rarefaction.",
                  style = "font-size:11pt")
            )
          )
        })
        
        output$text <- renderText({"You are all set! You can proceed to data analysis!"})
        
        qc_biom = infile$qc_biom
        output$downloadData2 <- downloadHandler(
          filename = function() {
            paste("biom.after.qc.Rdata")
          },
          content = function(file1) {
            save(qc_biom, file = file1)
          })
        
        rare_biom = infile$rare_biom
        output$downloadData3 <- downloadHandler(
          filename = function() {
            paste("rare.biom.after.qc.Rdata")
          },
          content = function(file1) {
            save(rare_biom, file = file1)
          })
        
        final_dat <<- infile$rare_biom 
        final_dat_qc <<- infile$qc_biom 
        
        incProgress(1/10, message = "Done")
        shinyjs::enable("run")
        shinyjs::enable("slider1")
        shinyjs::enable("slider2")
        shinyjs::enable("kingdom")
        shinyjs::enable("skip")
        shinyjs::enable("binwidth")
        shinyjs::enable("binwidth2")
      })
  })
  
  ##### DIVERSITY CALCULATION #####
  
  observeEvent(input$divCalcRun, {
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        shinyjs::disable("divCalcRun")
        
        rare.otu.tab <- otu_table(infile$rare_biom)
        
        incProgress(3/10, message = "Calculating Alpha Diversity")
        
        Tree <<- NULL
        
        try(Tree <<- phy_tree(infile$rare_biom), silent = TRUE)
        
        Is.tree <<- NULL
        
        if (is.null(Tree)){
          Is.tree <<- "withoutTree"
        } else {
          Is.tree <<- "withTree"
        }
        
        if (Is.tree == "withTree") {
          chooseData$alpha.div.rare = alpha.v1.func(infile$rare_biom)
        } else if (Is.tree == "withoutTree") {
          chooseData$alpha.div.rare = alpha.v1.func.no.tree(infile$rare_biom)
        }
        
        chooseData$alpha.div = chooseData$alpha.div.rare
        incProgress(3/10, message = "Calculating Beta Diversity")
        ds.Ks$res = Ds.Ks.func(infile$rare_biom, infile$qc_biom, Is.tree)
        
        output$divCalcDownload <- renderUI({
          tagList(
            box(title = strong("Download Data", style = "color:white"), width = NULL, status = "primary", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:15pt"),
                p("You can download alpha- and beta-diversity data.",
                  style = "font-size:11pt"), 
                h5("Alpha Diversity"),
                downloadButton("alphaDiv", "Download", width = '50%'),br(),
                h5("Beta Diversity"),
                downloadButton("betaDiv", "Download", width = '50%', style = "color:black; background-color: white"),
                br(), br(), 
                p("You can download the original count and normalized taxonomic data.",
                  style = "font-size:11pt"), 
                h5("Count (Original)"), 
                downloadButton("taxadataCount", "Download", width = '50%', style = "color:black; background-color: red3"), br(), 
                h5("CLR"),
                downloadButton("taxadataCLR", "Download", width = '50%', style = "color:black; background-color: red3"), br(), br()
            )
          )
        })
        
        alpha.div = chooseData$alpha.div
        
        output$alphaDiv <- downloadHandler(
          filename = function() {
            paste("Alpha.Diversity.txt")
          },
          content = function(alpha.file) {
            write.table(chooseData$alpha.div, file = alpha.file, row.names = TRUE, col.names = TRUE, sep = "\t")
          })
        
        output$betaDiv <- downloadHandler(
          filename = function() {
            paste("Beta.Diversity.zip")
          },
          content <- function(fname) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            
            if (Is.tree == "withTree") {
              
              dataFiles = c("Jaccard.txt", "Bray.Curtis.txt", "U.UniFrac.txt" ,"G.UniFrac.txt", "W.UniFrac.txt")
              
              write.table(as.data.frame(ds.Ks$res$Ks$Jaccard), file = "Jaccard.txt", 
                          row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(ds.Ks$res$Ks$Bray.Curtis), file = "Bray.Curtis.txt", 
                          row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(ds.Ks$res$Ks$U.UniFrac), file = "U.UniFrac.txt", 
                          row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(ds.Ks$res$DKs$G.UniFrac), file = "G.UniFrac.txt", 
                          row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(ds.Ks$res$Ks$W.UniFrac), file = "W.UniFrac.txt", 
                          row.names = TRUE, col.names = TRUE, sep = "\t")
              
              zip(zipfile=fname, files=dataFiles) 
            } 
            else if (Is.tree == "withoutTree") {
              
              dataFiles = c("Jaccard.txt", "Bray.Curtis.txt")
              
              write.table(as.data.frame(ds.Ks$res$Ks$Jaccard), file = "Jaccard.txt", 
                          row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(ds.Ks$res$Ks$Bray.Curtis), file = "Bray.Curtis.txt", 
                          row.names = TRUE, col.names = TRUE, sep = "\t")
              
              zip(zipfile=fname, files=dataFiles) 
            }
          }
        )
        
        incProgress(3/10, message = "Transformation in progress")
        
        rare.otu.tab <- otu_table(infile$rare_biom)
        rare.tax.tab <- tax_table(infile$rare_biom)
        no.rare.otu.tab <- otu_table(infile$qc_biom)
        no.rare.tax.tab <- tax_table(infile$qc_biom)
        
        rare.otu.tabNA <- otu_table(infile$rare_biomNA)
        rare.tax.tabNA <- tax_table(infile$rare_biomNA)
        no.rare.otu.tabNA <- otu_table(infile$qc_biomNA)
        no.rare.tax.tabNA <- tax_table(infile$qc_biomNA)
        
        chooseData$taxa.out <- tax.trans(no.rare.otu.tab, no.rare.tax.tab, rare.otu.tab, rare.tax.tab)
        # Subcompositional Invariant 
        original_taxa <<- chooseData$taxa.out 
        check_to_see_taxa <<- tax.trans.subcominv(no.rare.otu.tab, no.rare.tax.tab, rare.otu.tab, rare.tax.tab) 
        
        chooseData$taxa.outNA = tax.trans.na(no.rare.otu.tabNA, no.rare.tax.tabNA, rare.otu.tabNA, rare.tax.tabNA)
        chooseData$taxa.names.out <- taxa.names.rank(chooseData$taxa.out[[1]])
        chooseData$tax.tab <- rare.tax.tab
        chooseData$tax.tabNA = rare.tax.tabNA
        
        chooseData$NAadded <- add_NA(chooseData$taxa.outNA, chooseData$tax.tabNA)
        
        incProgress(3/10, message = "Transformation in progress")
        
        count_biom = chooseData$taxa.out$count
        rare_biom = chooseData$taxa.out$rare.count
        prop_biom = chooseData$taxa.out$prop
        clr_biom = chooseData$taxa.out$clr
        arc_biom = chooseData$taxa.out$arcsin
        
        output$taxadataCount <- downloadHandler(
          
          filename = function() {
            paste("Count.Data.zip")
          },
          content = function(count.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(count_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=count.file, files=dataFiles)
          }
        )
        
        incProgress(3/10, message = "Saving")
        
        output$taxadataCLR <- downloadHandler(
          filename = function() {
            paste("CLR.Transformed.Data.zip")
          },
          content = function(clr.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(clr_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=clr.file, files=dataFiles)
          }
        )
        
        output$taxadataArc <- downloadHandler(
          filename = function() {
            paste("Arcsine.Transformed.Data.zip")
          },
          content = function(clr.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(arc_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=clr.file, files=dataFiles)
          }
        )
        
        incProgress(1/10, message = "Done")
        shinyjs::enable("divCalcRun")
      })
  })
  
  ##### Non-Mediational Analysis ##### 
  
  observeEvent(input$runbtn_nonmed, {
    
    validate(
      if ((input$non.covariate == "Covariate(s)" ) & (is.null(input$non.covariateOptions))) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("runbtn_nonmed")
    
    withProgress({
      if(input$non.covariate == "Covariate(s)"){
        sam_dat0 <<- tryCatch(sam_no_miss_cov(chooseData$sam.dat, input$non.covariateOptions), error = function(e) {  
          message("No outcome is available!")
          showModal(modalDialog(div("No outcome is available!")))
          return(NULL)})
      }else{
        sam_dat0 <<-chooseData$sam.dat
      }
      
      incProgress(3/10, message = "Displaying Variables in Use")
      
      non.med.cat.treat <- c(non.cat.1, non.cat.2)
      non.med.bin.cat.treat <- c(input$non.med.rename_bin_var1, input$non.med.rename_bin_var2)
      
      sam_dat0 <<- beta.bin.recode.func(sam_dat0, input$non.treatvar, non.med.cat.treat, input$non.med.rename_bin_var1, input$non.med.rename_bin_var2)
    
      surv.dat0 <<- tryCatch(surv.events.dat.func(sam_dat0, 
                                                 input$non.treatvar, 
                                                 input$non.outvar, 
                                                 input$non.censorvar), error = function(e) {  
                                                   message("No outcome is available!")
                                                   showModal(modalDialog(div("No outcome is available!")))
                                                   return(NULL)})
      
      
      num.treat1_size <<- tryCatch(surv.num.treat1(surv.dat0, "treat"), error = function(e) {  
        message("No outcome is available!")
        showModal(modalDialog(div("No outcome is available!")))
        return(NULL)})
      num.treat2_size <<- tryCatch(surv.num.treat2(surv.dat0, "treat"), error = function(e) {  
        message("No outcome is available!")
        showModal(modalDialog(div("No outcome is available!")))
        return(NULL)})
      
      num.event <<- tryCatch(surv.num.uncen(surv.dat0), error = function(e) {  
        message("No outcome is available!")
        showModal(modalDialog(div("No outcome is available!")))
        return(NULL)})
      num.cen <<- tryCatch(surv.num.censr(surv.dat0), error = function(e) {  
        message("No outcome is available!")
        showModal(modalDialog(div("No outcome is available!")))
        return(NULL)})
      
      to_check <<-  unique(sam_dat0[[input$non.treatvar]] )
      
      output$treatment_size <- renderValueBox({
        treat.ls <- unique(sam_dat0[[input$non.treatvar]] )
        
        valueBox(
          value = tags$p(paste0(num.treat1_size), style = "font-size: 75%;"),
          paste("Sample Size (", treat.ls[1],")", sep="" ), icon = icon("user-circle"), color = "red")
      })
      
      output$control_size <- renderValueBox({
        treat.ls <- unique( sam_dat0[[input$non.treatvar]] )
        
        valueBox(
          value = tags$p(paste0(num.treat2_size), style = "font-size: 75%;"),
          paste("Sample Size (", treat.ls[2],")", sep="" ), icon = icon("user-circle"), color = "orange")
      })
      
      output$censored_size <- renderValueBox({
        
        valueBox(
          value = tags$p(paste0(num.cen), style = "font-size: 75%;"),
          "Sample Size (Censored)", icon = icon("user-circle"), color = "yellow")
      })
      
      output$uncensored_size <- renderValueBox({
        
        valueBox(
          value = tags$p(paste0(num.event), style = "font-size: 75%;"),
          "Sample Size (Event)", icon = icon("user-circle"), color = "teal")
      })
      
      validate(
        if (num.treat1_size < 2 & num.treat2_size < 2 & num.event < 2 ) {
          if( num.treat1_size < 2 & num.treat2_size < 2 ){
            showNotification("Error: Too few sample Treatment size. Select different Treatment variable.",
                             type = "error")
          }
          if( num.event < 2 ){
            showNotification("Error: Too few sample Event size. Select different Event variable.",
                             type = "error")
          }
        }
      )
      
      #########################################################
      incProgress(6/10, message = "Calculating")
      
      if (input$non.covariate == "None") {
        
        surv.dat <<- tryCatch(surv.events.dat.func(sam_dat0, 
                                                  input$non.treatvar, 
                                                  input$non.outvar, 
                                                  input$non.censorvar), error = function(e) {  
                                                    message("No outcome is available!")
                                                    showModal(modalDialog(div("No outcome is available!")))
                                                    return(NULL)})
        
        if (input$non.chooseMethod == "Cox model") {
          cox.fit <<- tryCatch(surv.Cox.fit.func1(surv.dat), error = function(e) {  
            message("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)})
          
          
          output$non.display_results = renderUI({
            tagList(
              box(title = strong( "Survival Curve", style = "color:white"), 
                  align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                  plotOutput("surv.plot2", height = 400, width = 500)
              )
            )
          })
          
          output$surv.plot2 = renderPlot({
            isolate(surv.Cox.plot.func(cox.fit, surv.dat))
          })
        }
        else if (input$non.chooseMethod == "Kaplan-Meier analysis (default)") {
          
          surv.KM.fit <<- surv.KM.fit.func( surv.dat )
          
          output$non.display_results = renderUI({
            tagList(
              box(title = strong( "Kaplan-Meier Curve", style = "color:white"), 
                  align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                  plotOutput("surv.plot3", height = 400, width = 500)
              )
            )
          })
          output$surv.plot3 = try(renderPlot({
            isolate(surv.KM.plot.func( surv.KM.fit, surv.dat, input$non.test ))
          }), silent = TRUE) 
        }} 
      else {
        surv.dat <<-tryCatch(surv.events.dat.func2(sam_dat0, 
                                                  input$non.treatvar, 
                                                  input$non.outvar, 
                                                  input$non.censorvar,
                                                  input$non.covariateOptions), error = function(e) {  
                                                    message("No outcome is available!")
                                                    showModal(modalDialog(div("No outcome is available!")))
                                                    return(NULL)}) 
        
        cox.fit <<- tryCatch(surv.Cox.fit.func2(surv.dat, input$non.covariateOptions), error = function(e) {  
          message("No outcome is available!")
          showModal(modalDialog(div("No outcome is available!")))
          return(NULL)})
        
        output$non.display_results = renderUI({
          tagList(
            box(title = strong( "Adjusted Survival Curve", style = "color:white"), 
                align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                plotOutput("surv.plot1", height = 400, width = 500)
            )
          )
        })

        dat_3 <<- input$non.covariateOptions
        
        output$surv.plot1 = renderPlot({
          isolate(surv.Cox.plot.func(cox.fit, surv.dat, input$non.covariateOptions))
        })
      }
      
      if(input$non.chooseMethod == "Cox model"){
        ref_string = REFERENCE_CHECK(method_name = isolate(input$non.chooseMethod) )
        
      } else if (input$non.chooseMethod == "Kaplan-Meier analysis (default)") {
        ref_string = REFERENCE_CHECK(method_name = isolate(input$non.chooseMethod), sig_test = isolate(input$non.test) )
      }
      
      if (is.null(ref_string)) {
        shinyjs::hide("non.references")
      } else {
        shinyjs::show("non.references")
        output$non.references = renderUI({
          tagList(
            box(title = strong("References"), width = NULL, status = "primary", solidHeader = TRUE,
                HTML(paste(ref_string, collapse="<br/>"))
            )
          )
        })
      }
      shinyjs::enable("runbtn_nonmed")
    }) 
  }) 
  ##### Mediational Analysis ##### 
  
  ##### ALPHA DIVERSITY #####
  
  observeEvent(input$runbtn_alpha,{
    
    shinyjs::disable("runbtn_alpha") 
    shinyjs::disable("treatvars")
    shinyjs::disable("outvars")
    shinyjs::disable("censorvar")
    shinyjs::disable("outcome_types")
    shinyjs::disable("covariates")
    shinyjs::disable("covariatesOptions")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        incProgress(0.05, message = "Data Processing")
        
        ## Variables
        Treatment = input$treatvar
        Event = input$censorvar
        Outcome = input$outvar
        Interaction = "FALSE"
        Is.covariates = input$covariate
        Covariates = input$covariateOptions
        Covariates.regression = ifelse(is.null(input$covariateOptions), NA, paste0(input$covariateOptions, collapse = " + "))
        Method = input$chooseMethod
        
        ## Mediator data
        Mediator.list <<- colnames(chooseData$alpha.div)
        
        ## Data
        Data <- cbind(chooseData$alpha.div, chooseData$sam.dat)
        
        if (sum(is.na(Data[[Outcome]])) >= 1) {
          ind = which(is.na(Data[[Outcome]]))
          Data.na = Data[-ind, ]
        } 
        else {Data.na = Data}
        
        Data.na.stand <- Data.na
        
        for (med in Mediator.list) {
          Data.na.stand[[med]] = (Data.na[[med]] - mean(Data.na[[med]])) / sd(Data.na[[med]])
        }
        
        cov_check <<- input$covariateOptions 
        sam.dat.check <<- chooseData$sam.dat
        
        if (Is.covariates == "Covariate(s)") {
          for (i in 1:length(input$covariateOptions)) {
            if (is.na(sum(as.numeric(unlist(na.omit(chooseData$sam.dat[, input$covariateOptions[i]])))))){
              Data.na.stand[[input$covariateOptions[i]]] <- as.factor(Data.na.stand[[input$covariateOptions[i]]])
            }else{
              Data.na.stand[[input$covariateOptions[i]]] <- as.numeric(Data.na.stand[[input$covariateOptions[i]]])
            }
        }
        }
        
        Treatment.type <- try(is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$treatvar), silent = TRUE) 

        # (1) Imai Method ---------------------------------------------
        
        if (input$chooseMethod %in% c("Imai Method (Default)", "Imai Method")) {
          
          Med.out <<- tryCatch(Mediation(Data.na.stand, Mediator.list,
                                        Treatment, Event,  Outcome, Is.covariates, Covariates.regression,
                                        Interaction, Treatment.type, 3000),
                              error = function(e) {
                                message("No outcome is available!")
                                showModal(modalDialog(div("No outcome is available!")))
                                return(NULL)})
          
          
          ## Visualization
          incProgress(0.05, message = "Displaying Results")
          
          if (Treatment.type == "Binary") {
            output$alpha.display_results = renderUI({
              fluidRow(style = "position:relative",
                       tabBox(width = 12, title = strong("Forest Plot", style = "color:black"),
                              tabPanel("ACME (Mediation)", align = "center",
                                       plotOutput("forestplot_ACME1", height = 500, width = 800)),
                              tabPanel("ADE", align = "center",
                                       plotOutput("forestplot_ADE", height = 500, width = 800)),
                              tabPanel("Total Effect", align = "center",
                                       plotOutput("forestplot_TOTAL", height = 500, width = 800))))
            })
          }
          else if (Treatment.type == "Continuous") {
            
            output$alpha.display_results = renderUI({
              fluidRow(style = "position:relative",
                       tabBox(width = 12, title = strong("Forest Plot", style = "color:black"),
                              tabPanel("ACME (Mediation)", align = "center",
                                       plotOutput("forestplot_ACME1", height = 500, width = 800)),
                              tabPanel("ADE", align = "center",
                                       plotOutput("forestplot_ADE", height = 500, width = 800)),
                              tabPanel("Total Effect", align = "center",
                                       plotOutput("forestplot_TOTAL", height = 500, width = 800))))
            })
          }
          
          output$forestplot_ACME1 <- renderPlot({
            tryCatch(Forestplot(Med.out, Is.tree, "Continuous", Treatment.type, Interaction, effect = "ACME"), 
                     error = function(e) {
                       message("No outcome is available!")
                       showModal(modalDialog(div("No outcome is available!")))
                       return(NULL)})
          })
      
          output$forestplot_ADE <- renderPlot({
            tryCatch(Forestplot(Med.out, Is.tree, "Continuous", Treatment.type, Interaction, effect = "ADE"), 
                     error = function(e) {
                       message("No outcome is available!")
                       showModal(modalDialog(div("No outcome is available!")))
                       return(NULL)})
          })
          
          output$forestplot_TOTAL <- renderPlot({
            tryCatch(Forestplot(Med.out, Is.tree, "Continuous", Treatment.type, Interaction, effect = "Total Effect"), 
                     error = function(e) {
                       message("No outcome is available!")
                       showModal(modalDialog(div("No outcome is available!")))
                       return(NULL)})
          })
        } 
      }
    )
    
    ref_string = REFERENCE_CHECK(method_name = isolate(input$chooseMethod))
    if (is.null(ref_string)) {
      shinyjs::hide("alpha.references")
    } else {
      shinyjs::show("alpha.references")
      output$alpha.references = renderUI({
        tagList(
          box(title = strong("Reference", style = "color:white"), width = NULL, status = "primary", solidHeader = TRUE,
              HTML(paste(ref_string, collapse="<br/>"))
          )
        )
      })
    }
    
    shinyjs::enable("runbtn_alpha")
    shinyjs::enable("treatvars")
    shinyjs::enable("outvars")
    shinyjs::enable("censorvar")
    shinyjs::enable("outcome_types")
    shinyjs::enable("covariates")
    shinyjs::enable("covariatesOptions")
    
  }, ignoreNULL = TRUE, ignoreInit = TRUE) 
  
  ##### BETA DIVERSITY #####
  
  observeEvent(input$runbtn_beta,{
    
    shinyjs::disable("runbtn_beta")
    shinyjs::disable("beta.outvars")
    shinyjs::disable("beta.censorvar")
    shinyjs::disable("beta.treatvar")
    shinyjs::disable("beta.outcome_types")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        incProgress(1/10, message = "Data Processing")
        
        ## Variables
        Treatment = input$beta.treatvar
        Outcome = input$beta.outvar
        Censor = input$beta.censorvar
        Method = input$beta.chooseMethod
        Covariates = input$beta.covariateOptions
        Ds.Ks = ds.Ks$res
        
        Treatment.type = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, Treatment)
        Outcome.type =  is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, Outcome)
        
        ## Data
        Data <- chooseData$sam.dat
        
        ## Remove rows with missing values on treatment and outcome
        ind1 <- c() 
        ind2 <- c()
        
        if (sum(is.na(Data[[Treatment]])) >= 1) {
          ind1 = which(is.na(Data[[Treatment]]))
        }
        if (sum(is.na(Data[[Outcome]])) >= 1) {
          ind2 = which(is.na(Data[[Outcome]]))
        }
        
        ind.missing <<- unique(c(ind1, ind2))
        
        if (length(ind.missing) == 0) {
          ind.missing = "NA"
          Data.na <<-  Data
        } else {
          Data.na <<- Data[-ind.missing, ]
        }
        
        ## Mediator data
        beta.Mediator <<- Ds.Ks$Ks
        beta.Mediator.list <<- names(beta.Mediator)
        
        if (Is.tree == "withTree") {
          if ("NA" %in% ind.missing) {
            beta.Mediator.na <<- beta.Mediator
          } else {
            beta.Mediator.na <<- list(JAC = beta.Mediator[[1]][-ind.missing, -ind.missing],
                                      BC = as.matrix(beta.Mediator[[2]])[-ind.missing, -ind.missing],
                                      UniFrac = beta.Mediator[[3]][-ind.missing, -ind.missing],
                                      GUniFrac = beta.Mediator[[4]][-ind.missing, -ind.missing],
                                      WUniFrac = beta.Mediator[[5]][-ind.missing, -ind.missing])
          }
        } 
        else if (Is.tree == "withoutTree") {
          if ("NA" %in% ind.missing) {
            beta.Mediator.na <<- beta.Mediator
          } else {
            beta.Mediator.na <<- list(JAC = beta.Mediator[[1]][-ind.missing, -ind.missing],
                                      BC = as.matrix(beta.Mediator[[2]])[-ind.missing, -ind.missing])
          }
        }
        
        incProgress(1/10, message = "Renaming")
        
        beta.categors_treat = c(beta.categos$cat1, beta.categos$cat2)
        beta.bin_categos_treat = c(input$beta.rename_bin_var1, input$beta.rename_bin_var2)
        rename.catsbin_ref1 = beta.bin_categos_treat[which(beta.categors_treat == beta.categos$cat1)]
        rename.catsbin_com1 = beta.bin_categos_treat[which(beta.categors_treat != beta.categos$cat1)]
        
        if (Treatment.type == "Binary") {
          beta.bin.ori.cat_treatvar <<- beta.bin.cat.ori.func(Data.na, Treatment)
          beta.Data <<- try(beta.bin.recode.func(Data.na, Treatment, beta.bin.ori.cat_treatvar,
                                                       rename.catsbin_ref1, rename.catsbin_com1), silent = TRUE) 
        } 
        else if (Treatment.type == "Continuous") {
          beta.Data <<- try(beta.con.recode.func(Data.na, Treatment, input$beta.rename_con_var1), silent = TRUE) 
        }
        
        incProgress(1/10, message = "Calculating DACT")
        
        beta_med_model_out <<- beta_med_model(beta.Data, beta.Mediator.na, Treatment, Treatment.type, Covariates)
        beta_out_model_out <<- beta_out_model(beta.Data, beta.Mediator.na, Outcome, Censor, Treatment, Covariates)
        
        dact.result <<- beta_dact(beta_med_model_out, beta_out_model_out, beta.Mediator.list)
        
        if (Treatment.type == "Binary") {
          beta.Treatvar.out <<- try(beta.bin.out.func(beta.Data, Ds.Ks, Treatment,
                                                      rename.catsbin_ref1, rename.catsbin_com1), silent = TRUE) 
          Treatvar.out1 <<- beta.Treatvar.out
        } else {
          beta.Treatvar.out <<- try(beta.con.out.func(beta.Data, Ds.Ks, Treatment,
                                                      input$beta.rename_con_var1, ind.missing), silent = TRUE) 
          Treatvar.out2 <<- beta.Treatvar.out
        }
        
        beta.Outvar.out <<- try(beta.con.out.func(beta.Data, Ds.Ks, Outcome,
                                                  Outcome, ind.missing, input$beta.censorvar), silent = TRUE) 
        Outvar.out2 <<- beta.Outvar.out
        
        incProgress(3/10, message = "Displaying Results")
        
        if (Is.tree == "withTree") {
          
          output$beta.display_results = renderUI({
            fluidRow(style = "position:relative",
                     tabBox(width = 12, title = strong("Principal Coordinate Analysis Plot", style = "color:black"),
                            tabPanel("DACT", align = "center",
                                     column(width = 5, offset = 1,
                                            plotOutput("beta_graph_plot1", height = 800)),
                                     column(width = 5,
                                            plotOutput("beta_graph_plot2", height = 800)))))
          })
          
          output$beta_graph_plot1 <- renderPlot({
            if (Treatment.type == "Binary"){
              try(MedTest.bin.con.plot1(dact.result, beta.Treatvar.out, beta.Outvar.out), silent = TRUE) 
            }else if (Treatment.type == "Continuous"){
              try(MedTest.con.con.plot1(dact.result, beta.Treatvar.out, beta.Outvar.out), silent = TRUE) 
            }
          })
          
          output$beta_graph_plot2 <- renderPlot({
            if (Treatment.type == "Binary"){
              try(MedTest.bin.con.plot2(dact.result, beta.Treatvar.out, beta.Outvar.out), silent = TRUE) 
            }else if (Treatment.type == "Continuous"){
              try(MedTest.con.con.plot2(dact.result, beta.Treatvar.out, beta.Outvar.out), silent = TRUE) 
            }
          })
        }
        else if (Is.tree == "withoutTree") {
          
          output$beta.display_results = renderUI({
            fluidRow(style = "position:relative",
                     tabBox(width = 12, title = strong("Principal Coordinate Analysis Plot", style = "color:black"),
                            tabPanel("MedTest", align = "center",
                                     plotOutput("beta_graph_plot", height = 800, width = 600))))
          })
          
          output$beta_graph_plot <- renderPlot({
            
            if (Treatment.type == "Binary") {
              try(MedTest.bin.con.plot(dact.result, beta.Treatvar.out, beta.Outvar.out), silent = TRUE) 
            } 
            else {
              try(MedTest.con.con.plot(dact.result, beta.Treatvar.out, beta.Outvar.out), silent = TRUE) 
            }
          })
        }
      }
    )
    
    ref_string = REFERENCE_CHECK(method_name = "Beta Analysis")
    if (is.null(ref_string)) {
      shinyjs::hide("beta.references")
    } else {
      shinyjs::show("beta.references")
      output$beta.references = renderUI({
        tagList(
          box(title = strong("Reference", style = "color:white"), width = NULL, status = "primary", solidHeader = TRUE,
              HTML(paste(ref_string, collapse="<br/>"))
          )
        )
      })
    }
    
    incProgress(7/10, message = "Displaying Results")
    
    shinyjs::enable("runbtn_beta")
    shinyjs::enable("beta.outvars")
    shinyjs::enable("beta.censorvar")
    shinyjs::enable("beta.treatvar")
    shinyjs::enable("beta.outcome_types")
  }) 
  
  ##### TAXONOMIC ANALYSIS #####
  
  observeEvent(input$runbtn_taxa, {
    
    shinyjs::disable("runbtn_taxa")
    shinyjs::disable("choose_taxa")
    shinyjs::disable("include_species_taxa")
    shinyjs::disable("reg_taxa")
    shinyjs::disable("sim_taxa")
    shinyjs::disable("outcome_taxa")
    shinyjs::disable("treat_taxa")
    shinyjs::disable("covariates_taxa")
    
    if (input$cov_taxa == "None"){
      covariate_taxa <- NULL
    }else{
      covariate_taxa <- input$covariates_taxa
    }
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        taxa_ori <<- chooseData$taxa.out[[taxa.types$dataType]]
        
        taxa <<- list()
        for (i in 1:length(taxa_ori)){
          taxa[[i]] <<- data.frame(scale(taxa_ori[[i]]))
          rownames(taxa[[i]]) <<- rownames(taxa_ori[[i]])
          colnames(taxa[[i]]) <<- colnames(taxa_ori[[i]])
        }
        
        if (input$include_species_taxa == "Phylum - Species") {
          include = TRUE
        } else {
          include = FALSE
        }
        
        interaction_taxa = FALSE 
        
        names(taxa) <- names(taxa_ori)
        
        if (input$choose_taxa %in% c("Imai Method", "Imai Method (Default)")){
          
          if (length(table(chooseData$sam.dat[, input$outcome_taxa])) == 2){
            taxa_reg_type <<- "logistic"
          }else{
            taxa_reg_type <<- "linear"
          }
          
          taxa_sim_type <<- "bootstrap"
          boot.iterations <<- "perc"
          
          taxa_unique <<- name_split_list(taxa)
          
          if (binary_cont(chooseData$sam.dat[,input$treat_taxa]) == "binary"){
            
            med_result <<- tryCatch(mediation.taxon.total(chooseData$sam.dat, taxa, taxa_unique, input$treat_taxa, covariate_taxa, input$outcome_taxa, input$censor_taxa, interaction_taxa, n.sim = 3000, inc = include), 
                                    error = function(e) {
                                      message("No outcome is available!")
                                      showModal(modalDialog(div("No outcome is available!")))
                                      return(NULL)})
            
            result_before_qc <<- med_result 
            
            result_origin <<- tryCatch(med_result_concate(taxa, med_result, taxa_unique, include), 
                                       error = function(e) {  
                                         message("No outcome is available!")
                                         showModal(modalDialog(div("No outcome is available!")))
                                         return(NULL)})
            
            result <<- tryCatch(q_convert_tax_med(result_origin), 
                                error = function(e) {  
                                  message("No outcome is available!")
                                  showModal(modalDialog(div("No outcome is available!")))
                                  return(NULL)})
            
            result_after_qc <<- result 
            
            result_te <<- tryCatch(taxa.med.sep(result, "total_effect", include), 
                                   error = function(e) {  
                                     message("No outcome is available!")
                                     showModal(modalDialog(div("No outcome is available!")))
                                     return(NULL)})
            result_aa <<- tryCatch(taxa.med.sep(result, "acme_average", include), 
                                   error = function(e) {  
                                     message("No outcome is available!")
                                     showModal(modalDialog(div("No outcome is available!")))
                                     return(NULL)})
            result_ac <<- tryCatch(taxa.med.sep(result, "acme_control", include), 
                                   error = function(e) {  
                                     message("No outcome is available!")
                                     showModal(modalDialog(div("No outcome is available!")))
                                     return(NULL)})
            result_at <<- tryCatch(taxa.med.sep(result, "acme_treated", include), 
                                   error = function(e) {  
                                     message("No outcome is available!")
                                     showModal(modalDialog(div("No outcome is available!")))
                                     return(NULL)})
            result_ad <<- tryCatch(taxa.med.sep(result, "ade_average", include), 
                                   error = function(e) {  
                                     message("No outcome is available!")
                                     showModal(modalDialog(div("No outcome is available!")))
                                     return(NULL)})
            
            nrow <- try(taxa.forest.plot.pages(result_aa, species.include = include), silent = TRUE)
            nrow_1 <- try(taxa.forest.plot.pages(result_ac, species.include = include), silent = TRUE)
            nrow_2 <- try(taxa.forest.plot.pages(result_at, species.include = include), silent = TRUE)
            nrow_3 <- try(taxa.forest.plot.pages(result_ad, species.include = include), silent = TRUE)
            nrow_4 <- try(taxa.forest.plot.pages(result_te, species.include = include), silent = TRUE)
            
            height_forest <- c() 
            
            result_total <- list(result_te, result_aa, result_ac, result_at, result_ad) 
            
            for (j in 1:5){
              
              sig.by.rank <- list() 
              result <- result_total[[j]]
              for(i in 1:5+include) {
                out <- result[[i]]
                ind.sig <- which(out[,"Q.val"] < 0.05) 
                
                sig.by.rank[[i]] <- as.numeric(ind.sig)
              }
              
              sig.num <- length(unlist(sig.by.rank))
              
              if(sig.num ==0){height_forest <- c(height_forest, 200)} else if (sig.num ==1){
                height_forest <- c(height_forest, 150)
              } else if (sig.num <5 & sig.num >1){
                height_forest <-  c(height_forest, 300)
              }else if (sig.num <= 20 & sig.num >=5){
                height_forest <-  c(height_forest, 500)
              } else {height_forest <-  c(height_forest, 800)}
            }
            
            taxa.name.rank <<- taxa.chat.rank.name(result_aa, taxa.names.rank(taxa, include), include, "Est", TRUE)
            rank.names <- names(table(taxa.name.rank[,1]))
            
            if(include){
              rank.order <- match(rank.names, c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
            }else{
              rank.order <- match(rank.names, c("Phylum", "Class", "Order", "Family", "Genus"))
            }
            
            rank.names <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")[sort(rank.order)]
            
            if (length(taxa.name.rank[,2]) == 0){output$taxa_display = renderUI({
              tagList(
                tabBox(title = strong(NULL, style = "color:black", side = "right"), width = NULL,
                       
                       tabPanel("ACME (Mediation)", align = "center", tagList(
                         box(title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                             do.call(tabsetPanel, lapply(1:nrow, function(i){
                               tabPanel(title = paste0("Page ", i), align = "center",
                                        plotOutput(paste0("aa_co_forest", i), height = as.numeric(height_forest[2]), width = 750))
                             }))
                         ),
                         box( title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                              tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_aa")))
                       )),
                       
                       
                       tabPanel("ADE", align = "center",
                                tagList(
                                  box( title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                                       do.call(tabsetPanel, lapply(1:nrow_3, function(i){
                                         tabPanel(title = paste0("Page ", i), align = "center",
                                                  plotOutput(paste0("ad_forest", i), height = as.numeric(height_forest[5]), width = 750))
                                       }))
                                  ),
                                  box( title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                                       tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_ad")))
                                )
                       ),
                       
                       tabPanel("Total Effect", align = "center",
                                tagList(
                                  box(title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                                      do.call(tabsetPanel, lapply(1:nrow_4, function(i){
                                        tabPanel(title = paste0("Page ", i), align = "center",
                                                 plotOutput(paste0("te_forest", i), height = as.numeric(height_forest[1]), width = 750))
                                      }))
                                  ),
                                  box(title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                                      tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_te")))
                                )
                       )
                )
              )
            })}else{
              output$taxa_display = renderUI({
                tagList(
                  tabBox(title = strong(NULL, style = "color:black", side = "right"), width = NULL,
                         
                         tabPanel("ACME (Mediation)", align = "center", 
                                  tagList(
                                    box(title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                                        do.call(tabsetPanel, lapply(1:nrow, function(i){
                                          tabPanel(title = paste0("Page ", i), align = "center",
                                                   plotOutput(paste0("aa_co_forest", i), height = as.numeric(height_forest[2]), width = 750))
                                        }))
                                    ),
                                    box(title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                                        tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_aa"))),
                                    
                                  )),
                         
                         tabPanel("ADE", align = "center",
                                  tagList(
                                    box( title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                                         do.call(tabsetPanel, lapply(1:nrow_3, function(i){
                                           tabPanel(title = paste0("Page ", i), align = "center",
                                                    plotOutput(paste0("ad_forest", i), height = as.numeric(height_forest[5]), width = 750))
                                         }))
                                    ),
                                    box( title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                                         tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_ad")))
                                  )
                         ),
                         
                         tabPanel("Total Effect", align = "center",
                                  tagList(
                                    box(title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                                        do.call(tabsetPanel, lapply(1:nrow_4, function(i){
                                          tabPanel(title = paste0("Page ", i), align = "center",
                                                   plotOutput(paste0("te_forest", i), height = as.numeric(height_forest[1]), width = 750))
                                        }))
                                    ),
                                    box(title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                                        tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_te")))
                                  )
                         )
                  )
                )
              })
            }
            
            incProgress(1/10, message = "Graph")
            
            forestplot.data_0 <<- tryCatch(taxa.forest.plot.pages1(result_aa, taxa.names.rank(taxa, include), include, "Est", TRUE), 
                                           error = function(e) {  
                                             message("No outcome is available!")
                                             showModal(modalDialog(div("No outcome is available!")))
                                             return(NULL)})
            nrow <- tryCatch(taxa.forest.plot.pages(result_aa, species.include = include), 
                             error = function(e) {  
                               message("No outcome is available!")
                               showModal(modalDialog(div("No outcome is available!")))
                               return(NULL)})
            
            taxa.name.rank.0 <<- taxa.chat.rank.name(result_aa, taxa.names.rank(taxa, include), include, "Est", TRUE)
            
            if(include == TRUE){
              ranks.names <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
            }else{
              ranks.names <- c("Phylum", "Class", "Order", "Family", "Genus")
            }
            
            
            lapply(1:nrow, function(j) {
              output[[paste0("aa_co_forest", j)]] <- try(renderPlot({
                taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data_0, page = j)
              }), silent = TRUE)
            })
            
            forestplot.data_1 <<- tryCatch(taxa.forest.plot.pages1(result_ac, taxa.names.rank(taxa, include), include, "Est", TRUE), 
                                           error = function(e) {  
                                             message("No outcome is available!")
                                             showModal(modalDialog(div("No outcome is available!")))
                                             return(NULL)})
            nrow_1 <- try(taxa.forest.plot.pages(result_ac, species.include = include), silent = TRUE)
            
            taxa.name.rank.1 <<- taxa.chat.rank.name(result_ac, taxa.names.rank(taxa, include), include, "Est", TRUE)
            
            
            lapply(1:nrow_1, function(j) {
              output[[paste0("ac_co_forest", j)]] <- try(renderPlot({
                taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data_1, page = j)
              }), silent = TRUE)
            })
            
            forestplot.data_2 <<- tryCatch(taxa.forest.plot.pages1(result_at, taxa.names.rank(taxa, include), include, "Est", TRUE), 
                                           error = function(e) {  
                                             message("No outcome is available!")
                                             showModal(modalDialog(div("No outcome is available!")))
                                             return(NULL)})
            nrow_2 <- try(taxa.forest.plot.pages(result_at, species.include = include), silent = TRUE)
            
            taxa.name.rank.2 <<- taxa.chat.rank.name(result_at, taxa.names.rank(taxa, include), include, "Est", TRUE)
            
            lapply(1:nrow_2, function(j) {
              output[[paste0("at_forest", j)]] <- renderPlot({
                taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data_2, page = j)
              })
            })
            
            forestplot.data_3 <<-  tryCatch(taxa.forest.plot.pages1(result_ad, taxa.names.rank(taxa, include), include, "Est", TRUE), 
                                            error = function(e) {  
                                              message("No outcome is available!")
                                              showModal(modalDialog(div("No outcome is available!")))
                                              return(NULL)})
            nrow_3 <- try(taxa.forest.plot.pages(result_ad, species.include = include), silent = TRUE)
            
            taxa.name.rank.3 <<- taxa.chat.rank.name(result_ad, taxa.names.rank(taxa, include), include, "Est", TRUE)
            
            lapply(1:nrow_3, function(j) {
              output[[paste0("ad_forest", j)]] <- try(renderPlot({
                taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data_3, page = j)
              }), silent = TRUE)
            })
            
            forestplot.data_4 <<- tryCatch(taxa.forest.plot.pages1(result_te, taxa.names.rank(taxa, include), include, "Est", TRUE), 
                                           error = function(e) {  
                                             message("No outcome is available!")
                                             showModal(modalDialog(div("No outcome is available!")))
                                             return(NULL)})
            nrow_4 <- try(taxa.forest.plot.pages(result_te, species.include = include), silent = TRUE)
            
            taxa.name.rank.4 <<- taxa.chat.rank.name(result_te, taxa.names.rank(taxa, include), include, "Est", TRUE)
            
            lapply(1:nrow_4, function(j) {
              output[[paste0("te_forest", j)]] <- try(renderPlot({
                taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data_4, page = j)
              }), silent = TRUE)
            })
            
            output$taxa_display_dend_aa = try(renderUI({
              box(title = NULL, width = NULL,  solidHeader = TRUE,
                  grVizOutput("dendrogram_0", height = 1000))
            }), silent = TRUE) 
            
            output$taxa_display_dend_ac = try(renderUI({
              box(title = NULL, width = NULL,  solidHeader = TRUE,
                  grVizOutput("dendrogram_1", height = 1000))
            }), silent = TRUE)
            
            output$taxa_display_dend_at = try(renderUI({
              box(title = NULL, width = NULL,  solidHeader = TRUE,
                  grVizOutput("dendrogram_2", height = 1000))
            }), silent = TRUE)
            
            output$taxa_display_dend_ad = try(renderUI({
              box(title = NULL, width = NULL,  solidHeader = TRUE,
                  grVizOutput("dendrogram_3", height = 1000))
            }), silent = TRUE)
            
            output$taxa_display_dend_te = try(renderUI({
              box(title = NULL, width = NULL, solidHeader = TRUE,
                  grVizOutput("dendrogram_4", height = 1000))
            }), silent = TRUE)
            
            if (include){
              output$dendrogram_0 = try(renderGrViz({
                taxa.sig.dend(result_aa, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_1 = try(renderGrViz({
                taxa.sig.dend(result_ac, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_2 = try(renderGrViz({
                taxa.sig.dend(result_at, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_3 = try(renderGrViz({
                taxa.sig.dend(result_ad, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_4 = try(renderGrViz({
                taxa.sig.dend(result_te, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
            }else{
              
              output$dendrogram_0 = try(renderGrViz({
                taxa.sig.dend.genus(result_aa, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_1 = try(renderGrViz({
                taxa.sig.dend.genus(result_ac, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_2 = try(renderGrViz({
                taxa.sig.dend.genus(result_at, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_3 = try(renderGrViz({
                taxa.sig.dend.genus(result_ad, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_4 = try(renderGrViz({
                taxa.sig.dend.genus(result_te, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
            }
            
            output$downloadTable_taxa = renderUI({
              tagList(
                p(" ", style = "margin-bottom: +15px;"),
                box(title = strong("Download Output Table", style ="color: #FFFFFF; border-color: #ff7f50"), width = NULL, status = "primary", solidHeader = TRUE,
                    p("You can download data analysis outputs.",
                      style = "font-size:11pt"),
                    h5("Data Analysis Outputs"),
                    downloadButton("tdownloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3"),
                    h5("Dendrogram"),
                    downloadButton("gdownload", "Download", width = '50%', style = "color:black; background-color: red3")
                )
              )
            })
            
            
            result_down <- list() 
            for (i in 1:5+include){
              result_down[[i]] <- list(result_te[[i]], result_aa[[i]], result_ac[[i]], result_at[[i]], result_ad[[i]])
              names(result_down[[i]]) <- c("Total Effect", "ACME (Overall)", "ACME (Control)", "ACME (Treated)", "ADE")
            }
            
            if (include){
              output$tdownloadTabl1 <- downloadHandler(
                filename = function() {
                  paste("Taxa.Analysis.Output.zip")
                }, 
                content = function(DA.file) {
                  temp <- setwd(tempdir())
                  on.exit(setwd(temp))
                  dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
                  
                  capture.output(result_down[[1]], file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[2]], file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[3]], file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[4]], file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[5]], file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[6]], file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  
                  zip(zipfile=DA.file, files=dataFiles)
                }
              )
            }else{
              output$tdownloadTabl1 <- downloadHandler(
                filename = function() {
                  paste("Taxa.Analysis.Output.zip")
                }, 
                content = function(DA.file) {
                  temp <- setwd(tempdir())
                  on.exit(setwd(temp))
                  dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
                  
                  capture.output(result_down[[1]], file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[2]], file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[3]], file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[4]], file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[5]], file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  
                  zip(zipfile=DA.file, files=dataFiles)
                }
              )
            }
            
            if (include){
              
              output$dendrogram_0 = try(renderGrViz({
                taxa.sig.dend(result_aa, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_1 = try(renderGrViz({
                taxa.sig.dend(result_ac, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_2 = try(renderGrViz({
                taxa.sig.dend(result_at, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_3 = try(renderGrViz({
                taxa.sig.dend(result_ad, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_4 = try(renderGrViz({
                taxa.sig.dend(result_te, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              
              output$gdownload <- downloadHandler(
                filename = function() {
                  paste("Taxa.Analysis.Graphical.Output.zip")
                },
                content = function(file) {
                  temp <- setwd(tempdir())
                  on.exit(setwd(temp))
                  dataFiles = c("Total Effect.html", "ACME (Overall).html", "ACME (Control).html", "ACME (Treated).html", "ADE.html")
                  
                  htmlwidgets::saveWidget(as_widget(taxa.sig.dend(result_te, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "Total Effect.html")
                  htmlwidgets::saveWidget(as_widget(taxa.sig.dend(result_aa, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "ACME (Overall).html")
                  htmlwidgets::saveWidget(as_widget(taxa.sig.dend(result_ac, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "ACME (Control).html")
                  htmlwidgets::saveWidget(as_widget(taxa.sig.dend(result_at, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "ACME (Treated).html")
                  htmlwidgets::saveWidget(as_widget(taxa.sig.dend(result_ad, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "ADE.html")
                  
                  zip(zipfile=file, files=dataFiles)
                }
              )
              
            }else{
              output$dendrogram_0 = try(renderGrViz({
                taxa.sig.dend.genus(result_aa, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_1 = try(renderGrViz({
                taxa.sig.dend.genus(result_ac, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_2 = try(renderGrViz({
                taxa.sig.dend.genus(result_at, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_3 = try(renderGrViz({
                taxa.sig.dend.genus(result_ad, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$dendrogram_4 = try(renderGrViz({
                taxa.sig.dend.genus(result_te, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
              
              output$gdownload <- downloadHandler(
                filename = function() {
                  paste("Taxa.Analysis.Graphical.Output.zip")
                },
                content = function(file) {
                  temp <- setwd(tempdir())
                  on.exit(setwd(temp))
                  dataFiles = c("Total Effect.html", "ACME (Overall).html", "ACME (Control).html", "ACME (Treated).html", "ADE.html")
                  
                  htmlwidgets::saveWidget(as_widget(taxa.sig.dend.genus(result_te, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "Total Effect.html")
                  htmlwidgets::saveWidget(as_widget(taxa.sig.dend.genus(result_aa, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "ACME (Overall).html")
                  htmlwidgets::saveWidget(as_widget(taxa.sig.dend.genus(result_ac, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "ACME (Control).html")
                  htmlwidgets::saveWidget(as_widget(taxa.sig.dend.genus(result_at, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "ACME (Treated).html")
                  htmlwidgets::saveWidget(as_widget(taxa.sig.dend.genus(result_ad, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "ADE.html")
                  
                  zip(zipfile=file, files=dataFiles)
                }
              )
            }
          }
          else {
            
            chooseData$sam.dat[[input$outcome_taxa]] <- as.numeric(chooseData$sam.dat[[input$outcome_taxa]])
            
            med_result <<- tryCatch(mediation.taxon.total(chooseData$sam.dat, taxa, taxa_unique, input$treat_taxa, covariate_taxa, input$outcome_taxa, input$censor_taxa, interaction_taxa, n.sim = 3000, inc = include),  #need to be adjusted (n.sim)  #here need change 
                                    error = function(e) {
                                      message("No outcome is available!")
                                      showModal(modalDialog(div("No outcome is available!")))
                                      return(NULL)})
            
            
            result <<- tryCatch(med_result_concate(taxa, med_result, taxa_unique, include), 
                                error = function(e) {  
                                  message("No outcome is available!")
                                  showModal(modalDialog(div("No outcome is available!")))
                                  return(NULL)})
            
            result<- tryCatch(q_convert_tax_med(result), 
                              error = function(e) {  
                                message("No outcome is available!")
                                showModal(modalDialog(div("No outcome is available!")))
                                return(NULL)})
            
            
            result_te <<- tryCatch(taxa.med.sep(result, "total_effect", include), error = function(e) {  
              message("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)})
            
            result_aa <<- tryCatch(taxa.med.sep(result, "acme_average", include), 
                                   error = function(e) {  
                                     message("No outcome is available!")
                                     showModal(modalDialog(div("No outcome is available!")))
                                     return(NULL)})
            result_ad <<- tryCatch(taxa.med.sep(result, "ade_average", include), 
                                   error = function(e) {  
                                     message("No outcome is available!")
                                     showModal(modalDialog(div("No outcome is available!")))
                                     return(NULL)})
            
            nrow_1 <- try(taxa.forest.plot.pages(result_aa, species.include = include), silent = TRUE)
            nrow_2 <- try(taxa.forest.plot.pages(result_ad, species.include = include), silent = TRUE)
            nrow_3 <- try(taxa.forest.plot.pages(result_te, species.include = include), silent = TRUE)
            
            height_forest <- c() 
            
            result_total <- list() 
            result_total[[1]] <- result_te 
            result_total[[2]] <- result_aa
            result_total[[3]] <- result_ad
            
            for (j in 1:3){
              
              sig.by.rank <- list() 
              result <- result_total[[j]]
              for(i in 1:5+include) {
                out <- result[[i]]
                ind.sig <- which(out[,"Q.val"] < 0.05) 
                
                sig.by.rank[[i]] <- as.numeric(ind.sig)
              }
              
              sig.num <- length(unlist(sig.by.rank))
              
              if(sig.num ==0){height_forest <- c(height_forest, 200)} else if (sig.num ==1){
                height_forest <- c(height_forest, 150)
              } else if (sig.num <5 & sig.num >1){
                height_forest <-  c(height_forest, 300)
              }else if (sig.num <= 20 & sig.num >=5){
                height_forest <-  c(height_forest, 500)
              } else {height_forest <-  c(height_forest, 800)}
            }
            
            taxa.name.rank <<- taxa.chat.rank.name(result_aa, taxa.names.rank(taxa, include), include, "Est", TRUE)
            rank.names <- names(table(taxa.name.rank[,1]))
            
            if(include){
              rank.order <- match(rank.names, c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
            }else{
              rank.order <- match(rank.names, c("Phylum", "Class", "Order", "Family", "Genus"))
            }
            
            rank.names <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")[sort(rank.order)]
            
            
            if (length(taxa.name.rank[,2]) == 0){
              output$taxa_display = try(renderUI({
                tagList(
                  tabBox(title = strong(NULL, style = "color:black", side = "right"), width = NULL,
                         tabPanel("ACME (Mediation)", align = "left",
                                  tagList(
                                    box(title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                                        do.call(tabsetPanel, lapply(1:nrow_1, function(i){
                                          tabPanel(title = paste0("Page ", i), align = "center",
                                                   plotOutput(paste0("aa_forest", i), height = as.numeric(height_forest[2]), width = 750))
                                        }))),
                                    box(title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                                        tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_ac")))
                                    
                                    
                                  )),
                         tabPanel("ADE", align = "center",
                                  tagList(
                                    box( title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                                         do.call(tabsetPanel, lapply(1:nrow_2, function(i){
                                           tabPanel(title = paste0("Page ", i), align = "center",
                                                    plotOutput(paste0("ad_forest", i), height = as.numeric(height_forest[3]), width = 750))
                                         }))),
                                    box( title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                                         tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_ad"))),
                                    
                                  )),
                         
                         tabPanel("Total Effect", align = "center",
                                  tagList(
                                    box( title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                                         do.call(tabsetPanel, lapply(1:nrow_3, function(i){
                                           tabPanel(title = paste0("Page ", i), align = "center",
                                                    plotOutput(paste0("te_forest", i), height = as.numeric(height_forest[1]), width = 750))
                                         }))),
                                    box( title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                                         tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_te")))
                                    
                                  )
                         )
                  )
                )
              }), silent = TRUE)
            }else{
              output$taxa_display = try(renderUI({
                tagList(
                  tabBox(title = strong(NULL, style = "color:black", side = "right"), width = NULL,
                         tabPanel("ACME (Mediation)", align = "left",
                                  tagList(
                                    box(title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                                        do.call(tabsetPanel, lapply(1:nrow_1, function(i){
                                          tabPanel(title = paste0("Page ", i), align = "center",
                                                   plotOutput(paste0("aa_forest", i), height = as.numeric(height_forest[2]), width = 750))
                                        }))),
                                    box(title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                                        tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_ac")))
                                  )),
                         tabPanel("ADE", align = "center",
                                  tagList(
                                    box( title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                                         do.call(tabsetPanel, lapply(1:nrow_2, function(i){
                                           tabPanel(title = paste0("Page ", i), align = "center",
                                                    plotOutput(paste0("ad_forest", i), height = as.numeric(height_forest[3]), width = 750))
                                         }))),
                                    box( title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                                         tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_ad"))),
                                    
                                  )),
                         
                         tabPanel("Total Effect", align = "center",
                                  tagList(
                                    box( title = strong("Forest Plot", style = "color:white", side = "left"), width = NULL, status = "primary", solidHeader = TRUE, 
                                         do.call(tabsetPanel, lapply(1:nrow_3, function(i){
                                           tabPanel(title = paste0("Page ", i), align = "center",
                                                    plotOutput(paste0("te_forest", i), height = as.numeric(height_forest[1]), width = 750))
                                         }))),
                                    box( title = strong("Dendrogram", style = "color:white", side = "left"), width = NULL, solidHeader = TRUE, status = "primary", 
                                         tabPanel(title = NULL, align = "center", uiOutput("taxa_display_dend_te")))
                                    
                                  )
                         )
                  )
                )
              }), silent = TRUE)
            }
            
            incProgress(1/10, message = "Graph")
            
            name.taxa.here <<- taxa.names.rank(taxa, include)
            
            forestplot.data_1 <<- tryCatch(taxa.forest.plot.pages1(result_aa, taxa.names.rank(taxa, include), include, "Est", TRUE), 
                                           error = function(e) {  
                                             message("No outcome is available!")
                                             showModal(modalDialog(div("No outcome is available!")))
                                             return(NULL)})
            
            
            taxa.name.rank.1 <<- taxa.chat.rank.name(result_aa, taxa.names.rank(taxa, include), include, "Est", TRUE)
            
            
            lapply(1:nrow_1, function(j) {
              output[[paste0("aa_forest", j)]] <- try(renderPlot({
                taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data_1, page = j)
              }), silent = TRUE)
            })
            
            forestplot.data_2 <- try(taxa.forest.plot.pages1(result_ad, taxa.names.rank(taxa, include), include, "Est", TRUE), silent = TRUE) 
            
            lapply(1:nrow_2, function(j) {
              output[[paste0("ad_forest", j)]] <- try(renderPlot({
                taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data_2, page = j)
              }), silent = TRUE)
            })
            
            taxa.name.rank.2 <<- taxa.chat.rank.name(result_ad, taxa.names.rank(taxa, include), include, "Est", TRUE)
            
            forestplot.data_3 <- tryCatch(taxa.forest.plot.pages1(result_te, taxa.names.rank(taxa, include), include, "Est", TRUE), 
                                          error = function(e) {  
                                            message("No outcome is available!")
                                            showModal(modalDialog(div("No outcome is available!")))
                                            return(NULL)})
            
            taxa.name.rank.3 <<- taxa.chat.rank.name(result_te, taxa.names.rank(taxa, include), include, "Est", TRUE)
            
            lapply(1:nrow_3, function(j){
              output[[paste0("te_forest", j)]] <- try(renderPlot({
                taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data_3, page = j)
              }), silent = TRUE)
            })
            
            output$taxa_display_dend_aa = try(renderUI({
              box(title = NULL, width = NULL,  solidHeader = TRUE,
                  grVizOutput("dendrogram_1", height = 1000))
            }), silent = TRUE)
            
            output$taxa_display_dend_ad = renderUI({
              box(title = NULL, width = NULL,  solidHeader = TRUE,
                  grVizOutput("dendrogram_2", height = 1000))
            })
            
            output$taxa_display_dend_te = renderUI({
              box(title = NULL, width = NULL, solidHeader = TRUE,
                  grVizOutput("dendrogram_3", height = 1000))
            })
            
            output$dendrogram_1 = try(renderGrViz({
              taxa.sig.dend(result_aa, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
            }), silent = TRUE)
            
            output$dendrogram_2 = try(renderGrViz({
              taxa.sig.dend(result_ad, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
            }), silent = TRUE)
            
            output$dendrogram_3 = try(renderGrViz({
              taxa.sig.dend(result_te, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
            }), silent = TRUE)
            
            
            output$downloadTable_taxa = renderUI({
              tagList(
                p(" ", style = "margin-bottom: +15px;"),
                box(title = strong("Download Output Table", style = "color: #FFFFFF; border-color: #ff7f50"), width = NULL, status = "primary", solidHeader = TRUE,
                    p("You can download data analysis outputs.",
                      style = "font-size:11pt"),
                    h5("Data Analysis Outputs"),
                    downloadButton("tdownloadTabl1", "Download", width = '50%')
                )
              )
            })
            
            result_down <- list() 
            for (i in 1:5+include){
              result_down[[i]] <- list(result_aa[[i]], result_te[[i]],  result_ad[[i]])
              names(result_down[[i]]) <- c("ACME (Overall)", "Total Effect", "ADE")
            }
            
            if (include){
              output$tdownloadTabl1 <- downloadHandler(
                filename = function() {
                  paste("Taxa.Analysis.Output.zip")
                }, 
                content = function(DA.file) {
                  temp <- setwd(tempdir())
                  on.exit(setwd(temp))
                  dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
                  
                  capture.output(result_down[[1]], file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[2]], file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[3]], file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[4]], file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[5]], file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[6]], file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  
                  zip(zipfile=DA.file, files=dataFiles)
                }
              )
            }else{
              output$tdownloadTabl1 <- downloadHandler(
                filename = function() {
                  paste("Taxa.Analysis.Output.zip")
                }, 
                content = function(DA.file) {
                  temp <- setwd(tempdir())
                  on.exit(setwd(temp))
                  dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
                  
                  capture.output(result_down[[1]], file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[2]], file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[3]], file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[4]], file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  capture.output(result_down[[5]], file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                  
                  zip(zipfile=DA.file, files=dataFiles)
                }
              )
            }
            
            output$gdownload <- downloadHandler(
              filename = function() {
                paste("Taxa.Analysis.Graphical.Output.zip")
              },
              content = function(file) {
                temp <- setwd(tempdir())
                on.exit(setwd(temp))
                dataFiles = c("ACME.html", "Total Effect.html",  "ADE.html")
                
                htmlwidgets::saveWidget(as_widget(taxa.sig.dend(result_aa, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "ACME.html")
                htmlwidgets::saveWidget(as_widget(taxa.sig.dend(result_te, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "Total Effect.html")
                htmlwidgets::saveWidget(as_widget(taxa.sig.dend(result_ad, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file = "ADE.html")
                
                zip(zipfile=file, files=dataFiles)
              }
            )
          }
          
          ref_string = REFERENCE_CHECK(method_name = isolate(input$choose_taxa), FDR = "Yes")
          if (is.null(ref_string)) {
            shinyjs::hide("taxa_references")
          } else {
            shinyjs::show("taxa_references")
            output$taxa_references = renderUI({
              tagList(
                box(title = strong("References", style = "color: #FFFFFF; border-color: #ff7f50"), width = NULL, status = "primary", solidHeader = TRUE,
                    HTML(paste(ref_string, collapse="<br/>"))
                )
              )
            })
          } 
          
          if(include == TRUE){
            ranks.names <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
          }else{
            ranks.names <- c("Phylum", "Class", "Order", "Family", "Genus")
          }
        } 
      }) 
    
    shinyjs::disable("runbtn_chat")
    shinyjs::enable("runbtn_taxa")
    shinyjs::enable("choose_taxa")
    shinyjs::enable("include_species_taxa")
    shinyjs::enable("reg_taxa")
    shinyjs::enable("sim_taxa")
    shinyjs::enable("outcome_taxa")
    shinyjs::enable("treat_taxa")
    shinyjs::enable("covariates_taxa")
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
}

shinyApp(ui = ui, server = server)
