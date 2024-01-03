# MiMedSurv

Title: MiMedSurv: A Unified Cloud Platform for Microbiome Causal Mediation Analysis with Survival Responses

Version: 1.0.0

Description: MiMedSurv (microbiome mediation analysis with survival responses) is a unified web cloud computing platform for comprehensive microbiome mediation analysis with survival (i.e., time-to-event) responses. MiMedSurv surveys the roles of the microbiome as a mediator between a treatment (e.g., medical intervention, environmental exposure) and a survival response on the host’s health or disease (e.g., time-to-disease, time-to-cure). The main features of MiMedSurv are as follows. First, MiMedSurv conducts basic exploratory non-mediational survival analysis, not involving microbiome, to survey the disparity in survival time between medical treatments (e.g., treatment vs. placebo, new treatment vs. old treatment) / environmental exposures (e.g., rural vs. urban, smoking vs. non-smoking). (see Non-Mediational Analysis). Second, MiMedSurv identifies the mediating roles of the microbiome in various aspects (see Mediational Analysis): (i) as a microbial ecosystem using ecological measures (e.g., alpha- and beta-diversity indices) (see Community-level Analysis) and (ii) as individual microbial taxa in various hierarchies (e.g., phyla, classes, orders, families, genera, species) (see Taxonomy-level Analysis). We also stress that MiMedSurv can conduct covariate-adjusted analysis to control for potential confounding factors (e.g., age, sex) to enhance the causality of the results especially for observational studies. MiMedSurv also provides user-friendly data preprocessing and analytic modules and makes nice visualizations.

NeedsCompilation: No

Depends: R(≥ 4.1.0)

Imports: Bioconductor ('BiocParallel', 'biomformat', 'phyloseq'); CRAN ('betareg', 'BiasedUrn', 'BiocManager', 'bios2mds', 'CompQuadForm', 'dashboardthemes', 'devtools', 'DiagrammeR', 'dirmult', 'dplyr', 'DT', 'ecodist', 'entropart', 'erer', 'fBasics', 'forestplot', 'fossil', 'ggplot2', 'ggthemes', 'googleVis', 'gridExtra', 'gridGraphics', 'compositions', 'GUniFrac', 'htmltools', 'ICSNP', 'MiRKAT', 'nlme', 'patchwork', 'phangorn', 'picante', 'plotly', 'PMCMRplus', 'quantreg', 'remotes', 'reticulate', 'rgl', 'rmarkdown', 'robCompositions', 'robustbase', 'seqinr', 'shiny', 'shinydashboard', 'shinyjs', 'shinyWidgets', 'stringr', 'tidyverse', 'vegan', 'xtable', 'zip', 'bda', 'mediation', 'survival', 'survminer', 'coin'); GitHub ('DACT')

License: GPL 1, GPL 2 

**URLs**: Web Server (http://mimedsurv.micloud.kr), GitHub (http://github.com/yj7599/mimedsurvgit) 

**Maintainer**: Hyojung Jang (hyojung.jang@northwestern.edu)

**Reference**: Jang, H., Koh, H. MiMedSurv: A Unified Web Cloud Computing Platform for Microbiome Causal Mediation Analysis with Survival Responses (In Review)


## Prerequites

#### Notice: For the local implementation, you do not need to install all the pre-requite R packages individually. You only need to install the 'shiny' package, and then run a simple command in 'Launch App' below. Then, all the pre-requisite R packages will be installed and imported automatically. 

shiny
```
install.packages('shiny')
```

## Launch App

```
library(shiny)

runGitHub('MiMedSurvGit', 'yj7599', ref = 'main')
```

## Troubleshooting Tips

If you have any problems for using MiMedSurv, please report in issues (https://github.com/YJ7599/mimedsurvgit/issues) or email Hyo Jung Jang (hyojung.jang@northwestern.edu).

