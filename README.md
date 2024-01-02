# MiMedSurv

Title: Comprehensive microbiome mediation analysis with survival (i.e., time-to-event) responses. 

Version: 1.0.0

Maintainer: Hyo Jung Jang <hyojung.jang@northwestern.edu>

Description: MiMedSurv (microbiome mediation analysis with survival responses) is a unified web cloud computing platform for comprehensive microbiome mediation analysis with survival (i.e., time-to-event) responses. MiMedSurv surveys the roles of the microbiome as a mediator between the treatment (e.g., environmental, behavioral or medical exposures) and the survival response on the host’s health or disease (e.g., time to the onset of a disease). The main features of MiMedSurv are as follows. First, MiMedSurv conducts basic non-mediational analysis, not involving the microbiome in the analysis, to survey the disparity in survival time between the treatment groups (e.g., 'placebo and treatment groups' or 'old treatment and new treatment groups') (see Non-Mediational Analysis). Second, MiMedSurv conducts microbiome mediational analysis (see Mediational Analysis) in various aspects (1) as a whole microbial ecosystem using different ecological measures (e.g., alpha- and beta-diversity indices) (see Community-level Analysis) or (2) as individual microbial taxa (e.g., phyla, classes, orders, families, genera, species) using different data normalization methods (see Taxonomy-level Analysis). Importantly, MiMedSurv enables covariate-adjusted analysis to control for potential confounding factors (e.g., age, gender) for both non-mediational and mediational analyses to enhance the causality of the results especially for observational studies. Finally, MiMedSurv provides flexible and user-friendly data processing and analytic modules and makes nice graphical representations.

NeedsCompilation: No

Depends: R(≥ 4.1.0)

Imports: Bioconductor ('BiocParallel', 'biomformat', 'phyloseq'); CRAN ('betareg', 'BiasedUrn', 'BiocManager', 'bios2mds', 'CompQuadForm', 'dashboardthemes', 'devtools', 'DiagrammeR', 'dirmult', 'dplyr', 'DT', 'ecodist', 'entropart', 'erer', 'fBasics', 'forestplot', 'fossil', 'ggplot2', 'ggthemes', 'googleVis', 'gridExtra', 'gridGraphics', 'compositions', 'GUniFrac', 'htmltools', 'ICSNP', 'lme4', 'lmerTest', 'MiRKAT', 'nlme', 'patchwork', 'phangorn', 'picante', 'plotly', 'PMCMRplus', 'quantreg', 'remotes', 'reticulate', 'rgl', 'rmarkdown', 'robCompositions', 'robustbase', 'seqinr', 'shiny', 'shinydashboard', 'shinyjs', 'shinyWidgets', 'stringr', 'tidyverse', 'vegan', 'xtable', 'zCompositions', 'zip', 'bda', 'mediation', 'survival', 'survminer', 'coin'); GitHub ('DACT', 'chatgpt')

License: GPL 1, GPL 2 

**URLs**: Web Server (http://mimedsurv.micloud.kr), GitHub (http://github.com/yj7599/mimedsurvgit) 

**Maintainer**: Hyojung Jang (hyojung.jang@northwestern.edu)

**Reference**: Jang, H. and Koh, H. MiMedSurv: A Unified Web Cloud Computing Service for Microbiome Mediation Analysis with Survival Responses (In review)

## Prerequites

#### Notice: For the local implementation, you do not need to install all the pre-requite R packages individually. You only need to install the 'shiny' package, and then run a simple command in 'Launch App' below. Then, all the pre-requisite R packages will be installed and imported automatically. 

shiny
```
install.packages('shiny')
```

## Launch App

```
library(shiny)

runGitHub('mimedsurvgit', 'yj7599', ref = 'main')
```

## Troubleshooting Tips

If you have any problems for using MiMedSurv, please report in issues (https://github.com/YJ7599/mimedsurvgit/issues) or email Hyo Jung Jang (hyojung.jang@northwestern.edu).

