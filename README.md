# RegimeSwitchesPantanal  

## Overview  
This repository contains the data, scripts, and outputs associated with the study on **hydrological regime shifts and drought classification** in the Paraguay River Basin (PRB). The research applies a combination of **Hidden Markov Models (HMM)** and **Multinomial Logistic Regression (MLR)** to analyze low-water levels, drought persistence, and climatological influences.  

The repository is structured to foster **transparency** and **reproducibility** of all analyses and visualizations presented in the study.  

---

## Repository Structure  

- **data/**: Contains raw input data, including hydrological and climate time series used in the study.  
- **outputTables/**: Includes tables summarizing the outputs of statistical models and simulations.  
- **plots/**: Folder with generated figures, including transition matrices, histograms, and teleconnection analysis.  
- **scripts/**: Contains R scripts used for data processing, modeling, and visualization.  
- **tables/**: Pre-processed tables used as intermediate or final results for the analysis.  
- **damages/**: Placeholder for additional files or context-specific outputs.  

**Key Files**:  
- `MainScript.R`: The primary R script that executes the workflow, from data loading to model fitting and plotting results.  
- `README.md`: This documentation file.  
- `RegimeSwitchesPantanal.Rproj`: R project file for setting up the working environment.  
- `66825000.csv` & `67100000.csv`: Example datasets used for stream gauge analysis.  

### Notes on Data Availability

The land use data utilized in this study are not included in this repository due to file size constraints. However, these data are available from Dias et al. (2016), which provides detailed land cover and land use classification for the study area. The dataset can be accessed at http://www.biosfera.dea.ufv.br/pt-BR/banco/uso-do-solo-agricola-no-brasil-1940-2012---dias-et-al-2016.

If you encounter any difficulties accessing these files or have questions about their use, please feel free to contact me for assistance.
---

## Dependencies  

The following R packages are required to run the analysis:  
- `depmixS4`: For Hidden Markov Model (HMM) implementation.  
- `nnet`: For multinomial logistic regression.  
- `ggplot2`: For creating visualizations.  
- `dplyr` and `tidyr`: For data manipulation.  
- `terra` or `raster`: For geospatial operations.  

To install these dependencies, run the following code in R:  
```r  
install.packages(c(
  "tidyverse", "gridExtra", "lubridate", "ForestFit", "caret", "nnet",
  "randomForest", "kernlab", "MASS", "cowplot", "markovchain", "e1071", 
  "rpart", "rpart.plot", "dismo", "mlbench", "gbm", "xgboost", 
  "varSelRF", "Gmisc", "pdp", "depmixS4", "plyr", "ggpubr", 
  "diagram", "zoo", "ggdist", "tidyquant"
))
```


## How to Reproduce the Results

1. Clone the repository:

```r  
bash
git clone https://github.com/msuassuna/RegimeSwitchesPantanal.git  
```

2. Set up the environment:

- Open the R project file RegimeSwitchesPantanal.Rproj in RStudio
- Ensure all required packages are installed

3. Run the Main Script:
- Execute MainScript.R step-by-step to reproduce the analyses and generate the outputs

## Outputs
- Figures: Generated plots for transition matrices, drought class analysis, and teleconnection patterns.
- Tables: Outputs summarizing model parameters, probabilities, and statistical summaries.

## Contact
For questions, feedback, or collaboration opportunities, please contact:
Marcus Suassuna

Email: marcus.santos@sgb.gov.br / msuassuna@gmail.com

GitHub: github.com/msuassuna
