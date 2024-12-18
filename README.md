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
install.packages(c("depmixS4", "nnet", "ggplot2", "dplyr", "tidyr", "terra"))  



How to Reproduce the Results
Clone the repository:

bash
Copiar código
git clone https://github.com/msuassuna/RegimeSwitchesPantanal.git  
Set up the environment:

Open the R project file RegimeSwitchesPantanal.Rproj in RStudio.
Ensure all required packages are installed.
Run the Main Script:
Execute MainScript.R to reproduce the analyses and generate the outputs.

Outputs
Figures: Generated plots for transition matrices, drought class analysis, and teleconnection patterns.
Tables: Outputs summarizing model parameters, probabilities, and statistical summaries.
License
This repository is licensed under the MIT License.

Contact
For questions, feedback, or collaboration opportunities, please contact:
Marcus Suassuna

Email: example@email.com
GitHub: github.com/msuassuna
javascript
Copiar código
