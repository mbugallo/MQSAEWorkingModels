# MQSAEWorkingModels
R code for the paper Generalized M-Quantile Small Area Estimation under Working Models, including simulations and real-data applications.

**Authors:** María Bugallo [1], Nicola Salvati [2] and Ray Chambers [3]

**Affiliations**:

[1] Department of Mathematics, University of A Coruña, Spain.

[2] Department of Economics and Management, University of Pisa, Italy.

[3] College of Business and Economics, Australian National University, Australia.


##  Repository Overview
This repository contains all the code and materials required to reproduce the analyses, simulations, and figures presented in the main document of the paper Generalized M-Quantile Small Area Estimation under Working Models and its appendices.
All datasets used in this repository are either: (1) generated directly in R, or (2) loaded from existing R libraries. No external raw data files are required beyond those available through standard R packages or the referenced sources.

## File and Folder Structure

### `binaryMQ/`
This folder contains the scripts for the binary M-quantile case. These scripts are based on the methodology developed in:

Chambers, R., Salvati, N., & Tzavidis, N. (2016). *Semiparametric small area estimation for binary outcomes, with application to unemployment estimation for local authorities in the UK*. Journal of the Royal Statistical Society: Series A (Statistics in Society), 179(2), 453–479.  

The original materials can be obtained from:  
[http://wileyonlinelibrary.com/journal/rss-datasets](http://wileyonlinelibrary.com/journal/rss-datasets)  

The code has been adapted and integrated into this repository for reproducibility and consistency with the rest of the project.

---

### `results/`
Running the scripts will automatically generate and store outputs in this directory.

---

### `AuxFunctions.R`
This script contains the auxiliary functions used throughout the repository. All main scripts source this file when needed.

---

### `lossfunction.R`
Graphical representation of the GALI mean, variance, etc.; and loss functions. This script produces: Figure 1 of the main document; Results included in Appendix A.

---

### `Simulation1Continuous.R`
Simulation study for continuous simulated data. Corresponds to Section 5 of the main document. Additional results are provided in Appendix E. This script generates fully simulated continuous data and evaluates the performance of the methods.

---

### `Simulation1Binary.R`
Simulation study for binary simulated data. Corresponds to Section 6 of the main document. Additional results are provided in Appendix G. This script generates fully simulated binary data and evaluates the proposed methodology under this framework.

---

### `Countydata.R`
Application to real data. Corresponds to Section 7 of the main document. This dataset contains agricultural information for 12 counties in Iowa, where US farm survey data were combined with LANDSAT satellite imagery to construct predictors of county-level mean hectares of corn and soybeans.  

From: Battese, G. E., Harter, R. M., & Fuller, W. A. (1988). *An error component model for prediction of county crop areas using survey and satellite data*. Journal of the American Statistical Association, 83, 28–36.

---

### `SimuCountydata.R`
Simulation study based on the continuous real-data structure. Corresponds to Section 7 of the main document. Additional results are provided in Appendix F. This script generates simulated continuous responses based on the structure of the real dataset.
