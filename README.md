# About the Package

The R package `CMAverse` provides a suite of functions conducting causal mediation analysis including `cmdag` for DAG visualization, `cmest` for statistical modeling and `cmsens` for sensitivity analysis.

## DAG Visualization

`cmdag` visualizes the scientific setting via a DAG.

## Statistical Modeling

`cmest` implements six causal mediation analysis approaches including the regression-based approach by [Valeri et al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3659198/) and [Vanderweele et al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287269/), 
the weighting-based approach by [Vanderweele et al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287269/), the inverse odd-ratio weighting 
approach by [Tchetgen Tchetgen et al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3954805/), the natural effect model by
[Vansteelandt et al. (2012)](https://www.degruyter.com/view/journals/em/1/1/article-p131.xml?language=en), the marginal structural model by [VanderWeele 
et al. (2009)](https://pubmed.ncbi.nlm.nih.gov/19234398), and the g-formula approach by [Lin et al. (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5285457/). 

## Multiple Imputation

`cmest` provides options to perform multiple imputation for the dataset with missing values via the `mice` package, estimate the causal effects by each of the imputed datasets and pool the results together.

## Sensitivity Analysis

`cmsens` conducts sensitivity analysis for unmeasured confounding via the E-value approach by [Vanderweele et al. (2017)](https://pubmed.ncbi.nlm.nih.gov/28693043/) and sensitivity analysis for measurement error via the regression calibration correction by [Carroll et al. (1995)](https://www.taylorfrancis.com/books/9780429139635), the SIMEX correction by [Cook et al. (1994)](https://www.jstor.org/stable/2290994?seq=1#metadata_info_tab_contents) and the MCSIMEX correction by [Küchenhoff et al. (2006)](https://pubmed.ncbi.nlm.nih.gov/16542233/). The sensitivity analysis for measurement error is currently available for the regression-based approach and the g-formulaa approach.


# Installation

The latest version can be installed via:

```R
devtools::install_github("LindaValeri/CMAverse")
```

Load CMAverse:

```R
library(CMAverse)
```

# References

Valeri L, Vanderweele TJ (2013). Mediation analysis allowing for exposure-mediator interactions and causal interpretation: theoretical assumptions and implementation with SAS and SPSS macros. Psychological Methods. 18(2): 137 - 150.
 
VanderWeele TJ, Vansteelandt S (2014). Mediation analysis with multiple mediators. Epidemiologic Methods. 2(1): 95 - 115.

Tchetgen Tchetgen EJ (2013). Inverse odds ratio-weighted estimation for causal mediation analysis. Statistics in medicine. 32: 4567 - 4580.

Nguyen QC, Osypuk TL, Schmidt NM, Glymour MM, Tchetgen Tchetgen EJ. Practical guidance for conducting mediation analysis with multiple mediators using inverse odds ratio weighting (2015). American Journal of Epidemiology. 181(5): 349 - 356.

VanderWeele TJ. Marginal structural models for the estimation of direct and indirect effects (2009). Epidemiology. 20(1): 18 - 26.

VanderWeele TJ, Tchetgen Tchetgen EJ (2017). Mediation analysis with time varying exposures and mediators. Journal of the Royal Statistical Society: Series B (Statistical Methodology). 79(3): 917 - 938.

Lin SH, Young J, Logan R, Tchetgen Tchetgen EJ, VanderWeele TJ (2017). Parametric mediational g-formula approach to mediation analysis with time-varying exposures, mediators, and confounders. Epidemiology. 28: 266 - 274.

Vansteelandt S, Bekaert M, Lange T. (2012). Imputation Strategies for the Estimation of Natural Direct and Indirect Effects. Epidemiologic Methods. 1(1): 131 - 158.

Steen J, Loeys T, Moerkerke B, Vansteelandt S (2017). Medflex: an R package for flexible mediation analysis using natural effect models. Journal of Statistical Software. 76(11).

Imai K, Keele L, Tingley D. A general approach to causal mediation analysis (2010). Psychological Methods. 15(4): 309 - 334.

VanderWeele TJ. A unification of mediation and interaction: a 4-way decomposition (2014). Epidemiology. 25(5): 749 - 61.

Schomaker M, Heumann C. Bootstrap inference when using multiple imputation (2018). Statistics in Medicine. 37(14): 2252 - 2266.

VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-Value (2017). Annals of Internal Medicine. 167(4): 268 - 274.

Smith LH, VanderWeele TJ. Mediational E-values: Approximate sensitivity analysis for unmeasured mediator-outcome confounding (2019). Epidemiology. 30(6): 835 - 837.

Carrol RJ, Ruppert D, Stefanski LA, Crainiceanu C. Measurement Error in Nonlinear Models: A Modern Perspective, Second Edition (2006). London: Chapman & Hall.

Cook JR, Stefanski LA. Simulation-extrapolation estimation in parametric measurement error models (1994). Journal of the American Statistical Association, 89(428): 1314 - 1328.

Küchenhoff H, Mwalili SM, Lesaffre E. A general method for dealing with misclassification in regression: the misclassification SIMEX (2006). Biometrics. 62(1): 85 - 96.

Stefanski LA, Cook JR. Simulation-extrapolation: the measurement error jackknife (1995). Journal of the American Statistical Association. 90(432): 1247 - 56.