About the Package
=================

The R package `CMAverse` provides a suite of functions conducting causal
mediation analysis including `cmdag` for DAG visualization, `cmest` for
statistical modeling and `cmsens` for sensitivity analysis.

DAG Visualization
-----------------

`cmdag` visualizes the scientific setting via a DAG.

Statistical Modeling
--------------------

`cmest` implements six causal mediation analysis approaches including
*the regression-based approach* by [Valeri et
al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3659198/) and
[Vanderweele et
al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287269/), *the
weighting-based approach* by [Vanderweele et
al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287269/), *the
inverse odd-ratio weighting approach* by [Tchetgen Tchetgen et
al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3954805/), *the
natural effect model* by [Vansteelandt et
al. (2012)](https://www.degruyter.com/view/journals/em/1/1/article-p131.xml?language=en),
*the marginal structural model* by [VanderWeele et
al. (2009)](https://pubmed.ncbi.nlm.nih.gov/19234398), and *the
g-formula approach* by [Lin et
al. (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5285457/).

Multiple Imputation
-------------------

`cmest` provides options to perform multiple imputation for a dataset
with missing values via the `mice` package, estimate the causal effects
with each of the imputed datasets and pool the results together.

Sensitivity Analysis
--------------------

`cmsens` conducts sensitivity analysis for unmeasured confounding via
the *E-value* approach by [Vanderweele et
al. (2017)](https://pubmed.ncbi.nlm.nih.gov/28693043/) and sensitivity
analysis for measurement error via *regression calibration* by [Carroll
et al. (1995)](https://www.taylorfrancis.com/books/9780429139635),
*SIMEX* by [Cook et
al. (1994)](https://www.jstor.org/stable/2290994?seq=1#metadata_info_tab_contents)
or *MCSIMEX* by [Küchenhoff et
al. (2006)](https://pubmed.ncbi.nlm.nih.gov/16542233/). The sensitivity
analysis for measurement error is currently available for *the
regression-based approach* and *the g-formulaa approach*.

Installation
============

The latest version can be installed via:

    devtools::install_github("LindaValeri/CMAverse")

Load CMAverse:

    library(CMAverse)

Quickstart Guide
================

We illustrate the general workflow of the `CMAverse` package by a quick
example. Firstly, let’s simulate some data and plot the DAG of the
scientific setting. The simulated data contains a binary exposure, a
binary mediator, a continuous mediator, a continuous outcome and two
pre-exposure confounders.

    n <- 500
    C1 <- rnorm(n, mean = 1, sd = 1)
    C2 <- rbinom(n, 1, 0.6)
    C2[which(C2 == 0)] <- "C2_0"
    C2[which(C2 == 1)] <- "C2_1"
    pa <- exp(0.2 - 0.5*C1 + 0.1*(C2 == "C2_1"))/(1 + exp(0.2 - 0.5*C1 + 0.1*(C2 == "C2_1")))
    A <- rbinom(n, 1, pa)
    A[which(A == 0)] <- "control"
    A[which(A == 1)] <- "treat"
    pm <- exp(1 + 0.5*(A == "treat") - 1.5*C1 + 0.5*(C2 == "C2_1"))/
      (1 + exp(1 + 0.5*(A == "treat") - 1.5*C1 + 0.5*(C2 == "C2_1")))
    M1 <- rbinom(n, 1, pm)
    M2 <- rnorm(n, 2 + 0.8*(A == "treat") - M1 + 0.5*C1 + 2*(C2 == "C2_1"), 1)
    Y <- rnorm(n, mean = 0.5 + 0.4*(A == "treat") + 0.5*M1 + 0.6*M2 + 0.3*(A == "treat")*M1 +
                 0.2*(A == "treat")*M2 - 0.3*C1 + 2*(C2=="C2_1"), sd = 1)
    data <- data.frame(A, M1, M2, Y, C1, C2)

The DAG can be plotted using the `cmdag` function.

    cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"), 
          prec = c("C1", "C2"), postc = NULL,
          node = FALSE, text_col = "black")

![](README_files/figure-markdown_strict/unnamed-chunk-4-1.png)

Then, we estimate the causal effects using the `cmest` function. We use
the regression-based approach for illustration. The reference values for
the exposure are set to be 0 and 1. The reference values for the two
mediators are set to be 0.

    est <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                    mediator = c("M1", "M2"), prec = c("C1", "C2"), EMint = TRUE,
                    mreg = list("logistic", "linear"), yreg = "linear",
                    astar = 0, a = 1, mval = list(0, 0),
                    estimation = "imputation", inference = "bootstrap", nboot = 200)

Summarizing and plotting the results:

    summary(est)

    ## 
    ## Causal Mediation Analysis Via the Regression-based Approach
    ##  
    ## Direct counterfactual imputation estimation with 
    ##  bootstrap standard errors, percentile confidence intervals and p-values 
    ##  
    ##               Estimate Std.error   95% CIL 95% CIU  P.val    
    ## cde           0.521487  0.268116 -0.013986   1.042   0.07 .  
    ## pnde          1.197079  0.102549  0.966198   1.390 <2e-16 ***
    ## tnde          1.338773  0.102742  1.134994   1.545 <2e-16 ***
    ## pnie          0.423595  0.070399  0.306427   0.571 <2e-16 ***
    ## tnie          0.565289  0.075643  0.416000   0.718 <2e-16 ***
    ## te            1.762368  0.106978  1.545383   1.955 <2e-16 ***
    ## pm            0.320755  0.039574  0.242353   0.397 <2e-16 ***
    ## intref        0.675592  0.228491  0.230643   1.093 <2e-16 ***
    ## intmed        0.141694  0.051340  0.044126   0.230 <2e-16 ***
    ## cde(prop)     0.295901  0.147598 -0.008183   0.564   0.07 .  
    ## intref(prop)  0.383343  0.136178  0.119061   0.633 <2e-16 ***
    ## intmed(prop)  0.080400  0.029396  0.023272   0.135 <2e-16 ***
    ## pnie(prop)    0.240356  0.036970  0.171181   0.317 <2e-16 ***
    ## pm(overall)   0.320755  0.039574  0.242353   0.397 <2e-16 ***
    ## int(overall)  0.463743  0.162601  0.154421   0.761 <2e-16 ***
    ## pe(overall)   0.704099  0.147598  0.435990   1.008 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    plot(est) +
      theme(axis.text.x = element_text(angle = 30, vjust = 0.8))

![](README_files/figure-markdown_strict/unnamed-chunk-7-1.png)

Lastly, let’s conduct sensitivity analysis for the results. Sensitivity
analysis for unmeasured confounding:

    cmsens(object = est, sens = "uc")

    ## Sensitivity Analysis For Unmeasured Confounding 
    ## 
    ## Evalues on the ratio scale: 
    ##         estRR  lowerRR  upperRR Evalue.estRR Evalue.lowerRR Evalue.upperRR
    ## cde  1.223388 0.998856 1.498393     1.746160       1.000000             NA
    ## pnde 1.588564 1.470019 1.716669     2.555503       2.301245             NA
    ## tnde 1.678019 1.552573 1.813601     2.744662       2.478806             NA
    ## pnie 1.177950 1.116875 1.242365     1.635789       1.478172             NA
    ## tnie 1.244283 1.175099 1.317540     1.795606       1.628706             NA
    ## te   1.976623 1.823003 2.143188     3.366016       3.047885             NA

Assume that the continuous pre-exposure confounder was measured with
error. Sensitivity analysis using regression calibration with a set of
assumed standard deviations of the measurement error 0.1, 0.2 and 0.3:

    me1 <- cmsens(object = est, sens = "me", MEmethod = "rc", 
                  MEvariable = "C1", MEvartype = "con", MEerror = c(0.1, 0.2, 0.3))

Summarizing and plotting the results:

    summary(me1)

    ## Sensitivity Analysis For Measurement Error 
    ##  
    ## The variable measured with error: C1
    ## Type of the variable measured with error: continuous
    ## 
    ## Measurement error 1: 
    ## 0.1
    ## Measurement error correction for measurement error 1: 
    ##              Estimate Std.error  95% CIL 95% CIU  P.val    
    ## cde           0.51964   0.26855 -0.03334   1.058   0.06 .  
    ## pnde          1.18710   0.09957  1.00923   1.397 <2e-16 ***
    ## tnde          1.33236   0.10642  1.14380   1.536 <2e-16 ***
    ## pnie          0.42287   0.07443  0.27507   0.560 <2e-16 ***
    ## tnie          0.56813   0.07559  0.43909   0.715 <2e-16 ***
    ## te            1.75523   0.10374  1.58755   1.938 <2e-16 ***
    ## pm            0.32368   0.03891  0.24456   0.396 <2e-16 ***
    ## intref        0.66746   0.24402  0.20852   1.175 <2e-16 ***
    ## intmed        0.14526   0.05419  0.04309   0.246   0.02 *  
    ## cde(prop)     0.29605   0.14848 -0.02002   0.579   0.06 .  
    ## intref(prop)  0.38027   0.14270  0.11083   0.692 <2e-16 ***
    ## intmed(prop)  0.08276   0.03030  0.02382   0.137   0.02 *  
    ## pnie(prop)    0.24092   0.04019  0.16061   0.315 <2e-16 ***
    ## pm(overall)   0.32368   0.03891  0.24456   0.396 <2e-16 ***
    ## int(overall)  0.46303   0.16970  0.13510   0.813 <2e-16 ***
    ## pe(overall)   0.70395   0.14848  0.42066   1.020 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------
    ## 
    ## Measurement error 2: 
    ## 0.2
    ## Measurement error correction for measurement error 2: 
    ##              Estimate Std.error 95% CIL 95% CIU  P.val    
    ## cde           0.51370   0.24178 0.04764   0.987   0.02 *  
    ## pnde          1.20586   0.09923 0.99641   1.370 <2e-16 ***
    ## tnde          1.33409   0.10586 1.13356   1.532 <2e-16 ***
    ## pnie          0.44440   0.06725 0.29954   0.549 <2e-16 ***
    ## tnie          0.57264   0.08255 0.41022   0.717 <2e-16 ***
    ## te            1.77850   0.10696 1.54511   1.958 <2e-16 ***
    ## pm            0.32198   0.04198 0.24938   0.407 <2e-16 ***
    ## intref        0.69216   0.21433 0.22117   1.067 <2e-16 ***
    ## intmed        0.12823   0.05112 0.05458   0.248 <2e-16 ***
    ## cde(prop)     0.28884   0.13703 0.02623   0.569   0.02 *  
    ## intref(prop)  0.38918   0.12273 0.12503   0.626 <2e-16 ***
    ## intmed(prop)  0.07210   0.02793 0.03132   0.135 <2e-16 ***
    ## pnie(prop)    0.24988   0.03640 0.17128   0.309 <2e-16 ***
    ## pm(overall)   0.32198   0.04198 0.24938   0.407 <2e-16 ***
    ## int(overall)  0.46129   0.14745 0.16730   0.737 <2e-16 ***
    ## pe(overall)   0.71116   0.13703 0.43099   0.974 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------
    ## 
    ## Measurement error 3: 
    ## 0.3
    ## Measurement error correction for measurement error 3: 
    ##               Estimate Std.error   95% CIL 95% CIU  P.val    
    ## cde           0.502278  0.286518 -0.014643   1.088   0.06 .  
    ## pnde          1.174673  0.101812  0.983094   1.349 <2e-16 ***
    ## tnde          1.319282  0.104802  1.119086   1.505 <2e-16 ***
    ## pnie          0.444247  0.074372  0.327600   0.604 <2e-16 ***
    ## tnie          0.588857  0.078207  0.446115   0.743 <2e-16 ***
    ## te            1.763530  0.101974  1.574002   1.939 <2e-16 ***
    ## pm            0.333908  0.041084  0.262401   0.413 <2e-16 ***
    ## intref        0.672395  0.253474  0.152410   1.080   0.01 ** 
    ## intmed        0.144610  0.054859  0.032333   0.247   0.01 ** 
    ## cde(prop)     0.284814  0.159225 -0.008935   0.620   0.06 .  
    ## intref(prop)  0.381278  0.148396  0.092493   0.645   0.01 ** 
    ## intmed(prop)  0.082000  0.030976  0.016967   0.137   0.01 ** 
    ## pnie(prop)    0.251908  0.040025  0.191593   0.346 <2e-16 ***
    ## pm(overall)   0.333908  0.041084  0.262401   0.413 <2e-16 ***
    ## int(overall)  0.463278  0.176766  0.109666   0.768   0.01 ** 
    ## pe(overall)   0.715186  0.159225  0.379946   1.009 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------

    plot(me1) +
      theme(axis.text.x = element_text(angle = 30, vjust = 0.8))

![](README_files/figure-markdown_strict/unnamed-chunk-11-1.png)

Then, assume that the exposure was measured with error. Sensitivity
analysis using MCSIMEX with two assumed misclassification matrices:

    me2 <- cmsens(object = est, sens = "me", MEmethod = "simex", MEvariable = "A", MEvartype = "cat", 
                  MEerror = list(matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2), 
                                 matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)))

Summarizing and plotting the results:

    summary(me2)

    ## Sensitivity Analysis For Measurement Error 
    ##  
    ## The variable measured with error: A
    ## Type of the variable measured with error: categorical
    ## 
    ## Measurement error 1: 
    ##      [,1] [,2]
    ## [1,] 0.95 0.05
    ## [2,] 0.05 0.95
    ## 
    ## Measurement error correction for measurement error 1: 
    ##              Estimate Std.error 95% CIL 95% CIU  P.val    
    ## cde           0.62807   0.28731 0.10604   1.255   0.02 *  
    ## pnde          1.37069   0.12126 1.10864   1.585 <2e-16 ***
    ## tnde          1.55006   0.13849 1.26308   1.786 <2e-16 ***
    ## pnie          0.41672   0.08066 0.27627   0.566 <2e-16 ***
    ## tnie          0.59609   0.08902 0.41907   0.771 <2e-16 ***
    ## te            1.96678   0.12625 1.69858   2.184 <2e-16 ***
    ## pm            0.30308   0.04189 0.22502   0.388 <2e-16 ***
    ## intref        0.74262   0.25963 0.25009   1.178   0.01 ** 
    ## intmed        0.17937   0.06935 0.04687   0.320   0.02 *  
    ## cde(prop)     0.31934   0.14362 0.05134   0.626   0.02 *  
    ## intref(prop)  0.37758   0.13206 0.12435   0.605   0.01 ** 
    ## intmed(prop)  0.09120   0.03368 0.02392   0.151   0.02 *  
    ## pnie(prop)    0.21188   0.04167 0.14156   0.284 <2e-16 ***
    ## pm(overall)   0.30308   0.04189 0.22502   0.388 <2e-16 ***
    ## int(overall)  0.46878   0.16201 0.16029   0.746   0.01 ** 
    ## pe(overall)   0.68066   0.14362 0.37429   0.949 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------
    ## 
    ## Measurement error 2: 
    ##      [,1] [,2]
    ## [1,]  0.9  0.1
    ## [2,]  0.1  0.9
    ## 
    ## Measurement error correction for measurement error 2: 
    ##              Estimate Std.error 95% CIL 95% CIU  P.val    
    ## cde           0.79094   0.37141 0.16541   1.549   0.03 *  
    ## pnde          1.57150   0.14836 1.26901   1.867 <2e-16 ***
    ## tnde          1.77735   0.15300 1.48966   2.065 <2e-16 ***
    ## pnie          0.38910   0.07948 0.24675   0.565 <2e-16 ***
    ## tnie          0.59495   0.08720 0.43769   0.754 <2e-16 ***
    ## te            2.16645   0.13708 1.89950   2.419 <2e-16 ***
    ## pm            0.27462   0.04094 0.20024   0.356 <2e-16 ***
    ## intref        0.78056   0.29852 0.15330   1.283   0.02 *  
    ## intmed        0.20585   0.08664 0.02830   0.358   0.03 *  
    ## cde(prop)     0.36508   0.16400 0.08034   0.691   0.03 *  
    ## intref(prop)  0.36030   0.14421 0.06418   0.599   0.02 *  
    ## intmed(prop)  0.09502   0.03974 0.01265   0.164   0.03 *  
    ## pnie(prop)    0.17960   0.03736 0.11354   0.259 <2e-16 ***
    ## pm(overall)   0.27462   0.04094 0.20024   0.356 <2e-16 ***
    ## int(overall)  0.45531   0.17933 0.08554   0.757   0.01 ** 
    ## pe(overall)   0.63492   0.16400 0.30928   0.920 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------

    plot(me2) +
      theme(axis.text.x = element_text(angle = 30, vjust = 0.8))

![](README_files/figure-markdown_strict/unnamed-chunk-14-1.png)

References
==========

Valeri L, Vanderweele TJ (2013). Mediation analysis allowing for
exposure-mediator interactions and causal interpretation: theoretical
assumptions and implementation with SAS and SPSS macros. Psychological
Methods. 18(2): 137 - 150.

VanderWeele TJ, Vansteelandt S (2014). Mediation analysis with multiple
mediators. Epidemiologic Methods. 2(1): 95 - 115.

Tchetgen Tchetgen EJ (2013). Inverse odds ratio-weighted estimation for
causal mediation analysis. Statistics in medicine. 32: 4567 - 4580.

Nguyen QC, Osypuk TL, Schmidt NM, Glymour MM, Tchetgen Tchetgen EJ.
Practical guidance for conducting mediation analysis with multiple
mediators using inverse odds ratio weighting (2015). American Journal of
Epidemiology. 181(5): 349 - 356.

VanderWeele TJ. Marginal structural models for the estimation of direct
and indirect effects (2009). Epidemiology. 20(1): 18 - 26.

VanderWeele TJ, Tchetgen Tchetgen EJ (2017). Mediation analysis with
time varying exposures and mediators. Journal of the Royal Statistical
Society: Series B (Statistical Methodology). 79(3): 917 - 938.

Lin SH, Young J, Logan R, Tchetgen Tchetgen EJ, VanderWeele TJ (2017).
Parametric mediational g-formula approach to mediation analysis with
time-varying exposures, mediators, and confounders. Epidemiology. 28:
266 - 274.

Vansteelandt S, Bekaert M, Lange T. (2012). Imputation Strategies for
the Estimation of Natural Direct and Indirect Effects. Epidemiologic
Methods. 1(1): 131 - 158.

Steen J, Loeys T, Moerkerke B, Vansteelandt S (2017). Medflex: an R
package for flexible mediation analysis using natural effect models.
Journal of Statistical Software. 76(11).

VanderWeele TJ. A unification of mediation and interaction: a 4-way
decomposition (2014). Epidemiology. 25(5): 749 - 61.

Imai K, Keele L, Tingley D. A general approach to causal mediation
analysis (2010). Psychological Methods. 15(4): 309 - 334.

Schomaker M, Heumann C. Bootstrap inference when using multiple
imputation (2018). Statistics in Medicine. 37(14): 2252 - 2266.

VanderWeele TJ, Ding P. Sensitivity analysis in observational research:
introducing the E-Value (2017). Annals of Internal Medicine. 167(4): 268
- 274.

Smith LH, VanderWeele TJ. Mediational E-values: Approximate sensitivity
analysis for unmeasured mediator-outcome confounding (2019).
Epidemiology. 30(6): 835 - 837.

Carrol RJ, Ruppert D, Stefanski LA, Crainiceanu C. Measurement Error in
Nonlinear Models: A Modern Perspective, Second Edition (2006). London:
Chapman & Hall.

Cook JR, Stefanski LA. Simulation-extrapolation estimation in parametric
measurement error models (1994). Journal of the American Statistical
Association, 89(428): 1314 - 1328.

Küchenhoff H, Mwalili SM, Lesaffre E. A general method for dealing with
misclassification in regression: the misclassification SIMEX (2006).
Biometrics. 62(1): 85 - 96.

Stefanski LA, Cook JR. Simulation-extrapolation: the measurement error
jackknife (1995). Journal of the American Statistical Association.
90(432): 1247 - 56.
