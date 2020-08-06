
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CMAverse <img src="man/figures/logo.png" align="right" width="300" />

## About the Package

The R package `CMAverse` provides a suite of functions conducting causal
mediation analysis including `cmdag` for DAG visualization, `cmest` for
statistical modeling and `cmsens` for sensitivity analysis.

### DAG Visualization

`cmdag` visualizes the scientific setting via a DAG.

### Statistical Modeling

`cmest` implements six causal mediation analysis approaches including
*the regression-based approach* by [Valeri et
al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3659198/) and
[Vanderweele et
al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287269/),
*the weighting-based approach* by [Vanderweele et
al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287269/),
*the inverse odd-ratio weighting approach* by [Tchetgen Tchetgen et
al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3954805/),
*the natural effect model* by [Vansteelandt et
al. (2012)](https://www.degruyter.com/view/journals/em/1/1/article-p131.xml?language=en),
*the marginal structural model* by [VanderWeele et
al. (2009)](https://pubmed.ncbi.nlm.nih.gov/19234398), and *the
g-formula approach* by [Lin et
al. (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5285457/).

|                                              |   rb    |   wb   | iorw | ne |  msm   | gformula\[1\] |
| -------------------------------------------- | :-----: | :----: | :--: | :-: | :----: | :-----------: |
| Continuous Y\[2\]                            |    √    |   √    |  √   | √  |   √    |       √       |
| Binary Y                                     |    √    |   √    |  √   | √  |   √    |       √       |
| Count Y                                      |    √    |   √    |  √   | √  |   √    |       √       |
| Nominal Y                                    |    √    |   √    |  √   | ×  |   √    |       √       |
| Ordinal Y                                    |    √    |   √    |  √   | ×  |   √    |       √       |
| Survival Y                                   |    √    |   √    |  √   | ×  |   √    |       √       |
| Continuous M                                 |    √    |   √    |  √   | √  |   ×    |       √       |
| Binary M                                     |    √    |   √    |  √   | √  |   √    |       √       |
| Nominal M                                    |    √    |   √    |  √   | √  |   √    |       √       |
| Ordinal M                                    |    √    |   √    |  √   | √  |   √    |       √       |
| Count M                                      |    √    |   √    |  √   | √  |   ×    |       √       |
| M of Any Type                                |    ×    |   √    |  √   | √  |   ×    |       ×       |
| Continuous A                                 |    √    | ×\[3\] |  ×   | √  | ×\[4\] |       √       |
| Binary A                                     |    √    |   √    |  √   | √  |   √    |       √       |
| Nominal A                                    |    √    |   √    |  √   | √  |   √    |       √       |
| Ordinal A                                    |    √    |   √    |  √   | √  |   √    |       √       |
| Count A                                      |    √    | ×\[5\] |  ×   | √  | ×\[6\] |       √       |
| Multiple Mediators                           |    √    |   √    |  √   | √  |   √    |       √       |
| Post-exposure Confounding                    |    ×    |   ×    |  ×   | ×  |   √    |       √       |
| 2-way Decomposition                          |    √    |   √    |  √   | √  |   √    |       √       |
| 4-way Decomposition                          |    √    |   √    |  ×   | √  |   √    |       √       |
| Estimation: Closed-form Parameter Function   | √\[7\]  |   ×    |  ×   | ×  |   ×    |       ×       |
| Estimation: Direct Counterfactual Imputation |    √    |   √    |  √   | √  |   √    |       √       |
| Inference: Delta Method                      | √\[8\]  |   ×    |  ×   | ×  |   ×    |       ×       |
| Inference: Bootstrapping                     |    √    |   √    |  √   | √  |   √    |       √       |
| Marginal Effects                             | √\[9\]  |   √    |  √   | √  |   √    |       √       |
| Effects Conditional On C                     | √\[10\] |   ×    |  ×   | ×  |   ×    |       ×       |

Table: Supported Data Types and Functionalities of `cmest`

### Multiple Imputation

`cmest` provides options to perform multiple imputation for a dataset
with missing values via the `mice` package, estimate the causal effects
with each of the imputed datasets and pool the results together.

### Sensitivity Analysis

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
regression-based approach* and *the g-formula approach*.

## Installation

The latest version can be installed via:

``` r
devtools::install_github("LindaValeri/CMAverse")
```

Load `CMAverse`:

``` r
library(CMAverse)
```

## Quickstart Guide

We illustrate the general workflow of the `CMAverse` package by a quick
example. Firstly, let’s simulate some data and plot the DAG of the
scientific setting. The simulated dataset contains a binary exposure, a
binary mediator, a continuous mediator, a continuous outcome and two
pre-exposure confounders.

``` r
n <- 30
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
```

The DAG can be plotted using the `cmdag` function.

``` r
cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"), 
      prec = c("C1", "C2"), postc = NULL,
      node = FALSE, text_col = "black")
```

![](man/figuresplot_dag-1.png)<!-- -->

Then, we estimate the causal effects using the `cmest` function. We use
the regression-based approach for illustration. The reference values for
the exposure are set to be 0 and 1. The reference values for the two
mediators are set to be 0.

``` r
est <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                mediator = c("M1", "M2"), prec = c("C1", "C2"), EMint = TRUE,
                mreg = list("logistic", "linear"), yreg = "linear",
                astar = 0, a = 1, mval = list(0, 0),
                estimation = "imputation", inference = "bootstrap", nboot = 5)
```

Summarizing and plotting the results:

``` r
summary(est)
#> 
#> Causal Mediation Analysis Via the Regression-based Approach
#>  
#> Direct counterfactual imputation estimation with 
#>  bootstrap standard errors, percentile confidence intervals and p-values 
#>  
#>              Estimate Std.error  95% CIL 95% CIU  P.val    
#> cde           1.79585   1.61455 -1.53279   2.100    0.4    
#> pnde          1.67472   0.42547  1.00895   2.059 <2e-16 ***
#> tnde          1.64392   0.49248  1.10432   2.248 <2e-16 ***
#> pnie          0.54145   0.31529  0.35501   1.013 <2e-16 ***
#> tnie          0.51065   0.48752  0.49978   1.560 <2e-16 ***
#> te            2.18537   0.28037  2.01052   2.614 <2e-16 ***
#> pm            0.23367   0.17623  0.19991   0.601 <2e-16 ***
#> intref       -0.12113   1.22188 -0.04504   2.542    0.8    
#> intmed       -0.03080   0.63256 -0.37638   1.134    0.8    
#> cde(prop)     0.82176   0.63819 -0.57693   0.816    0.4    
#> intref(prop) -0.05543   0.46847 -0.01730   0.983    0.8    
#> intmed(prop) -0.01410   0.25527 -0.18832   0.437    0.8    
#> pnie(prop)    0.24776   0.14997  0.13637   0.462 <2e-16 ***
#> pm(overall)   0.23367   0.17623  0.19991   0.601 <2e-16 ***
#> int(overall) -0.06952   0.69954 -0.05447   1.402    0.8    
#> pe(overall)   0.17824   0.63819  0.18425   1.577 <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
plot(est) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8))
```

![](man/figuresplot_cmest-1.png)<!-- -->

Lastly, let’s conduct sensitivity analysis for the results. Sensitivity
analysis for unmeasured confounding:

``` r
cmsens(object = est, sens = "uc")
#> Sensitivity Analysis For Unmeasured Confounding 
#> 
#> Evalues on the ratio scale: 
#>         estRR   lowerRR  upperRR Evalue.estRR Evalue.lowerRR Evalue.upperRR
#> cde  2.249041 0.5407336 9.354305     3.925092       1.000000             NA
#> pnde 2.129394 1.4626200 3.100134     3.680175       2.285200             NA
#> tnde 2.099994 1.3595681 3.243658     3.619856       2.058751             NA
#> pnie 1.276815 0.9665999 1.686589     1.871325       1.000000             NA
#> tnie 1.259187 0.8187949 1.936445     1.830470       1.000000             NA
#> te   2.681304 2.0934113 3.434295     4.804530       3.606342             NA
```

Assume that the continuous pre-exposure confounder was measured with
error. Sensitivity analysis using regression calibration with a set of
assumed standard deviations of the measurement error 0.1, 0.2 and 0.3:

``` r
me1 <- cmsens(object = est, sens = "me", MEmethod = "rc", 
              MEvariable = "C1", MEvartype = "con", MEerror = c(0.1, 0.2, 0.3))
```

Summarizing and plotting the results:

``` r
summary(me1)
#> Sensitivity Analysis For Measurement Error 
#>  
#> The variable measured with error: C1
#> Type of the variable measured with error: continuous
#> 
#> Measurement error 1: 
#> 0.1
#> Measurement error correction for measurement error 1: 
#>              Estimate Std.error  95% CIL 95% CIU  P.val    
#> cde           1.80927   2.34837 -0.97823   5.011    0.4    
#> pnde          1.71727   0.53662  1.17056   2.329 <2e-16 ***
#> tnde          1.68646   0.69584  0.87117   2.494 <2e-16 ***
#> pnie          0.53380   0.41713  0.05375   1.003 <2e-16 ***
#> tnie          0.50299   0.61884  0.03566   1.432 <2e-16 ***
#> te            2.22026   0.36509  1.86147   2.715 <2e-16 ***
#> pm            0.22655   0.24556  0.01681   0.549 <2e-16 ***
#> intref       -0.09200   2.20948 -3.28280   2.272    0.8    
#> intmed       -0.03081   0.71733 -0.87866   0.954    0.8    
#> cde(prop)     0.81489   1.22351 -0.36883   2.737    0.4    
#> intref(prop) -0.04144   1.06365 -1.79622   0.862    0.8    
#> intmed(prop) -0.01388   0.33629 -0.48456   0.362    0.8    
#> pnie(prop)    0.24042   0.21482  0.02055   0.535 <2e-16 ***
#> pm(overall)   0.22655   0.24556  0.01681   0.549 <2e-16 ***
#> int(overall) -0.05531   1.38729 -2.27170   1.224    0.8    
#> pe(overall)   0.18511   1.22351 -1.73676   1.369    0.8    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ----------------------------------------------------------------
#> 
#> Measurement error 2: 
#> 0.2
#> Measurement error correction for measurement error 2: 
#>              Estimate Std.error  95% CIL 95% CIU  P.val    
#> cde           1.85567   2.56296 -0.92920   5.628    0.4    
#> pnde          1.73427   0.65656  1.02717   2.595 <2e-16 ***
#> tnde          1.69365   0.43897  0.96058   1.887 <2e-16 ***
#> pnie          0.53945   0.61087  0.22719   1.724 <2e-16 ***
#> tnie          0.49883   0.43985  0.09402   1.075 <2e-16 ***
#> te            2.23310   0.35178  1.93156   2.762 <2e-16 ***
#> pm            0.22338   0.20873  0.04250   0.507 <2e-16 ***
#> intref       -0.12139   1.98068 -3.06329   2.001    0.4    
#> intmed       -0.04063   0.92740 -1.54648   0.815    0.4    
#> cde(prop)     0.83098   0.99070 -0.45356   2.046    0.4    
#> intref(prop) -0.05436   0.81552 -1.11261   0.949    0.4    
#> intmed(prop) -0.01819   0.37079 -0.55848   0.387    0.4    
#> pnie(prop)    0.24157   0.21929  0.10746   0.634 <2e-16 ***
#> pm(overall)   0.22338   0.20873  0.04250   0.507 <2e-16 ***
#> int(overall) -0.07255   1.18358 -1.66760   1.336    0.4    
#> pe(overall)   0.16902   0.99070 -1.04562   1.454    0.8    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ----------------------------------------------------------------
#> 
#> Measurement error 3: 
#> 0.3
#> Measurement error correction for measurement error 3: 
#>              Estimate Std.error  95% CIL 95% CIU  P.val    
#> cde           1.96205   2.54447  1.31162   7.215 <2e-16 ***
#> pnde          1.85066   1.38063  1.26316   4.534 <2e-16 ***
#> tnde          1.77078   0.72456  1.25796   3.044 <2e-16 ***
#> pnie          0.62185   0.10614  0.19210   0.427 <2e-16 ***
#> tnie          0.54197   0.76464 -1.27433   0.408    0.8    
#> te            2.39263   0.67340  1.65641   3.293 <2e-16 ***
#> pm            0.22652   0.26426 -0.37945   0.242    0.8    
#> intref       -0.11139   1.17818 -2.68128  -0.037 <2e-16 ***
#> intmed       -0.07988   0.69668 -1.49867   0.000    0.4    
#> cde(prop)     0.82004   0.60535  0.77631   2.187 <2e-16 ***
#> intref(prop) -0.04656   0.34143 -0.80764  -0.018 <2e-16 ***
#> intmed(prop) -0.03339   0.20352 -0.44878   0.000    0.4    
#> pnie(prop)    0.25990   0.07376  0.06933   0.242 <2e-16 ***
#> pm(overall)   0.22652   0.26426 -0.37945   0.242    0.8    
#> int(overall) -0.07994   0.54330 -1.25642  -0.018 <2e-16 ***
#> pe(overall)   0.17996   0.60535 -1.18709   0.224    0.8    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ----------------------------------------------------------------
```

``` r
plot(me1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8))
```

![](man/figuresplot_cmsens_me_con-1.png)<!-- -->

Then, assume that the exposure was measured with error. Sensitivity
analysis using MCSIMEX with two assumed misclassification matrices:

``` r
me2 <- cmsens(object = est, sens = "me", MEmethod = "simex", MEvariable = "A", 
              MEvartype = "cat", B = 5,
              MEerror = list(matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2), 
                             matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)))
```

Summarizing and plotting the results:

``` r
summary(me2)
#> Sensitivity Analysis For Measurement Error 
#>  
#> The variable measured with error: A
#> Type of the variable measured with error: categorical
#> 
#> Measurement error 1: 
#>      [,1] [,2]
#> [1,] 0.95 0.05
#> [2,] 0.05 0.95
#> 
#> Measurement error correction for measurement error 1: 
#>              Estimate Std.error  95% CIL 95% CIU  P.val    
#> cde           1.07926   3.23894 -0.22949   7.928    0.4    
#> pnde          2.06299   1.01671  1.63218   4.129 <2e-16 ***
#> tnde          2.52009   0.62224  1.93876   3.470 <2e-16 ***
#> pnie          0.25112   0.50955 -0.88360   0.341    0.8    
#> tnie          0.70822   0.86597 -1.73235   0.430    0.4    
#> te            2.77121   0.43913  1.99121   3.064 <2e-16 ***
#> pm            0.25557   0.37790 -0.72618   0.223    0.4    
#> intref        0.98373   2.35650 -3.82223   1.885    0.8    
#> intmed        0.45711   0.74052 -1.31114   0.386    0.8    
#> cde(prop)     0.38945   1.36094 -0.14853   3.303    0.4    
#> intref(prop)  0.35498   1.02916 -1.59513   0.955    0.8    
#> intmed(prop)  0.16495   0.31141 -0.54516   0.197    0.8    
#> pnie(prop)    0.09062   0.19737 -0.34397   0.121    0.8    
#> pm(overall)   0.25557   0.37790 -0.72618   0.223    0.4    
#> int(overall)  0.51993   1.32466 -2.13546   1.152    0.8    
#> pe(overall)   0.61055   1.36094 -2.30344   1.149    0.8    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ----------------------------------------------------------------
#> 
#> Measurement error 2: 
#>      [,1] [,2]
#> [1,]  0.9  0.1
#> [2,]  0.1  0.9
#> 
#> Measurement error correction for measurement error 2: 
#>              Estimate Std.error  95% CIL 95% CIU  P.val    
#> cde           2.21042   1.92113  0.78428   5.078 <2e-16 ***
#> pnde          1.68468   0.68298  1.46855   3.034 <2e-16 ***
#> tnde          1.27483   0.71259  1.68595   3.451 <2e-16 ***
#> pnie          1.14790   0.44120 -0.55081   0.498    0.4    
#> tnie          0.73805   0.40885 -0.45716   0.563    0.8    
#> te            2.42273   0.43257  1.99550   2.973 <2e-16 ***
#> pm            0.30464   0.18473 -0.18031   0.282    0.8    
#> intref       -0.52574   1.40320 -2.04317   0.933    0.8    
#> intmed       -0.40985   0.65335 -0.69858   0.820    0.8    
#> cde(prop)     0.91237   0.71146  0.37023   1.969 <2e-16 ***
#> intref(prop) -0.21700   0.54362 -0.79305   0.347    0.8    
#> intmed(prop) -0.16917   0.24188 -0.27161   0.285    0.8    
#> pnie(prop)    0.47380   0.15896 -0.18628   0.173    0.4    
#> pm(overall)   0.30464   0.18473 -0.18031   0.282    0.8    
#> int(overall) -0.38617   0.78017 -1.06466   0.615    0.8    
#> pe(overall)   0.08763   0.71146 -0.96896   0.630    0.8    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ----------------------------------------------------------------
```

``` r
plot(me2) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8))
```

![](man/figuresplot_cmsens_me_cat-1.png)<!-- -->

## References

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

Valeri L, Lin X, VanderWeele TJ. Mediation analysis when a continuous
mediator is measured with error and the outcome follows a generalized
linear model (2014). Statistics in medicine, 33(28): 4875–4890.

1.  rb: the regression-based approach; wb: the weighting-based approach;
    iorw: the inverse odds ratio weighting approach; ne: the natural
    effect model; msm: the marginal structural model; gformula: the
    g-formula approach.

2.  Y denotes the outcome, A denotes the exposure, M denotes the
    mediator(s) and C denotes the pre-exposure confounder(s).

3.  continuous A is not supported when C is not empty; otherwise, it is
    supported.

4.  continuous A is not supported when C is not empty; otherwise, it is
    supported.

5.  count A is not supported when C is not empty; otherwise, it is
    supported.

6.  count A is not supported when C is not empty; otherwise, it is
    supported.

7.  closed-form parameter function estimation only supports the
    regression-based approach and a single mediator.

8.  delta method inference is available only when closed-form parameter
    function estimation is used.

9.  marginal effects are estimated when direct counterfactual imputation
    estimation is used.

10. conditional effects are estimated when closed-form parameter
    function estimation is used.
