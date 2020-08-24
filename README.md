<!-- README.md is generated from README.Rmd. Please edit that file -->

CMAverse <img src="man/figures/logo.png" align="right" width="240" />
=====================================================================

[![Project Status:
Active](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Travis build
status](https://travis-ci.com/BS1125/CMAverse.svg?branch=master)](https://travis-ci.com/BS1125/CMAverse)
[![R build
status](https://github.com/BS1125/CMAverse/workflows/R-CMD-check/badge.svg)](https://github.com/BS1125/CMAverse/actions)
[![Codecov test
coverage](https://codecov.io/gh/BS1125/CMAverse/branch/master/graph/badge.svg)](https://codecov.io/gh/BS1125/CMAverse?branch=master)

About the Package
-----------------

The R package `CMAverse` provides a suite of functions for reproducible
causal mediation analysis including `cmdag` for DAG visualization,
`cmest` for statistical modeling and `cmsens` for sensitivity analysis.
See the [package website](https://bs1125.github.io/CMAverse/) for
details.

### DAG Visualization

`cmdag` visualizes the scientific setting via a DAG.

### Statistical Modeling

`cmest` implements six causal mediation analysis approaches including
*the regression-based approach* by [Valeri et
al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3659198/) and
[VanderWeele et
al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287269/), *the
weighting-based approach* by [VanderWeele et
al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287269/), *the
inverse odd-ratio weighting approach* by [Tchetgen Tchetgen
(2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3954805/), *the
natural effect model* by [Vansteelandt et
al. (2012)](https://www.degruyter.com/view/journals/em/1/1/article-p131.xml?language=en),
*the marginal structural model* by [VanderWeele et
al. (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5560424/), and
*the g-formula approach* by [Lin et
al. (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5285457/).

<table>
<caption>Table: Supported Data Types and Functionalities of <code>cmest</code></caption>
<thead>
<tr class="header">
<th></th>
<th style="text-align: center;">rb</th>
<th style="text-align: center;">wb</th>
<th style="text-align: center;">iorw</th>
<th style="text-align: center;">ne</th>
<th style="text-align: center;">msm</th>
<th style="text-align: center;">gformula<a href="#fn1" class="footnote-ref" id="fnref1" role="doc-noteref"><sup>1</sup></a></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Continuous Y<a href="#fn2" class="footnote-ref" id="fnref2" role="doc-noteref"><sup>2</sup></a></td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Binary Y</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Count Y</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Nominal Y</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Ordinal Y</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Survival Y</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Continuous M</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Binary M</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Nominal M</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Ordinal M</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Count M</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>M of Any Type</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
</tr>
<tr class="odd">
<td>Continuous A</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×<a href="#fn3" class="footnote-ref" id="fnref3" role="doc-noteref"><sup>3</sup></a></td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×<a href="#fn4" class="footnote-ref" id="fnref4" role="doc-noteref"><sup>4</sup></a></td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Binary A</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Nominal A</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Ordinal A</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Count A</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×<a href="#fn5" class="footnote-ref" id="fnref5" role="doc-noteref"><sup>5</sup></a></td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×<a href="#fn6" class="footnote-ref" id="fnref6" role="doc-noteref"><sup>6</sup></a></td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Multiple Mediators</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Post-exposure Confounding</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>2-way Decomposition</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>4-way Decomposition</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Estimation: Closed-form Parameter Function</td>
<td style="text-align: center;">√<a href="#fn7" class="footnote-ref" id="fnref7" role="doc-noteref"><sup>7</sup></a></td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
</tr>
<tr class="odd">
<td>Estimation: Direct Counterfactual Imputation</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Inference: Delta Method</td>
<td style="text-align: center;">√<a href="#fn8" class="footnote-ref" id="fnref8" role="doc-noteref"><sup>8</sup></a></td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
</tr>
<tr class="odd">
<td>Inference: Bootstrapping</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Marginal Effects</td>
<td style="text-align: center;">√<a href="#fn9" class="footnote-ref" id="fnref9" role="doc-noteref"><sup>9</sup></a></td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Effects Conditional On C</td>
<td style="text-align: center;">√<a href="#fn10" class="footnote-ref" id="fnref10" role="doc-noteref"><sup>10</sup></a></td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
</tr>
</tbody>
</table>
<section class="footnotes" role="doc-endnotes">
<ol>
<li id="fn1" role="doc-endnote"><p>rb: the regression-based approach; wb: the weighting-based approach; iorw: the inverse odds ratio weighting approach; ne: the natural effect model; msm: the marginal structural model; gformula: the g-formula approach.<a href="#fnref1" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn2" role="doc-endnote"><p>Y denotes the outcome, A denotes the exposure, M denotes the mediator(s) and C denotes the pre-exposure confounder(s).<a href="#fnref2" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn3" role="doc-endnote"><p>continuous A is not supported when C is not empty; otherwise, it is supported.<a href="#fnref3" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn4" role="doc-endnote"><p>continuous A is not supported when C is not empty; otherwise, it is supported.<a href="#fnref4" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn5" role="doc-endnote"><p>count A is not supported when C is not empty; otherwise, it is supported.<a href="#fnref5" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn6" role="doc-endnote"><p>count A is not supported when C is not empty; otherwise, it is supported.<a href="#fnref6" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn7" role="doc-endnote"><p>closed-form parameter function estimation only supports the regression-based approach and a single mediator.<a href="#fnref7" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn8" role="doc-endnote"><p>delta method inference is available only when closed-form parameter function estimation is used.<a href="#fnref8" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn9" role="doc-endnote"><p>marginal effects are estimated when direct counterfactual imputation estimation is used.<a href="#fnref9" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn10" role="doc-endnote"><p>conditional effects are estimated when closed-form parameter function estimation is used.<a href="#fnref10" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
</ol>
</section>

### Multiple Imputation

`cmest` provides options to perform multiple imputation for a dataset
with missing values via the `mice` package, estimate the causal effects
with each of the imputed datasets and pool the results together.

### Sensitivity Analysis

`cmsens` conducts sensitivity analysis for unmeasured confounding via
the *E-value* approach by [VanderWeele et
al. (2017)](https://pubmed.ncbi.nlm.nih.gov/28693043/) and [Smith et
al. (2019)](https://pubmed.ncbi.nlm.nih.gov/31348008/), and sensitivity
analysis for measurement error via *regression calibration* by [Carroll
et al. (1995)](https://www.taylorfrancis.com/books/9780429139635) and
*SIMEX* by [Cook et
al. (1994)](https://www.jstor.org/stable/2290994?seq=1#metadata_info_tab_contents)
and [Küchenhoff et
al. (2006)](https://pubmed.ncbi.nlm.nih.gov/16542233/). The sensitivity
analysis for measurement error is currently available for *the
regression-based approach* and *the g-formula approach*.

Installation
------------

The latest version can be installed via:

    devtools::install_github("BS1125/CMAverse")

Load `CMAverse`:

    library(CMAverse)

Quickstart Guide
----------------

We illustrate the general workflow of the `CMAverse` package by a quick
example. Firstly, let’s simulate some data and plot the DAG of the
scientific setting. The simulated dataset contains a binary exposure, a
binary mediator, a continuous mediator, a continuous outcome and two
pre-exposure confounders.

    n <- 100
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
          basec = c("C1", "C2", "C3"), postc = NULL,
          node = FALSE, text_col = "black")

![](man/figures/plot_dag-1.png)

Then, we estimate the causal effects using the `cmest` function. We use
the regression-based approach for illustration. The reference values for
the exposure are set to be 0 and 1. The reference values for the two
mediators are set to be 0.

    est <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                    mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                    mreg = list("logistic", "linear"), yreg = "linear",
                    astar = 0, a = 1, mval = list(0, 0),
                    estimation = "imputation", inference = "bootstrap", nboot = 20)

Summarizing and plotting the results:

    summary(est)

    ## # Outcome Regression: 
    ## 
    ## Call:
    ## glm(formula = Y ~ A + M1 + M2 + A * M1 + A * M2 + C1 + C2, family = gaussian(), 
    ##     data = getCall(x$regsumm$yreg)$data, weights = getCall(x$regsumm$yreg)$weights)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -2.04693  -0.53712   0.02205   0.58158   2.27407  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.26713    0.43215   2.932  0.00425 ** 
    ## Atreat      -0.86779    0.61501  -1.411  0.16161    
    ## M1          -0.04839    0.31170  -0.155  0.87697    
    ## M2           0.53167    0.11817   4.499 1.99e-05 ***
    ## C1          -0.30268    0.10606  -2.854  0.00534 ** 
    ## C2C2_1       1.83641    0.29854   6.151 1.97e-08 ***
    ## Atreat:M1    1.03168    0.39837   2.590  0.01117 *  
    ## Atreat:M2    0.35068    0.12888   2.721  0.00779 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.8228192)
    ## 
    ##     Null deviance: 440.631  on 99  degrees of freedom
    ## Residual deviance:  75.699  on 92  degrees of freedom
    ## AIC: 273.95
    ## 
    ## Number of Fisher Scoring iterations: 2
    ## 
    ## # Mediator Regression: 
    ## 
    ## Call:
    ## glm(formula = M1 ~ A + C1 + C2, family = binomial(), data = getCall(x$regsumm$mreg[[1L]])$data, 
    ##     weights = getCall(x$regsumm$mreg[[1L]])$weights)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.6894  -0.8718  -0.1670   0.7951   2.3486  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  1.08423    0.65565   1.654   0.0982 .  
    ## Atreat       0.05031    0.50534   0.100   0.9207    
    ## C1          -1.58663    0.35549  -4.463 8.07e-06 ***
    ## C2C2_1       0.76626    0.53989   1.419   0.1558    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 138.589  on 99  degrees of freedom
    ## Residual deviance:  99.309  on 96  degrees of freedom
    ## AIC: 107.31
    ## 
    ## Number of Fisher Scoring iterations: 5
    ## 
    ## 
    ## Call:
    ## glm(formula = M2 ~ A + M1 + C1 + C2, family = gaussian(), data = getCall(x$regsumm$mreg[[2L]])$data, 
    ##     weights = getCall(x$regsumm$mreg[[2L]])$weights)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -2.20635  -0.67665   0.01575   0.55827   2.50908  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   2.1163     0.2747   7.705 1.25e-11 ***
    ## Atreat        0.8060     0.1960   4.113 8.29e-05 ***
    ## M1           -1.2360     0.2287  -5.405 4.80e-07 ***
    ## C1            0.2852     0.1064   2.681  0.00866 ** 
    ## C2C2_1        2.2721     0.2049  11.088  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.8951861)
    ## 
    ##     Null deviance: 237.532  on 99  degrees of freedom
    ## Residual deviance:  85.043  on 95  degrees of freedom
    ## AIC: 279.59
    ## 
    ## Number of Fisher Scoring iterations: 2
    ## 
    ## # Causal Mediation Analysis via the Regression-based Approach
    ##  
    ## Direct counterfactual imputation estimation with 
    ##  bootstrap standard errors, percentile confidence intervals and p-values 
    ##  
    ##              Estimate Std.error  95% CIL 95% CIU  P.val    
    ## cde          -0.86779   0.65365 -1.87101   0.259    0.2    
    ## pnde          0.74058   0.22560  0.37923   1.111 <2e-16 ***
    ## tnde          1.00528   0.26771  0.70567   1.575 <2e-16 ***
    ## pnie          0.44970   0.14327  0.15184   0.630 <2e-16 ***
    ## tnie          0.71440   0.16437  0.41717   0.949 <2e-16 ***
    ## te            1.45498   0.24416  1.14724   1.929 <2e-16 ***
    ## intref        1.60837   0.64125  0.60196   2.720 <2e-16 ***
    ## intmed        0.26471   0.13269  0.12715   0.566 <2e-16 ***
    ## cde(prop)    -0.59643   0.44885 -1.39017   0.149    0.2    
    ## intref(prop)  1.10542   0.40454  0.39795   1.749 <2e-16 ***
    ## intmed(prop)  0.18193   0.07459  0.07879   0.334 <2e-16 ***
    ## pnie(prop)    0.30908   0.10013  0.10438   0.439 <2e-16 ***
    ## pm            0.49101   0.11292  0.29700   0.723 <2e-16 ***
    ## int           1.28735   0.46559  0.48987   2.030 <2e-16 ***
    ## pe            1.59643   0.44885  0.85062   2.390 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Reference values: 
    ## $a
    ## [1] "treat"
    ## 
    ## $astar
    ## [1] "control"
    ## 
    ## $mval
    ## $mval[[1]]
    ## [1] 0
    ## 
    ## $mval[[2]]
    ## [1] 0

    ggcmest(est) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.8))

![](man/figures/plot_cmest-1.png)

Lastly, let’s conduct sensitivity analysis for the results. Sensitivity
analysis for unmeasured confounding:

    cmsens(object = est, sens = "uc")

    ## Sensitivity Analysis For Unmeasured Confounding 
    ## 
    ## Evalues on the ratio scale: 
    ##          estRR   lowerRR  upperRR Evalue.estRR Evalue.lowerRR Evalue.upperRR
    ## cde  0.6877599 0.3962087 1.193850     2.266466             NA              1
    ## pnde 1.3763585 1.1378010 1.664933     2.096084       1.533768             NA
    ## tnde 1.5428319 1.2308983 1.933816     2.457981       1.764014             NA
    ## pnie 1.2140646 1.0758340 1.370056     1.723857       1.361465             NA
    ## tnie 1.3609082 1.1846767 1.563356     2.061738       1.652418             NA
    ## te   1.8730975 1.5243880 2.301576     3.151924       2.418463             NA

Assume that the continuous pre-exposure confounder was measured with
error. Sensitivity analysis for measurement error using regression
calibration with a set of assumed standard deviations of the measurement
error 0.1, 0.2 and 0.3:

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
    ## cde          -0.87002   0.61819 -1.89922   0.174    0.2    
    ## pnde          0.74054   0.22713  0.42010   1.193 <2e-16 ***
    ## tnde          0.98168   0.23526  0.71285   1.487 <2e-16 ***
    ## pnie          0.47960   0.19472  0.15297   0.818 <2e-16 ***
    ## tnie          0.72074   0.24985  0.39572   1.142 <2e-16 ***
    ## te            1.46128   0.22463  1.15663   1.931 <2e-16 ***
    ## intref        1.61056   0.49791  0.98416   2.566 <2e-16 ***
    ## intmed        0.24114   0.13630  0.12018   0.556 <2e-16 ***
    ## cde(prop)    -0.59538   0.42222 -1.32777   0.097    0.2    
    ## intref(prop)  1.10216   0.36833  0.55857   1.808 <2e-16 ***
    ## intmed(prop)  0.16502   0.07780  0.07732   0.339 <2e-16 ***
    ## pnie(prop)    0.32820   0.11522  0.11205   0.517 <2e-16 ***
    ## pm            0.49322   0.13704  0.27182   0.729 <2e-16 ***
    ## int           1.26718   0.41533  0.67220   2.097 <2e-16 ***
    ## pe            1.59538   0.42222  0.90270   2.328 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------
    ## 
    ## Measurement error 2: 
    ## 0.2
    ## Measurement error correction for measurement error 2: 
    ##              Estimate Std.error  95% CIL 95% CIU  P.val    
    ## cde          -0.87707   0.59069 -1.77290   0.183    0.2    
    ## pnde          0.79177   0.16210  0.58571   1.098 <2e-16 ***
    ## tnde          1.07635   0.23109  0.66568   1.467 <2e-16 ***
    ## pnie          0.43490   0.19209  0.12845   0.734 <2e-16 ***
    ## tnie          0.71948   0.19523  0.41077   0.983 <2e-16 ***
    ## te            1.51125   0.23278  1.10932   1.839 <2e-16 ***
    ## intref        1.66884   0.56056  0.77432   2.559 <2e-16 ***
    ## intmed        0.28458   0.19851 -0.04160   0.597    0.3    
    ## cde(prop)    -0.58036   0.38245 -1.17387   0.115    0.2    
    ## intref(prop)  1.10428   0.34723  0.50275   1.646 <2e-16 ***
    ## intmed(prop)  0.18831   0.12689 -0.03255   0.387    0.3    
    ## pnie(prop)    0.28778   0.11798  0.08485   0.488 <2e-16 ***
    ## pm            0.47608   0.09224  0.29155   0.614 <2e-16 ***
    ## int           1.29259   0.45522  0.54852   2.012 <2e-16 ***
    ## pe            1.58036   0.38245  0.88486   2.174 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------
    ## 
    ## Measurement error 3: 
    ## 0.3
    ## Measurement error correction for measurement error 3: 
    ##              Estimate Std.error  95% CIL 95% CIU  P.val    
    ## cde          -0.89027   0.71029 -2.06766   0.322    0.3    
    ## pnde          0.72419   0.21914  0.49316   1.177 <2e-16 ***
    ## tnde          1.06049   0.25086  0.77695   1.612 <2e-16 ***
    ## pnie          0.38582   0.18250  0.25844   0.860 <2e-16 ***
    ## tnie          0.72212   0.25880  0.33365   1.259 <2e-16 ***
    ## te            1.44631   0.25485  1.13055   1.978 <2e-16 ***
    ## intref        1.61446   0.67105  0.57277   2.870 <2e-16 ***
    ## intmed        0.33630   0.15915  0.05279   0.615 <2e-16 ***
    ## cde(prop)    -0.61554   0.44158 -1.22522   0.234    0.3    
    ## intref(prop)  1.11626   0.39383  0.43702   1.606 <2e-16 ***
    ## intmed(prop)  0.23252   0.08258  0.05130   0.335 <2e-16 ***
    ## pnie(prop)    0.26676   0.09689  0.15121   0.461 <2e-16 ***
    ## pm            0.49929   0.12455  0.24150   0.717 <2e-16 ***
    ## int           1.34878   0.45662  0.51018   1.920 <2e-16 ***
    ## pe            1.61554   0.44158  0.76626   2.225 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------

    ggcmsens(me1) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.8))

![](man/figures/plot_cmsens_me_con-1.png)

Then, assume that the exposure was measured with error. Sensitivity
analysis for measurement error using SIMEX with two assumed
misclassification matrices:

    me2 <- cmsens(object = est, sens = "me", MEmethod = "simex", MEvariable = "A", 
                  MEvartype = "cat", B = 20,
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
    ##              Estimate Std.error  95% CIL 95% CIU  P.val    
    ## cde          -0.69884   0.93427 -2.76621   0.499    0.2    
    ## pnde          0.92925   0.23365  0.42667   1.173 <2e-16 ***
    ## tnde          1.26685   0.25047  0.73463   1.670 <2e-16 ***
    ## pnie          0.34132   0.13041  0.17664   0.604 <2e-16 ***
    ## tnie          0.67891   0.28064  0.31394   1.255 <2e-16 ***
    ## te            1.60817   0.29137  0.99481   2.075 <2e-16 ***
    ## intref        1.62810   0.78348  0.66422   3.422 <2e-16 ***
    ## intmed        0.33759   0.24728  0.05100   0.925 <2e-16 ***
    ## cde(prop)    -0.43456   0.63850 -1.81496   0.333    0.2    
    ## intref(prop)  1.01239   0.54729  0.42472   2.245 <2e-16 ***
    ## intmed(prop)  0.20992   0.13635  0.03803   0.512 <2e-16 ***
    ## pnie(prop)    0.21224   0.07019  0.12534   0.357 <2e-16 ***
    ## pm            0.42217   0.13673  0.24027   0.690 <2e-16 ***
    ## int           1.22232   0.63518  0.46985   2.595 <2e-16 ***
    ## pe            1.43456   0.63850  0.66741   2.815 <2e-16 ***
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
    ##               Estimate Std.error   95% CIL 95% CIU  P.val    
    ## cde          -0.736215  0.954986 -2.840858   0.194    0.2    
    ## pnde          1.050472  0.320708  0.422315   1.521 <2e-16 ***
    ## tnde          1.488021  0.299204  0.959927   1.990 <2e-16 ***
    ## pnie          0.294209  0.238333  0.003238   0.825    0.1 .  
    ## tnie          0.731758  0.304975  0.572168   1.497 <2e-16 ***
    ## te            1.782230  0.293253  1.400750   2.284 <2e-16 ***
    ## intref        1.786687  0.871184  0.851292   3.928 <2e-16 ***
    ## intmed        0.437549  0.224634  0.209823   0.966 <2e-16 ***
    ## cde(prop)    -0.413086  0.523835 -1.718662   0.093    0.2    
    ## intref(prop)  1.002501  0.486884  0.422431   2.178 <2e-16 ***
    ## intmed(prop)  0.245506  0.112432  0.104197   0.465 <2e-16 ***
    ## pnie(prop)    0.165079  0.115842  0.002552   0.381    0.1 .  
    ## pm            0.410586  0.139007  0.278730   0.776 <2e-16 ***
    ## int           1.248007  0.574570  0.526628   2.584 <2e-16 ***
    ## pe            1.413086  0.523835  0.906988   2.719 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------

    ggcmsens(me2) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.8))

![](man/figures/plot_cmsens_me_cat-1.png)

References
----------

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
introducing the E-Value (2017). Annals of Internal Medicine. 167(4): 
268 - 274.

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
