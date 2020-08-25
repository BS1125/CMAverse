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

    n <- 1000
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
                    estimation = "imputation", inference = "bootstrap", nboot = 300)

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
    ## -3.11180  -0.70348  -0.05408   0.71066   2.92588  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.31968    0.13595   2.352 0.018891 *  
    ## Atreat       0.59792    0.21105   2.833 0.004703 ** 
    ## M1           0.54236    0.10530   5.151 3.13e-07 ***
    ## M2           0.59697    0.03847  15.518  < 2e-16 ***
    ## C1          -0.22520    0.04191  -5.373 9.66e-08 ***
    ## C2C2_1       2.11888    0.09159  23.134  < 2e-16 ***
    ## Atreat:M1    0.19515    0.14563   1.340 0.180528    
    ## Atreat:M2    0.15599    0.04410   3.537 0.000423 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 1.026529)
    ## 
    ##     Null deviance: 5142.0  on 999  degrees of freedom
    ## Residual deviance: 1018.3  on 992  degrees of freedom
    ## AIC: 2874
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
    ## -2.8091  -0.8498   0.1506   0.8511   2.6784  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   1.0168     0.1752   5.803  6.5e-09 ***
    ## Atreat        0.5640     0.1535   3.675 0.000238 ***
    ## C1           -1.4609     0.1065 -13.720  < 2e-16 ***
    ## C2C2_1        0.3799     0.1551   2.450 0.014293 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1386.2  on 999  degrees of freedom
    ## Residual deviance: 1036.5  on 996  degrees of freedom
    ## AIC: 1044.5
    ## 
    ## Number of Fisher Scoring iterations: 4
    ## 
    ## 
    ## Call:
    ## glm(formula = M2 ~ A + M1 + C1 + C2, family = gaussian(), data = getCall(x$regsumm$mreg[[2L]])$data, 
    ##     weights = getCall(x$regsumm$mreg[[2L]])$weights)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.2332  -0.6571  -0.0929   0.6352   4.0297  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.94294    0.08630   22.52   <2e-16 ***
    ## Atreat       0.71537    0.06631   10.79   <2e-16 ***
    ## M1          -0.96908    0.07521  -12.88   <2e-16 ***
    ## C1           0.52887    0.03780   13.99   <2e-16 ***
    ## C2C2_1       1.98339    0.06470   30.66   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.9997585)
    ## 
    ##     Null deviance: 2556.77  on 999  degrees of freedom
    ## Residual deviance:  994.76  on 995  degrees of freedom
    ## AIC: 2844.6
    ## 
    ## Number of Fisher Scoring iterations: 2
    ## 
    ## # Causal Mediation Analysis via the Regression-based Approach
    ##  
    ## Direct counterfactual imputation estimation with 
    ##  bootstrap standard errors, percentile confidence intervals and p-values 
    ##  
    ##              Estimate Std.error 95% CIL 95% CIU   P.val    
    ## cde           0.59792   0.24485 0.16374   1.106 0.00667 ** 
    ## pnde          1.18295   0.07499 1.04955   1.338 < 2e-16 ***
    ## tnde          1.30013   0.07619 1.16256   1.453 < 2e-16 ***
    ## pnie          0.42247   0.04903 0.33826   0.530 < 2e-16 ***
    ## tnie          0.53964   0.05610 0.43808   0.661 < 2e-16 ***
    ## te            1.72259   0.08284 1.57245   1.897 < 2e-16 ***
    ## intref        0.58503   0.21325 0.17925   1.002 0.00667 ** 
    ## intmed        0.11718   0.04407 0.03243   0.200 0.00667 ** 
    ## cde(prop)     0.34711   0.13902 0.09358   0.611 0.00667 ** 
    ## intref(prop)  0.33962   0.12526 0.10512   0.585 0.00667 ** 
    ## intmed(prop)  0.06802   0.02532 0.01883   0.118 0.00667 ** 
    ## pnie(prop)    0.24525   0.02567 0.20125   0.299 < 2e-16 ***
    ## pm            0.31327   0.02856 0.26117   0.371 < 2e-16 ***
    ## int           0.40764   0.14961 0.12697   0.695 0.00667 ** 
    ## pe            0.65289   0.13902 0.38887   0.906 < 2e-16 ***
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
    ##         estRR  lowerRR  upperRR Evalue.estRR Evalue.lowerRR Evalue.upperRR
    ## cde  1.271034 1.048883 1.540237     1.857971       1.275317             NA
    ## pnde 1.607193 1.515356 1.704595     2.595058       2.399069             NA
    ## tnde 1.684536 1.586790 1.788303     2.758373       2.551733             NA
    ## pnie 1.184657 1.139949 1.231120     1.652371       1.539366             NA
    ## tnie 1.241667 1.188204 1.297535     1.789453       1.661093             NA
    ## te   1.995598 1.870013 2.129618     3.405142       3.145526             NA

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
    ##              Estimate Std.error 95% CIL 95% CIU   P.val    
    ## cde           0.59615   0.23743 0.15934   1.072 0.00667 ** 
    ## pnde          1.17879   0.07968 1.02692   1.331 < 2e-16 ***
    ## tnde          1.29464   0.07577 1.15037   1.448 < 2e-16 ***
    ## pnie          0.42616   0.05025 0.33458   0.537 < 2e-16 ***
    ## tnie          0.54201   0.05840 0.43434   0.663 < 2e-16 ***
    ## te            1.72080   0.08077 1.56541   1.875 < 2e-16 ***
    ## intref        0.58264   0.19814 0.18732   0.951 < 2e-16 ***
    ## intmed        0.11585   0.04125 0.03730   0.194 < 2e-16 ***
    ## cde(prop)     0.34644   0.13329 0.09642   0.610 0.00667 ** 
    ## intref(prop)  0.33859   0.11736 0.11053   0.551 < 2e-16 ***
    ## intmed(prop)  0.06732   0.02389 0.02065   0.112 < 2e-16 ***
    ## pnie(prop)    0.24765   0.02645 0.19941   0.303 < 2e-16 ***
    ## pm            0.31497   0.03141 0.25857   0.378 < 2e-16 ***
    ## int           0.40591   0.14032 0.13119   0.660 < 2e-16 ***
    ## pe            0.65356   0.13329 0.39036   0.904 < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------
    ## 
    ## Measurement error 2: 
    ## 0.2
    ## Measurement error correction for measurement error 2: 
    ##              Estimate Std.error 95% CIL 95% CIU   P.val    
    ## cde           0.59045   0.23372 0.16016   1.053 0.00667 ** 
    ## pnde          1.18096   0.08027 1.01211   1.328 < 2e-16 ***
    ## tnde          1.29833   0.07961 1.13269   1.434 < 2e-16 ***
    ## pnie          0.43328   0.05429 0.33379   0.545 < 2e-16 ***
    ## tnie          0.55065   0.06674 0.43077   0.676 < 2e-16 ***
    ## te            1.73161   0.08576 1.55784   1.902 < 2e-16 ***
    ## intref        0.59052   0.19856 0.18005   0.943 < 2e-16 ***
    ## intmed        0.11737   0.04197 0.03654   0.197 < 2e-16 ***
    ## cde(prop)     0.34098   0.13394 0.09668   0.605 0.00667 ** 
    ## intref(prop)  0.34102   0.11655 0.11421   0.559 < 2e-16 ***
    ## intmed(prop)  0.06778   0.02401 0.02194   0.112 < 2e-16 ***
    ## pnie(prop)    0.25022   0.02858 0.20146   0.305 < 2e-16 ***
    ## pm            0.31800   0.03429 0.25709   0.388 < 2e-16 ***
    ## int           0.40880   0.13956 0.13579   0.668 < 2e-16 ***
    ## pe            0.65902   0.13394 0.39540   0.903 < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------
    ## 
    ## Measurement error 3: 
    ## 0.3
    ## Measurement error correction for measurement error 3: 
    ##              Estimate Std.error 95% CIL 95% CIU   P.val    
    ## cde           0.57940   0.23907 0.09092   0.985 0.02667 *  
    ## pnde          1.17241   0.08245 1.00005   1.317 < 2e-16 ***
    ## tnde          1.29099   0.07763 1.13812   1.435 < 2e-16 ***
    ## pnie          0.44777   0.05511 0.34188   0.548 < 2e-16 ***
    ## tnie          0.56635   0.06869 0.43270   0.700 < 2e-16 ***
    ## te            1.73876   0.08090 1.59069   1.894 < 2e-16 ***
    ## intref        0.59301   0.19890 0.22873   0.974 0.00667 ** 
    ## intmed        0.11858   0.04336 0.04509   0.220 0.00667 ** 
    ## cde(prop)     0.33323   0.13546 0.05492   0.562 0.02667 *  
    ## intref(prop)  0.34105   0.11620 0.13288   0.563 0.00667 ** 
    ## intmed(prop)  0.06820   0.02476 0.02563   0.123 0.00667 ** 
    ## pnie(prop)    0.25752   0.02924 0.20049   0.309 < 2e-16 ***
    ## pm            0.32572   0.03643 0.25840   0.400 < 2e-16 ***
    ## int           0.40925   0.13992 0.16280   0.674 0.00667 ** 
    ## pe            0.66677   0.13546 0.43767   0.945 < 2e-16 ***
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
                  MEvartype = "cat", B = 200,
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
    ##              Estimate Std.error 95% CIL 95% CIU   P.val    
    ## cde           0.76428   0.26985 0.20462   1.269 0.00667 ** 
    ## pnde          1.38883   0.09021 1.20685   1.550 < 2e-16 ***
    ## tnde          1.52720   0.09390 1.32879   1.697 < 2e-16 ***
    ## pnie          0.42948   0.05565 0.32518   0.540 < 2e-16 ***
    ## tnie          0.56784   0.06711 0.44379   0.701 < 2e-16 ***
    ## te            1.95668   0.09526 1.76508   2.127 < 2e-16 ***
    ## intref        0.62455   0.23180 0.20197   1.097 0.01333 *  
    ## intmed        0.13836   0.05647 0.04133   0.266 0.00667 ** 
    ## cde(prop)     0.39060   0.13600 0.10513   0.653 0.00667 ** 
    ## intref(prop)  0.31919   0.11888 0.10557   0.560 0.01333 *  
    ## intmed(prop)  0.07071   0.02840 0.02024   0.133 0.00667 ** 
    ## pnie(prop)    0.21949   0.02697 0.16870   0.269 < 2e-16 ***
    ## pm            0.29021   0.03100 0.23058   0.352 < 2e-16 ***
    ## int           0.38990   0.14631 0.12500   0.686 0.01333 *  
    ## pe            0.60940   0.13600 0.34741   0.895 < 2e-16 ***
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
    ## cde           0.92531   0.31322 0.38200   1.536 <2e-16 ***
    ## pnde          1.59300   0.10560 1.37862   1.788 <2e-16 ***
    ## tnde          1.76278   0.10093 1.57267   1.941 <2e-16 ***
    ## pnie          0.42056   0.05546 0.31616   0.534 <2e-16 ***
    ## tnie          0.59033   0.07415 0.44763   0.725 <2e-16 ***
    ## te            2.18333   0.10059 1.95738   2.370 <2e-16 ***
    ## intref        0.66769   0.25291 0.16121   1.106 <2e-16 ***
    ## intmed        0.16978   0.07145 0.03143   0.312 <2e-16 ***
    ## cde(prop)     0.42381   0.14055 0.18098   0.708 <2e-16 ***
    ## intref(prop)  0.30581   0.11821 0.07262   0.515 <2e-16 ***
    ## intmed(prop)  0.07776   0.03273 0.01461   0.143 <2e-16 ***
    ## pnie(prop)    0.19262   0.02451 0.14987   0.244 <2e-16 ***
    ## pm            0.27038   0.03273 0.20568   0.340 <2e-16 ***
    ## int           0.38357   0.14971 0.08752   0.659 <2e-16 ***
    ## pe            0.57619   0.14055 0.29173   0.819 <2e-16 ***
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
