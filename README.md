<!-- README.md is generated from README.Rmd. Please edit that file -->

CMAverse <img src="man/figures/logo.png" align="right" width="240" />
=====================================================================

[![Project Status:
Active](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

[![Travis build
status](https://travis-ci.com/BS1125/CMAverse.svg?branch=master)](https://travis-ci.com/BS1125/CMAverse)

About the Package
-----------------

The R package `CMAverse` provides a suite of functions for reproducible
causal mediation analysis including `cmdag` for DAG visualization,
`cmest` for statistical modeling and `cmsens` for sensitivity analysis.

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

    devtools::install_github("LindaValeri/CMAverse")

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
    ##     Min       1Q   Median       3Q      Max  
    ## -2.2973  -0.7771   0.0005   0.6298   3.5332  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.56627    0.46919   1.207   0.2306    
    ## Atreat       0.36856    0.77315   0.477   0.6347    
    ## M1           0.40059    0.37959   1.055   0.2940    
    ## M2           0.64642    0.12499   5.172 1.35e-06 ***
    ## C1          -0.32253    0.15443  -2.089   0.0395 *  
    ## C2C2_1       1.90660    0.31312   6.089 2.60e-08 ***
    ## Atreat:M1   -0.09295    0.50941  -0.182   0.8556    
    ## Atreat:M2    0.27940    0.16091   1.736   0.0858 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 1.216538)
    ## 
    ##     Null deviance: 534.43  on 99  degrees of freedom
    ## Residual deviance: 111.92  on 92  degrees of freedom
    ## AIC: 313.05
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
    ## -2.3778  -0.6126   0.3279   0.6663   1.8612  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   1.2423     0.5550   2.238   0.0252 *  
    ## Atreat        0.2813     0.5479   0.513   0.6077    
    ## C1           -1.7660     0.3732  -4.732 2.22e-06 ***
    ## C2C2_1        1.3757     0.5666   2.428   0.0152 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 131.791  on 99  degrees of freedom
    ## Residual deviance:  88.777  on 96  degrees of freedom
    ## AIC: 96.777
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
    ## -2.86482  -0.64760  -0.00601   0.67966   2.26602  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   2.0860     0.2944   7.085 2.41e-10 ***
    ## Atreat        0.8358     0.2219   3.766 0.000287 ***
    ## M1           -1.0277     0.2750  -3.737 0.000318 ***
    ## C1            0.4389     0.1407   3.119 0.002402 ** 
    ## C2C2_1        1.8047     0.2354   7.668 1.50e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 1.114778)
    ## 
    ##     Null deviance: 243.28  on 99  degrees of freedom
    ## Residual deviance: 105.90  on 95  degrees of freedom
    ## AIC: 301.52
    ## 
    ## Number of Fisher Scoring iterations: 2
    ## 
    ## # Causal Mediation Analysis via the Regression-based Approach
    ##  
    ## Direct counterfactual imputation estimation with 
    ##  bootstrap standard errors, percentile confidence intervals and p-values 
    ##  
    ##              Estimate Std.error  95% CIL 95% CIU  P.val    
    ## cde           0.36856   0.55259 -0.76140   1.150    0.7    
    ## pnde          1.15734   0.21379  0.78125   1.565 <2e-16 ***
    ## tnde          1.34145   0.18381  1.08131   1.733 <2e-16 ***
    ## pnie          0.50599   0.15557  0.20548   0.757 <2e-16 ***
    ## tnie          0.69010   0.30381  0.16081   1.308 <2e-16 ***
    ## te            1.84743   0.25482  1.55444   2.395 <2e-16 ***
    ## intref        0.78878   0.43032  0.14606   1.656 <2e-16 ***
    ## intmed        0.18411   0.18851 -0.07720   0.574    0.1 .  
    ## cde(prop)     0.19950   0.30667 -0.42774   0.678    0.7    
    ## intref(prop)  0.42696   0.22522  0.09010   0.884 <2e-16 ***
    ## intmed(prop)  0.09966   0.09165 -0.04724   0.287    0.1 .  
    ## pnie(prop)    0.27389   0.06687  0.12854   0.362 <2e-16 ***
    ## pm            0.37354   0.13233  0.10275   0.603 <2e-16 ***
    ## int           0.52661   0.29662  0.10572   1.111 <2e-16 ***
    ## pe            0.80050   0.30667  0.32219   1.428 <2e-16 ***
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
    ## cde  1.155292 0.756546 1.764202     1.578857       1.000000             NA
    ## pnde 1.573480 1.335766 1.853497     2.523406       2.005471             NA
    ## tnde 1.691131 1.468999 1.946854     2.772239       2.299034             NA
    ## pnie 1.219180 1.082199 1.373500     1.736113       1.380453             NA
    ## tnie 1.310340 1.038249 1.653737     1.948031       1.237526             NA
    ## te   2.061793 1.696139 2.506275     3.541387       2.782762             NA

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
    ## cde           0.36576   0.89672 -1.12581   1.832    0.6    
    ## pnde          1.16218   0.28097  0.68957   1.492 <2e-16 ***
    ## tnde          1.36234   0.21780  0.99732   1.678 <2e-16 ***
    ## pnie          0.51907   0.16686  0.29478   0.881 <2e-16 ***
    ## tnie          0.71923   0.23257  0.38485   1.218 <2e-16 ***
    ## te            1.88141   0.21993  1.59596   2.327 <2e-16 ***
    ## intref        0.79642   0.69557 -0.41297   1.902    0.4    
    ## intmed        0.20016   0.15216 -0.08423   0.428    0.1 .  
    ## cde(prop)     0.19441   0.47017 -0.65529   0.896    0.6    
    ## intref(prop)  0.42331   0.39134 -0.20372   1.106    0.4    
    ## intmed(prop)  0.10639   0.08233 -0.04093   0.239    0.1 .  
    ## pnie(prop)    0.27589   0.08062  0.15552   0.431 <2e-16 ***
    ## pm            0.38228   0.12313  0.21758   0.598 <2e-16 ***
    ## int           0.52970   0.46386 -0.24280   1.297    0.2    
    ## pe            0.80559   0.47017  0.10413   1.655    0.1 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------
    ## 
    ## Measurement error 2: 
    ## 0.2
    ## Measurement error correction for measurement error 2: 
    ##              Estimate Std.error  95% CIL 95% CIU  P.val    
    ## cde           0.35664   0.79873 -0.92101   2.166    0.4    
    ## pnde          1.13630   0.29115  0.73758   1.744 <2e-16 ***
    ## tnde          1.38360   0.25768  1.01982   1.857 <2e-16 ***
    ## pnie          0.55935   0.19116  0.29930   0.909 <2e-16 ***
    ## tnie          0.80665   0.32767  0.30956   1.453 <2e-16 ***
    ## te            1.94295   0.28051  1.45631   2.368 <2e-16 ***
    ## intref        0.77966   0.66137 -0.60875   1.659    0.2    
    ## intmed        0.24730   0.17665 -0.06557   0.565    0.3    
    ## cde(prop)     0.18356   0.43235 -0.49931   1.163    0.4    
    ## intref(prop)  0.40128   0.35814 -0.32986   0.997    0.2    
    ## intmed(prop)  0.12728   0.08373 -0.03589   0.241    0.3    
    ## pnie(prop)    0.28789   0.08397  0.16075   0.422 <2e-16 ***
    ## pm            0.41517   0.13482  0.16559   0.618 <2e-16 ***
    ## int           0.52856   0.42438 -0.36575   1.180    0.2    
    ## pe            0.81644   0.43235 -0.16278   1.499    0.2    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------
    ## 
    ## Measurement error 3: 
    ## 0.3
    ## Measurement error correction for measurement error 3: 
    ##              Estimate Std.error  95% CIL 95% CIU  P.val    
    ## cde           0.33860   1.03908 -1.41101   1.813    0.7    
    ## pnde          1.08980   0.28362  0.47505   1.411 <2e-16 ***
    ## tnde          1.34772   0.22211  0.87051   1.685 <2e-16 ***
    ## pnie          0.58171   0.16696  0.32956   0.915 <2e-16 ***
    ## tnie          0.83964   0.33934  0.39549   1.501 <2e-16 ***
    ## te            1.92944   0.23526  1.52038   2.300 <2e-16 ***
    ## intref        0.75120   0.84854 -0.47863   1.937    0.7    
    ## intmed        0.25792   0.23205 -0.11178   0.587    0.4    
    ## cde(prop)     0.17549   0.57983 -0.78330   0.981    0.7    
    ## intref(prop)  0.38934   0.46676 -0.26119   1.108    0.7    
    ## intmed(prop)  0.13368   0.11622 -0.06154   0.265    0.4    
    ## pnie(prop)    0.30149   0.08159  0.18644   0.468 <2e-16 ***
    ## pm            0.43517   0.15752  0.23556   0.724 <2e-16 ***
    ## int           0.52302   0.56678 -0.32010   1.316    0.7    
    ## pe            0.82451   0.57983  0.01906   1.783    0.1 .  
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
    ## cde           0.50876   0.92410 -1.48974   1.630    0.5    
    ## pnde          1.32927   0.27017  0.91237   1.778 <2e-16 ***
    ## tnde          1.64985   0.27399  1.23670   2.092 <2e-16 ***
    ## pnie          0.59994   0.19770  0.19636   0.814 <2e-16 ***
    ## tnie          0.92053   0.31338  0.19928   1.385 <2e-16 ***
    ## te            2.24979   0.27137  1.62592   2.516 <2e-16 ***
    ## intref        0.82051   0.83187 -0.24873   2.406    0.3    
    ## intmed        0.32058   0.21955 -0.02404   0.745    0.2    
    ## cde(prop)     0.22614   0.43861 -0.65199   0.837    0.5    
    ## intref(prop)  0.36470   0.37299 -0.12273   1.048    0.3    
    ## intmed(prop)  0.14250   0.09535 -0.01256   0.325    0.2    
    ## pnie(prop)    0.26667   0.08746  0.09704   0.362 <2e-16 ***
    ## pm            0.40916   0.12750  0.11281   0.604 <2e-16 ***
    ## int           0.50720   0.45705 -0.10468   1.373    0.2    
    ## pe            0.77386   0.43861  0.16331   1.652 <2e-16 ***
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
    ##              Estimate Std.error  95% CIL 95% CIU  P.val    
    ## cde          -0.05193   0.91528 -0.89695   2.457    0.4    
    ## pnde          1.44043   0.37912  0.97120   2.300 <2e-16 ***
    ## tnde          1.81998   0.29951  1.31589   2.339 <2e-16 ***
    ## pnie          0.39143   0.17756  0.27146   0.888 <2e-16 ***
    ## tnie          0.77098   0.33561  0.21134   1.443 <2e-16 ***
    ## te            2.21140   0.25755  2.00501   2.776 <2e-16 ***
    ## intref        1.49236   0.69827 -0.26180   2.096    0.2    
    ## intmed        0.37955   0.23651 -0.09377   0.715    0.2    
    ## cde(prop)    -0.02348   0.37905 -0.38989   0.974    0.4    
    ## intref(prop)  0.67485   0.29539 -0.10731   0.885    0.2    
    ## intmed(prop)  0.17163   0.09655 -0.04070   0.301    0.2    
    ## pnie(prop)    0.17700   0.07618  0.10524   0.378 <2e-16 ***
    ## pm            0.34864   0.13500  0.10051   0.580 <2e-16 ***
    ## int           0.84648   0.37499 -0.11402   1.186    0.2    
    ## pe            1.02348   0.37905  0.02615   1.390 <2e-16 ***
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
