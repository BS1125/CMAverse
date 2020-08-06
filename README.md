<!-- README.md is generated from README.Rmd. Please edit that file -->

CMAverse <img src="man/figures/logo.png" align="right" width="120" />
=====================================================================

About the Package
-----------------

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
<hr />
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

The DAG can be plotted using the `cmdag` function.

    cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"), 
          prec = c("C1", "C2"), postc = NULL,
          node = FALSE, text_col = "black")

![](man/figuresplot_dag-1.png)

Then, we estimate the causal effects using the `cmest` function. We use
the regression-based approach for illustration. The reference values for
the exposure are set to be 0 and 1. The reference values for the two
mediators are set to be 0.

    est <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                    mediator = c("M1", "M2"), prec = c("C1", "C2"), EMint = TRUE,
                    mreg = list("logistic", "linear"), yreg = "linear",
                    astar = 0, a = 1, mval = list(0, 0),
                    estimation = "imputation", inference = "bootstrap", nboot = 5)

Summarizing and plotting the results:

    summary(est)
    #> 
    #> Causal Mediation Analysis Via the Regression-based Approach
    #>  
    #> Direct counterfactual imputation estimation with 
    #>  bootstrap standard errors, percentile confidence intervals and p-values 
    #>  
    #>              Estimate Std.error  95% CIL 95% CIU  P.val    
    #> cde          -0.79603   1.23001 -2.58259   0.026    0.4    
    #> pnde          1.28025   0.48883 -0.01206   1.153    0.4    
    #> tnde          1.65581   0.49865  0.36531   1.480 <2e-16 ***
    #> pnie          0.17425   0.17082 -0.03401   0.376    0.4    
    #> tnie          0.54981   0.18089  0.25787   0.689 <2e-16 ***
    #> te            1.83006   0.54988  0.53216   1.827 <2e-16 ***
    #> pm            0.30044   0.35679  0.28307   1.063 <2e-16 ***
    #> intref        2.07627   1.00779  1.02890   3.341 <2e-16 ***
    #> intmed        0.37556   0.18627  0.13857   0.589 <2e-16 ***
    #> cde(prop)    -0.43497   1.82064 -4.17835   0.013    0.4    
    #> intref(prop)  1.13454   1.50025  0.62200   4.116 <2e-16 ***
    #> intmed(prop)  0.20522   0.27173  0.10317   0.745 <2e-16 ***
    #> pnie(prop)    0.09522   0.14888 -0.02149   0.341    0.4    
    #> pm(overall)   0.30044   0.35679  0.28307   1.063 <2e-16 ***
    #> int(overall)  1.33976   1.76771  0.78085   4.861 <2e-16 ***
    #> pe(overall)   1.43497   1.82064  0.98670   5.178 <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    plot(est) +
      theme(axis.text.x = element_text(angle = 30, vjust = 0.8))

![](man/figuresplot_cmest-1.png)

Lastly, let’s conduct sensitivity analysis for the results. Sensitivity
analysis for unmeasured confounding:

    cmsens(object = est, sens = "uc")
    #> Sensitivity Analysis For Unmeasured Confounding 
    #> 
    #> Evalues on the ratio scale: 
    #>          estRR   lowerRR  upperRR Evalue.estRR Evalue.lowerRR Evalue.upperRR
    #> cde  0.7363063 0.2919112 1.857233     2.055545             NA              1
    #> pnde 1.6360966 1.1327186 2.363175     2.656251       1.520446             NA
    #> tnde 1.8902995 1.2990795 2.750588     3.187578       1.922400             NA
    #> pnie 1.0693029 0.9403700 1.215914     1.341527       1.000000             NA
    #> tnie 1.2354421 1.0782746 1.415518     1.774771       1.368794             NA
    #> te   2.0213027 1.3365947 3.056772     3.458091       2.007334             NA

Assume that the continuous pre-exposure confounder was measured with
error. Sensitivity analysis using regression calibration with a set of
assumed standard deviations of the measurement error 0.1, 0.2 and 0.3:

    me1 <- cmsens(object = est, sens = "me", MEmethod = "rc", 
                  MEvariable = "C1", MEvartype = "con", MEerror = c(0.1, 0.2, 0.3))

Summarizing and plotting the results:

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
    #> cde          -0.79820   2.05942 -4.75515  -0.273 <2e-16 ***
    #> pnde          1.23070   0.97048 -0.37641   2.061    0.4    
    #> tnde          1.59187   0.54750  1.13974   2.436 <2e-16 ***
    #> pnie          0.11461   0.36516 -0.69258   0.204    0.8    
    #> tnie          0.47578   0.45751 -0.03185   1.085    0.4    
    #> te            1.70648   0.58479  0.70442   2.198 <2e-16 ***
    #> pm            0.27881   0.75418 -0.02771   1.668    0.4    
    #> intref        2.02890   1.25282  1.39297   4.406 <2e-16 ***
    #> intmed        0.36117   0.72871  0.12311   1.747 <2e-16 ***
    #> cde(prop)    -0.46775   3.45270 -7.38684  -0.133 <2e-16 ***
    #> intref(prop)  1.18894   2.71558  1.00033   6.719 <2e-16 ***
    #> intmed(prop)  0.21164   1.25728  0.09077   2.716 <2e-16 ***
    #> pnie(prop)    0.06716   0.53134 -1.06312   0.153    0.8    
    #> pm(overall)   0.27881   0.75418 -0.02771   1.668    0.4    
    #> int(overall)  1.40058   3.96922  1.16239   9.435 <2e-16 ***
    #> pe(overall)   1.46775   3.45270  1.13316   8.387 <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> ----------------------------------------------------------------
    #> 
    #> Measurement error 2: 
    #> 0.2
    #> Measurement error correction for measurement error 2: 
    #>              Estimate Std.error  95% CIL 95% CIU  P.val    
    #> cde          -0.80530   0.37662 -1.02077  -0.165 <2e-16 ***
    #> pnde          1.23266   0.27164  1.13750   1.805 <2e-16 ***
    #> tnde          1.59550   0.33676  1.37133   2.094 <2e-16 ***
    #> pnie          0.12076   0.29531 -0.12280   0.575    0.4    
    #> tnie          0.48361   0.13086  0.34616   0.677 <2e-16 ***
    #> te            1.71627   0.38668  1.49001   2.475 <2e-16 ***
    #> pm            0.28178   0.03110  0.22808   0.305 <2e-16 ***
    #> intref        2.03796   0.51477  1.33203   2.518 <2e-16 ***
    #> intmed        0.36284   0.23032  0.10179   0.631 <2e-16 ***
    #> cde(prop)    -0.46921   0.18842 -0.52098  -0.108 <2e-16 ***
    #> intref(prop)  1.18744   0.16689  0.88014   1.252 <2e-16 ***
    #> intmed(prop)  0.21141   0.11874  0.04194   0.322 <2e-16 ***
    #> pnie(prop)    0.07036   0.13411 -0.06096   0.249    0.4    
    #> pm(overall)   0.28178   0.03110  0.22808   0.305 <2e-16 ***
    #> int(overall)  1.39885   0.24205  1.03935   1.562 <2e-16 ***
    #> pe(overall)   1.46921   0.18842  1.10822   1.521 <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> ----------------------------------------------------------------
    #> 
    #> Measurement error 3: 
    #> 0.3
    #> Measurement error correction for measurement error 3: 
    #>              Estimate Std.error  95% CIL 95% CIU  P.val    
    #> cde          -0.81947   2.13735 -2.01148   3.312    0.8    
    #> pnde          1.31452   0.75082  0.19375   2.105 <2e-16 ***
    #> tnde          1.68671   0.40920  0.92934   1.798 <2e-16 ***
    #> pnie          0.16463   0.19096 -0.18778   0.242    0.8    
    #> tnie          0.53682   0.50781 -0.15871   0.950    0.8    
    #> te            1.85133   0.45301  0.96949   2.001 <2e-16 ***
    #> pm            0.28996   0.40127 -0.07873   0.853    0.8    
    #> intref        2.13399   1.82454 -1.31816   3.003    0.8    
    #> intmed        0.37219   0.48587 -0.35246   0.789    0.4    
    #> cde(prop)    -0.44264   1.48394 -2.07658   1.663    0.8    
    #> intref(prop)  1.15267   1.59956 -0.67414   3.065    0.8    
    #> intmed(prop)  0.20104   0.35827 -0.17234   0.712    0.4    
    #> pnie(prop)    0.08892   0.12400 -0.12629   0.141    0.8    
    #> pm(overall)   0.28996   0.40127 -0.07873   0.853    0.8    
    #> int(overall)  1.35371   1.55774 -0.77412   3.124    0.4    
    #> pe(overall)   1.44264   1.48394 -0.66263   3.077    0.4    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> ----------------------------------------------------------------

    plot(me1) +
      theme(axis.text.x = element_text(angle = 30, vjust = 0.8))

![](man/figuresplot_cmsens_me_con-1.png)

Then, assume that the exposure was measured with error. Sensitivity
analysis using MCSIMEX with two assumed misclassification matrices:

    me2 <- cmsens(object = est, sens = "me", MEmethod = "simex", MEvariable = "A", 
                  MEvartype = "cat", B = 5,
                  MEerror = list(matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2), 
                                 matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)))

Summarizing and plotting the results:

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
    #> cde           0.93447   2.26096 -3.06639   2.565    0.8    
    #> pnde          1.35819   0.93235 -0.05559   2.160    0.4    
    #> tnde          1.57364   0.78995  0.82240   2.604 <2e-16 ***
    #> pnie          0.14561   0.34534 -0.45514   0.418    0.8    
    #> tnie          0.36106   0.46025 -0.14848   1.012    0.4    
    #> te            1.71925   0.65997  0.93636   2.480 <2e-16 ***
    #> pm            0.21001   0.49298 -0.08477   1.074    0.4    
    #> intref        0.42372   1.78265 -1.09703   3.113    0.4    
    #> intmed        0.21545   0.39077 -0.04326   0.912    0.4    
    #> cde(prop)     0.54353   1.95489 -3.28545   1.368    0.8    
    #> intref(prop)  0.24646   1.56863 -0.57259   3.243    0.4    
    #> intmed(prop)  0.12532   0.42551 -0.01797   0.956    0.4    
    #> pnie(prop)    0.08469   0.20174 -0.26618   0.228    0.8    
    #> pm(overall)   0.21001   0.49298 -0.08477   1.074    0.4    
    #> int(overall)  0.37177   1.97333 -0.59056   4.199    0.4    
    #> pe(overall)   0.45647   1.95489 -0.36822   4.285    0.4    
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
    #> cde          -2.94931   1.76611 -4.51957  -0.287 <2e-16 ***
    #> pnde          1.05633   0.53605  0.65996   1.973 <2e-16 ***
    #> tnde          2.01633   0.33569  1.61175   2.456 <2e-16 ***
    #> pnie         -0.20615   0.23444 -0.60042  -0.014    0.4    
    #> tnie          0.75385   0.59776 -0.38049   0.919    0.8    
    #> te            1.81018   0.32588  1.21715   1.961 <2e-16 ***
    #> pm            0.41645   0.38401 -0.31700   0.579    0.8    
    #> intref        4.00564   1.41588  1.78985   5.218 <2e-16 ***
    #> intmed        0.95999   0.46489  0.04349   1.192    0.4    
    #> cde(prop)    -1.62929   1.15128 -2.95347  -0.178 <2e-16 ***
    #> intref(prop)  2.21284   0.94549  1.02183   3.403 <2e-16 ***
    #> intmed(prop)  0.53033   0.30004  0.02150   0.775    0.4    
    #> pnie(prop)   -0.11388   0.14756 -0.33850  -0.008    0.4    
    #> pm(overall)   0.41645   0.38401 -0.31700   0.579    0.8    
    #> int(overall)  2.74318   1.16479  1.30742   4.157 <2e-16 ***
    #> pe(overall)   2.62929   1.15128  1.17806   3.953 <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> ----------------------------------------------------------------

    plot(me2) +
      theme(axis.text.x = element_text(angle = 30, vjust = 0.8))

![](man/figuresplot_cmsens_me_cat-1.png)

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
