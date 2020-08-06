<!-- README.md is generated from README.Rmd. Please edit that file -->

CMAverse <img src="man/figures/logo.png" align="right" width="250" />
=====================================================================

[![Project Status:
Active](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

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
    #>               Estimate Std.error   95% CIL 95% CIU P.val
    #> cde           -0.74954  31.26003  -2.69855  63.512   0.8
    #> pnde           0.69782   5.69243 -10.77027   1.118   0.8
    #> tnde           0.81557  21.42268 -42.45167   1.494   0.8
    #> pnie           0.42889  16.10370  -0.14838  32.948   0.4
    #> tnie           0.54664   0.55096  -0.10711   1.267   0.4
    #> te             1.24446   5.32598  -9.50408   1.673   0.4
    #> pm             0.43926   0.48642  -0.17628   0.951   0.8
    #> intref         1.44736  36.94434 -74.28200   3.584   0.8
    #> intmed         0.11775  15.73744 -31.68140   0.621   0.8
    #> cde(prop)     -0.60230   4.16165  -6.09392   4.399   0.8
    #> intref(prop)   1.16304   4.55613  -4.33725   7.158   0.8
    #> intmed(prop)   0.09462   1.56398  -0.68669   3.008   0.4
    #> pnie(prop)     0.34464   1.90408  -3.09782   1.616   0.8
    #> pm(overall)    0.43926   0.48642  -0.17628   0.951   0.8
    #> int(overall)   1.25766   6.00771  -5.01528  10.163   0.8
    #> pe(overall)    1.60230   4.16165  -3.39922   7.094   0.8

    plot(est) +
      theme(axis.text.x = element_text(angle = 30, vjust = 0.8))

![](man/figuresplot_cmest-1.png)

Lastly, let’s conduct sensitivity analysis for the results. Sensitivity
analysis for unmeasured confounding:

    cmsens(object = est, sens = "uc")
    #> Sensitivity Analysis For Unmeasured Confounding 
    #> 
    #> Evalues on the ratio scale: 
    #>          estRR      lowerRR      upperRR Evalue.estRR Evalue.lowerRR
    #> cde  0.7748469 7.106482e-10 8.448453e+08     1.902960             NA
    #> pnde 1.2680597 2.866857e-02 5.608844e+01     1.851083              1
    #> tnde 1.3199069 8.453793e-07 2.060796e+06     1.969712              1
    #> pnie 1.1571549 2.556614e-05 5.237425e+04     1.583597              1
    #> tnie 1.2044675 8.346551e-01 1.738134e+00     1.700728              1
    #> te   1.5273368 4.407035e-02 5.293258e+01     2.424789              1
    #>      Evalue.upperRR
    #> cde               1
    #> pnde             NA
    #> tnde             NA
    #> pnie             NA
    #> tnie             NA
    #> te               NA

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
    #>                Estimate  Std.error    95% CIL 95% CIU  P.val    
    #> cde           -0.750607   2.717236  -4.874211   0.935    0.8    
    #> pnde           0.600007   0.407310  -0.183471   0.712    0.4    
    #> tnde           0.772065   0.552929  -0.145503   1.035    0.8    
    #> pnie           0.491369   0.446456   0.148428   1.191 <2e-16 ***
    #> tnie           0.663426   0.637768   0.192393   1.568 <2e-16 ***
    #> te             1.263434   0.961663   0.008923   2.143    0.4    
    #> pm             0.525098   1.943008  -3.284767   0.761    0.4    
    #> intref         1.350614   2.635971  -0.627730   5.356    0.8    
    #> intmed         0.172058   0.220462  -0.116376   0.378    0.8    
    #> cde(prop)     -0.594100  37.182657  -2.120807  74.868    0.8    
    #> intref(prop)   1.069003  35.245261 -70.619786   2.389    0.8    
    #> intmed(prop)   0.136183   0.581681  -1.124890   0.176    0.8    
    #> pnie(prop)     0.388915   1.402690  -2.185063   0.973    0.4    
    #> pm(overall)    0.525098   1.943008  -3.284767   0.761    0.4    
    #> int(overall)   1.205185  35.806482 -71.744677   2.565    0.8    
    #> pe(overall)    1.594100  37.182657 -73.867781   3.121    0.8    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> ----------------------------------------------------------------
    #> 
    #> Measurement error 2: 
    #> 0.2
    #> Measurement error correction for measurement error 2: 
    #>              Estimate Std.error  95% CIL 95% CIU  P.val    
    #> cde          -0.75412   5.22537 -7.76519   4.563    0.8    
    #> pnde          0.62702   0.30391  0.60213   1.339 <2e-16 ***
    #> tnde          0.78123   0.32689  0.62308   1.453 <2e-16 ***
    #> pnie          0.47492   0.26602 -0.10643   0.551    0.4    
    #> tnie          0.62913   0.38329 -0.05255   0.901    0.4    
    #> te            1.25615   0.42165  0.81349   1.711 <2e-16 ***
    #> pm            0.50084   0.24072 -0.05760   0.550    0.4    
    #> intref        1.38115   5.02381 -3.23321   8.399    0.8    
    #> intmed        0.15421   0.14675 -0.00979   0.351    0.4    
    #> cde(prop)    -0.60035   5.35333 -9.63146   2.672    0.8    
    #> intref(prop)  1.09951   5.44289 -1.89004  10.405    0.8    
    #> intmed(prop)  0.12276   0.08979 -0.01345   0.213    0.4    
    #> pnie(prop)    0.37807   0.19949 -0.12326   0.337    0.4    
    #> pm(overall)   0.50084   0.24072 -0.05760   0.550    0.4    
    #> int(overall)  1.22227   5.38419 -1.81768  10.391    0.8    
    #> pe(overall)   1.60034   5.35333 -1.67211  10.631    0.8    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> ----------------------------------------------------------------
    #> 
    #> Measurement error 3: 
    #> 0.3
    #> Measurement error correction for measurement error 3: 
    #>              Estimate Std.error  95% CIL 95% CIU  P.val    
    #> cde          -0.76129   1.33306 -3.20424  -0.430 <2e-16 ***
    #> pnde          0.60669   0.37375  0.23471   1.179 <2e-16 ***
    #> tnde          0.74254   0.28941  0.80468   1.483 <2e-16 ***
    #> pnie          0.46111   0.35651 -0.04244   0.760    0.4    
    #> tnie          0.59697   0.56361  0.07590   1.347    0.4    
    #> te            1.20365   0.33816  1.22784   2.086 <2e-16 ***
    #> pm            0.49596   0.33641  0.04104   0.849    0.4    
    #> intref        1.36798   1.14929  1.06884   3.719 <2e-16 ***
    #> intmed        0.13585   0.26965  0.08448   0.692 <2e-16 ***
    #> cde(prop)    -0.63248   0.86535 -2.05258  -0.252 <2e-16 ***
    #> intref(prop)  1.13652   0.83032  0.62796   2.382 <2e-16 ***
    #> intmed(prop)  0.11287   0.14903  0.06456   0.389 <2e-16 ***
    #> pnie(prop)    0.38309   0.22556 -0.04363   0.484    0.4    
    #> pm(overall)   0.49596   0.33641  0.04104   0.849    0.4    
    #> int(overall)  1.24939   0.87982  0.80807   2.771 <2e-16 ***
    #> pe(overall)   1.63248   0.86535  1.25234   3.053 <2e-16 ***
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
    #>               Estimate Std.error   95% CIL 95% CIU  P.val    
    #> cde            1.66250   5.38504  -9.96214   3.529    0.4    
    #> pnde           0.66439   1.03384  -0.82488   1.551    0.4    
    #> tnde           0.71316   0.88486   0.26286   2.479 <2e-16 ***
    #> pnie           0.68170   0.45187  -0.59191   0.451    0.4    
    #> tnie           0.73047   0.24284   0.47805   1.016 <2e-16 ***
    #> te             1.39486   1.16927  -0.28677   2.522    0.4    
    #> pm             0.52369   0.66993  -0.94646   0.489    0.4    
    #> intref        -0.99811   4.68259  -2.14747   9.219    0.4    
    #> intmed         0.04877   0.41312   0.26672   1.161 <2e-16 ***
    #> cde(prop)      1.19188  10.05963  -1.99766  19.802    0.8    
    #> intref(prop)  -0.71557   9.39357 -17.85534   2.593    0.8    
    #> intmed(prop)   0.03496   1.21130  -2.19237   0.372    0.4    
    #> pnie(prop)     0.48873   0.55790   0.01861   1.269 <2e-16 ***
    #> pm(overall)    0.52369   0.66993  -0.94646   0.489    0.4    
    #> int(overall)  -0.68060  10.59821 -20.04771   2.965    0.8    
    #> pe(overall)   -0.19188  10.05963 -18.80181   2.998    0.8    
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
    #> cde          -2.76063   3.55151 -6.38787   1.841    0.8    
    #> pnde          1.16821   0.84776 -0.04832   2.053    0.4    
    #> tnde          1.13573   1.60192 -0.19131   3.772    0.4    
    #> pnie          0.18603   0.51782 -0.56827   0.725    0.8    
    #> tnie          0.15355   0.46162  0.11120   1.181 <2e-16 ***
    #> te            1.32176   1.13656  0.51115   3.203 <2e-16 ***
    #> pm            0.11617   0.47594  0.14018   1.162 <2e-16 ***
    #> intref        3.92884   4.13804 -1.34167   8.364    0.8    
    #> intmed       -0.03248   0.82938 -0.15466   1.719    0.8    
    #> cde(prop)    -2.08860   2.93955 -4.57456   2.521    0.8    
    #> intref(prop)  2.97243   3.26657 -2.65981   5.411    0.8    
    #> intmed(prop) -0.02457   0.33753 -0.30968   0.525    0.8    
    #> pnie(prop)    0.14075   0.73407 -0.17485   1.449    0.8    
    #> pm(overall)   0.11617   0.47594  0.14018   1.162 <2e-16 ***
    #> int(overall)  2.94786   3.49111 -2.96949   5.547    0.8    
    #> pe(overall)   3.08860   2.93955 -1.52053   5.575    0.8    
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
