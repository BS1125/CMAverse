<!-- README.md is generated from README.Rmd. Please edit that file -->

CMAverse <img src="man/figures/logo.png" align="right" width="200" />
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
    #>               Estimate Std.error   95% CIL 95% CIU  P.val    
    #> cde           0.157820  1.043697 -0.509488   1.915    0.8    
    #> pnde          0.884084  0.391870  0.281059   1.175 <2e-16 ***
    #> tnde          0.881445  0.511871  0.034728   1.196    0.4    
    #> pnie          0.191833  0.270172 -0.055581   0.565    0.4    
    #> tnie          0.189194  0.285642 -0.199211   0.507    0.4    
    #> te            1.073278  0.332616  0.574208   1.377 <2e-16 ***
    #> pm            0.176277  0.319253 -0.226791   0.541    0.4    
    #> intref        0.726264  1.344480 -1.580572   1.641    0.4    
    #> intmed       -0.002639  0.223692 -0.273795   0.262    0.8    
    #> cde(prop)     0.147045  1.755498 -0.548080   3.466    0.8    
    #> intref(prop)  0.676679  1.984259 -2.977920   1.765    0.4    
    #> intmed(prop) -0.002459  0.294839 -0.496321   0.216    0.8    
    #> pnie(prop)    0.178736  0.461546 -0.061783   1.010    0.4    
    #> pm(overall)   0.176277  0.319253 -0.226791   0.541    0.4    
    #> int(overall)  0.674219  2.205041 -3.452921   1.610    0.4    
    #> pe(overall)   0.852955  1.755498 -2.466233   1.548    0.4    
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
    #>         estRR   lowerRR  upperRR Evalue.estRR Evalue.lowerRR Evalue.upperRR
    #> cde  1.080335 0.3976073 2.935369     1.374935       1.000000             NA
    #> pnde 1.541665 1.0592526 2.243781     2.455484       1.309779             NA
    #> tnde 1.539674 0.9430319 2.513802     2.451223       1.000000             NA
    #> pnie 1.098477 0.8480441 1.422865     1.427377       1.000000             NA
    #> tnie 1.097059 0.8344933 1.442238     1.423370       1.000000             NA
    #> te   1.691297 1.2299154 2.325757     2.772587       1.761683             NA

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
    #>               Estimate Std.error   95% CIL 95% CIU  P.val    
    #> cde            0.15648   3.82863  -5.35497   3.795    0.8    
    #> pnde           0.82160   0.99174  -0.08374   2.333    0.4    
    #> tnde           0.85828   0.82589   0.50644   2.385 <2e-16 ***
    #> pnie           0.23205   0.42894  -0.35021   0.613    0.8    
    #> tnie           0.26873   0.35972  -0.23138   0.597    0.4    
    #> te             1.09033   0.68715   0.47209   2.109 <2e-16 ***
    #> pm             0.24646   0.55379  -0.09910   1.260    0.4    
    #> intref         0.66513   2.96604  -1.55974   5.271    0.8    
    #> intmed         0.03668   0.41372  -0.04396   0.856    0.4    
    #> cde(prop)      0.14351   6.33758 -12.06050   1.799    0.8    
    #> intref(prop)   0.61002   5.81960  -0.93223  11.822    0.8    
    #> intmed(prop)   0.03364   0.94790  -0.03883   1.935    0.4    
    #> pnie(prop)     0.21282   0.51409  -0.73301   0.531    0.8    
    #> pm(overall)    0.24646   0.55379  -0.09910   1.260    0.4    
    #> int(overall)   0.64367   6.75585  -0.92019  13.757    0.8    
    #> pe(overall)    0.85649   6.33758  -0.79906  13.061    0.8    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> ----------------------------------------------------------------
    #> 
    #> Measurement error 2: 
    #> 0.2
    #> Measurement error correction for measurement error 2: 
    #>              Estimate Std.error  95% CIL 95% CIU  P.val    
    #> cde           0.15206   1.17683 -1.89738   1.047    0.8    
    #> pnde          0.72537   0.30383  0.16490   0.916 <2e-16 ***
    #> tnde          0.93943   0.24902  0.31661   0.905 <2e-16 ***
    #> pnie          0.43056   0.47459 -0.07958   0.984    0.8    
    #> tnie          0.64462   0.56698 -0.20874   0.985    0.8    
    #> te            1.36999   0.57943  0.25960   1.546 <2e-16 ***
    #> pm            0.47053   0.71971 -0.74129   0.838    0.8    
    #> intref        0.57331   0.90497 -0.15674   2.062    0.4    
    #> intmed        0.21406   0.16783 -0.12916   0.259    0.8    
    #> cde(prop)     0.11100   1.33958 -2.17299   0.679    0.8    
    #> intref(prop)  0.41848   1.57489 -0.09916   3.546    0.4    
    #> intmed(prop)  0.15625   0.29901 -0.45714   0.277    0.8    
    #> pnie(prop)    0.31428   0.44712 -0.28415   0.698    0.8    
    #> pm(overall)   0.47053   0.71971 -0.74129   0.838    0.8    
    #> int(overall)  0.57472   1.48749 -0.09868   3.148    0.4    
    #> pe(overall)   0.88900   1.33958  0.32141   3.173 <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> ----------------------------------------------------------------
    #> 
    #> Measurement error 3: 
    #> 0.3
    #> Measurement error correction for measurement error 3: 
    #>                Estimate  Std.error    95% CIL 95% CIU  P.val    
    #> cde           0.1431537  1.5983526 -1.3771712   2.370    0.8    
    #> pnde          0.7927200  0.4416151 -0.0038698   1.021    0.4    
    #> tnde          0.9078601  0.3867782  0.3032882   1.137 <2e-16 ***
    #> pnie          0.3154644  0.3688073 -0.0038862   0.812    0.4    
    #> tnie          0.4306045  0.4878730 -0.0007049   1.240    0.4    
    #> te            1.2233245  0.3520672  0.9135701   1.782 <2e-16 ***
    #> pm            0.3519953  0.3991242 -0.0228303   0.998    0.4    
    #> intref        0.6495663  1.4511080 -1.8054670   1.511    0.8    
    #> intmed        0.1151401  0.2662153 -0.2024188   0.454    0.8    
    #> cde(prop)     0.1170202  1.3534117 -1.1387930   2.113    0.8    
    #> intref(prop)  0.5309845  1.3600182 -1.6334465   1.648    0.8    
    #> intmed(prop)  0.0941206  0.2229648 -0.1862224   0.367    0.8    
    #> pnie(prop)    0.2578747  0.3183256 -0.0165466   0.726    0.4    
    #> pm(overall)   0.3519953  0.3991242 -0.0228303   0.998    0.4    
    #> int(overall)  0.6251051  1.4946480 -1.8103278   1.661    0.8    
    #> pe(overall)   0.8829798  1.3534117 -1.1133394   2.139    0.4    
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
    #> cde           -0.22215   5.19387 -11.89332   0.428    0.4    
    #> pnde           0.61920   1.85699  -2.70154   1.573    0.4    
    #> tnde           0.96955   0.69523   0.07939   1.677    0.4    
    #> pnie           0.28946   1.36567  -2.72336   0.425    0.4    
    #> tnie           0.63982   0.95929  -0.66396   1.564    0.8    
    #> te             1.25902   1.39427  -1.30747   2.075    0.8    
    #> pm             0.50819   1.54726  -1.05597   2.740    0.8    
    #> intref         0.84135   3.36210   1.14548   9.211 <2e-16 ***
    #> intmed         0.35035   2.09423  -0.49497   4.216    0.8    
    #> cde(prop)     -0.17645   7.05190  -3.71012  12.428    0.8    
    #> intref(prop)   0.66826   7.82695 -14.00734   4.385    0.8    
    #> intmed(prop)   0.27828   1.98955  -2.93998   2.049    0.8    
    #> pnie(prop)     0.22991   0.95491  -0.28435   1.966    0.8    
    #> pm(overall)    0.50819   1.54726  -1.05597   2.740    0.8    
    #> int(overall)   0.94654   7.79527 -12.34392   4.969    0.8    
    #> pe(overall)    1.17645   7.05190 -11.42807   4.710    0.8    
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
    #>               Estimate Std.error   95% CIL 95% CIU  P.val    
    #> cde            0.19647   5.49040  -8.47836   4.038    0.8    
    #> pnde           1.15026   1.84218  -1.39287   2.626    0.8    
    #> tnde           1.20758   1.26473  -0.54455   2.343    0.4    
    #> pnie           0.23829   0.60821  -0.93480   0.473    0.4    
    #> tnie           0.29561   0.32129  -0.04016   0.699    0.4    
    #> te             1.44587   1.57229  -0.69384   2.618    0.8    
    #> pm             0.20445   0.64083  -1.14324   0.174    0.8    
    #> intref         0.95379   3.75741  -1.41187   7.085    0.8    
    #> intmed         0.05732   0.87174  -0.37901   1.634    0.8    
    #> cde(prop)      0.13588   5.14102   1.16042  12.085 <2e-16 ***
    #> intref(prop)   0.65967   4.56784 -10.07289  -0.247 <2e-16 ***
    #> intmed(prop)   0.03964   1.04315  -2.30945   0.099    0.4    
    #> pnie(prop)     0.16481   0.67300  -0.25386   1.360    0.4    
    #> pm(overall)    0.20445   0.64083  -1.14324   0.174    0.8    
    #> int(overall)   0.69931   5.58770 -12.38233  -0.293 <2e-16 ***
    #> pe(overall)    0.86412   5.14102 -11.08497  -0.160 <2e-16 ***
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
