<!-- README.md is generated from README.Rmd. Please edit that file -->

CMAverse <img src="man/figures/logo.png" align="right" width="200" />
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
    #>               Estimate Std.error   95% CIL 95% CIU  P.val    
    #> cde          -0.422462  1.391651 -1.916640   1.197    0.4    
    #> pnde          0.497881  0.766442  0.208589   1.933 <2e-16 ***
    #> tnde          0.650642  0.197800  0.606790   1.064 <2e-16 ***
    #> pnie          0.683343  0.435277  0.009697   1.102    0.4    
    #> tnie          0.836104  0.741239 -0.274757   1.509    0.4    
    #> te            1.333985  0.424994  0.828100   1.753 <2e-16 ***
    #> pm            0.626771  0.429307 -0.140242   0.872    0.4    
    #> intref        0.920343  0.713249  0.736137   2.378 <2e-16 ***
    #> intmed        0.152761  0.591729 -0.873172   0.407    0.4    
    #> cde(prop)    -0.316692  0.882464 -1.417483   0.645    0.4    
    #> intref(prop)  0.689920  0.520817  0.468821   1.738 <2e-16 ***
    #> intmed(prop)  0.114515  0.372755 -0.491325   0.347    0.4    
    #> pnie(prop)    0.512257  0.269173 -0.017086   0.637    0.4    
    #> pm(overall)   0.626771  0.429307 -0.140242   0.872    0.4    
    #> int(overall)  0.804435  0.862602 -0.018877   1.955    0.4    
    #> pe(overall)   1.316692  0.882464  0.355133   2.417 <2e-16 ***
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
    #> cde  0.8567984 0.3165065 2.319395     1.608803             NA              1
    #> pnde 1.1997869 0.6932823 2.076338     1.689380       1.000000             NA
    #> tnde 1.2687469 1.1012933 1.461662     1.852675       1.435290             NA
    #> pnie 1.2840160 0.9403629 1.753256     1.887904       1.000000             NA
    #> tnie 1.3578172 0.7988769 2.307824     2.054847       1.000000             NA
    #> te   1.6290912 1.2018938 2.208130     2.641439       1.694494             NA

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
    #> cde          -0.42452   1.07011 -0.98809   1.566    0.4    
    #> pnde          0.61305   0.44173  0.35418   1.381 <2e-16 ***
    #> tnde          0.74492   0.31046  0.52294   1.210 <2e-16 ***
    #> pnie          0.75554   0.47145  0.19665   1.315 <2e-16 ***
    #> tnie          0.88741   0.57762  0.37075   1.725 <2e-16 ***
    #> te            1.50046   0.59493  1.02560   2.365 <2e-16 ***
    #> pm            0.59142   0.22880  0.27863   0.786 <2e-16 ***
    #> intref        1.03757   0.85150 -0.23647   1.906    0.4    
    #> intmed        0.13187   0.24581 -0.18320   0.414    0.4    
    #> cde(prop)    -0.28292   0.56654 -0.73799   0.691    0.4    
    #> intref(prop)  0.69150   0.60169 -0.08778   1.425    0.4    
    #> intmed(prop)  0.08789   0.11743 -0.07786   0.190    0.4    
    #> pnie(prop)    0.50354   0.18865  0.15938   0.596 <2e-16 ***
    #> pm(overall)   0.59142   0.22880  0.27863   0.786 <2e-16 ***
    #> int(overall)  0.77939   0.67530 -0.15300   1.568    0.4    
    #> pe(overall)   1.28292   0.56654  0.30858   1.738 <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> ----------------------------------------------------------------
    #> 
    #> Measurement error 2: 
    #> 0.2
    #> Measurement error correction for measurement error 2: 
    #>               Estimate Std.error   95% CIL 95% CIU  P.val    
    #> cde          -0.431083  0.623191 -1.056635   0.391    0.8    
    #> pnde          0.488798  0.449650  0.281434   1.282 <2e-16 ***
    #> tnde          0.659727  0.445833  0.321212   1.362 <2e-16 ***
    #> pnie          0.637035  0.527422 -0.169599   1.013    0.4    
    #> tnie          0.807964  0.561800  0.017178   1.247    0.4    
    #> te            1.296762  0.704427  0.660914   2.313 <2e-16 ***
    #> pm            0.623062  0.313882  0.005606   0.786    0.4    
    #> intref        0.919882  0.574976  0.366004   1.742 <2e-16 ***
    #> intmed        0.170929  0.131836 -0.056937   0.261    0.4    
    #> cde(prop)    -0.332430  0.404078 -0.692423   0.172    0.8    
    #> intref(prop)  0.709368  0.360804  0.380351   1.243 <2e-16 ***
    #> intmed(prop)  0.131812  0.144164 -0.091954   0.260    0.4    
    #> pnie(prop)    0.491251  0.373091 -0.240207   0.638    0.4    
    #> pm(overall)   0.623062  0.313882  0.005606   0.786    0.4    
    #> int(overall)  0.841180  0.436368  0.400675   1.303 <2e-16 ***
    #> pe(overall)   1.332430  0.404078  0.828462   1.692 <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> ----------------------------------------------------------------
    #> 
    #> Measurement error 3: 
    #> 0.3
    #> Measurement error correction for measurement error 3: 
    #>              Estimate Std.error  95% CIL 95% CIU  P.val    
    #> cde          -0.44357   0.96294 -1.23707   0.872    0.4    
    #> pnde          0.47112   0.58004 -0.34512   1.124    0.4    
    #> tnde          0.61979   0.40127 -0.18595   0.731    0.4    
    #> pnie          0.73597   0.49075  0.14903   1.361 <2e-16 ***
    #> tnie          0.88463   0.39428  0.29618   1.259 <2e-16 ***
    #> te            1.35576   0.88017 -0.03077   2.092    0.4    
    #> pm            0.65250   0.93417 -1.18949   0.873    0.4    
    #> intref        0.91470   0.62501  0.23363   1.719 <2e-16 ***
    #> intmed        0.14866   0.29608 -0.39491   0.283    0.4    
    #> cde(prop)    -0.32718   3.38938 -0.84133   6.532    0.8    
    #> intref(prop)  0.67468   2.45917 -4.34217   0.977    0.4    
    #> intmed(prop)  0.10965   0.44667 -0.79361   0.234    0.8    
    #> pnie(prop)    0.54285   0.50843 -0.40379   0.653    0.4    
    #> pm(overall)   0.65250   0.93417 -1.18949   0.873    0.4    
    #> int(overall)  0.78433   2.88944 -5.13578   1.199    0.8    
    #> pe(overall)   1.32718   3.38938 -5.53166   1.841    0.4    
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
    #> cde          -0.63979   0.82340 -1.15969   0.862    0.8    
    #> pnde          0.43371   0.65534 -0.16001   1.253    0.4    
    #> tnde          0.63960   0.49237  0.44614   1.370 <2e-16 ***
    #> pnie          0.54310   0.48613 -0.20714   0.990    0.4    
    #> tnie          0.74899   0.34027  0.32378   1.108 <2e-16 ***
    #> te            1.18270   0.85789  0.23900   2.320 <2e-16 ***
    #> pm            0.63329   0.68335  0.23164   1.858 <2e-16 ***
    #> intref        1.07350   0.43384  0.32890   1.338 <2e-16 ***
    #> intmed        0.20589   0.23076  0.09967   0.606 <2e-16 ***
    #> cde(prop)    -0.54096   1.64955 -3.39675   0.360    0.8    
    #> intref(prop)  0.90767   1.06496  0.18820   2.539 <2e-16 ***
    #> intmed(prop)  0.17409   1.37466  0.05066   2.953 <2e-16 ***
    #> pnie(prop)    0.45920   0.73095 -1.10669   0.555    0.4    
    #> pm(overall)   0.63329   0.68335  0.23164   1.858 <2e-16 ***
    #> int(overall)  1.08176   2.29539  0.23971   5.492 <2e-16 ***
    #> pe(overall)   1.54096   1.64955  0.63979   4.397 <2e-16 ***
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
    #>              Estimate Std.error 95% CIL 95% CIU  P.val    
    #> cde           -0.4913    1.2935 -3.6031  -0.687 <2e-16 ***
    #> pnde           0.7672    0.4312 -0.6431   0.367    0.8    
    #> tnde           1.2053    0.4048  0.3892   1.323 <2e-16 ***
    #> pnie           1.5440    0.4272  0.2996   1.349 <2e-16 ***
    #> tnie           1.9821    0.6015  1.0764   2.450 <2e-16 ***
    #> te             2.7493    0.3497  1.0073   1.840 <2e-16 ***
    #> pm             0.7209    0.2583  0.7493   1.355 <2e-16 ***
    #> intref         1.2585    1.3384  0.8385   3.753 <2e-16 ***
    #> intmed         0.4380    0.3084  0.5027   1.135 <2e-16 ***
    #> cde(prop)     -0.1787    1.2675 -3.4415  -0.450 <2e-16 ***
    #> intref(prop)   0.4578    1.2261  0.4672   3.313 <2e-16 ***
    #> intmed(prop)   0.1593    0.1362  0.3543   0.683 <2e-16 ***
    #> pnie(prop)     0.5616    0.2250  0.1880   0.747 <2e-16 ***
    #> pm(overall)    0.7209    0.2583  0.7493   1.355 <2e-16 ***
    #> int(overall)   0.6171    1.2399  1.0188   3.858 <2e-16 ***
    #> pe(overall)    1.1787    1.2675  1.4500   4.441 <2e-16 ***
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
