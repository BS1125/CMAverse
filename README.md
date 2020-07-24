<img src="man/figures/logo.png" width="50%" style="display: block; margin: auto 0 auto auto;" />

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

<table>
<caption>Supported Data Types and Functionalities of <code>cmest</code></caption>
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
regression-based approach* and *the g-formula approach*.

Installation
============

The latest version can be installed via:

    devtools::install_github("LindaValeri/CMAverse")

Load `CMAverse`:

    library(CMAverse)

Quickstart Guide
================

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
          prec = c("C1", "C2"), postc = NULL,
          node = FALSE, text_col = "black")

![](man/figures/plot_dag-1.png)

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
    ##              Estimate Std.error 95% CIL 95% CIU  P.val    
    ## cde           0.50192   0.21744 0.07152   0.911   0.03 *  
    ## pnde          1.23082   0.07086 1.10602   1.376 <2e-16 ***
    ## tnde          1.37577   0.07457 1.24440   1.522 <2e-16 ***
    ## pnie          0.49086   0.05204 0.38308   0.583 <2e-16 ***
    ## tnie          0.63581   0.05862 0.52685   0.748 <2e-16 ***
    ## te            1.86663   0.07704 1.70766   2.014 <2e-16 ***
    ## pm            0.34062   0.02734 0.28965   0.390 <2e-16 ***
    ## intref        0.72890   0.19583 0.36740   1.110 <2e-16 ***
    ## intmed        0.14495   0.04033 0.07272   0.230 <2e-16 ***
    ## cde(prop)     0.26889   0.11424 0.03848   0.497   0.03 *  
    ## intref(prop)  0.39049   0.10565 0.19740   0.612 <2e-16 ***
    ## intmed(prop)  0.07765   0.02106 0.03812   0.119 <2e-16 ***
    ## pnie(prop)    0.26296   0.02568 0.20630   0.308 <2e-16 ***
    ## pm(overall)   0.34062   0.02734 0.28965   0.390 <2e-16 ***
    ## int(overall)  0.46815   0.12493 0.23719   0.725 <2e-16 ***
    ## pe(overall)   0.73111   0.11424 0.50267   0.962 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    plot(est) +
      theme(axis.text.x = element_text(angle = 30, vjust = 0.8))

![](man/figures/plot_cmest-1.png)

Lastly, let’s conduct sensitivity analysis for the results. Sensitivity
analysis for unmeasured confounding:

    cmsens(object = est, sens = "uc")

    ## Sensitivity Analysis For Unmeasured Confounding 
    ## 
    ## Evalues on the ratio scale: 
    ##         estRR  lowerRR  upperRR Evalue.estRR Evalue.lowerRR Evalue.upperRR
    ## cde  1.217625 1.030504 1.438722     1.732392       1.207803             NA
    ## pnde 1.620686 1.534912 1.711254     2.623650       2.441026             NA
    ## tnde 1.715515 1.620100 1.816549     2.823430       2.622409             NA
    ## pnie 1.212352 1.164894 1.261743     1.719743       1.603169             NA
    ## tnie 1.283289 1.226843 1.342331     1.886233       1.754386             NA
    ## te   2.079808 1.960424 2.206463     3.578405       3.332589             NA

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
    ##              Estimate Std.error 95% CIL 95% CIU  P.val    
    ## cde           0.49960   0.21135 0.07527   0.967   0.04 *  
    ## pnde          1.23013   0.07929 1.09667   1.397 <2e-16 ***
    ## tnde          1.37450   0.07663 1.23527   1.541 <2e-16 ***
    ## pnie          0.49434   0.05377 0.40500   0.595 <2e-16 ***
    ## tnie          0.63872   0.06288 0.52792   0.758 <2e-16 ***
    ## te            1.86884   0.08220 1.72517   2.048 <2e-16 ***
    ## pm            0.34177   0.03034 0.28742   0.412 <2e-16 ***
    ## intref        0.73052   0.17646 0.39263   1.119 <2e-16 ***
    ## intmed        0.14437   0.03780 0.07347   0.228 <2e-16 ***
    ## cde(prop)     0.26733   0.10963 0.04101   0.499   0.04 *  
    ## intref(prop)  0.39090   0.09625 0.19737   0.603 <2e-16 ***
    ## intmed(prop)  0.07725   0.02004 0.03716   0.122 <2e-16 ***
    ## pnie(prop)    0.26452   0.02588 0.21814   0.316 <2e-16 ***
    ## pm(overall)   0.34177   0.03034 0.28742   0.412 <2e-16 ***
    ## int(overall)  0.46815   0.11494 0.23475   0.725 <2e-16 ***
    ## pe(overall)   0.73267   0.10963 0.50129   0.959 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------
    ## 
    ## Measurement error 2: 
    ## 0.2
    ## Measurement error correction for measurement error 2: 
    ##              Estimate Std.error 95% CIL 95% CIU  P.val    
    ## cde           0.49219   0.19219 0.11754   0.853   0.01 ** 
    ## pnde          1.21853   0.07155 1.09341   1.360 <2e-16 ***
    ## tnde          1.36615   0.07560 1.23568   1.518 <2e-16 ***
    ## pnie          0.49865   0.05854 0.40067   0.618 <2e-16 ***
    ## tnie          0.64627   0.06866 0.52092   0.786 <2e-16 ***
    ## te            1.86480   0.08819 1.70322   2.050 <2e-16 ***
    ## pm            0.34656   0.02975 0.28861   0.403 <2e-16 ***
    ## intref        0.72634   0.16973 0.42906   1.045 <2e-16 ***
    ## intmed        0.14762   0.03615 0.07958   0.222 <2e-16 ***
    ## cde(prop)     0.26394   0.10150 0.06426   0.450   0.01 ** 
    ## intref(prop)  0.38950   0.09215 0.22538   0.569 <2e-16 ***
    ## intmed(prop)  0.07916   0.01854 0.04289   0.113 <2e-16 ***
    ## pnie(prop)    0.26740   0.02672 0.22060   0.328 <2e-16 ***
    ## pm(overall)   0.34656   0.02975 0.28861   0.403 <2e-16 ***
    ## int(overall)  0.46866   0.10863 0.26379   0.668 <2e-16 ***
    ## pe(overall)   0.73606   0.10150 0.54965   0.936 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------
    ## 
    ## Measurement error 3: 
    ## 0.3
    ## Measurement error correction for measurement error 3: 
    ##              Estimate Std.error 95% CIL 95% CIU  P.val    
    ## cde           0.47807   0.18633 0.09341   0.830   0.02 *  
    ## pnde          1.20710   0.07228 1.06028   1.343 <2e-16 ***
    ## tnde          1.35412   0.07359 1.20259   1.489 <2e-16 ***
    ## pnie          0.51579   0.05598 0.41954   0.628 <2e-16 ***
    ## tnie          0.66281   0.06486 0.55487   0.780 <2e-16 ***
    ## te            1.86991   0.08485 1.69052   2.027 <2e-16 ***
    ## pm            0.35446   0.02911 0.30183   0.409 <2e-16 ***
    ## intref        0.72903   0.15875 0.43605   1.058 <2e-16 ***
    ## intmed        0.14702   0.03629 0.06704   0.225 <2e-16 ***
    ## cde(prop)     0.25566   0.09751 0.04960   0.429   0.02 *  
    ## intref(prop)  0.38988   0.08653 0.23359   0.559 <2e-16 ***
    ## intmed(prop)  0.07862   0.01895 0.03674   0.115 <2e-16 ***
    ## pnie(prop)    0.27584   0.02599 0.23033   0.327 <2e-16 ***
    ## pm(overall)   0.35446   0.02911 0.30183   0.409 <2e-16 ***
    ## int(overall)  0.46850   0.10381 0.26742   0.686 <2e-16 ***
    ## pe(overall)   0.74434   0.09751 0.57051   0.950 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------

    plot(me1) +
      theme(axis.text.x = element_text(angle = 30, vjust = 0.8))

![](man/figures/plot_cmsens_me_con-1.png)

Then, assume that the exposure was measured with error. Sensitivity
analysis using MCSIMEX with two assumed misclassification matrices:

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
    ##              Estimate Std.error 95% CIL 95% CIU  P.val    
    ## cde           0.67069   0.24454 0.21348   1.121   0.01 ** 
    ## pnde          1.44012   0.09013 1.26336   1.613 <2e-16 ***
    ## tnde          1.61731   0.08573 1.46060   1.765 <2e-16 ***
    ## pnie          0.49755   0.05453 0.38354   0.600 <2e-16 ***
    ## tnie          0.67473   0.06564 0.55013   0.797 <2e-16 ***
    ## te            2.11485   0.09240 1.91705   2.296 <2e-16 ***
    ## pm            0.31905   0.02892 0.27102   0.375 <2e-16 ***
    ## intref        0.76943   0.19887 0.39908   1.164 <2e-16 ***
    ## intmed        0.17719   0.04684 0.08191   0.272 <2e-16 ***
    ## cde(prop)     0.31713   0.11199 0.10436   0.519   0.01 ** 
    ## intref(prop)  0.36382   0.09796 0.18080   0.575 <2e-16 ***
    ## intmed(prop)  0.08378   0.02226 0.03916   0.127 <2e-16 ***
    ## pnie(prop)    0.23526   0.02383 0.18933   0.281 <2e-16 ***
    ## pm(overall)   0.31905   0.02892 0.27102   0.375 <2e-16 ***
    ## int(overall)  0.44761   0.11849 0.22828   0.690 <2e-16 ***
    ## pe(overall)   0.68287   0.11199 0.48118   0.896 <2e-16 ***
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
    ## cde           0.83047   0.28236 0.29720   1.375 <2e-16 ***
    ## pnde          1.65142   0.09630 1.47902   1.836 <2e-16 ***
    ## tnde          1.86738   0.09925 1.68033   2.066 <2e-16 ***
    ## pnie          0.49659   0.05870 0.37617   0.623 <2e-16 ***
    ## tnie          0.71255   0.07099 0.57810   0.861 <2e-16 ***
    ## te            2.36397   0.10177 2.19258   2.577 <2e-16 ***
    ## pm            0.30142   0.02709 0.24741   0.345 <2e-16 ***
    ## intref        0.82095   0.23563 0.34472   1.287 <2e-16 ***
    ## intmed        0.21596   0.06292 0.08951   0.330 <2e-16 ***
    ## cde(prop)     0.35130   0.11686 0.12748   0.589 <2e-16 ***
    ## intref(prop)  0.34727   0.10228 0.14611   0.550 <2e-16 ***
    ## intmed(prop)  0.09135   0.02628 0.03898   0.144 <2e-16 ***
    ## pnie(prop)    0.21007   0.02353 0.16158   0.253 <2e-16 ***
    ## pm(overall)   0.30142   0.02709 0.24741   0.345 <2e-16 ***
    ## int(overall)  0.43863   0.12687 0.18057   0.678 <2e-16 ***
    ## pe(overall)   0.64870   0.11686 0.41101   0.873 <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------------------

    plot(me2) +
      theme(axis.text.x = element_text(angle = 30, vjust = 0.8))

![](man/figures/plot_cmsens_me_cat-1.png)

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
