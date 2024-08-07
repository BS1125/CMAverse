<!-- README.md is generated from README.Rmd. Please edit that file -->

# CMAverse: a suite of functions for causal mediation analysis<img src="man/figures/logo.png" align="right" width="240" />

[![Project Status:
Active](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![R build
status](https://github.com/BS1125/CMAverse/workflows/R-CMD-check/badge.svg)](https://github.com/BS1125/CMAverse/actions)
[![Codecov test
coverage](https://codecov.io/gh/BS1125/CMAverse/branch/master/graph/badge.svg)](https://codecov.io/gh/BS1125/CMAverse)

## About the Package

The R package `CMAverse` provides a suite of functions for reproducible
causal mediation analysis including `cmdag` for DAG visualization,
`cmest` for statistical modeling and `cmsens` for sensitivity analysis.

See the package [website](https://bs1125.github.io/CMAverse/) for a
quickstart guide, an overview of statistical modeling approaches and
examples.

Cite the paper: [CMAverse a suite of functions for reproducible causal
mediation
analyses](https://journals.lww.com/epidem/citation/9000/cmaverse__a_suite_of_functions_for_reproducible.98264.aspx)

We welcome your feedback and questions:

-   Email <bs3141@columbia.edu> for general questions

-   Email <zw2899@cumc.columbia.edu> for questions related to
    `cmest_multistate`

### DAG Visualization

`cmdag` visualizes causal relationships via a directed acyclic graph
(DAG).

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
*the g-formula approach* by [Robins
(1986)](https://www.sciencedirect.com/science/article/pii/0270025586900886).

`cmest` currently supports a single exposure, multiple sequential
mediators and a single outcome. When multiple mediators are of interest,
`cmest` estimates the joint mediated effect through the set of
mediators. `cmest` also allows for time varying confounders preceding
mediators. The two causal scenarios supported are:

1.  There are no confounders affected by the exposure:

![](man/figures/unnamed-chunk-1-1.png)

1.  There are mediator-outcome confounders affected by the exposure and
    these confounders precede all of the mediators:

![](man/figures/unnamed-chunk-2-1.png)

<table>
<caption>Table: Supported Data Types and Functionalities of
<code>cmest</code></caption>
<thead>
<tr class="header">
<th></th>
<th style="text-align: center;">rb</th>
<th style="text-align: center;">wb</th>
<th style="text-align: center;">iorw</th>
<th style="text-align: center;">ne</th>
<th style="text-align: center;">msm</th>
<th style="text-align: center;">gformula<a href="#fn1"
class="footnote-ref" id="fnref1"
role="doc-noteref"><sup>1</sup></a></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Continuous Y<a href="#fn2" class="footnote-ref" id="fnref2"
role="doc-noteref"><sup>2</sup></a></td>
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
<td style="text-align: center;">×</td>
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
<td style="text-align: center;">×<a href="#fn3" class="footnote-ref"
id="fnref3" role="doc-noteref"><sup>3</sup></a></td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×<a href="#fn4" class="footnote-ref"
id="fnref4" role="doc-noteref"><sup>4</sup></a></td>
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
<td style="text-align: center;">×<a href="#fn5" class="footnote-ref"
id="fnref5" role="doc-noteref"><sup>5</sup></a></td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×<a href="#fn6" class="footnote-ref"
id="fnref6" role="doc-noteref"><sup>6</sup></a></td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Mediator-outcome Confounder(s) Affected by A</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>2-way Decomposition</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>4-way Decomposition</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Estimation: Closed-form Parameter Function</td>
<td style="text-align: center;">√<a href="#fn7" class="footnote-ref"
id="fnref7" role="doc-noteref"><sup>7</sup></a></td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
</tr>
<tr class="even">
<td>Estimation: Direct Counterfactual Imputation</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Inference: Delta Method</td>
<td style="text-align: center;">√<a href="#fn8" class="footnote-ref"
id="fnref8" role="doc-noteref"><sup>8</sup></a></td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
</tr>
<tr class="even">
<td>Inference: Bootstrapping</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="odd">
<td>Marginal Effects</td>
<td style="text-align: center;">√<a href="#fn9" class="footnote-ref"
id="fnref9" role="doc-noteref"><sup>9</sup></a></td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
<td style="text-align: center;">√</td>
</tr>
<tr class="even">
<td>Effects Conditional on C</td>
<td style="text-align: center;">√<a href="#fn10" class="footnote-ref"
id="fnref10" role="doc-noteref"><sup>10</sup></a></td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
<td style="text-align: center;">×</td>
</tr>
</tbody>
</table>
<section class="footnotes footnotes-end-of-document"
role="doc-endnotes">
<hr />
<ol>
<li id="fn1" role="doc-endnote"><p>rb: the regression-based approach;
wb: the weighting-based approach; iorw: the inverse odds ratio weighting
approach; ne: the natural effect model; msm: the marginal structural
model; gformula: the g-formula approach.<a href="#fnref1"
class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn2" role="doc-endnote"><p>Y denotes the outcome, A denotes the
exposure, M denotes the mediator(s) and C denotes the exposure-outcome
confounder(s), the exposure-mediator confounder(s) and the
mediator-outcome confounder(s) not affected by the exposure.<a
href="#fnref2" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn3" role="doc-endnote"><p>continuous A is not supported when C
is not empty; otherwise, it is supported.<a href="#fnref3"
class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn4" role="doc-endnote"><p>continuous A is not supported when C
is not empty; otherwise, it is supported.<a href="#fnref4"
class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn5" role="doc-endnote"><p>count A is not supported when C is
not empty; otherwise, it is supported.<a href="#fnref5"
class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn6" role="doc-endnote"><p>count A is not supported when C is
not empty; otherwise, it is supported.<a href="#fnref6"
class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn7" role="doc-endnote"><p>closed-form parameter function
estimation only supports a single mediator.<a href="#fnref7"
class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn8" role="doc-endnote"><p>delta method inference is available
only when closed-form parameter function estimation is used.<a
href="#fnref8" class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn9" role="doc-endnote"><p>marginal effects are estimated when
direct counterfactual imputation estimation is used.<a href="#fnref9"
class="footnote-back" role="doc-backlink">↩︎</a></p></li>
<li id="fn10" role="doc-endnote"><p>conditional effects are estimated
when closed-form parameter function estimation is used.<a
href="#fnref10" class="footnote-back"
role="doc-backlink">↩︎</a></p></li>
</ol>
</section>

Table: Supported Data Types and Functionalities of `cmest`

### Multiple Imputation

`cmest` provides the option `multimp = TRUE` to perform multiple
imputations for a dataset with missing values.

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
al. (2006)](https://pubmed.ncbi.nlm.nih.gov/16542233/). Sensitivity
analysis for measurement error is currently available for *the
regression-based approach* and *the g-formula approach*.

## Installation

The latest version can be installed via:

    devtools::install_github("BS1125/CMAverse")

Load `CMAverse`:

    library(CMAverse)

## References

Valeri L, Vanderweele TJ (2013). Mediation analysis allowing for
exposure-mediator interactions and causal interpretation: theoretical
assumptions and implementation with SAS and SPSS macros. Psychological
Methods. 18(2): 137 - 150.

VanderWeele TJ, Vansteelandt S (2014). Mediation analysis with multiple
mediators. Epidemiologic Methods. 2(1): 95 - 115.

Tchetgen Tchetgen EJ (2013). Inverse odds ratio-weighted estimation for
causal mediation analysis. Statistics in medicine. 32: 4567 - 4580.

Nguyen QC, Osypuk TL, Schmidt NM, Glymour MM, Tchetgen Tchetgen EJ
(2015). Practical guidance for conducting mediation analysis with
multiple mediators using inverse odds ratio weighting. American Journal
of Epidemiology. 181(5): 349 - 356.

VanderWeele TJ, Tchetgen Tchetgen EJ (2017). Mediation analysis with
time varying exposures and mediators. Journal of the Royal Statistical
Society: Series B (Statistical Methodology). 79(3): 917 - 938.

Robins JM (1986). A new approach to causal inference in mortality
studies with a sustained exposure period-Application to control of the
healthy worker survivor effect. Mathematical Modelling. 7: 1393 - 1512.

Vansteelandt S, Bekaert M, Lange T (2012). Imputation Strategies for the
Estimation of Natural Direct and Indirect Effects. Epidemiologic
Methods. 1(1): 131 - 158.

Steen J, Loeys T, Moerkerke B, Vansteelandt S (2017). Medflex: an R
package for flexible mediation analysis using natural effect models.
Journal of Statistical Software. 76(11).

VanderWeele TJ (2014). A unification of mediation and interaction: a
4-way decomposition. Epidemiology. 25(5): 749 - 61.

Imai K, Keele L, Tingley D (2010). A general approach to causal
mediation analysis. Psychological Methods. 15(4): 309 - 334.

Schomaker M, Heumann C (2018). Bootstrap inference when using multiple
imputation. Statistics in Medicine. 37(14): 2252 - 2266.

VanderWeele TJ, Ding P (2017). Sensitivity analysis in observational
research: introducing the E-Value. Annals of Internal Medicine. 167(4):
268 - 274.

Smith LH, VanderWeele TJ (2019). Mediational E-values: Approximate
sensitivity analysis for unmeasured mediator-outcome confounding.
Epidemiology. 30(6): 835 - 837.

Carrol RJ, Ruppert D, Stefanski LA, Crainiceanu C (2006). Measurement
Error in Nonlinear Models: A Modern Perspective, Second Edition. London:
Chapman & Hall.

Cook JR, Stefanski LA (1994). Simulation-extrapolation estimation in
parametric measurement error models. Journal of the American Statistical
Association, 89(428): 1314 - 1328.

Küchenhoff H, Mwalili SM, Lesaffre E (2006). A general method for
dealing with misclassification in regression: the misclassification
SIMEX. Biometrics. 62(1): 85 - 96.

Stefanski LA, Cook JR. Simulation-extrapolation: the measurement error
jackknife (1995). Journal of the American Statistical Association.
90(432): 1247 - 56.

Valeri L, Lin X, VanderWeele TJ (2014). Mediation analysis when a
continuous mediator is measured with error and the outcome follows a
generalized linear model. Statistics in medicine, 33(28): 4875–4890.

Efron B (1987). Better Bootstrap Confidence Intervals. Journal of the
American Statistical Association. 82(397): 171-185.
