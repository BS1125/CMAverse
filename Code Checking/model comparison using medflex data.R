# Case 5: Binary Outcome and Continuous Mediator Without Exposure-mediator Interaction

## Causal Effects and Standard Errors Estimated By the Natural Effect Model

library("medflex")
library("CMA")
data("UPBdata")

impData=neImpute(UPB~attbin+negaff+age,family = binomial,nMed=1,data=UPBdata)
neMod1=neModel(UPB~attbin0+attbin1+age,family = binomial(),expData=impData,se="robust")

res1=summary(neEffdecomp(neMod))$coefficients

res1[, "Estimate"]=exp(summary(neEffdecomp(neMod))$coefficients[, "Estimate"])

for (i in 1:length(res1[, "Std. Error"]))
  res1[, "Std. Error"][i] <- msm::deltamethod(as.formula("~exp(x1)"), res1[, "Estimate"][i],
                                              summary(neEffdecomp(neMod))$coefficients[, "Std. Error"][i])
res1[,c("Estimate","Std. Error")]
#                         Estimate Std. Error
# natural direct effect   1.473154   2.011951
# natural indirect effect 1.405780   1.155350
# total effect            2.070930   3.688186

expData=neWeight(UPB~attbin+negaff+age,family = binomial,nMed=1,data=UPBdata)
neMod2=neModel(UPB~attbin0+attbin1+age,family = binomial(),expData=expData,se="robust")

res2=summary(neEffdecomp(neMod2))$coefficients

res2[, "Estimate"]=exp(summary(neEffdecomp(neMod2))$coefficients[, "Estimate"])

for (i in 1:length(res2[, "Std. Error"]))
  res2[, "Std. Error"][i] <- msm::deltamethod(as.formula("~exp(x1)"), res2[, "Estimate"][i],
                                              summary(neEffdecomp(neMod2))$coefficients[, "Std. Error"][i])
res2[,c("Estimate","Std. Error")]
# Estimate Std. Error
# natural direct effect   1.387632   1.107489
# natural indirect effect 1.910029   3.062107
# total effect            2.650416   6.652369


causal_mediation(data = UPBdata, model = "ne",
                 outcome = "UPB", exposure = 'attbin', exposure.type = "binary",
                 mediator = c("negaff"), covariates = c("age"), EMint = FALSE,
                 yreg = "logistic")
#                         Estimate Std. Error
# natural direct effect   1.473154   2.011951
# natural indirect effect 1.405780   1.155350
# total effect            2.070930   3.688186

causal_mediation(data = UPBdata, model = "standard", est.method = "delta",
                 outcome = "UPB", exposure = 'attbin', exposure.type = "binary",
                 mediator = c("negaff"), covariates = c("age"), EMint = FALSE,
                 yreg = "logistic", mreg = "linear")
# $`decomp3way`
# cde_rr  cde_rr_se    pnde_rr pnde_rr_se    tnde_rr tnde_rr_se    pnie_rr pnie_rr_se    tnie_rr tnie_rr_se
# 1.5255287  0.3526752  1.5255287  0.3526752  1.5255287  0.3526752  1.4212528  0.1296386  1.4212528  0.1296386
# te_rr   te_rr_se         pm      pm_se
# 2.1681620  0.5101798  0.5501234  0.1420178
#
# $decomp4way
# cde_err         cde_err_se         intref_err      intref_err_se         intmed_err
# 0.51328355         0.35016113         0.01224514         0.03062008         0.22138045
# intmed_err_se            pie_err         pie_err_se          total_err       total_err_se
# 0.15376230         0.42125283         0.12963861         1.16816196         0.51017976
# cde_err_prop    cde_err_prop_se    intref_err_prop intref_err_prop_se    intmed_err_prop
# 0.43939417         0.14563562         0.01048240         0.02590011         0.18951177
# intmed_err_prop_se       pie_err_prop    pie_err_prop_se         overall_pm      overall_pm_se
# 0.04926173         0.36061166         0.16964168         0.55012344         0.14201776
# overall_int     overall_int_se         overall_pe      overall_pe_se
# 0.19999417         0.05538204         0.56060583         0.14563562

causal_mediation(data = UPBdata, model = "standard", est.method = "bootstrap",
                 outcome = "UPB", exposure = 'attbin', exposure.type = "binary",
                 mediator = c("negaff"), covariates = c("age"), EMint = FALSE,
                 yreg = "logistic", mreg = "linear")
# $`decomp3way`
# cde_rr  cde_rr_se    pnde_rr pnde_rr_se    tnde_rr tnde_rr_se    pnie_rr pnie_rr_se    tnie_rr tnie_rr_se
# 1.5255287  0.3628136  1.5255287  0.3628136  1.5255287  0.3628136  1.4212528  0.1331474  1.4212528  0.1331474
# te_rr   te_rr_se         pm      pm_se
# 2.1681620  0.5398181  0.5501234  0.2050216
#
# $decomp4way
# cde_err         cde_err_se         intref_err      intref_err_se         intmed_err
# 0.51328355         0.35642044         0.01224514         0.04347121         0.22138045
# intmed_err_se            pie_err         pie_err_se          total_err       total_err_se
# 0.16998658         0.42125283         0.13314739         1.16816196         0.53981812
# cde_err_prop    cde_err_prop_se    intref_err_prop intref_err_prop_se    intmed_err_prop
# 0.43939417         0.20239361         0.01048240         0.02851157         0.18951177
# intmed_err_prop_se       pie_err_prop    pie_err_prop_se         overall_pm      overall_pm_se
# 0.07895676         0.36061166         0.26976678         0.55012344         0.20502164
# overall_int     overall_int_se         overall_pe      overall_pe_se
# 0.19999417         0.08866213         0.56060583         0.20239361

causal_mediation(data = UPBdata, model = "standard", est.method = "simulation",
                 outcome = "UPB", exposure = 'attbin', exposure.type = "binary",
                 mediator = c("negaff"), covariates = c("age"), EMint = FALSE,
                 yreg = "logistic", mreg = "linear")
# $`decomp3way`
# cde_rr  cde_rr_se    pnde_rr pnde_rr_se    tnde_rr tnde_rr_se    pnie_rr pnie_rr_se    tnie_rr tnie_rr_se
# 1.32289510 0.20261692 1.34420452 0.21633124 1.29778233 0.18586176 1.27982703 0.08118826 1.23851690 0.07202260
# te_rr   te_rr_se         pm      pm_se
# 1.66132289 0.26381485 0.52194000 0.46371047
#
# $decomp4way
# cde_err         cde_err_se         intref_err      intref_err_se         intmed_err
# 0.36463497         0.22921285        -0.02043045         0.01577139         0.03729133
# intmed_err_se            pie_err         pie_err_se          total_err       total_err_se
# 0.02854533         0.27982703         0.08118826         0.66132289         0.26381485
# cde_err_prop    cde_err_prop_se    intref_err_prop intref_err_prop_se    intmed_err_prop
# 0.50511112         0.48091847        -0.02705112         0.02044808         0.04914968
# intmed_err_prop_se       pie_err_prop    pie_err_prop_se         overall_pm      overall_pm_se
# 0.04035837         0.47279032         0.50040351         0.52194000         0.46371047
# overall_int     overall_int_se         overall_pe      overall_pe_se
# 0.02209856         0.02146628         0.49488888         0.48091847

causal_mediation(data = UPBdata, model = "wb",
                 outcome = "UPB", exposure = 'attbin', exposure.type = "binary",
                 mediator = c("negaff"), covariates = c("age"), EMint = FALSE,
                 yreg = "logistic")
# RRcde     RRpnde     RRtnde     RRpnie     RRtnie       RRte         pm   RRcde_se  RRpnde_se  RRtnde_se
# 1.31242310 1.30340641 1.25469970 1.26654822 1.21921885 1.58913767 0.48499913 0.20345304 0.19482116 0.16188914
# RRpnie_se  RRtnie_se    RRte_se      pm_se
# 0.07511639 0.06478573 0.22876135 0.24215055

causal_mediation(data = UPBdata, model = "rb", est.method = "bootstrap",
                 outcome = "UPB", exposure = 'attbin', exposure.type = "binary",
                 mediator = c("negaff"), covariates = c("age"), EMint = FALSE,
                 yreg = "logistic", mreg = "linear")
# RRcde    RRpnde    RRtnde    RRpnie    RRtnie      RRte        pm  RRcde_se RRpnde_se RRtnde_se RRpnie_se
# 1.5255287 1.5255287 1.5255287 1.4212528 1.4212528 2.1681620 0.5501234 0.3854486 0.3854486 0.3854486 0.1335453
# RRtnie_se   RRte_se     pm_se
# 0.1335453 0.5614909 0.5967537

causal_mediation(data = UPBdata, model = "rb", est.method = "simulation",
                 outcome = "UPB", exposure = 'attbin', exposure.type = "binary",
                 mediator = c("negaff"), covariates = c("age"), EMint = FALSE,
                 yreg = "logistic", mreg = "linear")[[1:2]]
# $`decomp3way`
# cde_rr  cde_rr_se    pnde_rr pnde_rr_se    tnde_rr tnde_rr_se    pnie_rr pnie_rr_se    tnie_rr tnie_rr_se
# 1.3327660  0.1962687  1.3551020  0.2101584  1.3065619  0.1804419  1.2835601  0.0802581  1.2401888  0.0704087
# te_rr   te_rr_se         pm      pm_se
# 1.6773498  0.2570016  0.5234930  0.2153474
#
# $decomp4way
# cde_err         cde_err_se         intref_err      intref_err_se         intmed_err
# 0.37643356         0.22335135        -0.02133152         0.01600350         0.03868767
# intmed_err_se            pie_err         pie_err_se          total_err       total_err_se
# 0.02782670         0.28356005         0.08025810         0.67734976         0.25700165
# cde_err_prop    cde_err_prop_se    intref_err_prop intref_err_prop_se    intmed_err_prop
# 0.50393122         0.22615484        -0.02742418         0.01520308         0.04969281
# intmed_err_prop_se       pie_err_prop    pie_err_prop_se         overall_pm      overall_pm_se
# 0.02528033         0.47380016         0.23480283         0.52349296         0.21534742
# overall_int     overall_int_se         overall_pe      overall_pe_se
# 0.02226863         0.01227144         0.49606878         0.22615484


