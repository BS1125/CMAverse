rm(list=ls())

library("medflex")
data("UPBdata")
impFit=glm(UPB~factor(attbin)*negaff+age,family = binomial(),data=UPBdata)
expData=neImpute(impFit)
model1=neModel(UPB~attbin0*attbin1+age,family = binomial(),expData=expData,se="robust")

summary(neEffdecomp(model1))$coefficients[,"Estimate"]
causal_mediation(data = UPBdata, model = "sem", sem.method="delta",
                 outcome = "UPB", exposure = 'attbin', exposure.type = "binary",
                 mediator = 'negaff', covariates = c("age"), EMint = TRUE,
                 yreg = "logistic", mreg = "linear")$"decomp3way"

data1=read.csv("./R/sim_data_lung.csv")
data2=read.csv("./R/sim_data_bangladesh.csv")

delta1 <- causal_mediation(data1, outcome = "case", treatment = "snp", mediator = "smoking",
                 covariates = c("sex", "colgrad", "age"),interaction = TRUE, mreg = "logistic",
                 yreg = "logistic", method = "delta",vecc=c(1,1,1),
                 event = NULL, m_star = 1, a_star = 0, a = 1)

bootstrap1 <- causal_mediation(data1, outcome = "case", treatment = "snp", mediator = "smoking",
                           covariates = c("sex", "colgrad", "age"),interaction = TRUE, mreg = "logistic",
                           yreg = "logistic",method = "bootstrap", nboot = 200, vecc=c(1,1,1),
                           event = NULL, m_star = 1, a_star = 0, a = 1)

simulation1 <- causal_mediation(data1, outcome = "case", treatment = "snp", mediator = "smoking",
                               covariates = c("sex", "colgrad", "age"),interaction = TRUE, mreg = "logistic",
                               yreg = "logistic",method = "simulation", nsims=100, vecc=c(1,1,1),
                               event = NULL, m_star = 1, a_star = 0, a = 1)

delta1$decomp3way
bootstrap1$decomp3way
simulation1$decomp3way

delta1$decomp4way
bootstrap1$decomp4way
simulation1$decomp4way

delta2 <- causal_mediation(data2, outcome = "cognitive_raw", treatment = "ln_mn_c", mediator = "protein_c",
                           covariates = c("female", "approxage", "birthlength_c"),interaction = TRUE, mreg = "linear",
                           yreg = "linear",method = "delta",vecc=c(1,1,1),
                           event = NULL, m_star = 1, a_star = 0, a = 1)

bootstrap2 <- causal_mediation(data2, outcome = "cognitive_raw", treatment = "ln_mn_c", mediator = "protein_c",
                               covariates = c("female", "approxage", "birthlength_c"),interaction = TRUE, mreg = "linear",
                               yreg = "linear",method = "bootstrap", nboot = 100, vecc=c(1,1,1),
                               event = NULL, m_star = 1, a_star = 0, a = 1)

simulation2 <- causal_mediation(data2, outcome = "cognitive_raw", treatment = "ln_mn_c", mediator = "protein_c",
                           covariates = c("female", "approxage", "birthlength_c"),interaction = TRUE, mreg = "linear",
                           yreg = "linear",method = "simulation", nsims=100, vecc=c(1,1,1),
                           event = NULL, m_star = 1, a_star = 0, a = 1)

delta2$decomp3way
bootstrap2$decomp3way
simulation2$decomp3way

delta2$decomp4way
bootstrap2$decomp4way
simulation2$decomp4way
