## -----------------------------------------------------------------------------
set.seed(1)
expit <- function(x) exp(x)/(1+exp(x))
n <- 10000
C1 <- rnorm(n, mean = 1, sd = 0.1)
C2 <- rbinom(n, 1, 0.6)
A <- rbinom(n, 1, expit(0.2 + 0.5*C1 + 0.1*C2))
M1 <- rpois(n, exp(1 - 2*A + 0.5*C1 + 0.8*C2))
linpred1 <- 0.1 + 0.1*A + 0.4*M1 - 0.5*C1 + 0.1*C2
linpred2 <- 0.4 + 0.2*A - 0.1*M1 - C1 + 0.5*C2
probm0 = 1 / (1 + exp(linpred1) + exp(linpred2))
probm1 = exp(linpred1) / (1 + exp(linpred1) + exp(linpred2))
probm2 = exp(linpred2) / (1 + exp(linpred1) + exp(linpred2))
M2 = factor(sapply(1:n, FUN = function(x) sample(c(0, 1, 2), size = 1, replace = TRUE,
                                                 prob=c(probm0[x],
                                                        probm1[x],
                                                        probm2[x]))))
Y <- rbinom(n, 1, expit(1 + 0.8*A - 1.8*M1 + 0.5*(M2 == 1) + 0.8*(M2 == 2) + 
                          0.5*A*M1 - 0.4*A*(M2 == 1) - 1.4*A*(M2 == 2)  + 0.3*C1 - 0.6*C2))
data <- data.frame(A, M1, M2, Y, C1, C2)

## -----------------------------------------------------------------------------
library(CMAverse)
cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
      basec = c("C1", "C2"), postc = NULL, node = TRUE, text_col = "white")

## ----message=F,warning=F,results='hide'---------------------------------------
res_rb <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                mreg = list("poisson", "multinomial"), yreg = "logistic",
                astar = 0, a = 1, mval = list(0, 2), yref = 1,
                estimation = "imputation", inference = "bootstrap", nboot = 2)

## ----message=F,warning=F------------------------------------------------------
summary(res_rb)

## ----message=F,warning=F,results='hide'---------------------------------------
res_wb <- cmest(data = data, model = "wb", outcome = "Y", exposure = "A",
                mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                ereg = "logistic", yreg = "logistic",
                astar = 0, a = 1, mval = list(0, 2), yref = 1,
                estimation = "imputation", inference = "bootstrap", nboot = 2)

## ----message=F,warning=F------------------------------------------------------
summary(res_wb)

## ----message=F,warning=F,results='hide'---------------------------------------
res_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                  mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                  ereg = "logistic", yreg = "logistic",
                  astar = 0, a = 1, mval = list(0, 2), yref = 1,
                  estimation = "imputation", inference = "bootstrap", nboot = 2)

## ----message=F,warning=F------------------------------------------------------
summary(res_iorw)

## ----message=F,warning=F,results='hide'---------------------------------------
res_ne <- cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                yreg = "logistic",
                astar = 0, a = 1, mval = list(0, 2), yref = 1,
                estimation = "imputation", inference = "bootstrap", nboot = 2)

## ----message=F,warning=F------------------------------------------------------
summary(res_ne)

## ----message=F,warning=F,results='hide'---------------------------------------
res_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                      mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                      mreg = list("poisson", "multinomial"), yreg = "logistic",
                      astar = 0, a = 1, mval = list(0, 2), yref = 1,
                      estimation = "imputation", inference = "bootstrap", nboot = 2)

## ----message=F,warning=F------------------------------------------------------
summary(res_gformula)

