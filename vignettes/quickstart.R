## -----------------------------------------------------------------------------
library(CMAverse)

## ----message=F,warning=F------------------------------------------------------
set.seed(1)
n <- 100
C1 <- rnorm(n, mean = 1, sd = 1)
C2 <- rbinom(n, 1, 0.6)
pa <- exp(0.2 - 0.5*C1 + 0.1*C2)/(1 + exp(0.2 - 0.5*C1 + 0.1*C2))
A <- rbinom(n, 1, pa)
pm <- exp(1 + 0.5*A - 1.5*C1 + 0.5*C2)/ (1 + exp(1 + 0.5*A - 1.5*C1 + 0.5*C2))
M1 <- rbinom(n, 1, pm)
M2 <- rnorm(n, 2 + 0.8*A - M1 + 0.5*C1 + 2*C2, 1)
Y <- rnorm(n, mean = 0.5 + 0.4*A + 0.5*M1 + 0.6*M2 + 0.3*A*M1 + 0.2*A*M2 - 0.3*C1 + 2*C2, sd = 1)
data <- data.frame(A, M1, M2, Y, C1, C2)

## ----plot_dag,message=F,warning=F---------------------------------------------
cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"), 
      basec = c("C1", "C2"), postc = NULL, node = FALSE, text_col = "black")

## ----message=F,warning=F,results='hide'---------------------------------------
est <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                mreg = list("logistic", "linear"), yreg = "linear",
                astar = 0, a = 1, mval = list(1, 1),
                estimation = "imputation", inference = "bootstrap", nboot = 20)

## ----message=F,warning=F------------------------------------------------------
summary(est)

## ----plot_cmest,message=F,warning=F,fig.width=8,fig.height=5------------------
ggcmest(est) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.8))

## ----message=F,warning=F------------------------------------------------------
cmsens(object = est, sens = "uc")

## ----message=F,warning=F,results='hide'---------------------------------------
me1 <- cmsens(object = est, sens = "me", MEmethod = "rc", 
              MEvariable = "C1", MEvartype = "con", MEerror = c(0.1, 0.2, 0.3))

## ----message=F,warning=F------------------------------------------------------
summary(me1)

## ----plot_cmsens_me_con,message=F,warning=F,fig.width=8,fig.height=5----------
ggcmsens(me1) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.8))

## ----message=F,warning=F,results='hide'---------------------------------------
me2 <- cmsens(object = est, sens = "me", MEmethod = "simex", MEvariable = "A", 
              MEvartype = "cat", B = 5,
              MEerror = list(matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2), 
                             matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)))

## ----message=F,warning=F------------------------------------------------------
summary(me2)

## ----plot_cmsens_me_cat,message=F,warning=F,fig.width=8,fig.height=6----------
ggcmsens(me2) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.8))

