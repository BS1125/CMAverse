rm(list=ls())

n <- 1000
ID <- 1:n
# simulate pre-exposure confounders
C1 = rnorm(n, mean = 1, sd = 1)
C2 = sample(c("C2_0", "C2_1"), n, 1, c(0.6, 0.4))
# simulate a binary exposure
linpred <- 0.2 - 0.5*C1 + 0.1*(C2=="C2_1")
pa <- exp(linpred)/(1 + exp(linpred))
A <- rbinom(n, 1, pa)
# simulate a binary mediator
pm <- exp(linpred)/(1 + exp(linpred))
M1 <- rbinom(n, 1, pm)
# simulate a categorical mediator
linpred1 <- 0.1 + 0.1*A - 0.5*C1 + 0.1*(C2 == "C2_1")
linpred2 <- 0.4 + 0.2*A - C1 + 0.5*(C2 == "C2_1")
probm0 = 1 / (1 + exp(linpred1) + exp(linpred2))
probm1 = exp(linpred1) / (1 + exp(linpred1) + exp(linpred2))
probm2 = exp(linpred2) / (1 + exp(linpred1) + exp(linpred2))
M2 = sapply(1:n, FUN = function(x) sample(c("M2_0","M2_1","M2_2"), size = 1, replace = TRUE,
                                         prob=c(probm0[x],
                                                probm1[x],
                                                probm2[x])))
# simulate a continuous outcome
linpred <- 0.5 + 0.5*A + 0.2*M1 + 0.5*(M2 == "M2_1") - 0.4*(M2 == "M2_2") - 0.3*C1 - 2*(C2=="C2_1")
contY <- rnorm(n, mean = linpred, sd = 1)
# simulate a binary outcome
linpred <- -5 + 0.5*A + 0.2*M1 + 0.5*(M2 == "M2_1") - 0.4*(M2 == "M2_2") - 0.3*C1 - 2*(C2=="C2_1")
py <- exp(linpred)/(1 + exp(linpred))
binY <- rbinom(n, 1, py)

cma2020 <- data.frame(ID, C1, C2, A, M1, M2, contY, binY)

usethis::use_data(cma2020, internal = TRUE, overwrite = TRUE)

