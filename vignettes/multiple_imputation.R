## -----------------------------------------------------------------------------
set.seed(1)
expit <- function(x) exp(x)/(1+exp(x))
n <- 10000
C1 <- rnorm(n, mean = 1, sd = 0.1)
C2 <- rbinom(n, 1, 0.6)
A <- rbinom(n, 1, expit(0.2 + 0.5*C1 + 0.1*C2))
M <- rbinom(n, 1, expit(1 + 2*A + 1.5*C1 + 0.8*C2))
Y <- rbinom(n, 1, expit(-3 - 0.4*A - 1.2*M + 0.5*A*M + 0.3*C1 - 0.6*C2))
data_noNA <- data.frame(A, M, Y, C1, C2)
missing <- sample(1:(5*n), n*0.1, replace = FALSE)
C1[missing[which(missing <= n)]] <- NA
C2[missing[which((missing > n)*(missing <= 2*n) == 1)] - n] <- NA
A[missing[which((missing > 2*n)*(missing <= 3*n) == 1)] - 2*n] <- NA
M[missing[which((missing > 3*n)*(missing <= 4*n) == 1)] - 3*n] <- NA
Y[missing[which((missing > 4*n)*(missing <= 5*n) == 1)] - 4*n] <- NA
data <- data.frame(A, M, Y, C1, C2)

## -----------------------------------------------------------------------------
library(CMAverse)
cmdag(outcome = "Y", exposure = "A", mediator = "M",
      basec = c("C1", "C2"), postc = NULL, node = TRUE, text_col = "white")

