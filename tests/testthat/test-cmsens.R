context("cmsens corrects causal effects correctly")

test_that("sensitivity analysis for measurement error works correctly for binary Y and binary M", {

  # a continuous variable measured with error
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 1000
  C1 <- rnorm(n, mean = 1, sd = 1)
  error <- rnorm(n, sd = 0.05)
  C1_error <- C1 + error
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M <- rbinom(n, 1, pm)
  py <- expit(-3 + 0.8*A + 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  data <- data.frame(A, M, Y, C1, C1_error, C2)

  # naive results
  res_binbin_naive <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                      mediator = "M", basec = c("C1_error", "C2"), EMint = TRUE,
                      mreg = list("logistic"), yreg = "logistic",
                      astar = 0, a = 1, mval = list(1),
                      estimation = "paramfunc", inference = "delta")

  # cmsens-corrected results
  res_binbin_cmsens_rc <- cmsens(object = res_binbin_naive, sens = "me", MEmethod = "rc",
                              MEvariable = "C1_error", MEvartype = "con", MEerror = 0.1)
  res_binbin_cmsens_simex <- cmsens(object = res_binbin_naive, sens = "me", MEmethod = "simex",
                                 MEvariable = "C1_error", MEvartype = "con", MEerror = 0.1)

  # reference results
  yreg <- glm(Y ~ A*M + C1 + C2, family = binomial, data = data)
  mreg <- glm(M ~ A + C1 + C2, family = binomial, data = data)
  thetas <- unname(coef(yreg))
  betas <- unname(coef(mreg))
  beta0 <- betas[1]
  beta1 <- betas[2]
  beta2 <- betas[3]
  beta3 <- betas[4]
  theta0 <- thetas[1]
  theta1 <- thetas[2]
  theta2 <- thetas[3]
  theta3 <- thetas[6]
  theta4 <- thetas[4]
  theta5 <- thetas[5]
  m <- 1
  a <- 1
  astar <- 0
  meanc1 <- mean(C1)
  meanc2 <- mean(C2)
  cde_binbin <- exp((theta1+theta3*m)*(a-astar))
  pnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))
  pnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  te_binbin <- pnde_binbin*tnie_binbin
  cde_err_binbin <- exp(theta2*m)*(exp(theta1*(a-astar)+theta3*a*m)-exp(theta3*astar*m))*
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))/
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2+theta2+theta3*astar))
  intref_err_binbin <- pnde_binbin - 1 - cde_err_binbin
  intmed_err_binbin <- tnie_binbin * pnde_binbin - pnde_binbin - pnie_binbin + 1
  pnie_err_binbin <- pnie_binbin - 1
  total_err_binbin <- te_binbin - 1
  cde_err_prop_binbin <- cde_err_binbin/total_err_binbin
  intmed_err_prop_binbin <- intmed_err_binbin/total_err_binbin
  intref_err_prop_binbin <- intref_err_binbin/total_err_binbin
  pnie_err_prop_binbin <- pnie_err_binbin/total_err_binbin
  pm_binbin <- (pnie_err_binbin+intmed_err_binbin)/total_err_binbin
  int_binbin <- (intref_err_binbin+intmed_err_binbin)/total_err_binbin
  pe_binbin <- (intref_err_binbin+intmed_err_binbin+pnie_err_binbin)/total_err_binbin

  ref <- c(cde_binbin, pnde_binbin, tnde_binbin, pnie_binbin, tnie_binbin, te_binbin, cde_err_binbin,
           intref_err_binbin, intmed_err_binbin, pnie_err_binbin, cde_err_prop_binbin,
           intref_err_prop_binbin, intmed_err_prop_binbin, pnie_err_prop_binbin,
           pm_binbin, int_binbin, pe_binbin)
  # test
  expect_equal(summary(res_binbin_cmsens_rc)$summarydf[[1]]$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_cmsens_simex)$summarydf[[1]]$Estimate, ref, tolerance = 0.1)
  expect_equal(class(ggcmsens(res_binbin_cmsens_rc)), c("gg", "ggplot"))
  expect_equal(class(ggcmsens(res_binbin_cmsens_simex)), c("gg", "ggplot"))
  expect_equal(print(res_binbin_cmsens_rc), NULL)
  expect_equal(print(res_binbin_cmsens_simex), NULL)
  expect_equal(print(summary(res_binbin_cmsens_rc)), NULL)
  expect_equal(print(summary(res_binbin_cmsens_simex)), NULL)

  # a categorical variable measured with error
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 1000
  C1 <- rnorm(n, mean = 1, sd = 1)
  C2 <- rbinom(n, 1, 0.6)
  MEerror <- matrix(c(0.9,0.1,0.1,0.9), nrow = 2)
  C2_error <- C2
  for (j in 1:2) {
    C2_error[which(C2_error == c(0,1)[j])] <-
       sample(x = c(0,1), size = length(which(C2_error == c(0,1)[j])),
              prob = MEerror[, j], replace = TRUE)
   }
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M <- rbinom(n, 1, pm)
  py <- expit(-3 + 0.8*A + 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  data <- data.frame(A, M, Y, C1, C2, C2_error)

  # naive results
  res_binbin_naive <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                            mediator = "M", basec = c("C1", "C2_error"), EMint = TRUE,
                            mreg = list("logistic"), yreg = "logistic",
                            astar = 0, a = 1, mval = list(1),
                            estimation = "paramfunc", inference = "delta")

  # cmsens-corrected results
  res_binbin_cmsens_simex <- cmsens(object = res_binbin_naive, sens = "me", MEmethod = "simex",
                                    MEvariable = "C2_error", MEvartype = "cat", MEerror = list(MEerror))

  # reference results
  yreg <- glm(Y ~ A*M + C1 + C2, family = binomial, data = data)
  mreg <- glm(M ~ A + C1 + C2, family = binomial, data = data)
  thetas <- unname(coef(yreg))
  betas <- unname(coef(mreg))
  beta0 <- betas[1]
  beta1 <- betas[2]
  beta2 <- betas[3]
  beta3 <- betas[4]
  theta0 <- thetas[1]
  theta1 <- thetas[2]
  theta2 <- thetas[3]
  theta3 <- thetas[6]
  theta4 <- thetas[4]
  theta5 <- thetas[5]
  m <- 1
  a <- 1
  astar <- 0
  meanc1 <- mean(C1)
  meanc2 <- mean(C2)
  cde_binbin <- exp((theta1+theta3*m)*(a-astar))
  pnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))
  pnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  te_binbin <- pnde_binbin*tnie_binbin
  cde_err_binbin <- exp(theta2*m)*(exp(theta1*(a-astar)+theta3*a*m)-exp(theta3*astar*m))*
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))/
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2+theta2+theta3*astar))
  intref_err_binbin <- pnde_binbin - 1 - cde_err_binbin
  intmed_err_binbin <- tnie_binbin * pnde_binbin - pnde_binbin - pnie_binbin + 1
  pnie_err_binbin <- pnie_binbin - 1
  total_err_binbin <- te_binbin - 1
  cde_err_prop_binbin <- cde_err_binbin/total_err_binbin
  intmed_err_prop_binbin <- intmed_err_binbin/total_err_binbin
  intref_err_prop_binbin <- intref_err_binbin/total_err_binbin
  pnie_err_prop_binbin <- pnie_err_binbin/total_err_binbin
  pm_binbin <- (pnie_err_binbin+intmed_err_binbin)/total_err_binbin
  int_binbin <- (intref_err_binbin+intmed_err_binbin)/total_err_binbin
  pe_binbin <- (intref_err_binbin+intmed_err_binbin+pnie_err_binbin)/total_err_binbin

  ref <- c(cde_binbin, pnde_binbin, tnde_binbin, pnie_binbin, tnie_binbin, te_binbin, cde_err_binbin,
           intref_err_binbin, intmed_err_binbin, pnie_err_binbin, cde_err_prop_binbin,
           intref_err_prop_binbin, intmed_err_prop_binbin, pnie_err_prop_binbin,
           pm_binbin, int_binbin, pe_binbin)
  # test
  expect_equal(summary(res_binbin_cmsens_simex)$summarydf[[1]]$Estimate, ref, tolerance = 0.1)
  expect_equal(class(ggcmsens(res_binbin_cmsens_simex)), c("gg", "ggplot"))
  expect_equal(print(res_binbin_cmsens_simex), NULL)

})

test_that("sensitivity analysis for unmeasured confounding works correctly for binary Y and binary M", {

  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 1000
  C1 <- rnorm(n, mean = 1, sd = 1)
  error <- rnorm(n, sd = 0.05)
  C1_error <- C1 + error
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M <- rbinom(n, 1, pm)
  py <- expit(-3 + 0.8*A + 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  data <- data.frame(A, M, Y, C1, C1_error, C2)

  # naive results
  res_binbin_naive <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                            mediator = "M", basec = c("C1_error", "C2"), EMint = TRUE,
                            mreg = list("logistic"), yreg = "logistic",
                            astar = 0, a = 1, mval = list(1),
                            estimation = "paramfunc", inference = "delta")

  # cmsens results
  res_binbin_cmsens_uc <- cmsens(object = res_binbin_naive, sens = "uc")

  # evalue
  effect.pe <- res_binbin_naive$effect.pe[1:6]
  ci.lo <- res_binbin_naive$effect.ci.low[1:6]
  ci.up <- res_binbin_naive$effect.ci.high[1:6]
  evalues <- c()
  for (i in 1:6) {
    evalues <- rbind(evalues, evalues.RR(est = effect.pe[i], lo = ci.lo[i], hi = ci.up[i])[2, ])
  }

  # test
  expect_equal(as.vector(res_binbin_cmsens_uc$evalues[,4:6]), as.vector(evalues))

})

test_that("sensitivity analysis for unmeasured confounding works correctly for continuous Y and binary M", {

  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M <- rbinom(n, 1, pm)
  Y <- rnorm(n, -1 + 0.8*A + 0.5*M + 0.5*A*M + 0.3*C1 - 0.6*C2, 1)
  data <- data.frame(A, M, Y, C1, C2)

  # naive results
  res_contbin_naive <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                      mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                      mreg = list("logistic"), yreg = "linear",
                                      astar = 0, a = 1, mval = list(1),
                                      estimation = "paramfunc", inference = "delta")

  # cmsens results
  res_contbin_cmsens_uc <- cmsens(object = res_contbin_naive, sens = "uc")

  # evalue
  effect.pe <- res_contbin_naive$effect.pe[1:6]
  effect.se <- res_contbin_naive$effect.se[1:6]
  ci.lo <- res_contbin_naive$effect.ci.low[1:6]
  ci.up <- res_contbin_naive$effect.ci.high[1:6]
  evalues <- c()
  outcome.se <- sd(data$Y, na.rm = TRUE)
  d <- unname(effect.pe / outcome.se)
  sd <- unname(effect.se / outcome.se)
  for (i in 1:6) {
    evalues <- rbind(evalues, EValue::evalues.RR(est = exp(0.91 * d[i]), lo = exp(0.91 * d[i] - 1.78 * sd[i]),
                              hi = exp(0.91 * d[i] + 1.78 * sd[i]))[2, ])
  }

  # test
  expect_equal(as.vector(res_contbin_cmsens_uc$evalues[,4:6]), as.vector(evalues))

})


test_that("sensitivity analysis for measurement error works correctly for binary Y and binary M with multiple imputation", {
  
  # a continuous variable measured with error
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 1000
  C1 <- rnorm(n, mean = 1, sd = 1)
  error <- rnorm(n, sd = 0.05)
  C1_error <- C1 + error
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M <- rbinom(n, 1, pm)
  py <- expit(-3 + 0.8*A + 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  data_noNA <- data.frame(A, M, Y, C1, C1_error, C2)
  missing <- sample(1:(5*n), n*0.1, replace = FALSE)
  C1[missing[which(missing <= n)]] <- NA
  C1_error[missing[which(missing <= n)]] <- NA
  C2[missing[which((missing > n)*(missing <= 2*n) == 1)] - n] <- NA
  A[missing[which((missing > 2*n)*(missing <= 3*n) == 1)] - 2*n] <- NA
  M[missing[which((missing > 3*n)*(missing <= 4*n) == 1)] - 3*n] <- NA
  Y[missing[which((missing > 4*n)*(missing <= 5*n) == 1)] - 4*n] <- NA
  data <- data.frame(A, M, Y, C1, C1_error, C2)
  
  # naive results
  res_binbin_naive <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                            mediator = "M", basec = c("C1_error", "C2"), EMint = TRUE,
                            mreg = list("logistic"), yreg = "logistic",
                            astar = 0, a = 1, mval = list(1),
                            estimation = "paramfunc", inference = "delta", multimp = TRUE, m = 5)
  
  # cmsens-corrected results
  res_binbin_cmsens_rc <- cmsens(object = res_binbin_naive, sens = "me", MEmethod = "rc",
                                 MEvariable = "C1_error", MEvartype = "con", MEerror = 0.1)
  res_binbin_cmsens_simex <- cmsens(object = res_binbin_naive, sens = "me", MEmethod = "simex",
                                    MEvariable = "C1_error", MEvartype = "con", MEerror = 0.1)
  
  # reference results
  yreg <- glm(Y ~ A*M + C1 + C2, family = binomial, data = data_noNA)
  mreg <- glm(M ~ A + C1 + C2, family = binomial, data = data_noNA)
  thetas <- unname(coef(yreg))
  betas <- unname(coef(mreg))
  beta0 <- betas[1]
  beta1 <- betas[2]
  beta2 <- betas[3]
  beta3 <- betas[4]
  theta0 <- thetas[1]
  theta1 <- thetas[2]
  theta2 <- thetas[3]
  theta3 <- thetas[6]
  theta4 <- thetas[4]
  theta5 <- thetas[5]
  m <- 1
  a <- 1
  astar <- 0
  meanc1 <- mean(data_noNA$C1)
  meanc2 <- mean(data_noNA$C2)
  cde_binbin <- exp((theta1+theta3*m)*(a-astar))
  pnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))
  pnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  te_binbin <- pnde_binbin*tnie_binbin
  cde_err_binbin <- exp(theta2*m)*(exp(theta1*(a-astar)+theta3*a*m)-exp(theta3*astar*m))*
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))/
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2+theta2+theta3*astar))
  intref_err_binbin <- pnde_binbin - 1 - cde_err_binbin
  intmed_err_binbin <- tnie_binbin * pnde_binbin - pnde_binbin - pnie_binbin + 1
  pnie_err_binbin <- pnie_binbin - 1
  total_err_binbin <- te_binbin - 1
  cde_err_prop_binbin <- cde_err_binbin/total_err_binbin
  intmed_err_prop_binbin <- intmed_err_binbin/total_err_binbin
  intref_err_prop_binbin <- intref_err_binbin/total_err_binbin
  pnie_err_prop_binbin <- pnie_err_binbin/total_err_binbin
  pm_binbin <- (pnie_err_binbin+intmed_err_binbin)/total_err_binbin
  int_binbin <- (intref_err_binbin+intmed_err_binbin)/total_err_binbin
  pe_binbin <- (intref_err_binbin+intmed_err_binbin+pnie_err_binbin)/total_err_binbin
  
  ref <- c(cde_binbin, pnde_binbin, tnde_binbin, pnie_binbin, tnie_binbin, te_binbin, cde_err_binbin,
           intref_err_binbin, intmed_err_binbin, pnie_err_binbin, cde_err_prop_binbin,
           intref_err_prop_binbin, intmed_err_prop_binbin, pnie_err_prop_binbin,
           pm_binbin, int_binbin, pe_binbin)
  # test
  expect_equal(summary(res_binbin_cmsens_rc)$summarydf[[1]]$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_cmsens_simex)$summarydf[[1]]$Estimate, ref, tolerance = 0.1)
  expect_equal(class(ggcmsens(res_binbin_cmsens_rc)), c("gg", "ggplot"))
  expect_equal(class(ggcmsens(res_binbin_cmsens_simex)), c("gg", "ggplot"))
  expect_equal(print(res_binbin_cmsens_rc), NULL)
  expect_equal(print(res_binbin_cmsens_simex), NULL)
  expect_equal(print(summary(res_binbin_cmsens_rc)), NULL)
  expect_equal(print(summary(res_binbin_cmsens_simex)), NULL)
  
})