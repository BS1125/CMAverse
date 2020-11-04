context("cmest estimates causal effects correctly")

test_that("cmest works correctly for binary Y and binary M", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 100000
  C1 <- rnorm(n, mean = 1, sd = 0.1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M <- rbinom(n, 1, pm)
  py <- expit(-5 + 0.8*A - 0.6*M - 0.5*A*M + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  data <- data.frame(A, M, Y, C1, C2)
  yreg <- glm(Y ~ A*M + C1 + C2, family = binomial(), data = data)
  mreg <- glm(M ~ A + C1 + C2, family = binomial(), data = data)
  
  # results of cmest
  res_binbin_rb_param_delta <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                     mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                     mreg = list("logistic"), yreg = "logistic",
                                     astar = 0, a = 1, mval = list(1),
                                     estimation = "para", inference = "delt")
  res_binbin_rb_param_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                         mreg = list("logistic"), yreg = "logistic",
                                         astar = 0, a = 1, mval = list(1),
                                         estimation = "paramfunc", inference = "bootstrap")
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                        mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list("logistic"), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1),
                                        estimation = "impu", inference = "boot")
  res_binbin_wb <- cmest(data = data, model = "wb", outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                           mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1),
                           estimation = "imputation", inference = "bootstrap")
  res_binbin_ne <- cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap")
  
  # reference results
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
  expect_equal(summary(res_binbin_rb_param_delta)$summarydf$Estimate, ref)
  expect_equal(class(ggcmest(res_binbin_rb_param_delta)), c("gg", "ggplot"))
  expect_equal(summary(res_binbin_rb_param_bootstrap)$summarydf$Estimate, ref)
  expect_equal(summary(res_binbin_rb_impu_bootstrap)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_wb)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_ne)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_iorw)$summarydf$Estimate, ref[c(6,2,5,15)], tolerance = 0.1)
  expect_equal(summary(res_binbin_gformula)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(class(print(res_binbin_rb_param_delta)), "list")
  expect_equal(class(print(res_binbin_rb_param_bootstrap)), "list")
  expect_equal(class(print(res_binbin_rb_impu_bootstrap)), "list")
  expect_equal(class(print(res_binbin_wb)), "list")
  expect_equal(class(print(res_binbin_ne)), "list")
  expect_equal(class(print(res_binbin_iorw)), "list")
  expect_equal(class(print(res_binbin_gformula)), "list")
  expect_equal(class(print(summary(res_binbin_rb_param_delta))), "list")
  expect_equal(class(print(summary(res_binbin_rb_param_bootstrap))), "list")
  expect_equal(class(print(summary(res_binbin_rb_impu_bootstrap))), "list")
  expect_equal(class(print(summary(res_binbin_wb))), "list")
  expect_equal(class(print(summary(res_binbin_ne))), "list")
  expect_equal(class(print(summary(res_binbin_iorw))), "list")
  expect_equal(class(print(summary(res_binbin_gformula))), "list")
  
})


test_that("cmest works correctly for binary Y and binary M in a case control study", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 1000000
  C1 <- rnorm(n, mean = 1, sd = 0.1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  M <- rbinom(n, 1, expit(1 + 2*A + 1.5*C1 + 0.8*C2))
  py <- expit(-5 + 0.8*A - 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  yprevalence <- sum(Y)/n
  data <- data.frame(A, M, Y, C1, C2)
  case_indice <- sample(which(data$Y == 1), 2000, replace = FALSE)
  control_indice <- sample(which(data$Y == 0), 2000, replace = FALSE)
  data <- data[c(case_indice, control_indice), ]
  
  # yrare = TRUE
  # results of cmest
  res_binbin_rb_param_delta <- cmest(data = data, model = "rb", casecontrol = TRUE, yrare = TRUE,
                                     outcome = "Y", exposure = "A",
                                     mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                     mreg = list("logistic"), yreg = "loglinear",
                                     astar = 0, a = 1, mval = list(1),
                                     estimation = "paramfunc", inference = "delta")
  res_binbin_rb_param_bootstrap <- cmest(data = data, model = "rb", casecontrol = TRUE, yrare = TRUE,
                                         outcome = "Y", exposure = "A",
                                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                         mreg = list("logistic"), yreg = "loglinear",
                                         astar = 0, a = 1, mval = list(1),
                                         estimation = "paramfunc", inference = "bootstrap")
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb", casecontrol = TRUE, yrare = TRUE,
                                        outcome = "Y", exposure = "A",
                                        mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list("logistic"), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1),
                                        estimation = "imputation", inference = "bootstrap")
  res_binbin_wb <- cmest(data = data, model = "wb", casecontrol = TRUE, yrare = TRUE,
                         outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_iorw <- cmest(data = data, model = "iorw", casecontrol = TRUE, yrare = TRUE,
                           outcome = "Y", exposure = "A",
                           mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1),
                           estimation = "imputation", inference = "bootstrap")
  res_binbin_ne <- cmest(data = data, model = "ne", casecontrol = TRUE, yrare = TRUE,
                         outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_gformula <- cmest(data = data, model = "gformula", casecontrol = TRUE, yrare = TRUE,
                               outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap")
  
  # reference results
  thetas <- unname(coef(res_binbin_rb_param_delta$reg.output$yreg))
  betas <- unname(coef(res_binbin_rb_param_delta$reg.output$mreg[[1]]))
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
  meanc1 <- mean(data$C1)
  meanc2 <- mean(data$C2)
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
  expect_equal(unname(res_binbin_rb_param_delta$effect.pe), ref)
  expect_equal(unname(res_binbin_rb_param_bootstrap$effect.pe), ref)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_wb$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_ne$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_iorw$effect.pe), ref[c(6,2,5,15)], tolerance = 0.1)
  expect_equal(unname(res_binbin_gformula$effect.pe), ref, tolerance = 0.1)
  
  # yprevalence = yprevalence
  # results of cmest
  res_binbin_rb_param_delta <- cmest(data = data, model = "rb",
                                     casecontrol = TRUE, yprevalence = yprevalence,
                                     outcome = "Y", exposure = "A",
                                     mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                     mreg = list("logistic"), yreg = "loglinear",
                                     astar = 0, a = 1, mval = list(1),
                                     estimation = "paramfunc", inference = "delta")
  res_binbin_rb_param_delta_multinom <- cmest(data = data, model = "rb",
                                              casecontrol = TRUE, yprevalence = yprevalence,
                                              outcome = "Y", exposure = "A",
                                              mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                              mreg = list("multinomial"), yreg = "loglinear",
                                              astar = 0, a = 1, mval = list(1),
                                              estimation = "paramfunc", inference = "delta", 
                                              multimp = TRUE, m = 1)
  res_binbin_rb_param_bootstrap <- cmest(data = data, model = "rb",
                                         casecontrol = TRUE, yprevalence = yprevalence,
                                         outcome = "Y", exposure = "A",
                                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                         mreg = list("logistic"), yreg = "loglinear",
                                         astar = 0, a = 1, mval = list(1),
                                         estimation = "paramfunc", inference = "bootstrap")
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb",
                                        casecontrol = TRUE, yprevalence = yprevalence,
                                        outcome = "Y", exposure = "A",
                                        mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list("logistic"), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1),
                                        estimation = "imputation", inference = "bootstrap")
  res_binbin_wb <- cmest(data = data, model = "wb", casecontrol = TRUE, yprevalence = yprevalence,
                         outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_iorw <- cmest(data = data, model = "iorw", casecontrol = TRUE, yprevalence = yprevalence,
                           outcome = "Y", exposure = "A",
                           mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1),
                           estimation = "imputation", inference = "bootstrap")
  res_binbin_ne <- cmest(data = data, model = "ne", casecontrol = TRUE, yprevalence = yprevalence,
                         outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_gformula <- cmest(data = data, model = "gformula", casecontrol = TRUE, yprevalence = yprevalence,
                               outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap")
  
  # reference results
  thetas <- unname(coef(res_binbin_rb_param_delta$reg.output$yreg))
  betas <- unname(coef(res_binbin_rb_param_delta$reg.output$mreg[[1]]))
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
  meanc1 <- mean(data$C1)
  meanc2 <- mean(data$C2)
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
  expect_equal(unname(res_binbin_rb_param_delta$effect.pe), ref)
  expect_equal(unname(res_binbin_rb_param_delta_multinom$effect.pe), ref, tolerance = 0.001)
  expect_equal(class(print(res_binbin_rb_param_delta)), "list")
  expect_equal(class(print(res_binbin_rb_param_delta_multinom)), "list")
  expect_equal(unname(res_binbin_rb_param_bootstrap$effect.pe), ref)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_wb$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_ne$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_iorw$effect.pe), ref[c(6,2,5,15)], tolerance = 0.1)
  expect_equal(unname(res_binbin_gformula$effect.pe), ref, tolerance = 0.1)
  expect_equal(class(print(summary(res_binbin_rb_param_delta_multinom))), "list")
  
})


test_that("multiple imputation works correctly for binary Y and binary M ", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 0.1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M <- rbinom(n, 1, pm)
  py <- expit(-3 + 0.8*A - 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  data <- data.frame(A, M, Y, C1, C2)
  yreg <- glm(Y ~ A*M + C1 + C2, family = binomial(), data = data)
  mreg <- glm(M ~ A + C1 + C2, family = binomial(), data = data)
  
  # results of cmest
  res_binbin_rb_param_delta <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                     mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                     mreg = list("logistic"), yreg = "logistic",
                                     astar = 0, a = 1, mval = list(1),
                                     estimation = "paramfunc", inference = "delta", multimp = TRUE)
  res_binbin_rb_param_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                         mreg = list("logistic"), yreg = "logistic",
                                         astar = 0, a = 1, mval = list(1),
                                         estimation = "paramfunc", inference = "bootstrap", multimp = TRUE)
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                        mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list("logistic"), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1),
                                        estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  res_binbin_wb <- cmest(data = data, model = "wb", outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  res_binbin_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                           mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1),
                           estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  res_binbin_ne <- cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  res_binbin_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  
  # reference results
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
  expect_equal(summary(res_binbin_rb_param_delta)$summarydf$Estimate, ref)
  expect_equal(class(ggcmest(res_binbin_rb_param_delta)), c("gg", "ggplot"))
  expect_equal(summary(res_binbin_rb_param_bootstrap)$summarydf$Estimate, ref)
  expect_equal(summary(res_binbin_rb_impu_bootstrap)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_wb)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_ne)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_iorw)$summarydf$Estimate, ref[c(6,2,5,15)], tolerance = 0.1)
  expect_equal(summary(res_binbin_gformula)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(class(print(res_binbin_rb_param_delta)), "list")
  expect_equal(class(print(res_binbin_rb_param_bootstrap)), "list")
  expect_equal(class(print(res_binbin_rb_impu_bootstrap)), "list")
  expect_equal(class(print(res_binbin_wb)), "list")
  expect_equal(class(print(res_binbin_ne)), "list")
  expect_equal(class(print(res_binbin_iorw)), "list")
  expect_equal(class(print(res_binbin_gformula)), "list")
  expect_equal(class(print(summary(res_binbin_rb_param_delta))), "list")
  expect_equal(class(print(summary(res_binbin_rb_param_bootstrap))), "list")
  expect_equal(class(print(summary(res_binbin_rb_impu_bootstrap))), "list")
  expect_equal(class(print(summary(res_binbin_wb))), "list")
  expect_equal(class(print(summary(res_binbin_ne))), "list")
  expect_equal(class(print(summary(res_binbin_iorw))), "list")
  expect_equal(class(print(summary(res_binbin_gformula))), "list")
  
  # results of cmest
  res_binbin_rb_param_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                         mreg = list("logistic"), yreg = "logistic",
                                         astar = 0, a = 1, mval = list(1),
                                         estimation = "paramfunc", inference = "bootstrap", 
                                         boot.ci.type = "bca", multimp = TRUE)
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                        mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list("logistic"), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1),
                                        estimation = "imputation", inference = "bootstrap", 
                                        boot.ci.type = "bca", multimp = TRUE)
  res_binbin_wb <- cmest(data = data, model = "wb", outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap", 
                         boot.ci.type = "bca", multimp = TRUE)
  res_binbin_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                           mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1),
                           estimation = "imputation", inference = "bootstrap", 
                           boot.ci.type = "bca", multimp = TRUE)
  res_binbin_ne <- cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap", 
                         boot.ci.type = "bca", multimp = TRUE)
  res_binbin_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap", 
                               boot.ci.type = "bca", multimp = TRUE)
  
  expect_equal(summary(res_binbin_rb_param_bootstrap)$summarydf$Estimate, ref)
  expect_equal(summary(res_binbin_rb_impu_bootstrap)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_wb)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_ne)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_iorw)$summarydf$Estimate, ref[c(6,2,5,15)], tolerance = 0.1)
  expect_equal(summary(res_binbin_gformula)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(class(print(res_binbin_rb_param_bootstrap)), "list")
  expect_equal(class(print(res_binbin_rb_impu_bootstrap)), "list")
  expect_equal(class(print(res_binbin_wb)), "list")
  expect_equal(class(print(res_binbin_ne)), "list")
  expect_equal(class(print(res_binbin_iorw)), "list")
  expect_equal(class(print(res_binbin_gformula)), "list")
  expect_equal(class(print(summary(res_binbin_rb_param_bootstrap))), "list")
  expect_equal(class(print(summary(res_binbin_rb_impu_bootstrap))), "list")
  expect_equal(class(print(summary(res_binbin_wb))), "list")
  expect_equal(class(print(summary(res_binbin_ne))), "list")
  expect_equal(class(print(summary(res_binbin_iorw))), "list")
  expect_equal(class(print(summary(res_binbin_gformula))), "list")
})

