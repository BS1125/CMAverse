test_that("cmest works correctly for regressions with prior weights", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 1000000
  C1 <- rnorm(n, mean = 1, sd = 0.1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  M1 <- rbinom(n, 1, expit(0.5 + 0.2*A + 0.5*C1 - 0.8*C2))
  sum(M1)
  M2 <- rbinom(n, 1, expit(-1 + A + M1 - 0.5*C1 + 0.5*C2))
  sum(M2)
  py <- expit(-6 + 0.8*A - 1.2*M1 + 0.5*M2 + 0.3*A*M1 + 0.2*A*M2 + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  yprevalence <- sum(Y)/n
  data <- data.frame(A, M1, M2, Y, C1, C2)
  case_indice <- sample(which(data$Y == 1), 2000, replace = FALSE)
  control_indice <- sample(which(data$Y == 0), 2000, replace = FALSE)
  data <- data[c(case_indice, control_indice), ]
  mreg1 <- glm(M1 ~ A + C1 + C2, data = data, family = binomial, weights = rep(2, 4000))
  mreg2 <- glm(M2 ~ A + C1 + C2, data = data, family = binomial, weights = rep(2, 4000))
  mreg1_msm <- glm(M1 ~ A, data = data, family = binomial, weights = rep(2, 4000))
  mreg2_msm <- glm(M2 ~ A, data = data, family = binomial, weights = rep(2, 4000))
  
  # yrare = TRUE
  # results of cmest
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb", casecontrol = TRUE, yrare = TRUE,
                                        outcome = "Y", exposure = "A",
                                        mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list(mreg1, mreg2), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1, 1),
                                        estimation = "imputation", inference = "bootstrap")
  res_binbin_wb <- cmest(data = data, model = "wb", casecontrol = TRUE, yrare = TRUE,
                         outcome = "Y", exposure = "A",
                         mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1, 1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_iorw <- cmest(data = data, model = "iorw", casecontrol = TRUE, yrare = TRUE,
                           outcome = "Y", exposure = "A",
                           mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1, 1),
                           estimation = "imputation", inference = "bootstrap")
  res_binbin_msm <- cmest(data = data, model = "msm", casecontrol = TRUE, yrare = TRUE,
                          outcome = "Y", exposure = "A",
                          mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                          ereg = "logistic", yreg = "logistic", mreg = list(mreg1_msm, mreg2_msm),
                          wmnomreg = list("logistic", "logistic"), wmdenomreg = list("logistic", "logistic"),
                          astar = 0, a = 1, mval = list(1, 1),
                          estimation = "imputation", inference = "bootstrap")
  res_binbin_ne <- cmest(data = data, model = "ne", casecontrol = TRUE, yrare = TRUE,
                         outcome = "Y", exposure = "A",
                         mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1, 1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_gformula <- cmest(data = data, model = "gformula", casecontrol = TRUE, yrare = TRUE,
                               outcome = "Y", exposure = "A",
                               mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list(mreg1, mreg2), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1, 1),
                               estimation = "imputation", inference = "bootstrap")
  
  # test
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_wb$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe)[c(6,2,5,15)], 
               unname(res_binbin_iorw$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_ne$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_msm$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_gformula$effect.pe), tolerance = 0.1)
  
  # yprevalence = yprevalence
  # results of cmest
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb",
                                        casecontrol = TRUE, yprevalence = yprevalence,
                                        outcome = "Y", exposure = "A",
                                        mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list(mreg1, mreg2), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1, 1),
                                        estimation = "imputation", inference = "bootstrap")
  res_binbin_wb <- cmest(data = data, model = "wb", casecontrol = TRUE, yprevalence = yprevalence,
                         outcome = "Y", exposure = "A",
                         mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1, 1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_iorw <- cmest(data = data, model = "iorw", casecontrol = TRUE, yprevalence = yprevalence,
                           outcome = "Y", exposure = "A",
                           mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1, 1),
                           estimation = "imputation", inference = "bootstrap")
  res_binbin_msm <- cmest(data = data, model = "msm", casecontrol = TRUE, yprevalence = yprevalence,
                          outcome = "Y", exposure = "A",
                          mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                          ereg = "logistic", yreg = "logistic", mreg = list(mreg1_msm, mreg2_msm),
                          wmnomreg = list("logistic", "logistic"), wmdenomreg = list("logistic", "logistic"),
                          astar = 0, a = 1, mval = list(1, 1),
                          estimation = "imputation", inference = "bootstrap")
  res_binbin_ne <- cmest(data = data, model = "ne", casecontrol = TRUE, yprevalence = yprevalence,
                         outcome = "Y", exposure = "A",
                         mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1, 1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_gformula <- cmest(data = data, model = "gformula", casecontrol = TRUE, yprevalence = yprevalence,
                               outcome = "Y", exposure = "A",
                               mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list(mreg1, mreg2), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1, 1),
                               estimation = "imputation", inference = "bootstrap")
  
  # test
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_wb$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe)[c(6,2,5,15)], 
               unname(res_binbin_iorw$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_ne$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_msm$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_gformula$effect.pe), tolerance = 0.1)
  
})