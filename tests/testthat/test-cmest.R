test_that("cmest works correctly for binary Y and binary M with postc", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 0.1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  L <- rnorm(n, mean = 0.5 - A + 0.2*C1 + 0.6*C2)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2 + L)
  M <- rbinom(n, 1, pm)
  py <- expit(-3 + 0.8*A - 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2 + L)
  Y <- rbinom(n, 1, py)
  data <- data.frame(A, M, Y, C1, C2, L)
  
  # results of cmest
  expect_error(cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     mreg = list("logistic"), yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "paramfunc", inference = "delta"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  expect_error(cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     mreg = list("logistic"), yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "paramfunc", inference = "delta"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  expect_error(cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     mreg = list("logistic"), yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "paramfunc", inference = "bootstrap"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  expect_error(cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     mreg = list("logistic"), yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "imputation", inference = "bootstrap"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  expect_error(cmest(data = data, model = "wb", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     ereg = "logistic", yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "imputation", inference = "bootstrap"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  expect_error(cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     ereg = "logistic", yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "imputation", inference = "bootstrap"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  res_binbin_msm <- cmest(data = data, model = "msm", outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                          ereg = "logistic", yreg = "logistic", mreg = list("logistic"),
                          wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap")
  expect_error(cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "imputation", inference = "bootstrap"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  res_binbin_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic", postcreg = list("linear"),
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap")
  
  # test
  expect_equal(summary(res_binbin_gformula)$summarydf$Estimate,
               summary(res_binbin_msm)$summarydf$Estimate, tolerance = 0.1)
  expect_equal(class(print(res_binbin_msm)), "list")
  expect_equal(class(print(res_binbin_gformula)), "list")
  expect_equal(class(print(summary(res_binbin_msm))), "list")
  expect_equal(class(print(summary(res_binbin_gformula))), "list")
  
  res_binbin_msm <- cmest(data = data, model = "msm", outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                          ereg = "logistic", yreg = "logistic", mreg = list("logistic"),
                          wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  res_binbin_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic", postcreg = list("linear"),
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  
  # test
  expect_equal(summary(res_binbin_gformula)$summarydf$Estimate,
               summary(res_binbin_msm)$summarydf$Estimate, tolerance = 0.1)
  expect_equal(class(print(res_binbin_msm)), "list")
  expect_equal(class(print(res_binbin_gformula)), "list")
  expect_equal(class(print(summary(res_binbin_msm))), "list")
  expect_equal(class(print(summary(res_binbin_gformula))), "list")
  
})