context("svymultinom works correctly")

test_that("svymultinom works correctly", {
  
  set.seed(1)
  n <- 1000
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2 <- rnorm(n, mean = 1, sd = 1)
  x3 <- rbinom(n, size = 1, prob = 0.4)
  linearpred1 <- 1 + 0.3 * x1 - 0.5 * x2 - 0.2 * x3
  linearpred2 <- 2 + 1 * x1 - 2 * x2 - 1 * x3
  py2 <- exp(linearpred1) / (1 + exp(linearpred1) + exp(linearpred2))
  py3 <- exp(linearpred2) / (1 + exp(linearpred1) + exp(linearpred2))
  py1 <- 1 - py2 - py3
  y <- sapply(1:n, function(x) sample(size = 1, c(1:3), prob = c(py1[x], py2[x], py3[x])))
  w <- ifelse(x3 == 0, 0.4, 0.6)
  data <- data.frame(x1 = x1, x2 = x2, x3 = x3, y = y)
  
  naive_reg <- nnet::multinom(y ~ x1 + x2 + x3, weights = w, data = data, trace = FALSE)
  svydes <- survey::as.svrepdesign(survey::svydesign(~1, weights = ~w, data = data), type = "JK1")
  vcov_svy <- vcov(survey::withReplicates(svydes, function(w, data) {
    environment(formula) <- environment()
    as.vector(t(coef(eval(bquote(nnet::multinom(
      formula = y ~ x1 + x2 + x3, weights = w, data = data, trace = FALSE))))))
  }))
  attributes(vcov_svy)$means <- NULL
  
  reg <- svymultinom(y ~ x1 + x2 + x3, weights = w, data = data)
  reg_coef <- coef(reg)
  reg_vcov <- vcov(reg)
  dimnames(reg_vcov) <- list(NULL, NULL)
  reg_formula <- formula(reg)
  reg_pred <- predict(reg, newdata = data[1, ], type = "probs")
  reg_model <-  model.frame(reg)
  reg_summ <- summary(reg)
  reg_update <- update(reg, weights = w, data = data)
  
  # test
  expect_equal(reg_coef, coef(naive_reg))
  expect_equal(reg_vcov, vcov_svy)
  expect_equal(reg_formula, as.formula(y ~ x1 + x2 + x3))
  expect_equal(reg_pred, predict(naive_reg, newdata = data[1, ], type = "probs"))
  expect_equal(reg_model, model.frame(naive_reg))
  expect_equal(reg_summ$summarydf[, 1], as.vector(t(reg_coef)))
  expect_equal(reg_summ$summarydf[, 2], sqrt(diag(reg_vcov)))
  expect_equal(reg, reg_update)
  expect_equal(length(print(reg)), 64)
  expect_equal(class(print(summary(reg))), "data.frame")
  
})

