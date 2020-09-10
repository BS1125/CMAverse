context("rcreg works correctly")

test_that("rcreg works correctly for lm", {
  
  set.seed(1)
  n <- 10000
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2_true <- rnorm(n, mean = 1, sd = 1)
  error1 <- rnorm(n, mean = 0, sd = 0.5)
  x2_error <- x2_true + error1
  x3 <- rbinom(n, size = 1, prob = 0.4)
  y <- 1 + 2 * x1 + 4 * x2_true + 2 * x3  + rnorm(n, mean = 0, sd = 2)
  data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
                     x3 = x3, y = y)
  reg_naive <- lm(y ~ x1 + x2_error + x3, data = data)
  reg_true <- lm(y ~ x1 + x2_true + x3, data = data)
  reg_rc <- rcreg(reg = reg_naive, data = data, MEvariable = "x2_error",
                  MEerror = 0.5, variance = TRUE, nboot = 2)
  
  # test
  expect_equal(unname(coef(reg_rc)), unname(coef(reg_true)), tolerance = 0.1)
  expect_equal(formula(reg_rc), as.formula(y ~ x1 + x2_error + x3))
  expect_equal(sigma(reg_rc), 2, tolerance = 0.1)
  expect_equal(family(reg_rc), family(reg_true))
  expect_equal(as.numeric(predict(reg_rc, newdata = data[1, ], type = "response")), 
               as.numeric(unname(coef(reg_rc)) %*% c(1, as.numeric(data[1, c(1, 3, 4)]))))
  expect_equal(model.frame(reg_rc), model.frame(reg_naive))
  expect_equal(update(reg_rc, data = data)$RCcoef, coef(reg_rc))
  expect_equal(summary(reg_rc)$summarydf$Estimate, as.numeric(coef(reg_rc)))
  expect_equal(class(print(reg_rc)), "numeric")
  expect_equal(class(print(summary(reg_rc))), "data.frame")
})

test_that("rcreg works correctly for glm", {
  
  set.seed(1)
  n <- 10000
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2_true <- rnorm(n, mean = 1, sd = 1)
  error1 <- rnorm(n, mean = 0, sd = 0.5)
  x2_error <- x2_true + error1
  x3 <- rbinom(n, size = 1, prob = 0.4)
  linearpred <- 1 + 0.3 * x1 - 0.5 * x2_true - 0.2 * x3
  py <- exp(linearpred) / (1 + exp(linearpred))
  y <- rbinom(n, size = 1, prob = py)
  data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
                     x3 = x3, y = y)
  reg_naive <- glm(y ~ x1 + x2_error + x3, data = data, family = binomial("logit"))
  reg_true <- glm(y ~ x1 + x2_true + x3, data = data, family = binomial("logit"))
  reg_rc <- rcreg(reg = reg_naive, data = data, MEvariable = "x2_error",
                  MEerror = 0.5, variance = TRUE, nboot = 2)

  # test
  expect_equal(unname(coef(reg_rc)), unname(coef(reg_true)), tolerance = 0.1)
  expect_equal(formula(reg_rc), as.formula(y ~ x1 + x2_error + x3))
  expect_equal(family(reg_rc), family(reg_true))
  expect_equal(as.numeric(predict(reg_rc, newdata = data[1, ], type = "response")), 
               as.numeric(exp(t(unname(coef(reg_rc))) %*% c(1, as.numeric(data[1, c(1, 3, 4)])))/
                 (1+exp(t(unname(coef(reg_rc))) %*% c(1, as.numeric(data[1, c(1, 3, 4)]))))))
  expect_equal(model.frame(reg_rc), model.frame(reg_naive))
  expect_equal(update(reg_rc, data = data)$RCcoef, coef(reg_rc))
  expect_equal(summary(reg_rc)$summarydf$Estimate, as.numeric(coef(reg_rc)))
  
})

test_that("rcreg works correctly for multinom", {
  
  set.seed(1)
  n <- 10000
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2_true <- rnorm(n, mean = 1, sd = 1)
  error1 <- rnorm(n, mean = 0, sd = 0.5)
  x2_error <- x2_true + error1
  x3 <- rbinom(n, size = 1, prob = 0.4)
  linearpred1 <- 1 + 0.3 * x1 - 0.5 * x2_true - 0.2 * x3
  linearpred2 <- 2 + 1 * x1 - 2 * x2_true - 1 * x3
  py2 <- exp(linearpred1) / (1 + exp(linearpred1) + exp(linearpred2))
  py3 <- exp(linearpred2) / (1 + exp(linearpred1) + exp(linearpred2))
  py1 <- 1 - py2 - py3
  y <- sapply(1:n, function(x) sample(size = 1, c(1:3), prob = c(py1[x], py2[x], py3[x])))
  data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
                     x3 = x3, y = y)
  reg_naive <- nnet::multinom(factor(y) ~ x1 + x2_error + x3, data = data)
  reg_true <- nnet::multinom(factor(y) ~ x1 + x2_true + x3, data = data)
  reg_rc <- rcreg(reg = reg_naive, data = data, MEvariable = "x2_error",
                  MEerror = 0.5, variance = TRUE, nboot = 2)
  
  # test
  expect_equal(unname(coef(reg_rc)), as.numeric(t(unname(coef(reg_true)))), tolerance = 0.1)
  expect_equal(formula(reg_rc), as.formula(factor(y) ~ x1 + x2_error + x3))
  expect_equal(as.numeric(predict(reg_rc, newdata = data[1, ], type = "probs")),
               c(py1[1], py2[1], py3[1]), tolerance = 0.1)
  expect_equal(model.frame(reg_rc), model.frame(reg_naive))
  expect_equal(update(reg_rc, data = data)$RCcoef, coef(reg_rc))
  expect_equal(summary(reg_rc)$summarydf$Estimate, as.numeric(coef(reg_rc)))
  
})

test_that("rcreg works correctly for polr", {
  
  set.seed(1)
  n <- 10000
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2_true <- rnorm(n, mean = 1, sd = 1)
  error1 <- rnorm(n, mean = 0, sd = 0.5)
  x2_error <- x2_true + error1
  x3 <- rbinom(n, size = 1, prob = 0.4)
  linearpred1 <- 1 + 0.3 * x1 - 0.5 * x2_true - 0.2 * x3
  linearpred2 <- 2 + 0.3 * x1 - 0.5 * x2_true - 0.2 * x3
  py1 <- exp(linearpred1) / (1 + exp(linearpred1))
  py2 <- exp(linearpred2) / (1 + exp(linearpred2)) - py1
  py3 <- 1 - py1 - py2
  y <- sapply(1:n, function(x) sample(size = 1, c(1:3), prob = c(py1[x], py2[x], py3[x])))
  data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
                     x3 = x3, y = y)
  reg_naive <- MASS::polr(factor(y) ~ x1 + x2_error + x3, data = data)
  reg_true <- MASS::polr(factor(y) ~ x1 + x2_true + x3, data = data)
  reg_rc <- rcreg(reg = reg_naive, data = data, MEvariable = "x2_error",
                  MEerror = 0.5, variance = TRUE, nboot = 2)
  
  # test
  expect_equal(unname(coef(reg_rc))[1:3], as.numeric(t(unname(coef(reg_true)))), tolerance = 0.1)
  expect_equal(formula(reg_rc), as.formula(factor(y) ~ x1 + x2_error + x3))
  expect_equal(as.numeric(predict(reg_rc, newdata = data[1, ], type = "probs")),
               c(py1[1], py2[1], py3[1]), tolerance = 0.1)
  expect_equal(model.frame(reg_rc), model.frame(reg_naive))
  expect_equal(update(reg_rc, data = data)$RCcoef, coef(reg_rc))
  expect_equal(summary(reg_rc)$summarydf$Estimate, as.numeric(coef(reg_rc)))
  
})

test_that("rcreg works correctly for coxph", {
  
  set.seed(1)
  n <- 10000
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2_true <- rnorm(n, mean = 1, sd = 1)
  error1 <- rnorm(n, mean = 0, sd = 0.5)
  x2_error <- x2_true + error1
  x3 <- rbinom(n, size = 1, prob = 0.4)
  linearpred <- 1 + 0.3 * x1 - 0.5 * x2_true - 0.2 * x3
  y <- rexp(n, exp(linearpred))
  cen <- quantile(y, 0.75)
  ycen <- pmin(y, cen)
  delta <- as.numeric(y < cen)
  data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
                     x3 = x3, y = y, ycen = ycen, delta = delta)
  reg_naive <- survival::coxph(survival::Surv(ycen, delta) ~ x1 + x2_error + x3, data = data)
  reg_true <- survival::coxph(survival::Surv(ycen, delta) ~ x1 + x2_true + x3, data = data)
  reg_rc <- rcreg(reg = reg_naive, data = data, MEvariable = "x2_error",
                  MEerror = 0.5, variance = TRUE, nboot = 2)
  
  # test
  expect_equal(unname(coef(reg_rc)), as.numeric(t(unname(coef(reg_true)))), tolerance = 0.1)
  expect_equal(formula(reg_rc), as.formula(survival::Surv(ycen, delta) ~ x1 + x2_error + x3))
  expect_equal(model.frame(reg_rc), model.frame(reg_naive))
  expect_equal(update(reg_rc, data = data)$RCcoef, coef(reg_rc))
  expect_equal(summary(reg_rc)$summarydf$Estimate, as.numeric(coef(reg_rc)))
  
})