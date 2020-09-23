context("simexreg works correctly")

test_that("simexreg works correctly for lm", {
  
  # a continuous variable measured with error
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
  reg_simex <- simexreg(reg = reg_naive, data = data, MEvariable = "x2_error", MEvartype = "con",
                  MEerror = 0.5, variance = TRUE)
  
  # test
  expect_equal(unname(coef(reg_simex)), unname(coef(reg_true)), tolerance = 0.1)
  expect_equal(formula(reg_simex), as.formula(y ~ x1 + x2_error + x3))
  expect_equal(sigma(reg_simex), 2, tolerance = 0.1)
  expect_equal(family(reg_simex), family(reg_true))
  expect_equal(as.numeric(predict(reg_simex, newdata = data[1, ], type = "response")), 
               as.numeric(unname(coef(reg_simex)) %*% c(1, as.numeric(data[1, c(1, 3, 4)]))))
  expect_equal(model.frame(reg_simex), model.frame(reg_naive))
  expect_equal(update(reg_simex, data = data)$SIMEXcoef, coef(reg_simex), tolerance = 0.1)
  expect_equal(summary(reg_simex)$summarydf$Estimate, as.numeric(coef(reg_simex)))
  expect_equal(class(print(reg_simex)), "numeric")
  expect_equal(class(print(summary(reg_simex))), "data.frame")
  
  # a continuous dependent variable measured with error
  set.seed(1)
  n <- 10000
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2 <- rnorm(n, mean = 1, sd = 1)
  x3 <- rbinom(n, size = 1, prob = 0.4)
  error1 <- rnorm(n, mean = 0, sd = 0.5)
  y_true <- 1 + 2 * x1 + 4 * x2 + 2 * x3  + rnorm(n, mean = 0, sd = 2)
  y_error <- y_true + error1
  data <- data.frame(x1 = x1, x2 = x2,
                     x3 = x3, y_true = y_true, y_error = y_error)
  reg_naive <- lm(y_error ~ x1 + x2 + x3, data = data)
  reg_true <- lm(y_true ~ x1 + x2 + x3, data = data)
  reg_simex <- simexreg(reg = reg_naive, data = data, MEvariable = "y_error", MEvartype = "con",
                        MEerror = 0.5, variance = TRUE)
  
  # test
  expect_equal(unname(coef(reg_simex)), unname(coef(reg_true)), tolerance = 0.1)
  expect_equal(sigma(reg_simex), 2, tolerance = 0.1)
  
  # a categorical variable measured with error
  set.seed(1)
  n <- 10000
  x1 <- rnorm(n, mean = 5, sd = 3)
  x2_true <- sample(x = c(1:3), size = n, prob = c(0.2,0.3,0.5), replace = TRUE)
  MEerror <- matrix(c(0.9,0.05,0.05,0.06,0.92,0.02,0.02,0.03,0.95), nrow = 3)
  x2_error <- x2_true
  for (j in 1:3) {
    x2_error[which(x2_error == c(1:3)[j])] <-
      sample(x = c(1:3), size = length(which(x2_error == c(1:3)[j])),
             prob = MEerror[, j], replace = TRUE)
  }
  x2_true <- as.factor(x2_true)
  x2_error <- as.factor(x2_error)
  x3 <- rnorm(n, mean = 2, sd = 1)
  y <- 1 + 2 * x1 + 1.5*(x2_true == 2) + 2.5*(x2_true == 3) + 2 * x3  +
    rnorm(n, mean = 0, sd = 2)
  data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
                     x3 = x3, y = y)
  reg_naive <- lm(y ~ x1 + x2_error + x3, data = data)
  reg_true <- lm(y ~ x1 + x2_true + x3, data = data)
  reg_simex <- simexreg(reg = reg_naive, data = data, MEvariable = "x2_error", MEvartype = "cat",
                        MEerror = MEerror, variance = TRUE)
  
  # test
  expect_equal(unname(coef(reg_simex)), unname(coef(reg_true)), tolerance = 0.1)
  expect_equal(formula(reg_simex), as.formula(y ~ x1 + x2_error + x3))
  expect_equal(sigma(reg_simex), 2, tolerance = 0.1)
  expect_equal(family(reg_simex), family(reg_true))
  expect_equal(as.numeric(predict(reg_simex, newdata = data[1, ], type = "response")), 
               as.numeric(unname(coef(reg_simex)) %*% c(1, data[1,1], as.numeric(data[1,3] == 2), 
                                                        as.numeric(data[1,3] == 3), data[1,4])))
  expect_equal(model.frame(reg_simex), model.frame(reg_naive))
  expect_equal(update(reg_simex, data = data)$SIMEXcoef, coef(reg_simex), tolerance = 0.1)
  expect_equal(summary(reg_simex)$summarydf$Estimate, as.numeric(coef(reg_simex)))
  
})

test_that("simexreg works correctly for glm", {
  
  # a continuous variable measured with error
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
  reg_simex <- simexreg(reg = reg_naive, data = data, MEvariable = "x2_error", MEvartype = "con",
                  MEerror = 0.5, variance = TRUE)

  # test
  expect_equal(unname(coef(reg_simex)), unname(coef(reg_true)), tolerance = 0.1)
  expect_equal(formula(reg_simex), as.formula(y ~ x1 + x2_error + x3))
  expect_equal(family(reg_simex), family(reg_true))
  expect_equal(as.numeric(predict(reg_simex, newdata = data[1, ], type = "response")), 
               as.numeric(exp(t(unname(coef(reg_simex))) %*% c(1, as.numeric(data[1, c(1, 3, 4)])))/
                 (1+exp(t(unname(coef(reg_simex))) %*% c(1, as.numeric(data[1, c(1, 3, 4)]))))))
  expect_equal(model.frame(reg_simex), model.frame(reg_naive))
  expect_equal(update(reg_simex, data = data)$SIMEXcoef, coef(reg_simex), tolerance = 0.01)
  expect_equal(summary(reg_simex)$summarydf$Estimate, as.numeric(coef(reg_simex)))
  
  # a categorical variable measured with error
  set.seed(1)
  n <- 10000
  x1 <- rnorm(n, mean = 5, sd = 3)
  x2_true <- sample(x = c(1:3), size = n, prob = c(0.2,0.3,0.5), replace = TRUE)
  MEerror <- matrix(c(0.9,0.05,0.05,0.06,0.92,0.02,0.02,0.03,0.95), nrow = 3)
  x2_error <- x2_true
  for (j in 1:3) {
    x2_error[which(x2_error == c(1:3)[j])] <-
      sample(x = c(1:3), size = length(which(x2_error == c(1:3)[j])),
             prob = MEerror[, j], replace = TRUE)
  }
  x2_true <- as.factor(x2_true)
  x2_error <- as.factor(x2_error)
  x3 <- rnorm(n, mean = 2, sd = 1)
  linearpred <- 1 + 0.3 * x1 - 1.5*(x2_true == 2) - 2.5*(x2_true == 3) - 0.2 * x3
  py <- exp(linearpred) / (1 + exp(linearpred))
  y <- rbinom(n, size = 1, prob = py)
  data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
                     x3 = x3, y = y)
  reg_naive <- glm(y ~ x1 + x2_error + x3, data = data, family = binomial("logit"))
  reg_true <- glm(y ~ x1 + x2_true + x3, data = data, family = binomial("logit"))
  reg_simex <- simexreg(reg = reg_naive, data = data, MEvariable = "x2_error",
                        MEerror = MEerror, variance = TRUE, MEvartype = "cat")
  
  # test
  expect_equal(unname(coef(reg_simex)), unname(coef(reg_true)), tolerance = 0.1)
  expect_equal(formula(reg_simex), as.formula(y ~ x1 + x2_error + x3))
  expect_equal(family(reg_simex), family(reg_true))
  expect_equal(as.numeric(predict(reg_simex, newdata = data[1, ], type = "response")), 
               as.numeric(exp(t(unname(coef(reg_simex))) %*% c(1, data[1,1], as.numeric(data[1,3] == 2), 
                                                               as.numeric(data[1,3] == 3), data[1,4]))/
                            (1+exp(t(unname(coef(reg_simex))) %*% c(1, data[1,1], as.numeric(data[1,3] == 2), 
                                                                    as.numeric(data[1,3] == 3), data[1,4])))))
  expect_equal(model.frame(reg_simex), model.frame(reg_naive))
  expect_equal(update(reg_simex, data = data)$SIMEXcoef, coef(reg_simex), tolerance = 0.01)
  expect_equal(summary(reg_simex)$summarydf$Estimate, as.numeric(coef(reg_simex)))
  
})

test_that("simexreg works correctly for multinom", {
  
  # a continuous variable measured with error
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
  reg_simex <- simexreg(reg = reg_naive, data = data, MEvariable = "x2_error", MEvartype = "con",
                  MEerror = 0.5, variance = TRUE)
  
  # test
  expect_equal(unname(coef(reg_simex)), as.numeric(t(unname(coef(reg_true)))), tolerance = 0.1)
  expect_equal(formula(reg_simex), as.formula(factor(y) ~ x1 + x2_error + x3))
  expect_equal(as.numeric(predict(reg_simex, newdata = data[1, ], type = "probs")),
               c(py1[1], py2[1], py3[1]), tolerance = 0.1)
  expect_equal(model.frame(reg_simex), model.frame(reg_naive))
  expect_equal(update(reg_simex, data = data)$SIMEXcoef, coef(reg_simex), tolerance = 0.01)
  expect_equal(summary(reg_simex)$summarydf$Estimate, as.numeric(coef(reg_simex)))
  
  # a categorical variable measured with error
  set.seed(1)
  n <- 10000
  x1 <- rnorm(n, mean = 5, sd = 3)
  x2_true <- sample(x = c(1:3), size = n, prob = c(0.2,0.3,0.5), replace = TRUE)
  MEerror <- matrix(c(0.9,0.05,0.05,0.06,0.92,0.02,0.02,0.03,0.95), nrow = 3)
  x2_error <- x2_true
  for (j in 1:3) {
    x2_error[which(x2_error == c(1:3)[j])] <-
      sample(x = c(1:3), size = length(which(x2_error == c(1:3)[j])),
             prob = MEerror[, j], replace = TRUE)
  }
  x2_true <- as.factor(x2_true)
  x2_error <- as.factor(x2_error)
  x3 <- rnorm(n, mean = 2, sd = 1)
  linearpred1 <- 1 + 0.3 * x1 - 1.5*(x2_true == 2) - 2.5*(x2_true == 3) - 0.2 * x3
  linearpred2 <- 2 + 1 * x1 - 2*(x2_true == 2) - 0.5*(x2_true == 3) - 1 * x3
  py2 <- exp(linearpred1) / (1 + exp(linearpred1) + exp(linearpred2))
  py3 <- exp(linearpred2) / (1 + exp(linearpred1) + exp(linearpred2))
  py1 <- 1 - py2 - py3
  y <- sapply(1:n, function(x) sample(size = 1, c(1:3), prob = c(py1[x], py2[x], py3[x])))
  data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
                     x3 = x3, y = y)
  reg_naive <- nnet::multinom(factor(y) ~ x1 + x2_error + x3, data = data)
  reg_true <- nnet::multinom(factor(y) ~ x1 + x2_true + x3, data = data)
  reg_simex <- simexreg(reg = reg_naive, data = data, MEvariable = "x2_error",
                        MEerror = MEerror, variance = TRUE, MEvartype = "cat")
  
  # test
  expect_equal(unname(coef(reg_simex)), as.numeric(t(unname(coef(reg_true)))), tolerance = 0.1)
  expect_equal(formula(reg_simex), as.formula(factor(y) ~ x1 + x2_error + x3))
  expect_equal(as.numeric(predict(reg_simex, newdata = data[1, ], type = "probs")),
               c(py1[1], py2[1], py3[1]), tolerance = 0.1)
  expect_equal(model.frame(reg_simex), model.frame(reg_naive))
  expect_equal(update(reg_simex, data = data)$SIMEXcoef, coef(reg_simex), tolerance = 0.01)
  expect_equal(summary(reg_simex)$summarydf$Estimate, as.numeric(coef(reg_simex)))
  
})

test_that("simexreg works correctly for polr", {
  
  # a continuous variable measured with error
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
  reg_simex <- simexreg(reg = reg_naive, data = data, MEvariable = "x2_error", MEvartype = "con",
                  MEerror = 0.5, variance = TRUE)
  
  # test
  expect_equal(unname(coef(reg_simex))[1:3], as.numeric(t(unname(coef(reg_true)))), tolerance = 0.1)
  expect_equal(formula(reg_simex), as.formula(factor(y) ~ x1 + x2_error + x3))
  expect_equal(as.numeric(predict(reg_simex, newdata = data[1, ], type = "probs")),
               c(py1[1], py2[1], py3[1]), tolerance = 0.1)
  expect_equal(model.frame(reg_simex), model.frame(reg_naive))
  expect_equal(update(reg_simex, data = data)$SIMEXcoef, coef(reg_simex), tolerance = 0.01)
  expect_equal(summary(reg_simex)$summarydf$Estimate, as.numeric(coef(reg_simex)))
  
  # a categorical variable measured with error
  set.seed(1)
  n <- 10000
  x1 <- rnorm(n, mean = 5, sd = 1)
  x2_true <- sample(x = c(1:3), size = n, prob = c(0.2,0.3,0.5), replace = TRUE)
  MEerror <- matrix(c(0.9,0.05,0.05,0.06,0.92,0.02,0.02,0.03,0.95), nrow = 3)
  x2_error <- x2_true
  for (j in 1:3) {
    x2_error[which(x2_error == c(1:3)[j])] <-
      sample(x = c(1:3), size = length(which(x2_error == c(1:3)[j])),
             prob = MEerror[, j], replace = TRUE)
  }
  x2_true <- as.factor(x2_true)
  x2_error <- as.factor(x2_error)
  x3 <- rnorm(n, mean = 2, sd = 1)
  linearpred1 <- 0.3 * x1 - 0.5*(x2_true == 2) - 0.2*(x2_true == 3) - 0.2 * x3
  linearpred2 <- -1 + 0.2 * x1 - 0.8*(x2_true == 2) - 0.5*(x2_true == 3) - 0.1 * x3
  py2 <- exp(linearpred1) / (1 + exp(linearpred1) + exp(linearpred2))
  py3 <- exp(linearpred2) / (1 + exp(linearpred1) + exp(linearpred2))
  py1 <- 1 - py2 - py3
  y <- sapply(1:n, function(x) sample(size = 1, c(1:3), prob = c(py1[x], py2[x], py3[x])))
  data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
                     x3 = x3, y = y)
  reg_naive <- MASS::polr(factor(y) ~ x1 + x2_error + x3, data = data)
  reg_true <- MASS::polr(factor(y) ~ x1 + x2_true + x3, data = data)
  reg_simex <- simexreg(reg = reg_naive, data = data, MEvariable = "x2_error",
                        MEerror = MEerror, variance = TRUE, MEvartype = "cat")
  
  # test
  expect_equal(unname(coef(reg_simex))[1:4], as.numeric(t(unname(coef(reg_true)))), tolerance = 0.1)
  expect_equal(formula(reg_simex), as.formula(factor(y) ~ x1 + x2_error + x3))
  expect_equal(as.numeric(predict(reg_simex, newdata = data[1, ], type = "probs")),
               c(py1[1], py2[1], py3[1]), tolerance = 0.1)
  expect_equal(model.frame(reg_simex), model.frame(reg_naive))
  expect_equal(update(reg_simex, data = data)$SIMEXcoef, coef(reg_simex), tolerance = 0.01)
  expect_equal(summary(reg_simex)$summarydf$Estimate, as.numeric(coef(reg_simex)))

  })

test_that("simexreg works correctly for coxph", {
  
  # a continuous variable measured with error
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
  reg_simex <- simexreg(reg = reg_naive, data = data, MEvariable = "x2_error", MEvartype = "con",
                  MEerror = 0.5, variance = TRUE)
  
  # test
  expect_equal(unname(coef(reg_simex)), as.numeric(t(unname(coef(reg_true)))), tolerance = 0.1)
  expect_equal(formula(reg_simex), as.formula(survival::Surv(ycen, delta) ~ x1 + x2_error + x3))
  expect_equal(model.frame(reg_simex), model.frame(reg_naive))
  expect_equal(update(reg_simex, data = data)$SIMEXcoef, coef(reg_simex), tolerance = 0.01)
  expect_equal(summary(reg_simex)$summarydf$Estimate, as.numeric(coef(reg_simex)))
  
  # a categorical variable measured with error
  set.seed(1)
  n <- 10000
  x1 <- rnorm(n, mean = 5, sd = 3)
  x2_true <- sample(x = c(1:3), size = n, prob = c(0.2,0.3,0.5), replace = TRUE)
  MEerror <- matrix(c(0.9,0.05,0.05,0.06,0.92,0.02,0.02,0.03,0.95), nrow = 3)
  x2_error <- x2_true
  for (j in 1:3) {
    x2_error[which(x2_error == c(1:3)[j])] <-
      sample(x = c(1:3), size = length(which(x2_error == c(1:3)[j])),
             prob = MEerror[, j], replace = TRUE)
  }
  x2_true <- as.factor(x2_true)
  x2_error <- as.factor(x2_error)
  x3 <- rnorm(n, mean = 2, sd = 1)
  linearpred <- 1 + 0.3 * x1 - 1.5*(x2_true == 2) - 2.5*(x2_true == 3) - 0.2 * x3
  y <- rexp(n, exp(linearpred))
  cen <- quantile(y, 0.75)
  ycen <- pmin(y, cen)
  delta <- as.numeric(y < cen)
  data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
                     x3 = x3, y = y, ycen = ycen, delta = delta)
  reg_naive <- survival::coxph(survival::Surv(ycen, delta) ~ x1 + x2_error + x3, data = data)
  reg_true <- survival::coxph(survival::Surv(ycen, delta) ~ x1 + x2_true + x3, data = data)
  reg_simex <- simexreg(reg = reg_naive, data = data, MEvariable = "x2_error",
                        MEerror = MEerror, variance = TRUE, MEvartype = "cat")
  
  # test
  expect_equal(unname(coef(reg_simex)), as.numeric(t(unname(coef(reg_true)))), tolerance = 0.1)
  expect_equal(formula(reg_simex), as.formula(survival::Surv(ycen, delta) ~ x1 + x2_error + x3))
  expect_equal(model.frame(reg_simex), model.frame(reg_naive))
  expect_equal(update(reg_simex, data = data)$SIMEXcoef, coef(reg_simex), tolerance = 0.01)
  expect_equal(summary(reg_simex)$summarydf$Estimate, as.numeric(coef(reg_simex)))
  
})


test_that("simexreg works correctly for svyglm", {
  
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
  reg_naive <- survey::svyglm(y ~ x1 + x2_error + x3, family = binomial("logit"),
                              design = survey::svydesign(ids = ~1, data = data, w = rep(1,n)))
  reg_true <- survey::svyglm(y ~ x1 + x2_true + x3, family = binomial("logit"),
                             design = survey::svydesign(ids = ~1, data = data, w = rep(1,n)))
  reg_simex <- simexreg(reg = reg_naive, data = data, weights = rep(1,n), MEvariable = "x2_error",
                  MEerror = 0.5, MEvartype = "con", variance = FALSE)
  
  # test
  expect_equal(unname(coef(reg_simex)), unname(coef(reg_true)), tolerance = 0.1)
  expect_equal(formula(reg_simex), as.formula(y ~ x1 + x2_error + x3))
  expect_equal(family(reg_simex), family(reg_true))
  expect_equal(as.numeric(predict(reg_simex, newdata = data[1, ], type = "response")), 
               as.numeric(exp(t(unname(coef(reg_simex))) %*% c(1, as.numeric(data[1, c(1, 3, 4)])))/
                            (1+exp(t(unname(coef(reg_simex))) %*% c(1, as.numeric(data[1, c(1, 3, 4)]))))))
  expect_equal(model.frame(reg_simex), model.frame(reg_naive))
  
})
