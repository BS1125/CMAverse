#devtools::install_github("gerkelab/mediator")
#install.packages("devtools")
#install.packages("regmedint")
#install.packages("usethis")
#devtools::install_github("gerkelab/mediator")
#devtools::install_github("kaz-yos/regmedint")
library(mediator)
library(regmedint)
library(CMA)

data(vv2015)

# data(vv2015) + EM no int

summary(cmest(data = vv2015, model = "rb", outcome = "y", exposure = "x",
              mediator = "m", prec = "c", EMint = FALSE, mreg = "logistic", yreg = "linear",
              astar = 0, a = 1, mval = 1, estimation = "paramfunc", inference = "delta"))

coef(summary(regmedint(data = vv2015, yvar = "y", avar = "x", mvar = "m", cvar = c("c"),
                       a0 = 0, a1 = 1, m_cde = 1, c_cond = mean(vv2015$c),
                       mreg = "logistic", yreg = "linear", interaction = FALSE), exponentiate = TRUE))

mediator::mediator(data = vv2015, a = 1, a_star = 0, m = 1,
                   treat = "x",
                   out.model = glm(y ~ x + m + c,
                                   family = "gaussian",
                                   data = vv2015),
                   med.model = glm(m ~ x + c,
                                   family = "binomial",
                                   data = vv2015)) #nde CI not correct

# data(vv2015) + EM int

summary(cmest(data = vv2015, model = "rb", outcome = "y", exposure = "x",
              mediator = "m", prec = "c", EMint = TRUE, mreg = "logistic", yreg = "linear",
              astar = 0, a = 1, mval = 1, estimation = "paramfunc", inference = "delta"))

coef(summary(regmedint(data = vv2015, yvar = "y", avar = "x", mvar = "m", cvar = c("c"),
                       a0 = 0, a1 = 1, m_cde = 1, c_cond = mean(vv2015$c),
                       mreg = "logistic", yreg = "linear", interaction = TRUE), exponentiate = TRUE))

mediator::mediator(data = vv2015, a = 1, a_star = 0, m = 1,
                   treat = "x",
                   out.model = glm(y ~ x + m + c + x * m,
                                   family = "gaussian",
                                   data = vv2015),
                   med.model = glm(m ~ x + c,
                                   family = "binomial",
                                   data = vv2015)) #nde CI not correct


# mediation_example + binaryM + EM int

summary(cmest(data = mediation_example, model = "rb", outcome = "y", exposure = "x",
              mediator = "m_01", prec = "c", EMint = TRUE, mreg = "logistic", yreg = "logistic",
              astar = 0, a = 1, mval = 1, estimation = "paramfunc", inference = "delta"))

coef(summary(regmedint(data = mediation_example, yvar = "y", avar = "x", mvar = "m_01", cvar = c("c"),
                       a0 = 0, a1 = 1, m_cde = 1, c_cond = mean(mediation_example$c),
                       mreg = "logistic", yreg = "logistic", interaction = TRUE), exponentiate = TRUE))

mediator::mediator(data = mediation_example, treat = "x", a = 1, a_star = 0, m = 1,
                   out.model = glm(y ~ x + m_01 + c + x*m_01,
                                   family = "binomial",
                                   data = mediation_example),
                   med.model = glm(m_01 ~ x + c,
                                   family = "binomial",
                                   data = mediation_example))


# mediation_example + binaryM + EM no int

summary(cmest(data = mediation_example, model = "rb", outcome = "y", exposure = "x",
              mediator = "m_01", prec = "c", EMint = FALSE, mreg = "logistic", yreg = "logistic",
              astar = 0, a = 1, mval = 1, estimation = "paramfunc", inference = "delta"))

coef(summary(regmedint(data = mediation_example, yvar = "y", avar = "x", mvar = "m_01", cvar = c("c"),
                       a0 = 0, a1 = 1, m_cde = 1, c_cond = mean(mediation_example$c),
                       mreg = "logistic", yreg = "logistic", interaction = FALSE), exponentiate = TRUE))

mediator::mediator(data = mediation_example, treat = "x", a = 1, a_star = 0, m = 1,
                   out.model = glm(y ~ x + m_01 + c,
                                   family = "binomial",
                                   data = mediation_example),
                   med.model = glm(m_01 ~ x + c,
                                   family = "binomial",
                                   data = mediation_example))
