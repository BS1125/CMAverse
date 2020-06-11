rm(list=ls())

df_int <- read.csv("./Code Checking/Data Simulation and Code Checking/int_data_10000.txt", sep = " ")

df_int$M_cat_catA=as.factor(df_int$M_cat_catA)

library(CMA)

data=df_int
exposure="A_cont"
ereg="linear"
mediator=c("M_cont_contA","M_bin_contA")
outcome="contY_contM_contA_int"
prec=c("C_cont", "C_bin")
postc=NULL
postcreg=NULL
EMint = TRUE
event=NULL
mreg=list("linear", "logistic")
yreg="linear"
ereg=NULL
wmnomreg=NULL
wmdenomreg=NULL
model="gformula"
a=1
astar=0
mval=list(0)
precval=NULL
yref=NULL
indices=1:nrow(data)
estimation="imputation"
vecc=NULL
inference="bootstrap"
n=nrow(data)
nboot=10

formulas <- create_formulas(outcome=outcome, event=event, exposure=exposure, mediator=mediator,
                            prec=prec, postc=postc, EMint=EMint,
                             mreg=mreg, yreg=yreg, ereg=ereg, postcreg=postcreg,
                           wmreg=wmreg, model=model)

regressions <- run_regressions(formulas=formulas, exposure=exposure, mediator=mediator,
                              postc=postc, mreg=mreg, yreg=yreg, ereg=ereg, postcreg=postcreg,
                              wmreg=wmreg, model=model, data=data)

data=df_int
reg=regressions$outcome_regression
MEvariable=c("A_cont","M_cont_contA")

est <- est_step(indices=1:n, data=data, model=model,
         outcome=outcome, event=event, exposure=exposure, mediator=mediator,
         prec=prec, postc=postc, EMint=EMint,
         mreg=mreg, yreg=yreg, ereg=ereg, postcreg=postcreg,wmreg=wmreg,
                     astar=astar, a=a, mval=mval, yref=yref, vecc=vecc,
                     estimation=estimation)

inf <- inf_step(nboot=10, data=data, model=model,
         outcome=outcome, event=event, exposure=exposure, mediator=mediator,
         prec=prec, postc=postc, EMint=EMint,
         regressions = regressions,
         astar=astar, a=a, mval=mval, yref=yref, vecc=vecc,
         estimation=estimation, inference = inference)

cmest_out <- cmest(data=data, model=model,
      outcome=outcome, event=event, exposure=exposure, mediator=mediator,
      prec=prec, postc=postc, EMint=EMint,
      mreg=mreg, yreg=yreg, ereg=ereg, postcreg=postcreg,wmreg=wmreg,
      astar=astar, a=a, mval=mval, yref=yref, precval=precval,
      estimation=estimation, inference = inference,
      nboot = nboot)

summary(cmest_out)

cmsens.uc <- cmsens(cmest_out, sens = "uc")

MEvariable="A_cont"
MEvariable.type="continuous"
measurement.error=c(0.1,0.3,0.5)
lambda = c(0.5, 1, 1.5, 2)
B=5

Sys.time()
sens_out <- cmsens(cmest_out, sens = "me", MEvariable = MEvariable,
                   MEvariable.type = MEvariable.type, measurement.error = measurement.error,
                   B=B)
Sys.time()

plot(sens_out)


MEvariable="M_cont_contA"

sens_out <- cmsens(cmest_out, sens = "me", MEvariable = MEvariable,
                   MEvariable.type = MEvariable.type, measurement.error = measurement.error,
                   B=B)

plot(sens_out)


MEvariable="contY_contM_contA_int"

sens_out <- cmsens(cmest_out, sens = "me", MEvariable = MEvariable,
                   MEvariable.type = MEvariable.type, measurement.error = measurement.error,
                   B=B)

plot(sens_out)

MEvariable="C_cont"

sens_out <- cmsens(cmest_out, sens = "me", MEvariable = MEvariable,
                   MEvariable.type = MEvariable.type, measurement.error = measurement.error,
                   B=B)

plot(sens_out)

measurement.error = list(matrix(c(c(0.8, 0.1, 0.1, 0.2, 0.7, 0.1, 0.1, 0.2, 0.7)), nrow = 3),
                         matrix(c(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8)), nrow = 3))

paste("Misclassification Matrix = matrix(c(",
      paste(as.vector(measurement.error[[1]]), sep = "", collapse = ","), "), nrow = ",
      3, sep="")

cat("--------\n matrix=")
print(measurement.error[[1]])



reg=regressions$outcome_regression
measurement.error=0.2
MEvariable="contY_contM_contA_int"
MEvariable.type="continuous"
lambda = c(0.5, 1, 1.5, 2)
B=100
variance=TRUE

coef(reg)

param <- param_simex(reg, data, MEvariable, MEvariable.type, measurement.error,
                         lambda = c(0.5, 1, 1.5, 2), B = 1000,
                         variance = TRUE)
param$SIMEXcoef
param$SIMEXvar

try1=simex(reg,SIMEXvariable = MEvariable, B=1000, measurement.error = measurement.error, asymptotic = FALSE)
try2=simex(reg,SIMEXvariable = MEvariable,  B=1000, measurement.error = measurement.error, asymptotic = FALSE)
try3=simex(reg,SIMEXvariable = MEvariable,  B=1000, measurement.error = measurement.error, asymptotic = FALSE)

try1$coefficients
try2$coefficients
try3$coefficients


reg=regressions$outcome_regression
measurement.error=matrix(c(0.8,0.2,0.3,0.7),nrow=2)
MEvariable="C_bin"
MEvariable.type="categorical"
lambda = c(0.5, 1, 1.5, 2)
B=100
variance=TRUE

data=df_int
reg=regressions$exposure_regression
MEvariable=exposure
MEvariable.type="categorical"
measurement.error=matrix(c(0.8,0.1,0.1,0.1,0.6,0.3,0.2,0.1,0.7),nrow=3)
simexreg <- simex_reg(reg, data, MEvariable, MEvariable.type, measurement.error,
                     lambda = c(0.5, 1, 1.5, 2), B = 100)

predict(simexreg,type="probs")


options(contrasts = c("contr.treatment", "contr.poly"))
house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
house.multinom <- nnet::multinom(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
class(house.plr) #[1] "polr"
class(house.multinom) #[1] "multinom" "nnet"
summary(house.plr, digits = 3)
attributes(house.plr$terms)$dataClasses
attributes(house.multinom$terms)$dataClasses
predict(house.plr,type="prob")
predict(house.multinom,type="prob")

y=sample(c("a","b","c"),100,prob = c(0.2,0.3,0.5),replace = T)
x=rnorm(100)
z=rnorm(100)
try=nnet::multinom(y~x)
update(try,formula. = as.formula("y~z"))
predict(try,type="prob")[,"0"]
attributes(try$terms)$dataClasses

library(survival)
y <- rweibull(1000, shape=2, scale=5)
x <- rnorm(1000)
s <- survival::survreg(Surv(y)~x, dist="weibull")
attributes(s$terms)$dataClasses

test1 <- list(time=c(4,3,1,1,2,2,3),
              status=c(1,1,1,0,1,1,0),
              x=c(0,2,1,1,1,0,0),
              sex=c(0,0,0,0,1,1,1))

c <- coxph(Surv(time, status) ~ x + strata(sex), test1)
attributes(c$terms)$dataClasses

set.seed(20160227)
x<-seq(0,50,1)
y<-((runif(1,10,20)*x)/(runif(1,0,10)+x))+rnorm(51,0,1)
m<-nls(y~a1*x/(b+x))
m$dataClasses
class(m)#[1] "nls"

x=seq(0, pi * 2, 0.1)
Sample_data <- data.frame(y=sin(x) + rnorm(n = length(x), mean = 0, sd = sd(sin_x / 2)),
                          x=x)
gam_y <- gam(y ~ s(x), method = "REML",data=Sample_data)
summary(gam_y)
class(gam_y)#[1] "gam" "glm" "lm"
predict(gam_y, type = "response")
attributes(gam_y$terms)$dataClasses




First <- function() {
  message("When there exists at least one continuous mediator, we recommend you to use model = 'rb' or 'wb' or 'g-formula' or 'ne'")
  model <- readline(prompt = "Choose model:")
  x=1+2
  print(x)

}

if (interactive()){  if (exposure.type == "continuous" && !(model %in% c("rb", "ne", "g-formula"))) {
  message("When the exposure is continuous, we recommend you to use model = 'rb' or 'g-formula' or 'ne'")
  model <- readline(prompt = "Choose model:")
}}

if (interactive()){  if (((is.character(mreg)&&mreg == "linear")|
                          (!is.character(mreg)&&((inherits(mreg, "glm")|inherits(mreg, "lm"))&&
                                                 family(mreg)$family == "gaussian"))) &&
                         !(model %in% c("rb", "wb", "g-formula", "ne"))) {
  message("When there exists at least one continuous mediator, we recommend you to use model = 'rb' or 'wb' or 'g-formula' or 'ne'")
  model <- readline(prompt = "Choose model:")
}}

if (interactive()){  if ((EMMint|MMint) && !(model %in% c("wb", "iorw", "ne"))) {
  message("When EMMint = TRUE or MMint = TRUE, we recommend you to use model = 'wb' or'iorw' or 'ne'")
  model <- readline(prompt = "Choose model:")
}}

if (interactive()){  if (is.null(covariates.post) == FALSE && !(model %in% c("msm", "g-formula"))) {
  message("When covariates.post != NULL, we recommend you to use model = 'g-formula' or 'msm'")
  model <- readline(prompt = "Choose model:")
}}



library(mgcv)
set.seed(6)
## simulate some data from a three class model
n <- 1000
f1 <- function(x) sin(3*pi*x)*exp(-x)
f2 <- function(x) x^3
f3 <- function(x) .5*exp(-x^2)-.2
f4 <- function(x) 1
x1 <- runif(n);x2 <- runif(n)
eta1 <- 2*(f1(x1) + f2(x2))-.5
eta2 <- 2*(f3(x1) + f4(x2))-1
p <- exp(cbind(0,eta1,eta2))
p <- p/rowSums(p) ## prob. of each category
cp <- t(apply(p,1,cumsum)) ## cumulative prob.
## simulate multinomial response with these probabilities
## see also ?rmultinom
y <- apply(cp,1,function(x) min(which(x>runif(1))))-1
## plot simulated data...
plot(x1,x2,col=y+3)

## now fit the model...
b <- gam(list(y~s(x1)+s(x2),~s(x1)+s(x2)),family=multinom(K=2))
plot(b,pages=1)
gam.check(b)

## now a simple classification plot...
expand.grid(x1=seq(0,1,length=40),x2=seq(0,1,length=40)) -> gr
pp <- predict(b,newdata=gr,type="response")
#输出结果列名不是levels


quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
family(quine.nb1)$family

attributes(exposure_regression$terms)

data$A_bin=as.factor(data$A_bin)
try=update(exposure_regression, data=data)
attributes(try$terms)
