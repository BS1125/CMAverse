#class

#lm(): lm
#glm(): lm, glm
#multinom(): multinom, nnet
#gam(): lm, glm, gam
#polr(): polr
#coxph(): coxph
#survreg(): survreg
#glm.nb(): negbin, glm, lm

#predict.glm,predict.gam:type="response

###########################################Continuous###########################################

x=1:20
y=2*x+1+rnorm(20)

data <- data.frame(x=x,y=y)

# 1. lm()

reg <- lm(y ~ x, data = data, x=TRUE)
summary(reg)
reg$call$data
names(reg$call)
match(c("formula", "data", "subset", "weights", "na.action", "offset"),
      names(reg$call), 0L)
reg$model
reg$coefficients=c(20,60)
sum((fitted(reg)-y)^2)/18
sigma(reg)^2
reg$x
lm(y ~ x, data = data)
class(reg)
family(reg)$family
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response")
formula(reg)
reg$call$formula
coef(reg)
reg$coefficients
reg$model
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg, weights=c(100,rep(0,nrow(data)-1)))


# 2. glm(family = gaussian())

reg <- glm(y ~ x, data = data)
class(reg)
family(reg)$family
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response")
formula(reg)
coef(reg)
reg$coefficients
reg$model
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 3. glm(family = Gamma())
x <- runif(100, -1, 1)
a <- 0.5
b <- 1.2
y_true <- exp(0.5 + 1.2 * x)
y <- rgamma(100, rate = 10 / y_true, shape = 10)
data <- data.frame(x = x, y = y)
reg <- glm(y ~ x, family = Gamma(link = "log"), data = data)
family(reg)$family #[1] "Gamma"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response")
formula(reg)
coef(reg)
reg$coefficients
reg$model
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 4. glm(family = inverse.gaussian())

reg <- glm(y ~ x, family = inverse.gaussian(), data = data[data$x > 0, ])
family(reg)$family #[1] "inverse.gaussian"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response")
formula(reg)
coef(reg)
reg$coefficients
reg$model
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data[data$x > 0, ])-1))
update(reg)

# 5. glm(family = quasi())

reg <- glm(y ~ x, family = quasi(), data = data)
family(reg)$family #[1] "quasi"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response")
formula(reg)
coef(reg)
reg$coefficients
reg$model
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 6. gam(family = gaussian())

library(mgcv)
data <- gamSim(1,n=400,dist="normal",scale=2)
reg <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=data)
family(reg)$family #[1] "gaussian"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response")
formula(reg)
reg$model
coef(reg)
reg$coefficients
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
coef(update(reg))

# 7. gam(family = Gamma())

reg <- gam(y~s(x0)+s(x1)+s(x2)+s(x3), data=data[(data$y > 0)&(data$x0 > 0)&(data$x1 > 0)&(data$x2 > 0)&(data$x3 > 0), ], family = Gamma())
family(reg)$family #[1] "Gamma"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response")
formula(reg)
reg$model
coef(reg)
reg$coefficients
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data[(data$y > 0)&(data$x0 > 0)&(data$x1 > 0)&(data$x2 > 0)&(data$x3 > 0), ])-1))
update(reg)

# 8. gam(family = inverse.gaussian())

reg <- gam(y~s(x0)+s(x1)+s(x2)+s(x3), data=data[(data$y > 0)&(data$x0 > 0)&(data$x1 > 0)&(data$x2 > 0)&(data$x3 > 0), ], family = inverse.gaussian())
family(reg)$family #[1] "inverse.gaussian"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response")
formula(reg)
reg$model
coef(reg)
reg$coefficients
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data[(data$y > 0)&(data$x0 > 0)&(data$x1 > 0)&(data$x2 > 0)&(data$x3 > 0), ])-1))
update(reg)

# 9. gam(family = quasi())

reg <- gam(y~s(x0)+s(x1)+s(x2)+s(x3), data=data, family = quasi())
family(reg)$family #[1] "quasi"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response")
formula(reg)
reg$model
coef(reg)
reg$coefficients
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 10. gam(family=Tweedie())

library(mgcv)

dat <- gamSim(1,n=400,dist="poisson",scale=.2)
dat$y <- rTweedie(exp(dat$f),p=1.3,phi=.5) ## Tweedie response

## Fit a fixed p Tweedie, with wrong link ...
reg<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=Tweedie(1.25,power(.1)),
         data=dat)
family(reg)$family # "Tweedie(1.25)"
attributes(reg$terms)$dataClasses
predict(reg, newdata = dat, type = "response")
formula(reg)
reg$model
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(dat)-1))
update(reg)

# 11. gam(family=betar())

library(mgcv)
dat <- gamSim(1,n=400)
mu <- binomial()$linkinv(dat$f/4-2)
phi <- .5
a <- mu*phi;b <- phi - a;
dat$y <- rbeta(400,a,b)

reg <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=betar(link="logit"),data=dat)
family(reg)$family # "Beta regression(0.547)"
attributes(reg$terms)$dataClasses
predict(reg, newdata = dat, type = "response")
formula(reg)
reg$model
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(dat)-1))
update(reg)

# 12. gam(family=scat())

library(mgcv)

dat <- gamSim(1,n=400)
dat$y <- dat$f + rt(400,df=4)*2

reg <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=scat(link="identity"),data=dat)
family(reg)$family # "Scaled t(4.366,1.971)"
attributes(reg$terms)$dataClasses
predict(reg, newdata = dat, type = "response")
formula(reg)
reg$model
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(dat)-1))
update(reg)

# 13. gam(family=gaulss())
library(mgcv);library(MASS)
reg <- gam(list(accel~s(times,k=20,bs="ad"),~s(times)),
         data=mcycle,family=gaulss(b=0.02))
family(reg)$family # "gaulss"
attributes(reg$terms)$dataClasses
names(attributes(reg$terms)$dataClasses)
predict(reg, newdata = mcycle, type = "response")
formula(reg)
reg$model
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100000,rep(0,nrow(mcycle)-1))
update(reg)

# 14. gam(family=gevlss())
library(mgcv)
Fi.gev <- function(z,mu,sigma,xi) {
  ## GEV inverse cdf.
  xi[abs(xi)<1e-8] <- 1e-8 ## approximate xi=0, by small xi
  x <- mu + ((-log(z))^-xi-1)*sigma/xi
}

## simulate test data...
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 *
  (10 * x)^3 * (1 - x)^10
set.seed(1)
n <- 500
x0 <- runif(n);x1 <- runif(n);x2 <- runif(n)
mu <- f2(x2)
rho <- f0(x0)
xi <- (f1(x1)-4)/9
y <- Fi.gev(runif(n),mu,exp(rho),xi)
dat <- data.frame(y,x0,x1,x2);pairs(dat)

## fit model....
reg <- gam(list(y~s(x2),~s(x0),~s(x1)),family=gevlss,data=dat)
family(reg)$family # "gevlss"
attributes(reg$terms)$dataClasses
names(attributes(reg$terms)$dataClasses)
predict(reg, newdata = dat, type = "response")
formula(reg)[[1]][[2]]
all.vars(formula(reg)[[3]])
reg$model
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(dat)-1))
update(reg)

# 15. gam(family=mvn()) actually doesn't allow
library(mgcv)
## simulate some data...
V <- matrix(c(2,1,1,2),2,2)
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 *
  (10 * x)^3 * (1 - x)^10
n <- 300
x0 <- runif(n);x1 <- runif(n);
x2 <- runif(n);x3 <- runif(n)
y <- matrix(0,n,2)
for (i in 1:n) {
  mu <- c(f0(x0[i])+f1(x1[i]),f2(x2[i]))
  y[i,] <- rmvn(1,mu,V)
}
dat <- data.frame(y0=y[,1],y1=y[,2],x0=x0,x1=x1,x2=x2,x3=x3)

## fit model...

reg <- gam(list(y0~s(x0)+s(x1),y1~s(x2)+s(x3)),family=mvn(d=2),data=dat)
family(reg)$family # "Multivariate normal"
attributes(reg$terms)$dataClasses
predict(reg, newdata = dat, type = "response")
formula(reg)
reg$model
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(dat)-1))
update(reg)

###############################################Count##############################################

# 1. glm(family = poisson())

y <- c(18,17,15,20,10,20,25,13,12)
x <- gl(3,1,9)
data <- data.frame(y, x)
reg <- glm(y ~ x, family = poisson(), data = data)
family(reg)$family #[1] "poisson"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response") # expected counts, not necessarily integer
formula(reg)
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 2. glm(family = quasipoisson())

reg <- glm(y ~ x, family = quasipoisson(), data = data)
family(reg)$family #[1] "quasipoisson"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response") # expected counts, not necessarily integer
formula(reg)
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 3. glm.nb()

reg <- MASS::glm.nb(y ~ x, data = data)
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response") # expected counts, not necessarily integer
formula(reg)
class(reg)
coef(reg)
reg$coefficients
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 4. gam(family = poisson())

data <- gamSim(1,n=100,dist="poisson",scale=.1)
reg<-gam(y~s(x0,bs="cr")+s(x1,bs="cr")+s(x2,bs="cr")+
           s(x3,bs="cr"),family=poisson,data=data,method="REML")
family(reg)$family #[1] "poisson"
attributes(reg$terms)$dataClasses
reg_simex=simex(reg,SIMEXvariable = "x1",measurement.error = 0.1,asymptotic = FALSE)
predict(reg, newdata = data, type = "response") # expected counts, not necessarily integer
formula(reg)
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 5. gam(family = quasipoisson())

reg<-gam(y~s(x0,bs="cr")+s(x1,bs="cr")+s(x2,bs="cr")+
           s(x3,bs="cr"),family=quasipoisson,data=data,method="REML")
family(reg)$family #[1] "quasipoisson"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response") # expected counts, not necessarily integer
formula(reg)
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 6. gam(family = negbin())

reg<-gam(y~s(x0,bs="cr")+s(x1,bs="cr")+s(x2,bs="cr")+
           s(x3,bs="cr"),family=negbin(3),data=data,method="REML")
family(reg)$family #[1] "Negative Binomial(3)"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response") # expected counts, not necessarily integer
formula(reg)
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 7. gam(family = nb())

reg<-gam(y~s(x0,bs="cr")+s(x1,bs="cr")+s(x2,bs="cr")+
           s(x3,bs="cr"),family=nb,data=data,method="REML")
family(reg)$family #[1] "Negative Binomial(28.964)"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response") # expected counts, not necessarily integer
formula(reg)
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 8. gam(family = ziP())

rzip <- function(gamma,theta= c(-2,.3)) {
  ## generate zero inflated Poisson random variables, where
  ## lambda = exp(gamma), eta = theta[1] + exp(theta[2])*gamma
  ## and 1-p = exp(-exp(eta)).
  y <- gamma; n <- length(y)
  lambda <- exp(gamma)
  eta <- theta[1] + exp(theta[2])*gamma
  p <- 1- exp(-exp(eta))
  ind <- p > runif(n)
  y[!ind] <- 0
  np <- sum(ind)
  ## generate from zero truncated Poisson, given presence...
  y[ind] <- qpois(runif(np,dpois(0,lambda[ind]),1),lambda[ind])
  y
}

library(mgcv)
## Simulate some ziP data...
set.seed(1);n<-400
dat <- gamSim(1,n=n)
dat$y <- rzip(dat$f/4-1)

reg <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=ziP(),data=dat)
family(reg)$family #[1] "Zero inflated Poisson(-1.77,1.187)"
attributes(reg$terms)$dataClasses
predict(reg, newdata = dat, type = "response") # # expected counts, not necessarily integer
formula(reg)
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(dat)-1))
update(reg)

# 9. gam(family = ziplss())

library(mgcv)
## simulate some data...
f0 <- function(x) 2 * sin(pi * x); f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 *
  (10 * x)^3 * (1 - x)^10
n <- 500;set.seed(5)
x0 <- runif(n); x1 <- runif(n)
x2 <- runif(n); x3 <- runif(n)

## Simulate probability of potential presence...
eta1 <- f0(x0) + f1(x1) - 3
p <- binomial()$linkinv(eta1)
y <- as.numeric(runif(n)<p) ## 1 for presence, 0 for absence

## Simulate y given potentially present (not exactly model fitted!)...
ind <- y>0
eta2 <- f2(x2[ind])/3
y[ind] <- rpois(exp(eta2),exp(eta2))

## Fit ZIP model...
reg <- gam(list(y~s(x2)+s(x3),~s(x0)+s(x1)),family=ziplss())
family(reg)$family #[1] "ziplss"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data.frame(y=y,x2=x2,x3=x3,x0=x0,x1=x1), type = "response") # # expected counts, not necessarily integer
formula(reg)
formula(reg)[[1]][[1]]
formula(reg)[[1]][[2]]
formula(reg)[[1]][[3]]
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

###########################################Binary###########################################

# 1. glm(family = binomial())

x = rnorm(100)
y = rbinom(n,1,exp(x)/(1+exp(x)))
data <- data.frame(x = x, y = y)
reg <- glm(y~x, data = data,family = binomial()) # y should be 0 or 1 or factor
reg <- glm(y~x, data = data,family = binomial("probit"))
family(reg)$family #[1] "binomial"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response") # probability of success
formula(reg)
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 2. glm(family = quasibinomial())

reg <- glm(y~x, data = data,family = quasibinomial())
family(reg)$family #[1] "quasibinomial"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response") # probability of success
formula(reg)
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 3. gam(family = binomial())

reg <- gam(y~s(x), data = data,family = binomial())
family(reg)$family #[1] "binomial"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response") # probability of success
formula(reg)
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 4. gam(family = quasibinomial())

reg <- gam(y~s(x), data = data, family = quasibinomial())
family(reg)$family #[1] "quasibinomial"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response") # probability of success
formula(reg)
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

#############################################Nominal#############################################

# 1. multinom()

y <- sample(c("a","b","c"),100,prob = c(0.2,0.3,0.5),replace = T)
x <- rnorm(100)
data <- data.frame(y,x)
reg=nnet::multinom(y~x,data=data)
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "probs")
# probability of each category, with colnames. if 2 levels, only probability of success
formula(reg)
class(reg)
coef(reg) #matrix
as.vector(t(coef(reg)))
vcov(reg)
rownames(coef(reg))
colnames(coef(reg))
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)


# 2. gam(family = multinom()) outcome can't be factor

library(mgcv)

n <- 1000
f1 <- function(x) sin(3*pi*x)*exp(-x)
f2 <- function(x) x^3
f3 <- function(x) .5*exp(-x^2)-.2
f4 <- function(x) 1
x1 <- runif(n);x2 <- runif(n)
eta1 <- 2*(f1(x1) + f2(x2))-.5
p <- exp(cbind(0,eta1))
p <- p/rowSums(p) ## prob. of each category
cp <- t(apply(p,1,cumsum))
y <- apply(cp,1,function(x) min(which(x>runif(1))))-1

## now fit the model...
reg <- gam(list(y~s(x1)+s(x2)),family=multinom(K=1))

family(reg)$family #[1] "multinom"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data.frame(x1=x1,x2=x2), type = "response")
# probability of each category, no colnames.##columns=#levels
formula(reg)
coef(reg) #vector
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

############################################Ordinal#################################################

# 1. gam(family = ocat()) #outcome values start from 1, not factor. nlevel>2
library(mgcv)
n<-400
data <- gamSim(1,n=n)
data$f <- data$f - mean(data$f)

alpha <- c(-Inf,-1,0,5,Inf)
R <- length(alpha)-1
y <- data$f
u <- runif(n)
u <- data$f + log(u/(1-u))
for (i in 1:R) {
  y[u > alpha[i]&u <= alpha[i+1]] <- i
}
data$y <- y
reg <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=ocat(R=R),data=data)
family(reg)$family #[1] "Ordered Categorical(-1,0.1,6.24)"
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "response") # probability of each category, no colnames
formula(reg)
formula(reg)[[1]]
formula(reg)[[2]]
formula(reg)[[3]]
coef(reg)
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 2. polr() only can be used when nlevel>2. response must be a factor

library(MASS)
options(contrasts = c("contr.treatment", "contr.poly"))
reg <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
attributes(reg$terms)$dataClasses
predict(reg, newdata = housing, type = "probs") # probability of each category, with colnames
formula(reg)
coef(reg)
reg$coefficients
reg$zeta
c(coef(reg),reg$zeta)
vcov(reg)
dim(vcov(reg))
reg$call[["weights"]]=c(100,rep(0,nrow(housing)-1))
update(reg)

#########################################survival#################################################

# 1. coxph()

data <- data.frame(time=c(4,3,1,1,2,2,3),
              status=c(1,1,1,0,1,1,0),
              x=c(0,2,1,1,1,0,0),
              sex=c(0,0,0,0,1,1,1))
# Fit a stratified model
reg=coxph(Surv(time, status) ~ x + strata(sex), data)
attributes(reg$terms)$dataClasses
predict(reg, newdata = data, type = "risk")
formula(reg)
coef(reg)
reg$coefficients
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)

# 2. survreg()
library(survival)
reg=survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='weibull',
        scale=1)
attributes(reg$terms)$dataClasses
predict(reg, newdata = ovarian, type = "response")
formula(reg)[[1]]
formula(reg)[[2]]
formula(reg)[[3]]
coef(reg)
reg$coefficients
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(ovarian)-1))
update(reg)

# 3. gam(family=cox.ph())
library(mgcv)
library(survival) ## for data
col1 <- colon[colon$etype==1,] ## concentrate on single event
col1$differ <- as.factor(col1$differ)
col1$sex <- as.factor(col1$sex)

reg <- gam(time~s(age,by=sex)+sex+s(nodes)+perfor+rx+obstruct+adhere,
         family=cox.ph(),data=col1,weights=status)
family(reg)$family # "Cox PH"
reg$model
attributes(reg$terms)$dataClasses
predict(reg, newdata = ovarian, type = "response")
formula(reg)
coef(reg)
reg$coefficients
vcov(reg)
reg$call[["weights"]]=c(100,rep(0,nrow(data)-1))
update(reg)
