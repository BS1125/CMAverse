library(survival)

beta0 = -0.25
beta1 = 0.5
beta2 = 0.2
theta0 = -5
theta1 =0.8
theta2 = 1.8
theta3 = 0.2
theta4 = 0.1
meanc<-1
expit<-function(x){
  exp(x)/(1+exp(x))
}
n=10000

set.seed(100)
#generate binary exposure and confounder
A = rbinom(n,1,0.4)
C = rnorm(n,mean=1,sd=sqrt(1))

#generate Bin M and Cont M
linpred = (beta0+beta1*A+beta2*C)
probm = expit(linpred)
M_bin = rbinom(n,1,probm)
M_cont = rnorm(n,linpred,0.1)

##################################################no int###################################################

#generate survival data
linpred = exp(-(theta0+theta1*A+theta2*M_bin+theta4*C))
Y_surv_noint = rexp(n,linpred)
cen<-quantile(Y_surv_noint,probs = 0.75)
#censor observations after censoring time
Ycen_noint<-pmin(Y_surv_noint,cen)
#generate event  indicator
cens<-as.numeric(Ycen_noint<cen)
delta<-as.numeric(Y_surv_noint>cen)

#generate Cont Y, Bin Y, Count Y without int
linpred = theta0+theta1*A+theta2*M_bin+theta4*C
contY_binM_noint = rnorm(n,linpred,0.2)
binY_binM_noint = rbinom(n,1,expit(linpred))
countY_binM_noint = rpois(n,abs(linpred))
linpred = theta0+theta1*A+theta2*M_cont+theta4*C
contY_contM_noint = rnorm(n,linpred,0.2)
binY_contM_noint = rbinom(n,1,expit(linpred))
countY_contM_noint = rpois(n,abs(linpred))

noint_data = cbind(A,M_cont,M_bin,Ycen_noint,contY_binM_noint,binY_binM_noint,countY_binM_noint,
                   contY_contM_noint,binY_contM_noint,countY_contM_noint,cen,C,delta)

write.table(noint_data,file="noint_data_10000.txt")

summary(glm(contY_contM_noint~A+M_cont+C,family = gaussian()))
summary(glm(contY_binM_noint~A+M_bin+C,family = gaussian()))
summary(glm(binY_contM_noint~A+M_cont+C,family = binomial()))
summary(glm(binY_binM_noint~A+M_bin+C,family = binomial()))

######################################################int################################################

set.seed(101)
#generate binary exposure and confounder
A = rbinom(n,1,0.4)
C = rnorm(n,mean=1,sd=sqrt(1))

#generate Bin M and Cont M
linpred = (beta0+beta1*A+beta2*C)
probm = expit(linpred)
M_bin = rbinom(n,1,probm)
M_cont = rnorm(n,linpred,0.1)

#generate survival data
linpred = exp(-(theta0+theta1*A+theta2*M_bin+theta3*(A*M_bin)+theta4*C))
Y_surv_int = rexp(n,linpred)
cen<-quantile(Y_surv_int,probs = 0.75)
#censor observations after censoring time
Ycen_int<-pmin(Y_surv_int,cen)
#generate event  indicator
cens<-as.numeric(Ycen_int<cen)
delta<-as.numeric(Y_surv_int>cen)

#generate Cont Y, Bin Y, Count Y with int
linpred = theta0+theta1*A+theta2*M_bin+theta3*(A*M_bin)+theta4*C
contY_binM_int = rnorm(n,linpred,0.2)
binY_binM_int = rbinom(n,1,expit(linpred))
countY_binM_int = rpois(n,abs(linpred))
linpred = theta0+theta1*A+theta2*M_cont+theta3*(A*M_cont)+theta4*C
contY_contM_int = rnorm(n,linpred,0.2)
binY_contM_int = rbinom(n,1,expit(linpred))
countY_contM_int = rpois(n,abs(linpred))

int_data = cbind(A,M_cont,M_bin,Ycen_int,contY_binM_int,binY_binM_int,countY_binM_int,
                   contY_contM_int,binY_contM_int,countY_contM_int,cen,C,delta)
write.table(int_data,file="int_data_10000.txt")

sum(binY_binM_int)
sum(binY_contM_int)

summary(glm(contY_contM_int~A+M_cont+A*M_cont+C,family = gaussian()))
summary(glm(contY_binM_int~A+M_bin+A*M_bin+C,family = gaussian()))
summary(glm(binY_contM_int~A+M_cont+A*M_cont+C,family = binomial()))
summary(glm(binY_binM_int~A+M_bin+A*M_bin+C,family = binomial()))




x1 = rbinom(1000,1,0.4)          # some continuous variables
x2 = rnorm(1000,(0.5+1.5*x1))
z = -6 + 2*x1 + 3*x2 +0.5*x1*x2       # linear combination with a bias
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = rbinom(1000,1,pr)
sum(y)
df = data.frame(y=y,x1=x1,x2=x2)
glm(y~x1*x2,data=df,family="binomial")


##############################################multiple M without int####################################################

beta01 = -0.25
beta11 = 0.5
beta21 = 0.2

beta02 = -0.3
beta12 = 0.4
beta22 = 0.3

theta0 = -5
theta1 =0.8
theta2 = 1.8
theta3 = 1.2
theta4 = 0.1
meanc<-1
expit<-function(x){
  exp(x)/(1+exp(x))
}
n=10000

set.seed(100)
#generate binary exposure and confounder
A = rbinom(n,1,0.4)
C = rnorm(n,mean=1,sd=sqrt(1))

#generate Bin M and Cont M
linpred1 = (beta01+beta11*A+beta21*C)
linpred2 = (beta02+beta12*A+beta22*C)
probm1 = expit(linpred1)
probm2 = expit(linpred2)
M_bin1 = rbinom(n,1,probm1)
M_bin2 = rbinom(n,1,probm2)
M_cont1 = rnorm(n,linpred1,0.1)
M_cont2 = rnorm(n,linpred2,0.1)

#generate Cont Y, Bin Y, Count Y without int
linpred1 = theta0+theta1*A+theta2*M_cont1+theta3*M_cont2+theta4*C
linpred2 = theta0+theta1*A+theta2*M_bin1+theta3*M_cont1+theta4*C
linpred3 = theta0+theta1*A+theta2*M_bin1+theta3*M_bin2+theta4*C
contY_2contM_noint = rnorm(n,linpred1,0.2)
contY_binMcontM_noint = rnorm(n,linpred2,0.2)
contY_2binM_noint = rnorm(n,linpred3,0.2)
binY_2contM_noint = rbinom(n,1,expit(linpred1))
binY_binMcontM_noint = rbinom(n,1,expit(linpred2))
binY_2binM_noint = rbinom(n,1,expit(linpred3))

multipleM_noint_data = cbind(A,M_cont1,M_cont2,M_bin1,M_bin2,contY_2contM_noint,
                             contY_binMcontM_noint,contY_2binM_noint,binY_2contM_noint,
                             binY_binMcontM_noint,binY_binMcontM_noint,binY_2binM_noint,C)

write.table(multipleM_noint_data,file="multipleM_noint_data_10000.txt")


##############################################multiple M with exposure-mediator int##########################

beta01 = -0.25
beta11 = 0.5
beta21 = 0.2

beta02 = -0.3
beta12 = 0.4
beta22 = 0.3

theta0 = -5
theta1 =0.8
theta2 = 1.8
theta3 = 1.2
theta4 = 0.6
theta5 = 0.4
theta6 = 0.1
meanc<-1
expit<-function(x){
  exp(x)/(1+exp(x))
}
n=10000

set.seed(100)
#generate binary exposure and confounder
A = rbinom(n,1,0.4)
C = rnorm(n,mean=1,sd=sqrt(1))

#generate Bin M and Cont M
linpred1 = (beta01+beta11*A+beta21*C)
linpred2 = (beta02+beta12*A+beta22*C)
probm1 = expit(linpred1)
probm2 = expit(linpred2)
M_bin1 = rbinom(n,1,probm1)
M_bin2 = rbinom(n,1,probm2)
M_cont1 = rnorm(n,linpred1,0.1)
M_cont2 = rnorm(n,linpred2,0.1)

#generate Cont Y, Bin Y, Count Y without int
linpred1 = theta0+theta1*A+theta2*M_cont1+theta3*M_cont2+theta4*A*M_cont1+theta5*A*M_cont2+theta6*C
linpred2 = theta0+theta1*A+theta2*M_bin1+theta3*M_cont1+theta4*A*M_bin1+theta5*A*M_cont1+theta4*C
linpred3 = theta0+theta1*A+theta2*M_bin1+theta3*M_bin2+theta4*A*M_bin1+theta5*A*M_bin2+theta4*C
contY_2contM_EMint = rnorm(n,linpred1,0.2)
contY_binMcontM_EMint = rnorm(n,linpred2,0.2)
contY_2binM_EMint = rnorm(n,linpred3,0.2)
binY_2contM_EMint = rbinom(n,1,expit(linpred1))
binY_binMcontM_EMint = rbinom(n,1,expit(linpred2))
binY_2binM_EMint = rbinom(n,1,expit(linpred3))

multipleM_EMint_data = cbind(A,M_cont1,M_cont2,M_bin1,M_bin2,contY_2contM_EMint,
                             contY_binMcontM_EMint,contY_2binM_EMint,binY_2contM_EMint,
                             binY_binMcontM_EMint,binY_binMcontM_EMint,binY_2binM_EMint,C)

write.table(multipleM_EMint_data,file="multipleM_EMint_data_10000.txt")


##############################multiple M with mediator-mediator int##########################

beta01 = -0.25
beta11 = 0.5
beta21 = 0.2

beta02 = -0.3
beta12 = 0.4
beta22 = 0.3

theta0 = -5
theta1 =0.8
theta2 = 1.8
theta3 = 1.2
theta4 = 0.6
theta5 = 0.1
meanc<-1
expit<-function(x){
  exp(x)/(1+exp(x))
}
n=10000

set.seed(100)
#generate binary exposure and confounder
A = rbinom(n,1,0.4)
C = rnorm(n,mean=1,sd=sqrt(1))

#generate Bin M and Cont M
linpred1 = (beta01+beta11*A+beta21*C)
linpred2 = (beta02+beta12*A+beta22*C)
probm1 = expit(linpred1)
probm2 = expit(linpred2)
M_bin1 = rbinom(n,1,probm1)
M_bin2 = rbinom(n,1,probm2)
M_cont1 = rnorm(n,linpred1,0.1)
M_cont2 = rnorm(n,linpred2,0.1)

#generate Cont Y, Bin Y, Count Y without int
linpred1 = theta0+theta1*A+theta2*M_cont1+theta3*M_cont2+theta4*M_cont1*M_cont2+theta5*C
linpred2 = theta0+theta1*A+theta2*M_bin1+theta3*M_cont1+theta4*M_bin1*M_cont1+theta5*C
linpred3 = theta0+theta1*A+theta2*M_bin1+theta3*M_bin2+theta4*M_bin1*M_bin2+theta5*C
contY_2contM_MMint = rnorm(n,linpred1,0.2)
contY_binMcontM_MMint = rnorm(n,linpred2,0.2)
contY_2binM_MMint = rnorm(n,linpred3,0.2)
binY_2contM_MMint = rbinom(n,1,expit(linpred1))
binY_binMcontM_MMint = rbinom(n,1,expit(linpred2))
binY_2binM_MMint = rbinom(n,1,expit(linpred3))

multipleM_MMint_data = cbind(A,M_cont1,M_cont2,M_bin1,M_bin2,contY_2contM_MMint,
                             contY_binMcontM_MMint,contY_2binM_MMint,binY_2contM_MMint,
                             binY_binMcontM_MMint,binY_binMcontM_MMint,binY_2binM_MMint,C)

write.table(multipleM_MMint_data,file="multipleM_MMint_data_10000.txt")

##############################################multiple M with multiway int##########################

beta01 = -0.25
beta11 = 0.5
beta21 = 0.2

beta02 = -0.3
beta12 = 0.4
beta22 = 0.3

theta0 = -5
theta1 =0.8
theta2 = 1.8
theta3 = 1.2
theta4 = 0.6
theta5 = 0.4
theta6 = 0.3
theta7 = 0.5
theta8 = 0.1
meanc<-1
expit<-function(x){
  exp(x)/(1+exp(x))
}
n=10000

set.seed(100)
#generate binary exposure and confounder
A = rbinom(n,1,0.4)
C = rnorm(n,mean=1,sd=sqrt(1))

#generate Bin M and Cont M
linpred1 = (beta01+beta11*A+beta21*C)
linpred2 = (beta02+beta12*A+beta22*C)
probm1 = expit(linpred1)
probm2 = expit(linpred2)
M_bin1 = rbinom(n,1,probm1)
M_bin2 = rbinom(n,1,probm2)
M_cont1 = rnorm(n,linpred1,0.1)
M_cont2 = rnorm(n,linpred2,0.1)

#generate Cont Y, Bin Y, Count Y without int
linpred1 = theta0+theta1*A+theta2*M_cont1+theta3*M_cont2+theta4*A*M_cont1+theta5*A*M_cont2+
  theta6*M_cont1*M_cont2+theta7*A*M_cont1*M_cont2+theta8*C
linpred2 = theta0+theta1*A+theta2*M_bin1+theta3*M_cont1+theta4*A*M_bin1+theta5*A*M_cont1+
  theta6*M_bin1*M_cont1+theta7*A*M_bin1*M_cont1+theta8*C
linpred3 = theta0+theta1*A+theta2*M_bin1+theta3*M_bin2+theta4*A*M_bin1+theta5*A*M_bin2+
  theta6*M_bin1*M_bin2+theta7*A*M_bin1*M_bin2+theta8*C
contY_2contM_mint = rnorm(n,linpred1,0.2)
contY_binMcontM_mint = rnorm(n,linpred2,0.2)
contY_2binM_mint = rnorm(n,linpred3,0.2)
binY_2contM_mint = rbinom(n,1,expit(linpred1))
binY_binMcontM_mint = rbinom(n,1,expit(linpred2))
binY_2binM_mint = rbinom(n,1,expit(linpred3))

multipleM_mint_data = cbind(A,M_cont1,M_cont2,M_bin1,M_bin2,contY_2contM_mint,
                             contY_binMcontM_mint,contY_2binM_mint,binY_2contM_mint,
                             binY_binMcontM_mint,binY_binMcontM_mint,binY_2binM_mint,C)

write.table(multipleM_mint_data,file="multipleM_mint_data_10000.txt")
