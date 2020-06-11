rm(list=ls())

df_int <- read.csv("./Data Simulation and Code Checking/int_data_10000.txt", sep = " ")

df_int$M_cat_catA=as.factor(df_int$M_cat_catA)
df_int$A_cat=as.factor(df_int$A_cat)

# multinom()
reg <- nnet::multinom(A_cat~C_cont, data = df_int)
predict(reg, newdata = df_int, type = "probs") ##columns=#levels

reg <- nnet::multinom(A_bin~C_cont, data = df_int)
predict(reg, newdata = df_int, type = "probs") #vector of prob of success

# gam(family = multinom())

reg <- gam(list(A_cat~s(C_cont),~s(C_cont)),family=multinom(K=2),data=df_int)
predict(reg, newdata = df_int, type = "response") ##columns=#levels

reg <- gam(list(A_bin~s(C_cont)),family=multinom(K=1),data=df_int)
predict(reg, newdata = df_int, type = "response") ##columns=#levels

# gam(family = ocat())
library(mgcv)
data <- gamSim(1,n=10000)
data$y=df_int$A_cat+1
reg <- gam(y~s(x1),family=ocat(R=3),data=data)
predict(reg, newdata = data, type = "response") ##columns=#levels

# polr()
reg <- polr(factor(A_cat)~C_cont,data=df_int)
predict(reg, newdata = df_int, type = "probs") ##columns=#levels
