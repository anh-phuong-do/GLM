######################################################
# Work: Project GLM
# Name: Phuong Do
# Date: 28/4/23
# Question 2: Model Baseline and Adjacent
#######################################################

# Not cluster model, proportional odd and continuation
install.packages("tidyverse")
library(tidyverse )
#Import data
getwd()
setwd("/Users/macbookair/Documents/UHasselt:Master1:Sem2/GLM/Part 2/")
EG_data <- read.table("EG.dat", sep = ' ',header = T)

#Count the number of alive, malformed and dead for each litter
alive <- EG_data %>% count(id, response) %>% filter(response == 1) # all alive in each clusters
malformed <- EG_data %>% count(id, response) %>% filter(response == 2)
dead <-  EG_data %>% count(id, response) %>% filter(response == 3)
dose <- EG_data %>% count(id, dose)

#Merge the data, remove columns, remane and substitute the NA values with 0 values
data1 <- merge(dose, alive[, -c(2)], by = 'id', all.x = T)
data1 <- merge(data1, malformed[, -c(2)], by = 'id', all.x = T)
data1<- merge(data1, dead[, -c(2)], by = 'id', all.x = T)
colnames(data1) <- c('id', 'dose', 'size','alive', 'malformed', 'dead')

data1[is.na(data1)] <- 0

EG <- data1

# Continuation ratio model
install.packages("VGAM")
library(VGAM)
fit1 = vglm(cbind(alive, malformed, dead)~dose, cratio(parallel = TRUE), EG)
summary(fit1)
# proportional odd model
fit2 = vglm(cbind(alive, malformed, dead)~dose, propodds(reverse = F), EG)
summary(fit2)
install.packages("VGLM")





# Ordered Multivariate GEE

## with multgee
library(multgee)
#fit the Multivariate extension of GEE with independence working assumption --> not clustering then 
PO.gee.ind <- ordLORgee(response ~ as.factor(dose), data = EG_data, id = id, LORstr = 'independence', link = 'logit')
summary(PO.gee.ind) 
#Fit the MUltivariate extension of the GEE with the exchangable working assumption --> all littermates with the same correlation: why should it be different?
PO.gee.unif <- ordLORgee((response) ~ as.factor(dose), data = EG_data, id = id, LORstr = 'uniform', link = 'logit')
summary(PO.gee.unif)

##QIC tets for GEE --> different objects needed 
library(MuMIn)
QIC(PO.gee.unif) #QIC based on independence assumption comaprison --> a fitted model object of class "gee", "geepack", "geem", "wgee", or "yags".

library(wgeesel)
QIC.gee(PO.gee.unif) #no., it needs a wgee model 


## with geepack --> different resutls? --> DO NOT USE GEEPACK 
library(geepack)

# Independence working assumption
model_ind <- ordgee(ordered(response) ~ dose, data = EG_data, mean.link = 'logit', id = id, corstr = "independence")
summary(model_ind)

# Echangeable working assumption 
model_exch <- ordgee(ordered(response) ~ dose, data = EG_data, mean.link = 'logit', id = id, corstr = "exchangeable")
summary(model_exch) #if dose put as factor then NaN in the analysis --> too much complicated variance structure??



# Ordered Multivariate GLMM 

## 'No linear' dose effect 
EG_GLMM <- clmm(ordered(response) ~ as.factor(dose) + (1|id), link = 'logit', data = EG_data)
summary(EG_GLMM)

## 'Linear' dose effect 
Equ_GLMM <- clmm(ordered(response) ~ dose + (1|id), data = EG_data, link = 'logit')
summary(Equ_GLMM)

## LIkelihood ratio test --> dof = 2 (2 parms in reduced, 4 in full)
LRT.GLMM <- 2*(EG_GLMM$logLik - Equ_GLMM$logLik)
c(LRT.GLMM, 1 - pchisq(LRT.GLMM, 2)) #H0: linear increase due to dose  H1: not linear increase --> significant 
