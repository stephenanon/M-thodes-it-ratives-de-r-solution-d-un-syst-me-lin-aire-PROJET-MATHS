################### DEVOIR ANON NDA GUY SIMON STEPHANE ################################

######### PACKAGES

library("psych")
library("DescTools")
library("questionr")
library("labelled")
library("mfx")
library("stargazer")
library("robustbase")
library(ResourceSelection)
library("pROC")
library("caret")
library("FactoMineR") 
library("factoextra")
library("questionr")
library("forcats")
library("labelled")
library("dplyr")
library("plyr")
library(MASS) 
library(ggplot2) 

#### Description bréve de la base
BASEGLM = BASE_ind

names(BASEGLM)

lookfor(BASEGLM,details = TRUE)

######### TRAITEMENT DES VARIABLES

## Labeliser toutes les varibales

var_label(BASEGLM$milieu) = "Milieu de résidence"
var_label(BASEGLM$sexe) = "Sexe"
var_label(BASEGLM$age) = "Age"
var_label(BASEGLM$nbhandic) = "Nb handicaps"


## Labeliser toutes les variables qualitatives et choisir la modalité de référence

BASEGLM$sexe = as.factor(BASEGLM$sexe)
levels(BASEGLM$sexe) = c("Masculin","Féminin")

BASEGLM$milieu = as.factor(BASEGLM$milieu)
levels(BASEGLM$milieu) = c("Urbain","Rural")

BASEGLM$education = as.factor(BASEGLM$education)
levels(BASEGLM$education) = c("Aucun","Primaire","Secondaire","Supérieur")


######### ESTIMATION

## Choix de la loi
## Choix de la fonction de lien

# Poisson
model1 = glm(nbhandic ~ sexe + milieu + age + education, data = BASEGLM, family = poisson)

summary(model1)

# Binomiale négative
model2 = glm.nb(nbhandic ~ sexe + milieu + age + education, data = BASEGLM)

summary(model2)

######### DIAGNOSTIC DES RESIDUS ET RE-ESTIMATION

BASEGLM$residP1 = residuals.glm(model1, type = "pearson")
BASEGLM$residD1  = residuals.glm(model1, type = "deviance")
nobs = 1:length(BASEGLM$residD1 )

plot(nobs,BASEGLM$residP1,type = "p")
abline(a=5,b=0,col="red")
abline(a=-5,b=0,col="red")

model1m = glm(nbhandic ~ sexe + milieu + age + education, data = subset(BASEGLM, abs(residP1) < 5), family = poisson)

summary(model1m)

######### QUALITE DU MODELE

## Khi 2 de Pearson généralisé
model1m0 = glm(nbhandic ~ 1, data = subset(BASEGLM, abs(residP1) < 5), family = poisson)

anova(model1m,model1m0,test="Chisq")

## Pseudo R2

PseudoR2(model1m)

######### TESTS D'HYPOTHESE

## Significativité individuelle

model1m = glm(nbhandic ~ sexe + milieu + age + education, data = subset(BASEGLM, abs(residP1) < 5), family = poisson)
model1mr = glm(nbhandic ~ sexe + milieu + age, data = subset(BASEGLM, abs(residP1) < 5), family = poisson)
anova(model1m,model1mr,test="Chisq")

