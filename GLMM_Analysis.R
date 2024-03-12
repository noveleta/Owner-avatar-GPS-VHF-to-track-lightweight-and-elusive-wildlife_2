library(plyr)
library(rptR)
library(lattice)
library(lmerTest)
library(effects)
library(ggplot2)
library(sjPlot)
library(glmmTMB)
library(car)
library(performance)
library(usdm) # Collinearity
library(randomForestSRC) # RF
library(ggRandomForests) # RF
library(gbm) # BRT
library(dismo) # BRT
library(MuMIn) # Multi-model inference
library(DHARMa)
library(ggeffects)
library(splines)
library (MASS)

### establecer DIRECTORIO PARA LEER DATOS

setwd("C:/Users/Perico/Documents/RICOTÍ/ADEMUZ/ARTÍCULOS/Testing GPS")


### LO LEO DEL CSV
Datos_Gamma <- read.csv("Distancias_Gamma2.csv", sep=";")

### RESUMEN DE LOS DATOS E HISTOGRAMAS PARA VER LA DISTRIBUCION. PARECE POISSON O BINOMIAL NEGATIVA
str(Datos_Gamma)


### TRANSFORMO 'TAG', POSITION', 'SEPARATION' Y 'DATE' DE CARACTER A FACTOR
Datos_Gamma$Tag <- as.factor(Datos_Gamma$Tag)
Datos_Gamma$Position <- as.factor(Datos_Gamma$Position)
Datos_Gamma$Separation <- as.factor(Datos_Gamma$Separation)
Datos_Gamma$factor_aleatorio <- as.factor(Datos_Gamma$factor_aleatorio)


str(Datos_Gamma)
hist(Datos_Gamma$Distance, nclass = max(Datos_Gamma$Distance))
hist(Datos_Gamma$satellites, nclass = max(Datos_Gamma$satellites))

Datos_Gamma.sel <- subset(Datos_Gamma, Datos_Gamma$satellites > 4)
str(Datos_Gamma.sel)
summary (Datos_Gamma.sel)

### Exploro la correlación entre 'distance' y el valor'accuracy' del output del GPS 
cor (Datos_Gamma.sel[,13:14]) ##Para Pearson
cor (Datos_Gamma.sel[,13:14], method = "spearman") ##Para Spearman
pairs (Datos_Gamma.sel[,13:14], method = "spearman")
cor.test(Datos_Gamma.sel$Distance, Datos_Gamma.sel$Accuracy, method = "spearman")
cor.test(Datos_Gamma.sel$Distance, Datos_Gamma.sel$Accuracy, method = "pearson")



### ANALISIS DEL NUMERO DE SATELITES EN FUNCION DE LA POSICION (GAUSS DISTRIBUTION) ###############
modelD <- glmmTMB(satellites~Position + (1|Tag) + (1|factor_aleatorio), data=Datos_Gamma.sel)
summary(modelD)
confint(modelD)
hist(residuals(modelD))
shapiro.test(residuals(modelD))
plot(residuals(modelD))
plot(Datos_Gamma.sel$factor_aleatorio,residuals(modelD))
plot(Datos_Gamma.sel$satellites,residuals(modelD))
plot(Datos_Gamma.sel$Position,residuals(modelD))

r.squaredGLMM(modelD)


simulationOutput_D <- simulateResiduals(fittedModel = modelD, plot = T)
residuals(simulationOutput_D)
hist(simulationOutput_D)
testResiduals(simulationOutput_D)




plot_model(modelD, type = "pred", terms = c("Position"),
           ci.lvl = 0.95, 
           line.size = 1,
           title = "",
           axis.title=c("Type of arrangement","n? satelites (m)"),
           show.data=FALSE, dot.size = 2, alpha = 0.3)

pred.sat <- ggeffect(
  modelD,
  terms = c("satellites"),
  ci.lvl = 0.95,
  type = "re",
  typical = "mean",
  condition = NULL,
  back.transform = TRUE,
  interval = "confidence")

plot(pred.sat)

write.table(pred.sat,"pred.satellites.csv", sep=";" )

######### ANALISIS DE LA DISTANCIA PARA CADA NUMERO FIJO DE SATELITES ######### 

## prueba previa con todos los satelites juntos - no se incluira en el paper ####
FULL_Gamma <- glmmTMB(Distance ~ Position +  (1|Tag) + (1|factor_aleatorio), data=Datos_Gamma.sel, family=Gamma(link="log"))
summary(FULL_Gamma)

### PARA 5 SATELITES ###
Datos_Gamma.sat5 <- subset(Datos_Gamma, Datos_Gamma$satellites == 5)
FULL_Gamma.s5 <- glmmTMB(Distance ~ Position +  (1|Tag) + (1|factor_aleatorio), data=Datos_Gamma.sat5, family=Gamma(link="log"))
summary(FULL_Gamma.s5)
r.squaredGLMM(FULL_Gamma.s5)
pred.s5e <- ggeffect(
  FULL_Gamma.s5,
  terms = c("Position"),
  ci.lvl = 0.95,
  type = "re",
  typical = "mean",
  condition = NULL,
  back.transform = TRUE,
  interval = "confidence")

plot(pred.s5e)

write.table(pred.s5e,"5_SATELITES.csv", sep=";" )

### PARA 6 SATELITES ###
Datos_Gamma.sat6 <- subset(Datos_Gamma, Datos_Gamma$satellites == 6)
FULL_Gamma.s6 <- glmmTMB(Distance ~ Position +  (1|Tag) + (1|factor_aleatorio), data=Datos_Gamma.sat6, family=Gamma(link="log"))
summary(FULL_Gamma.s6)
r.squaredGLMM(FULL_Gamma.s6)
pred.s6e <- ggeffect(
  FULL_Gamma.s6,
  terms = c("Position"),
  ci.lvl = 0.95,
  type = "re",
  typical = "mean",
  condition = NULL,
  back.transform = TRUE,
  interval = "confidence")

plot(pred.s6e)

write.table(pred.s6e,"6_SATELITES.csv", sep=";" )


### PARA 7 SATELITES ###
Datos_Gamma.sat7 <- subset(Datos_Gamma, Datos_Gamma$satellites == 7)
FULL_Gamma.s7 <- glmmTMB(Distance ~ Position +  (1|Tag) + (1|factor_aleatorio), data=Datos_Gamma.sat7, family=Gamma(link="log"))
summary(FULL_Gamma.s7)
r.squaredGLMM(FULL_Gamma.s7)
pred.s7e <- ggeffect(
  FULL_Gamma.s7,
  terms = c("Position"),
  ci.lvl = 0.95,
  type = "re",
  typical = "mean",
  condition = NULL,
  back.transform = TRUE,
  interval = "confidence")

plot(pred.s7e)

write.table(pred.s7e,"7_SATELITES.csv", sep=";" )


### PARA 8 SATELITES ###
Datos_Gamma.sat8 <- subset(Datos_Gamma, Datos_Gamma$satellites == 8)
FULL_Gamma.s8 <- glmmTMB(Distance ~ Position +  (1|Tag) + (1|factor_aleatorio), data=Datos_Gamma.sat8, family=Gamma(link="log"))
summary(FULL_Gamma.s8)
r.squaredGLMM(FULL_Gamma.s8)
pred.s8e <- ggeffect(
  FULL_Gamma.s8,
  terms = c("Position"),
  ci.lvl = 0.95,
  type = "re",
  typical = "mean",
  condition = NULL,
  back.transform = TRUE,
  interval = "confidence")

plot(pred.s8e)

write.table(pred.s8e,"8_SATELITES.csv", sep=";" )


### PARA 9 SATELITES ###
Datos_Gamma.sat9 <- subset(Datos_Gamma, Datos_Gamma$satellites == 9)
FULL_Gamma.s9 <- glmmTMB(Distance ~ Position +  (1|Tag) + (1|factor_aleatorio), data=Datos_Gamma.sat9, family=Gamma(link="log"))
summary(FULL_Gamma.s9)
r.squaredGLMM(FULL_Gamma.s9)
pred.s9e <- ggeffect(
  FULL_Gamma.s9,
  terms = c("Position"),
  ci.lvl = 0.95,
  type = "re",
  typical = "mean",
  condition = NULL,
  back.transform = TRUE,
  interval = "confidence")

plot(pred.s9e)

write.table(pred.s9e,"9_SATELITES.csv", sep=";" )


### PARA 10 SATELITES ###
Datos_Gamma.sat10 <- subset(Datos_Gamma, Datos_Gamma$satellites == 10)
FULL_Gamma.s10 <- glmmTMB(Distance ~ Position +  (1|Tag) + (1|factor_aleatorio), data=Datos_Gamma.sat10, family=Gamma(link="log"))
summary(FULL_Gamma.s10)
r.squaredGLMM(FULL_Gamma.s10)
pred.s10e <- ggeffect(
  FULL_Gamma.s10,
  terms = c("Position"),
  ci.lvl = 0.95,
  type = "re",
  typical = "mean",
  condition = NULL,
  back.transform = TRUE,
  interval = "confidence")

plot(pred.s10e)

write.table(pred.s10e,"10_SATELITES.csv", sep=";" )

view (pred.s10e)

View(pred.s5e)

