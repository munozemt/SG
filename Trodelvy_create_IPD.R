### Algorithm to create a raw dataset from DigizeIt readings from a Kaplan-Meier curve ###

# Paquetes
pacman::p_load(MASS, splines, survival, IPDfromKM)

rm(list=ls()) # Limpiamos el entorno

# Cargar puntos Kaplan-Meier (tiempo, probabilidad)
OS_C <- read.csv("https://raw.githubusercontent.com/munozemt/SG/refs/heads/main/OS_chemo_no_met_A.csv", header = TRUE)
OS_SG <- read.csv("https://raw.githubusercontent.com/munozemt/SG/refs/heads/main/OS_SG_no_met_A.csv", header = TRUE)
PFS_C <- read.csv("https://raw.githubusercontent.com/munozemt/SG/refs/heads/main/PFS_chemo_no_met_A.csv", header = TRUE)
PFS_SG <- read.csv("https://raw.githubusercontent.com/munozemt/SG/refs/heads/main/PFS_SG_no_met_A.csv", header = TRUE)

# Cargar tablas resumen number at risk
atrisk_OS_C <- read.csv("https://raw.githubusercontent.com/munozemt/SG/refs/heads/main/OS_chemo_at_risk_A.csv", header=TRUE)
atrisk_OS_SG <- read.csv("https://raw.githubusercontent.com/munozemt/SG/refs/heads/main/OS_SG_at_risk_A.csv", header=TRUE)
atrisk_PFS_C <- read.csv("https://raw.githubusercontent.com/munozemt/SG/refs/heads/main/PFS_chemo_at_risk_A.csv", header=TRUE)
atrisk_PFS_SG <- read.csv("https://raw.githubusercontent.com/munozemt/SG/refs/heads/main/PFS_SG_at_risk_A.csv", header=TRUE)

# Objetos para alimentar a los comandos IPD
dat_OS_C <- OS_C[,2:3] # Solo columnas 2 y 3
trisk_OS_C <- atrisk_OS_C[1:8,2] # Times at risk
nrisk_OS_C <- atrisk_OS_C[1:8,5] # Number at risk

dat_OS_SG <- OS_SG[,2:3] # Solo columnas 2 y 3
trisk_OS_SG <- atrisk_OS_SG[1:8,2] # Times at risk
nrisk_OS_SG <- atrisk_OS_SG[1:8,5] # Number at risk

dat_PFS_C <- PFS_C[,2:3] # Solo columnas 2 y 3
trisk_PFS_C <- atrisk_PFS_C[1:8,2] # Times at risk
nrisk_PFS_C<- atrisk_PFS_C[1:8,5] # Number at risk

dat_PFS_SG <- PFS_SG[,2:3] # Solo columnas 2 y 3
trisk_PFS_SG <- atrisk_PFS_SG[1:8,2] # Times at risk
nrisk_PFS_SG <- atrisk_PFS_SG[1:8,5] # Number at risk

# Pre-procesar datos (maxy = 1, indica probabilidad expresada en el intervalo 0,1)
pre_OS_C <- preprocess(dat = dat_OS_C, trisk = trisk_OS_C, nrisk = nrisk_OS_C, maxy = 1)
pre_OS_SG <- preprocess(dat = dat_OS_SG, trisk = trisk_OS_SG, nrisk = nrisk_OS_SG, maxy = 1)
pre_PFS_C <- preprocess(dat = dat_PFS_C, trisk = trisk_PFS_C, nrisk = nrisk_PFS_C, maxy = 1)
pre_PFS_SG <- preprocess(dat = dat_PFS_SG, trisk = trisk_PFS_SG, nrisk = nrisk_PFS_SG, maxy = 1)

# Genera pseudo datos individuales, armID = 0 Chemo, armID = 1 Sacituzumab
ipd_OS_C <- getIPD(prep = pre_OS_C, armID = 1, tot.events = NULL)
ipd_OS_SG <- getIPD(prep = pre_OS_SG, armID = 1, tot.events = NULL)
ipd_PFS_C <- getIPD(prep = pre_PFS_C, armID = 1, tot.events = NULL)
ipd_PFS_SG <- getIPD(prep = pre_PFS_SG, armID = 1, tot.events = NULL)

# Checar si se rechaza o no la prueba Kolmogorov-Smirnov
summary(ipd_OS_C)
summary(ipd_OS_SG)
summary(ipd_PFS_C)
summary(ipd_PFS_SG)

# Inspecciona resultados
plot(ipd_OS_C)
plot(ipd_OS_SG)
plot(ipd_PFS_C)
plot(ipd_PFS_SG)

# Regresión de supervivencia, especificando distribución Weibull, debemos hacerlo para otras distribuciones (log-normal, log-logistic, gamma, exponential)
fit.weibull_OS_C <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_OS_C$IPD, dist = "weibull")
fit.lognormal_OS_C <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_OS_C$IPD, dist = "lnorm")
fit.loglogistic_OS_C <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_OS_C$IPD, dist = "llogis")
fit.gamma_OS_C <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_OS_C$IPD, dist = "gamma")
fit.exponential_OS_C <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_OS_C$IPD, dist = "exponential")

fit.weibull_OS_SG <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_OS_SG$IPD, dist = "weibull")
fit.lognormal_OS_SG <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_OS_SG$IPD, dist = "lnorm")
fit.loglogistic_OS_SG <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_OS_SG$IPD, dist = "llogis")
fit.gamma_OS_SG <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_OS_SG$IPD, dist = "gamma")
fit.exponential_OS_SG <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_OS_SG$IPD, dist = "exponential")

fit.weibull_PFS_C <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_PFS_C$IPD, dist = "weibull")
fit.lognormal_PFS_C <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_PFS_C$IPD, dist = "lnorm")
fit.loglogistic_PFS_C <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_PFS_C$IPD, dist = "llogis")
fit.gamma_PFS_C <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_PFS_C$IPD, dist = "gamma")
fit.exponential_PFS_C <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_PFS_C$IPD, dist = "exponential")

fit.weibull_PFS_SG <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_PFS_SG$IPD, dist = "weibull")
fit.lognormal_PFS_SG <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_PFS_SG$IPD, dist = "lnorm")
fit.loglogistic_PFS_SG <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_PFS_SG$IPD, dist = "llogis")
fit.gamma_PFS_SG <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_PFS_SG$IPD, dist = "gamma")
fit.exponential_PFS_SG <- flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd_PFS_SG$IPD, dist = "exponential")

aic_values_OS_C <- c(
  weibull = AIC(fit.weibull_OS_C),
  lognormal = AIC(fit.lognormal_OS_C),
  loglogistic = AIC(fit.loglogistic_OS_C),
  gamma = AIC(fit.gamma_OS_C),
  exponential = AIC(fit.exponential_OS_C)
)

aic_values_OS_SG <- c(
  weibull = AIC(fit.weibull_OS_SG),
  lognormal = AIC(fit.lognormal_OS_SG),
  loglogistic = AIC(fit.loglogistic_OS_SG),
  gamma = AIC(fit.gamma_OS_SG),
  exponential = AIC(fit.exponential_OS_SG)
)

aic_values_PFS_C <- c(
  weibull = AIC(fit.weibull_PFS_C),
  lognormal = AIC(fit.lognormal_PFS_C),
  loglogistic = AIC(fit.loglogistic_PFS_C),
  gamma = AIC(fit.gamma_PFS_C),
  exponential = AIC(fit.exponential_PFS_C)
)

aic_values_PFS_SG <- c(
  weibull = AIC(fit.weibull_PFS_SG),
  lognormal = AIC(fit.lognormal_PFS_SG),
  loglogistic = AIC(fit.loglogistic_PFS_SG),
  gamma = AIC(fit.gamma_PFS_SG),
  exponential = AIC(fit.exponential_PFS_SG)
)

# Imprimir AIC para cada modelo

print(aic_values_OS_C)
print(aic_values_OS_SG)
print(aic_values_PFS_C)
print(aic_values_PFS_SG)



