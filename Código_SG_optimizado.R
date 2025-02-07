### Cargar paquetes
pacman::p_load(MASS, splines, survival, flexsurv, IPDfromKM)

# Limpiar el entorno
rm(list = ls())

# Base URL para los archivos en GitHub
base_url <- "https://raw.githubusercontent.com/munozemt/SG/refs/heads/main/"

#### Función para cargar y preprocesar datos a partir de archivos del KM ####
processKM <- function(curve, arm) {
  # curve: "OS" o "PFS"
  # arm: "C" (chemo) o "SG" (Sacituzumab)
  
  # Determinar nombres de archivo basados en la combinación
  file_km <- paste0(curve, "_", ifelse(arm == "C", "chemo", "SG"), "_no_met_A.csv")
  file_atrisk <- paste0(curve, "_", ifelse(arm == "C", "chemo", "SG"), "_at_risk_A.csv")
  
  # Leer datos
  km_data <- read.csv(paste0(base_url, file_km), header = TRUE)
  atrisk_data <- read.csv(paste0(base_url, file_atrisk), header = TRUE)
  
  # Seleccionar columnas y filas correspondientes
  dat <- km_data[, 2:3]               # Se toman las columnas 2 y 3
  trisk <- atrisk_data[1:8, 2]          # Tiempo en riesgo (primeras 8 filas, columna 2)
  nrisk <- atrisk_data[1:8, 5]          # Número en riesgo (primeras 8 filas, columna 5)
  
  # Preprocesamiento de datos: maxy = 1 indica que la probabilidad está en el rango 0 a 1
  prep <- preprocess(dat = dat, trisk = trisk, nrisk = nrisk, maxy = 1)
  ipd <- getIPD(prep = prep, armID = 1, tot.events = NULL)
  
  return(ipd)
}

#### Función para ajustar modelos de supervivencia y obtener valores de AIC ####
fit_models <- function(ipd) {
  # Definir distribuciones a ajustar
  dists <- c("weibull", "lnorm", "llogis", "gamma", "exponential")
  # Ajusta modelos usando lapply
  fits <- lapply(dists, function(d) {
    flexsurvreg(formula = Surv(time, status) ~ 1, data = ipd$IPD, dist = d)
  })
  names(fits) <- dists
  # Extraer AICs
  aic_values <- sapply(fits, AIC)
  return(aic_values)
}

#### Lista con las combinaciones de curvas y brazos ####
curves <- list(
  "OS_C"    = list(curve = "OS",  arm = "C"),
  "OS_SG"   = list(curve = "OS",  arm = "SG"),
  "PFS_C"   = list(curve = "PFS", arm = "C"),
  "PFS_SG"  = list(curve = "PFS", arm = "SG")
)

#### Procesar datos y ajustar modelos para cada conjunto ####
# Inicializar listas para almacenar los IPD y los AIC
ipd_list <- list()
aic_list <- list()

for (name in names(curves)) {
  # Procesar datos
  ipd_list[[name]] <- processKM(curve = curves[[name]]$curve, arm = curves[[name]]$arm)
  # Mostrar sumario y gráfico para inspección
  print(summary(ipd_list[[name]]))
  plot(ipd_list[[name]], main = paste("Curva", name))
  
  # Ajuste de modelos y extracción de AIC
  aic_list[[name]] <- fit_models(ipd_list[[name]])
}

#### Imprimir valores de AIC para cada conjunto de datos ####
for (name in names(aic_list)) {
  cat("AIC para", name, ":\n")
  print(aic_list[[name]])
  cat("\n")
}

