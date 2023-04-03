# Introducción a los modelos de mezcla finitos
# Juan Sosa
# 2023

# Ejemplo: Datos sintéticos

# ajustes
rm(list = ls())
setwd("C:/Users/Juan Camilo/Dropbox/PAPERS/projects/BAYES")
source("finite mixture models - muestreador.R")
source("finite mixture models - em.R")

# parámetros de la simulación
n <- 100
H <- 2
theta <- c(0, 4)
sig2  <- c(1, 0.75) 
omega <- c(2/3, 1/3)

# simulación de datos
set.seed(123)
xi <- sample(x = 1:H, size = n, replace = TRUE, prob = omega)
y  <- rnorm(n = n, mean = theta[xi], sd = sqrt(sig2[xi]))

# previa
H      <- 2
mu0    <- 0
gam02  <- 1
nu0    <- 1
sig02  <- 1
alpha0 <- rep(1/H, H)

# ajuste del modelo usando el muestreador de Gibbs
# n_sams <- 5000
# n_burn <- 2500
# n_skip <- 10
# set.seed(1)
# tictoc::tic()
# muestras <- mcmc(y, H, mu0, gam02, nu0, sig02, alpha0, n_sams, n_burn, n_skip)
# tictoc::toc()
# save(muestras, file = "finite mixture models muestras ejemplo 1.RData")

# ajuste del modelo usando el algoritmo EM
# tol <- 0.0001
# set.seed(1)
# modas <- em(y, H, mu0, gam02, nu0, sig02, alpha0, tol, verbose = FALSE)
# save(modas, file = "finite mixture models modas ejemplo 1.RData")