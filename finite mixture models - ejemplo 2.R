# Introducción a los modelos de mezcla finitos
# Juan Sosa
# 2023

# Ejemplo: Datos sintéticos

# ajustes
rm(list = ls())
setwd("C:/Users/Juan Camilo/Dropbox/PAPERS/projects/BAYES")
source("finite mixture models - muestreador.R")

# datos
data(galaxy, package = "rebmix")
y <- galaxy$Velocity
n <- length(y)

# previa
H      <- 5
mu0    <- 0
gam02  <- 1
nu0    <- 1
sig02  <- 1
alpha0 <- rep(1/H, H)

# ajuste del modelo usando el muestreador de Gibbs
n_sams <- 5000
n_burn <- 2500
n_skip <- 10
set.seed(1)
tictoc::tic()
muestras <- mcmc(y, H, mu0, gam02, nu0, sig02, alpha0, n_sams, n_burn, n_skip)
tictoc::toc()
save(muestras, file = "finite mixture models muestras ejemplo 2.RData")