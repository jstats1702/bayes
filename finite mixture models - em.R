# Introducción a los modelos de mezcla finitos
# Juan Sosa
# 2023

# Algoritmo EM
# Mezcla univariada de Normales con parámetros de ubicación y escala específicos

# funciones auxiliares
#
dinvgamma <- function(x, a, b, log = FALSE) 
{
  # función de densidad distribución gamma inversa
  out <- a*log(b) - lgamma(a) - (a+1)*log(x) - b/x 
  if (log == FALSE) out <- exp(out)
  return(out)
}
#
logposphi <- function(y, theta, sig2, omega, mu0, gam02, nu0, sig02)
{
  # log-verosimilitud sin variables latentes
  out <- 0
  for (i in 1:length(y))
    out <- out + log(sum(omega*dnorm(x = y[i], mean = theta, sd = sqrt(sig2))))
  for (h in 1:length(theta))
    out <- out + dnorm(x = theta[h], mean = mu0, sd = sqrt(gam02), log = TRUE) + dinvgamma(x = sig2[h], a = nu0/2, b = nu0*sig02/2, log = TRUE)
  return(out)
}






# algoritmo EM
em <- function(y, H, mu0, gam02, nu0, sig02, alpha0, tol, verbose = TRUE)
{
  # ajustes
  y  <- as.numeric(y)
  yb <- mean(y)
  sy <- sd(y) 
  y  <- scale(y)
  n  <- length(y)
  # valores iniciales
  # omega <- c(gtools::rdirichlet(n = 1, alpha = alpha0))
  # theta <- rnorm(n = H, mean = mu0, sd = sqrt(gam02))
  # sig2  <- 1/rgamma(n = H, shape = nu0/2, rate = nu0*sig02/2)
  tmp   <- stats::kmeans(x = y, centers = H)
  omega <- table(tmp$cluster)/n
  theta <- rep(NA, H)
  sig2  <- rep(NA, H)
  nh    <- as.numeric(table(factor(x = tmp$cluster, levels = 1:H)))
  for (h in 1:H) {
    if (nh[h] > 0) {
      theta[h] <- mean(y[tmp$cluster == h])
      sig2 [h] <- var (y[tmp$cluster == h])
    }
  }
  # algoritmo
  lp   <- NULL
  b    <- 1
  cond <- TRUE
  while (cond) {
    # actualizar probabilidades    
    pi <- NULL
    xi <- NULL
    for (i in 1:n) {
      tmp <- NULL
      for (h in 1:H)
        tmp[h] <- omega[h]*dnorm(x = y[i], mean = theta[h], sd = sqrt(sig2[h]))
      tmp <- tmp/sum(tmp)
      pi <- rbind(pi, tmp)
      xi <- c(xi, which.max(tmp))
      rm(tmp)
    }
    # actualizar parámetros de la mezcla
    for (h in 1:H) {
      theta[h] <- (mu0/gam02 + sum(pi[,h]*y)/sig2[h])/(1/gam02 + sum(pi[,h])/sig2[h])
      sig2 [h] <- (nu0*sig02 + sum(pi[,h]*(y - theta[h])^2))/(nu0 + sum(pi[,h]) + 2)
      omega[h] <- (alpha0[h] + sum(pi[,h]) - 1)/(sum(alpha0) + sum(pi) - H)
    }
    # log-verosimilitud sin variables latentes
    lp[b] <- logposphi(y, theta, sig2, omega, mu0, gam02, nu0, sig02)
    # progreso
    if (verbose) cat("Iteración ", b, " -> log p = ", round(lp[b], 5), "\n", sep = "")
    # convergencia
    if(b > 1) 
      if (abs(lp[b] - lp[b-1]) < tol) 
        cond <- FALSE
    b <- b+1
  }
  return(list(pi = pi, xi = xi, theta = sy*theta + yb, sig2 = sy^2*sig2, omega = omega, lp = lp))
}