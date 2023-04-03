# Introducción a los modelos de mezcla finitos
# Juan Sosa
# 2023

# Muestreador de Gibbs
# Mezcla univariada de Normales con parámetros de ubicación y escala específicos

# distribuciones condicionales completas
#
sample_theta <- function (nh, ybh, H, mu0, gam02, theta, sig2) 
{
  for (h in 1:H) {
    if (nh[h] == 0) {
      theta[h] <- rnorm(n = 1, mean = mu0, sd = sqrt(gam02))
    } else {
      v2 <- 1/(1/gam02 + nh[h]/sig2[h])
      m  <- v2*(mu0/gam02 + nh[h]*ybh[h]/sig2[h])
      theta[h] <- rnorm(n = 1, mean = m, sd = sqrt(v2))
    }
  }
  return(theta)
}
#
sample_sig2 <- function (nh, ybh, ssh, H, nu0, sig02, theta, sig2)
{
  for (h in 1:H) {
    if (nh[h] == 0) {
      sig2[h] <- 1/rgamma(n = 1, shape = 0.5*nu0, rate = 0.5*nu0*sig02)
    } else {
      a <- 0.5*(nu0 + nh[h])
      b <- 0.5*(nu0*sig02 + ssh[h] + nh[h]*(ybh[h] - theta[h])^2)
      sig2[h] <- 1/rgamma(n = 1, shape = a, rate = b)
    }
  }
  return(sig2)
}
#
sample_xi <- function (H, omega, theta, sig2, xi, n, y) 
{
  for (i in 1:n) {
    lp <- rep(NA, H)
    for (h in 1:H) {
      lp[h] <- log(omega[h]) + dnorm(x = y[i], mean = theta[h], sd = sqrt(sig2[h]), log = TRUE)
    }
    xi[i] <- sample(x = 1:H, size = 1, replace = FALSE, prob = exp(lp - max(lp)))
  }
  return(xi)
}
#
sample_omega <- function (nh, alpha0)
{
  return(c(gtools::rdirichlet(n = 1, alpha = alpha0 + nh)))
}

# muestreador de Gibbs
mcmc <- function (y, H, mu0, gam02, nu0, sig02, alpha0, n_sams, n_burn, n_skip, verbose = TRUE) 
{
  # ajustes
  y  <- as.numeric(y)
  yb <- mean(y)
  sy <- sd(y) 
  y  <- scale(y)
  n  <- length(y)
  # número de iteraciones
  B <- n_burn + n_sams*n_skip
  ncat <- floor(0.1*B)
  # valores iniciales
  tmp   <- stats::kmeans(x = y, centers = H)
  xi    <- tmp$cluster
  omega <- as.numeric(table(factor(x = xi, levels = 1:H)))/H
  theta <- rnorm(n = H, mean = mu0, sd = sqrt(gam02))
  sig2  <- 1/rgamma(n = H, shape = nu0/2, rate = nu0*sig02/2)
  # almacenamiento
  THETA <- matrix(data = NA, nrow = n_sams, ncol = H) 
  SIG2  <- matrix(data = NA, nrow = n_sams, ncol = H)
  OMEGA <- matrix(data = NA, nrow = n_sams, ncol = H)
  XI    <- matrix(data = NA, nrow = n_sams, ncol = n)
  LP    <- matrix(data = NA, nrow = n_sams, ncol = 1)
  # cadena
  for (i in 1:B) {
    # actualizar estadísticos suficientes
    nh <- as.numeric(table(factor(x = xi, levels = 1:H)))
    ybh <- ssh <- rep(NA, H)
    for (h in 1:H) {
      if (nh[h] > 0) {
        indexh <- xi == h
        ybh[h] <- mean(y[indexh])
        ssh[h] <- sum((y[indexh] - ybh[h])^2)
      }
    }
    # actualizar parámetros
    sig2  <- sample_sig2  (nh, ybh, ssh, H, nu0, sig02, theta, sig2)
    theta <- sample_theta (nh, ybh, H, mu0, gam02, theta, sig2)
    omega <- sample_omega (nh, alpha0)
    xi    <- sample_xi    (H, omega, theta, sig2, xi, n, y)
    # almacenar y log-verosimilitud
    if (i > n_burn) {
      if (i%%n_skip == 0) {
        k <- (i - n_burn)/n_skip
        THETA[k,] <- sy*theta + yb
        SIG2 [k,] <- sy^2*sig2
        OMEGA[k,] <- omega
        XI   [k,] <- xi
        LP   [k,] <- sum(dnorm(x = y, mean = theta[xi], sd = sqrt(sig2[xi]), log = TRUE))
      }
    }
    # progreso
    if (verbose) if (i%%ncat == 0) cat(100*round(i/B, 1), "% completado \n", sep = "")
  }
  # salida
  return(list(THETA = THETA, SIG2 = SIG2, OMEGA = OMEGA, XI = XI, LP = LP))
}