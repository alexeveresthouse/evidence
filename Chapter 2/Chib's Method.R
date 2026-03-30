library(invgamma)
library(rstan)
library(matrixStats)
library(tidyverse)

set.seed(2026)

data <- c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)
datamean <- mean(data)

n <- length(data)

kappa <- 0.001

gibbs_complex<-function(N, a, b, kappa = 0.001, pmean = 0) 
{
  mat <- matrix(ncol = 2, nrow = N)
  mu <- 0
  sigma_sq <- 1
  mat[1, ] <- c(mu, sigma_sq)
  for (i in 2:N) 
  {
    mu <- rnorm(1, (n*datamean)/(n + kappa), sqrt(sigma_sq/(n + kappa)))
    sigma_sq <- rinvgamma(1, a + n/2 + 0.5, b + (sum((data - mu)^2) + (kappa*(mu - pmean)^2))/2)
    mat[i, ] <- c(mu, sigma_sq)
  }
  return(mat)
}

logMeanExp <- function(logx) {
  m <- max(logx)
  m + log(mean(exp(logx - m)))
}

mupost <- (n * datamean)/(n + kappa)

alpha <- a + n/2
beta <- b + (sum((data - mupost)^2) + kappa*mupost^2)/2

sigsqpost <- beta/(alpha + 1)

gibbsoutsimple <- gibbs_complex(50000, a = 1, b = 1)[-(1:(0.1*50000)), ] #  Eliminate "burn-in"
mugibbs <- gibbsoutsimple[, 1]
mugibbsmean <- mean(mugibbs)
sigsqgibbsmean <- mean(gibbsoutsimple[, 2])

logprior <- dnorm(mugibbsmean, 0, sqrt(sigsqgibbsmean/kappa), log = TRUE) + dinvgamma(sigsqgibbsmean, 1, 1, log = TRUE)
loglike <- sum(dnorm(data, mugibbsmean, sqrt(sigsqgibbsmean), log = TRUE))
logpost1 <- dnorm(mugibbsmean, (n*datamean)/(n + kappa), sqrt((sigsqgibbsmean)/(n + kappa)), log = TRUE)

sq_resids <- sapply(mugibbs,
  function(mu) sum((data - mu)^2))

log_post2_vals <- dinvgamma(sigsqgibbsmean, a + n/2 + 0.5, b + 0.5 * (sq_resids + kappa * mugibbs^2), log = TRUE)

logpost2 <- logMeanExp(log_post2_vals)

logpost <- logpost1 + logpost2

chib1 <- loglike + logprior - logpost # Log-scale calculation

# Analytic result with kappa term

kappa <- 0.001
biglogkappa <- sum(data^2)/2 - (n*mean(data))^2/(2*n + 2*kappa) + b
logevcomplextruekappa <- -0.5*n*log(2*pi) + 0.5*log(kappa) - 0.5*log(n + kappa) + a*log(b) - log(gamma(a)) + log(gamma(0.5*n + a)) - (0.5*n + a)*log(biglogkappa)

# Systematic analysis

data <- c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)
datamean <- mean(data)

n <- length(data)
pmean <- 0

niter <- c(100, 1000, 10000, 100000)
set.seed(2026)
a <- 1
b <- 1
kappa <- 0.001
nruns <- 100
chib <- matrix(nrow = nruns, ncol = length(niter))

for (k in 1:nruns){
  for (i in 1:length(niter)){
    gibbsoutsimple <- gibbs_complex(110000, a = 1, b = 1)[(110000 - niter[i]):110000, ]
    mugibbs <- gibbsoutsimple[, 1]
    mugibbsmean <- mean(mugibbs)
    sigsqgibbsmean <- mean(gibbsoutsimple[, 2])
    
    mupost <- (n * datamean)/(n + kappa)
    
    alpha <- a + n/2
    beta <- b + (sum((data - mupost)^2) + kappa*mupost^2)/2
    
    sigsqpost <- beta/(alpha + 1)
    
    logprior <- dnorm(mugibbsmean, 0, sqrt(sigsqgibbsmean/kappa), log = TRUE) + dinvgamma(sigsqgibbsmean, 1, 1, log = TRUE)
    loglike <- sum(dnorm(data, mugibbsmean, sqrt(sigsqgibbsmean), log = TRUE))
    logpost1 <- dnorm(mugibbsmean, (n*datamean)/(n + kappa), sqrt((sigsqgibbsmean)/(n + kappa)), log = TRUE)
    
    sq_resids <- sapply(mugibbs,
                        function(mu) sum((data - mu)^2))
    
    log_post2_vals <- dinvgamma(sigsqgibbsmean, a + n/2 + 0.5, b + 0.5 * (sq_resids + kappa * mugibbs^2), log = TRUE)
    
    logpost2 <- logMeanExp(log_post2_vals)
    
    logpost <- logpost1 + logpost2
    
    chib[k, i] <- loglike + logprior - logpost # Log-scale calculation
  }
}

colnames(chib) <- c("100", "1000", "10000", "100000")

df.chib <- as.data.frame(chib)
df.chib <- pivot_longer(df.chib, cols = everything(), names_to = "Number of samples", values_to = "Evidence")

palatinate <- rgb(104, 36, 109, maxColorValue = 255)

ggplot(df.chib, aes(x = `Number of samples`, y = Evidence, fill = "Evidence")) +
  geom_boxplot() +
  ylim(c(-12.3, -12.1)) +
  ylab("Log evidence") +
  stat_summary(
    geom = "crossbar", 
    width = 0.75, 
    fatten = 0.5, 
    color = "white", 
    fun.data = function(x){ 
      return(c(y = median(x), ymin = median(x), ymax = median(x))) 
    }
  ) +
  scale_fill_manual(values = palatinate) +
  theme_minimal() +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 16)) +
  geom_abline(intercept = logevcomplextruekappa, slope = 0, col = "red")
