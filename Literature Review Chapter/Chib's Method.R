library(invgamma)
library(rstan)

set.seed(2026)

data <- c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)
datamean <- mean(data)

n <- length(data)

gibbs_complex<-function(N, a, b) 
{
  mat <- matrix(ncol = 2, nrow = N)
  mu <- 0
  sigma_sq <- 1
  mat[1, ] <- c(mu, sigma_sq)
  for (i in 2:N) 
  {
    mu <- rnorm(1, (n*datamean)/(n + 1), sqrt(sigma_sq/(n + 1)))
    sigma_sq <- rinvgamma(1, (a + n + 1)/2, (b + sum((data - mu)^2) + mu^2)/2)
    mat[i, ] <- c(mu, sigma_sq)
  }
  return(mat)
}

postmeancheck <- stan("new complex model.stan",
                      data = list(n = 7,
                                  x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12),
                                  pmean = 0,
                                  a = 1,
                                  b = 1),
                      iter = 10000,
                      warmup = 5000)

mupost <- mean(extract(postmeancheck)$mu)
mugibbs <- extract(postmeancheck)$mu
sigsqpost <- mean(extract(postmeancheck)$sigma_sq) # posterior means here ok as distributions are unimodal

gibbsoutsimple <- gibbs_complex(50000, a = 1, b = 1)[-(1:0.1*50000), ] #  Eliminate "burn-in"

plot(gibbsoutsimple[,1], type = "l")
plot(gibbsoutsimple[,2], type = "l")

logprior <- dnorm(mupost, 0, sqrt(sigsqpost), log = TRUE) + dinvgamma(sigsqpost, 1, 1, log = TRUE)
loglike <- sum(dnorm(data, mupost, sqrt(sigsqpost), log = TRUE))
logpost1 <- dnorm(mupost, (n*datamean)/(n + 1), sqrt((sigsqpost)/(n + 1)), log = TRUE)
post2avg <- 0

for (i in 1:length(mugibbs)){
  post2avg <- post2avg + dinvgamma(sigsqpost, (1 + n + 1)/2, (1 + sum((data - mugibbs[i])^2) + mugibbs[i]^2)/2)
}

post2 <- post2avg / length(mugibbs)

logpost <- logpost1 + log(post2)

chib <- loglike + logprior - logpost # Log-scale calculation

exp(chib)
