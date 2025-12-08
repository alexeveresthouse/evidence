stanfit2 <- stan("complex model.stan",
                 data = list(n = 7,
                             x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)))

log_lik_sample <- as.matrix(stanfit2, pars = "log_lik")
likelihoods <- exp(log_lik_sample)

recips <-  1/likelihoods

harmonicmean <- 1/mean(recips)



stanfit2
library(durhamSLR)

diagnostics(as.array(stanfit2))

data <- c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)
nprior <- 10000

mu_prior <- rnorm(nprior, 0, 1)
sigma_prior <- 1/rgamma(nprior, 1, 2)

likelihood <- numeric(nprior)

for (i in 1:nprior){
  logliki <- sum(dnorm(data, mu_prior[i], sigma_prior[i], log = TRUE))
  likelihood[i] <- exp(logliki)
}

evidenceMC <- mean(likelihood)
