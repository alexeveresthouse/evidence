library(rstan)

stanfit <- stan("simple model.stan",
                data = list(n = 7,
                            sigma = 2,
                            x = c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)))

stanfit



log_lik_sample <- as.matrix(stanfit, pars = "log_lik")
likelihoods <- exp(log_lik_sample)

n <- 7
m <- min(likelihoods)

recips <-  1/likelihoods

harmonicmean <- 1/mean(recips)

logharmonicmean <- log(n) - m - log(sum(exp(-log_lik_sample - m)))
exp(logharmonicmean)


data <- c(0.4, 0.6, 1, -0.3, -0.9, 0.77, 0.12)
nprior <- 10000

mu_prior <- rnorm(nprior, 0, 1)

likelihood <- numeric(nprior)

for (i in 1:nprior){
  logliki <- sum(dnorm(data, mu_prior[i], 2, log = TRUE))
  likelihood[i] <- exp(logliki)
}

evidenceMC <- mean(likelihood)


# Analytic calculation

n <- 7
sigma <- 2
evidenceanalytic <- (2*pi*((sigma^2)/(n+sigma^2))^2)^(-1/2)*(2*pi*sigma^2)^(-n/2)*(2*pi)^(-0.5)*exp(-1/(2*sigma^2)*sum(data) + ((n)/(2*sigma^2) + 0.5)*((1/(2*sigma^2))/((n/(2*sigma^2) + 0.5))*sum(data))^2)
