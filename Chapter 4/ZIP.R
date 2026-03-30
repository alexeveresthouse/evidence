library(rstan)
library(tidyverse)
library(ggplot2)

##### Poisson Setup #####

set.seed(2026)
n <- 340  
m <- 4    

x1 <- rnorm(n, mean = 0, sd = 1)
x2 <- rnorm(n, mean = 0, sd = 1)
x3 <- rnorm(n, mean = 0, sd = 1)

x_matrix <- cbind(1, x1, x2, x3)

true_beta <- c(1.2, 0.5, -0.4, 0.3) 
true_p <- 0.30                      

log_eta <- x_matrix %*% true_beta
eta <- exp(log_eta)

y_synthetic <- numeric(n)

for (i in 1:n) {
  is_perfect_zero <- rbinom(1, size = 1, prob = true_p)
  if (is_perfect_zero == 1) {
    y_synthetic[i] <- 0
  } else {
    mu
    y_synthetic[i] <- rpois(1, lambda = eta[i])
  }
}

poisson_model <- stan_model("Poisson Power Posterior.stan")

stan_data <- list(
  n = n,
  m = m,
  x = x_matrix,
  y = y_synthetic,
  t = 1.0  # Set temperature to 1 to sample the full posterior!
)

# Run the model
fit_synthetic <- sampling(
  poisson_model, 
  data = stan_data,
  iter = 10000,
  warmup = 5000,
  chains = 4
)

# Check if Stan recovered the true parameters!
print(fit_synthetic, pars = c("beta"))
fitted_beta <- colMeans(rstan::extract(fit_synthetic)$beta)
y_fitted <- exp(x_matrix %*% fitted_beta)
sim_counts <- rpois(length(y_fitted), lambda = y_fitted)

plot_data <- data.frame(
  Count = c(y_synthetic, sim_counts),
  Source = rep(c("Actual Data", "Simulated (Poisson)"), each = length(y_synthetic))
)

ggplot(plot_data, aes(x = Count, fill = Source)) +
  geom_bar(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Actual Data" = "paleturquoise4", 
                               "Simulated (Poisson)" = "coral")) +
  scale_x_continuous(breaks = 0:max(plot_data$Count)) +
  labs(x = "Count Value",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 14))

##### ZIP setup #####

set.seed(2026)
n <- 340  
m <- 4    


x1 <- rnorm(n, mean = 0, sd = 1)
x2 <- rnorm(n, mean = 0, sd = 1)
x3 <- rnorm(n, mean = 0, sd = 1)

x_matrix <- cbind(1, x1, x2, x3)

true_beta <- c(1.2, 0.5, -0.4, 0.3) 
true_p <- 0.30                      

log_eta <- x_matrix %*% true_beta
eta <- exp(log_eta)

y_synthetic <- numeric(n)

for (i in 1:n) {
  is_perfect_zero <- rbinom(1, size = 1, prob = true_p)
  if (is_perfect_zero == 1) {
    y_synthetic[i] <- 0
  } else {
    y_synthetic[i] <- rpois(1, lambda = eta[i])
  }
}

zip_model <- stan_model("ZIP Power Posterior.stan")

stan_data <- list(n = n,
                  m = m,
                  x = x_matrix,
                  y = y_synthetic,
                  t = 1.0
)

# Run the model
fit_synthetic_zip <- sampling(zip_model, 
                              data = stan_data,
                              iter = 10000,
                              warmup = 5000,
                              chains = 4
)

# Check if Stan recovered the true parameters!
print(fit_synthetic_zip, pars = c("beta", "logit_p"))
fitted_beta <- colMeans(rstan::extract(fit_synthetic_zip)$beta)
fitted_p <- mean(rstan::extract(fit_synthetic_zip)$p)
y_fitted <- exp(x_matrix %*% fitted_beta)
fitted_is_perfect_zero <- rbinom(length(y_fitted), 1, prob = fitted_p)

sim_counts <- ifelse(fitted_is_perfect_zero == 1,
                     0,
                     rpois(length(y_fitted), lambda = y_fitted))

plot_data <- data.frame(
  Count = c(y_synthetic, sim_counts),
  Source = rep(c("Actual Data", "Simulated (ZIP)"), each = length(y_synthetic))
)

ggplot(plot_data, aes(x = Count, fill = Source)) +
  geom_bar(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Actual Data" = "paleturquoise4", 
                               "Simulated (ZIP)" = "coral")) +
  scale_x_continuous(breaks = 0:max(plot_data$Count)) +
  labs(x = "Count Value",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 14))

### Power posterior evidence estimate - base Poisson ###

# Create vector of temperatures; preferably denser closer to prior
ntemps <- 25
c <- 5
eps <- 0.01
t <- (seq(0, 1, length = ntemps)^c)
t <- eps + (1 - eps) * t
# Create vector to store log-likelihood evaluations in
logliks <- rep(NA, length(t))
# Repeat Stan evaluations at different temperatures
for (i in 1:length(t)){
  posti <- sampling(poisson_model,
                    data = list(n = n,
                                m = m,
                                t = t[i],
                                x = x_matrix,
                                y = y_synthetic),
                    iter = 10000,
                    warmup = 5000,
                    chains = 1,
                    refresh = 0,
                    show_messages = FALSE)
  
  logliks[i] <- mean(rstan::extract(posti)$log_lik)
  cat(sprintf("Completed Run %d\n", i))
}

# Update the log evidence using a trapezoidal approximation to the relevant integral
logevcomplex <- 0
for (i in 1:(length(t) - 1)){
  logevcomplex <- logevcomplex + (t[i + 1] - t[i])*(0.5)*(logliks[i + 1] + logliks[i])
}

# Final evidence value
logevcomplex


### Power Posterior evidence estimate - ZIP ###

# Create vector of temperatures; preferably denser closer to prior
ntemps <- 25
c <- 5
eps <- 0.01
t <- (seq(0, 1, length = ntemps)^c)
t <- eps + (1 - eps) * t
# Create vector to store log-likelihood evaluations in
logliks <- rep(NA, length(t))
# Repeat Stan evaluations at different temperatures
for (i in 1:length(t)){
  posti <- sampling(zip_model,
                data = list(n = n,
                            m = m,
                            t = t[i],
                            x = x_matrix,
                            y = y_synthetic),
                iter = 10000,
                warmup = 5000,
                chains = 1,
                refresh = 0,
                show_messages = FALSE)
  
  logliks[i] <- mean(rstan::extract(posti)$log_lik)
  cat(sprintf("Completed Run %d\n", i))
}

# Update the log evidence using a trapezoidal approximation to the relevant integral
logevcomplex <- 0
for (i in 1:(length(t) - 1)){
  logevcomplex <- logevcomplex + (t[i + 1] - t[i])*(0.5)*(logliks[i + 1] + logliks[i])
}

# Final evidence value
logevcomplex
