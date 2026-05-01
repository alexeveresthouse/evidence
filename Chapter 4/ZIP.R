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
    y_synthetic[i] <- rpois(1, lambda = eta[i])
  }
}

ggplot(as.data.frame(y_synthetic), aes(x = y_synthetic)) +
  geom_histogram(fill = "#26316f", binwidth = 1) +
  scale_x_continuous(breaks = 0:max(y_synthetic)) +
  labs(x = "Count Value",
       y = "Frequency") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 14),
        panel.background = element_rect(fill = "#f5f5ff", colour = NA))

poisson_model <- stan_model("Poisson Power Posterior.stan")

stan_data <- list(
  n = n,
  m = m,
  x = x_matrix,
  y = y_synthetic,
  t = 1.0 
)

fit_synthetic <- sampling(
  poisson_model, 
  data = stan_data,
  iter = 10000,
  warmup = 5000,
  chains = 4
)

print(fit_synthetic, pars = c("beta"))
fitted_beta <- colMeans(rstan::extract(fit_synthetic)$beta)
y_fitted <- exp(x_matrix %*% fitted_beta)
sim_counts <- rpois(length(y_fitted), lambda = y_fitted)

plot_data <- data.frame(
  Count = c(y_synthetic, sim_counts),
  Source = rep(c("Original Data", "Simulated Data (Poisson Model)"), each = length(y_synthetic))
)

ggplot(plot_data, aes(x = Count, fill = Source)) +
  geom_bar(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Original Data" = "#26316f", 
                               "Simulated Data (Poisson Model)" = "#cacbff")) +
  scale_x_continuous(breaks = 0:max(plot_data$Count)) +
  labs(x = "Count Value",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 14),
        panel.background = element_rect(fill = "#f5f5ff", colour = NA))

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

fit_synthetic_zip <- sampling(zip_model, 
                              data = stan_data,
                              iter = 10000,
                              warmup = 5000,
                              chains = 4
)

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
  Source = rep(c("Original Data", "Simulated Data (ZIP Model)"), each = length(y_synthetic))
)

ggplot(plot_data, aes(x = Count, fill = Source)) +
  geom_bar(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Original Data" = "#26316f", 
                               "Simulated Data (ZIP Model)" = "#cacbff")) +
  scale_x_continuous(breaks = 0:max(plot_data$Count)) +
  labs(x = "Count Value",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 14),
        panel.background = element_rect(fill = "#f5f5ff", colour = NA))

## See ZINB.R for all evidence calculations, including for these ZIP models.