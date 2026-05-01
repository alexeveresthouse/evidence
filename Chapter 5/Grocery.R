library(rstan)
library(loo)
library(tidyverse)
library(ggplot2)

### Data setup ###

grocerydata <- read.csv("perv2pub.csv")

grocery <- grocerydata[, c("DELIV_GROC", "R_AGE", "R_SEX", "EDUC", "WORKER",
                     "MEDCOND", "DRIVER")]

num_cols <- sapply(grocery, is.numeric)
grocery[num_cols][grocery[num_cols] < 0] <- NA

grocery$R_SEX <- factor(grocery$R_SEX, levels = c("1", "2"))
grocery$EDUC <- as.factor(grocery$EDUC)
grocery$WORKER <- factor(grocery$WORKER, levels = c("1", "2"))
grocery$MEDCOND <- factor(grocery$MEDCOND, levels = c("1", "2"))
grocery$DRIVER <- factor(grocery$DRIVER, levels = c("1", "2"))

grocery <- na.omit(grocery)

x_matrix <- model.matrix(~ R_AGE + R_SEX + EDUC + WORKER + MEDCOND + DRIVER,
                         data = grocery)

ggplot(grocery, aes(x = DELIV_GROC)) +
  geom_histogram(fill = "#26316f") +
  scale_x_continuous(breaks = 0:max(grocery$DELIV_GROC)) +
  labs(x = "Count Value",
       y = "Frequency") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 14),
        panel.background = element_rect(fill = "#f5f5ff", colour = NA))

### Fits ###

poisson_model <- stan_model("Poisson Power Posterior.stan")

grocery_data <- list(
  n = nrow(grocery),
  m = ncol(x_matrix),
  x = x_matrix,
  y = grocery$DELIV_GROC,
  t = 1.0 
)

fit_poisson <- sampling(poisson_model,
                        data = grocery_data,
                        iter = 10000,
                        warmup = 5000,
                        chains = 1)


zip_model <- stan_model("ZIP Power Posterior.stan")

fit_zip <- sampling(zip_model,
                    data = grocery_data,
                    iter = 10000,
                    warmup = 5000,
                    chains = 1)

negbin_model <- stan_model("Negative Binomial Power Posterior.stan")

fit_negbin <- sampling(negbin_model,
                       data = grocery_data,
                       iter = 10000,
                       warmup = 5000,
                       chains = 1)

zinb_model <- stan_model("ZINB Power Posterior.stan")

fit_zinb <- sampling(zinb_model,
                     data = grocery_data,
                     iter = 10000,
                     warmup = 5000,
                     chains = 1)

# Save models for future use

saveRDS(fit_poisson, file = "fit_poisson.rds")
saveRDS(fit_zip, file = "fit_zip.rds")
saveRDS(fit_negbin, file = "fit_negbin.rds")
saveRDS(fit_zinb, file = "fit_zinb.rds")

# Read models in

fit_poisson <- readRDS("fit_poisson.rds")
fit_zip <- readRDS("fit_zip.rds")
fit_negbin <- readRDS("fit_negbin.rds")
fit_zinb <- readRDS("fit_zinb.rds")

### ASIDE: POISSON GLM FIT ###

# poisson_glm <- glm(DELIV_GROC ~ R_AGE + R_SEX + EDUC + WORKER + MEDCOND + DRIVER,
#                    data = grocery)

# summary(poisson_glm)

# Poisson plot

fitted_beta <- colMeans(rstan::extract(fit_poisson)$beta)
y_fitted <- exp(x_matrix %*% fitted_beta)
sim_counts <- rpois(length(y_fitted), lambda = y_fitted)

plot_data <- data.frame(
  Count = c(grocery$DELIV_GROC, sim_counts),
  Source = rep(c("Grocery Data", "Simulated Data (Poisson Model)"), each = length(grocery$DELIV_GROC))
)

ggplot(plot_data, aes(x = Count, fill = Source)) +
  geom_bar(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Grocery Data" = "#26316f", 
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

# ZIP plot

fitted_beta <- colMeans(rstan::extract(fit_zip)$beta)
fitted_p <- mean(rstan::extract(fit_zip)$p)
y_fitted <- exp(x_matrix %*% fitted_beta)
fitted_is_perfect_zero <- rbinom(length(y_fitted), 1, prob = fitted_p)

sim_counts <- ifelse(fitted_is_perfect_zero == 1,
                     0,
                     rpois(length(y_fitted), lambda = y_fitted))

plot_data <- data.frame(
  Count = c(grocery$DELIV_GROC, sim_counts),
  Source = rep(c("Grocery Data", "Simulated Data (ZIP Model)"), each = length(grocery$DELIV_GROC))
)

ggplot(plot_data, aes(x = Count, fill = Source)) +
  geom_bar(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Grocery Data" = "#26316f", 
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

# Negative Binomial plot

fitted_beta <- colMeans(rstan::extract(fit_negbin)$beta)
fitted_phi <- mean(rstan::extract(fit_negbin)$phi)
y_fitted <- exp(x_matrix %*% fitted_beta)
sim_mu <- rgamma(length(y_fitted), fitted_phi, rate = fitted_phi / y_fitted)
sim_counts <- rpois(length(y_fitted), lambda = sim_mu)

plot_data <- data.frame(
  Count = c(grocery$DELIV_GROC, sim_counts),
  Source = rep(c("Grocery Data", "Simulated Data (Negative Binomial Model)"), each = length(grocery$DELIV_GROC))
)

ggplot(plot_data, aes(x = Count, fill = Source)) +
  geom_bar(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Grocery Data" = "#26316f", 
                               "Simulated Data (Negative Binomial Model)" = "#cacbff")) +
  scale_x_continuous(breaks = 0:max(plot_data$Count)) +
  labs(x = "Count Value",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 14),
        panel.background = element_rect(fill = "#f5f5ff", colour = NA))

# ZINB plot

fitted_beta <- colMeans(rstan::extract(fit_zinb)$beta)
fitted_phi <- mean(rstan::extract(fit_zinb)$phi)
fitted_p <- mean(rstan::extract(fit_zinb)$p)
y_fitted <- exp(x_matrix %*% fitted_beta)
sim_mu <- rgamma(length(y_fitted), fitted_phi, rate = fitted_phi / y_fitted)
fitted_is_perfect_zero <- rbinom(length(y_fitted), 1, prob = fitted_p)
sim_counts <- ifelse(fitted_is_perfect_zero == 1,
                     0,
                     rpois(length(y_fitted), lambda = sim_mu))

plot_data <- data.frame(
  Count = c(grocery$DELIV_GROC, sim_counts),
  Source = rep(c("Grocery Data", "Simulated Data (ZINB Model)"), each = length(grocery$DELIV_GROC))
)

ggplot(plot_data, aes(x = Count, fill = Source)) +
  geom_bar(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Grocery Data" = "#26316f", 
                               "Simulated Data (ZINB Model)" = "#cacbff")) +
  scale_x_continuous(breaks = 0:max(plot_data$Count)) +
  labs(x = "Count Value",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 14),
        panel.background = element_rect(fill = "#f5f5ff", colour = NA))

### Computing the evidence ###

n <- nrow(grocery)
m <- ncol(x_matrix)

ntemps <- 20
c <- 5
eps <- 0.01
t <- (seq(0, 1, length = ntemps)^c)
t <- eps + (1 - eps) * t

# Poisson model

# Create vector to store log-likelihood evaluations in
logliks_poisson <- rep(NA, length(t))
# Repeat Stan evaluations at different temperatures
for (i in 1:length(t)){
  posti <- sampling(poisson_model,
                    data = list(n = n,
                                m = m,
                                t = t[i],
                                x = x_matrix,
                                y = grocery$DELIV_GROC),
                    iter = 10000,
                    warmup = 5000,
                    chains = 1,
                    refresh = 0,
                    show_messages = FALSE)
  
  logliks_poisson[i] <- mean(rstan::extract(posti)$log_lik)
  cat(sprintf("Completed Run %d\n", i))
}

# Update the log evidence using a trapezoidal approximation to the relevant integral
logevcomplexpoisson <- 0
for (i in 1:(length(t) - 1)){
  logevcomplexpoisson <- logevcomplexpoisson +
    (t[i + 1] - t[i])*(0.5)*(logliks_poisson[i + 1] + logliks_poisson[i])
}

# ZIP model

# Create vector to store log-likelihood evaluations in
logliks_zip <- rep(NA, length(t))
# Repeat Stan evaluations at different temperatures
for (i in 1:length(t)){
  posti <- sampling(zip_model,
                    data = list(n = n,
                                m = m,
                                t = t[i],
                                x = x_matrix,
                                y = grocery$DELIV_GROC),
                    iter = 10000,
                    warmup = 5000,
                    chains = 1,
                    refresh = 0,
                    show_messages = FALSE)
  
  logliks_zip[i] <- mean(rstan::extract(posti)$log_lik)
  cat(sprintf("Completed Run %d\n", i))
}

# Update the log evidence using a trapezoidal approximation to the relevant integral
logevzip <- 0
for (i in 1:(length(t) - 1)){
  logevzip <- logevzip +
    (t[i + 1] - t[i])*(0.5)*(logliks_zip[i + 1] + logliks_zip[i])
}

# Negative binomial model

# Create vector to store log-likelihood evaluations in
logliks_negbin <- rep(NA, length(t))
# Repeat Stan evaluations at different temperatures
for (i in 1:length(t)){
  posti <- sampling(negbin_model,
                    data = list(n = n,
                                m = m,
                                t = t[i],
                                x = x_matrix,
                                y = grocery$DELIV_GROC),
                    iter = 10000,
                    warmup = 5000,
                    chains = 1,
                    refresh = 0,
                    show_messages = FALSE)
  
  logliks_negbin[i] <- mean(rstan::extract(posti)$log_lik)
  cat(sprintf("Completed Run %d\n", i))
}

# Update the log evidence using a trapezoidal approximation to the relevant integral
logevnegbin <- 0
for (i in 1:(length(t) - 1)){
  logevnegbin <- logevnegbin +
    (t[i + 1] - t[i])*(0.5)*(logliks_negbin[i + 1] + logliks_negbin[i])
}

# Create vector to store log-likelihood evaluations in
logliks_negbin <- rep(NA, length(t))
# Repeat Stan evaluations at different temperatures
for (i in 1:length(t)){
  posti <- sampling(negbin_model,
                    data = list(n = n,
                                m = m,
                                t = t[i],
                                x = x_matrix,
                                y = grocery$DELIV_GROC),
                    iter = 10000,
                    warmup = 5000,
                    chains = 1,
                    refresh = 0,
                    show_messages = FALSE)
  
  logliks_negbin[i] <- mean(rstan::extract(posti)$log_lik)
  cat(sprintf("Completed Run %d\n", i))
}

# Update the log evidence using a trapezoidal approximation to the relevant integral
logevnegbin <- 0
for (i in 1:(length(t) - 1)){
  logevnegbin <- logevnegbin +
    (t[i + 1] - t[i])*(0.5)*(logliks_negbin[i + 1] + logliks_negbin[i])
}

# ZINB model

# Create vector to store log-likelihood evaluations in
logliks_zinb <- rep(NA, length(t))
# Repeat Stan evaluations at different temperatures
for (i in 1:length(t)){
  posti <- sampling(zinb_model,
                    data = list(n = n,
                                m = m,
                                t = t[i],
                                x = x_matrix,
                                y = grocery$DELIV_GROC),
                    iter = 10000,
                    warmup = 5000,
                    chains = 1,
                    refresh = 0,
                    show_messages = FALSE)
  
  logliks_zinb[i] <- mean(rstan::extract(posti)$log_lik)
  cat(sprintf("Completed Run %d\n", i))
}

# Update the log evidence using a trapezoidal approximation to the relevant integral
logevzinb <- 0
for (i in 1:(length(t) - 1)){
  logevzinb <- logevzinb +
    (t[i + 1] - t[i])*(0.5)*(logliks_zinb[i + 1] + logliks_zinb[i])
}

### Posterior model probabilities ###

logSumExp <- function(logx){
  c <- max(logx)
  return(c + log(sum(exp(logx - c))))
}

evidence1 <- -8971.715
evidence2 <- -5900.029
evidence3 <- -5465.750
evidence4 <- -5447.424

logevidences <- c(evidence1, evidence2, evidence3, evidence4)

probs <- c(0.97, rep(0.01, 3))
logprob <- log(probs)

numerator <- logprob + logevidences

denominator <- logSumExp(numerator)

logposts <- numerator - denominator
posts <- exp(logposts)

### DIC ####

log_lik_poisson <- extract_log_lik(fit_poisson)
log_lik_zip <- extract_log_lik(fit_zip)
log_lik_negbin <- extract_log_lik(fit_negbin)
log_lik_zinb <- extract_log_lik(fit_zinb)

deviance_poisson <- -2 * rowSums(log_lik_poisson)
Dbar_poisson <- mean(deviance_poisson)
pD_poisson <- var(deviance_poisson) / 2

DIC_poisson <- Dbar_poisson + pD_poisson

deviance_zip <- -2 * rowSums(log_lik_zip)
Dbar_zip <- mean(deviance_zip)
pD_zip <- var(deviance_zip) / 2

DIC_zip <- Dbar_zip + pD_zip

deviance_negbin <- -2 * rowSums(log_lik_negbin)
Dbar_negbin <- mean(deviance_negbin)
pD_negbin <- var(deviance_negbin) / 2

DIC_negbin <- Dbar_negbin + pD_zip

deviance_zinb <- -2 * rowSums(log_lik_zinb)
Dbar_zinb <- mean(deviance_zinb)
pD_zinb <- var(deviance_zinb) / 2

DIC_zinb <- Dbar_zinb + pD_zinb
