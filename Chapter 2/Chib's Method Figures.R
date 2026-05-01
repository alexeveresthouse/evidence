set.seed(2026)

x <- seq(-5, 5, 10/10000)

y1 <- dnorm(x, mean = -2, sd = 1)
y2 <- dnorm(x, mean = 2, sd = 0.5)

E1 <- cbind.data.frame(x, y1)
E2 <- cbind.data.frame(x, y2)

ggplot(E1) +
  geom_line(aes(x = x, y = y1), linewidth = 1.5, colour = "#26316f") +
  xlab("x") +
  ylab("Density") +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "#f5f5ff"))

ggplot(E2) +
  geom_line(aes(x = x, y = y2), linewidth = 1.5, colour = "#26316f") +
  xlab("x") +
  ylab("Density") +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "#f5f5ff"))

y <- y1 + y2

E3 <- cbind.data.frame(x, y)

ggplot(E3) +
  geom_line(aes(x = x, y = y), linewidth = 1.5, colour = "#e5bd0a") +
  xlab("x") +
  ylab("Density") +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "#f5f5ff"))
