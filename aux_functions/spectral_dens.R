spde_spectral <- function(order, kappa, sigma, alpha, w) {
  # Get the coefficient table based on order
  table_name <- sprintf("m_brasil%dt", order)
  coef_table <- get(table_name, envir = asNamespace("rSPDE"))
  
  # Get fractional part of alpha for coefficient lookup
  alpha_frac <- alpha - floor(alpha)
  
  # Find the row corresponding to the fractional part of alpha
  row_idx <- which.min(abs(coef_table$alpha - alpha_frac))
  coeffs <- coef_table[row_idx, ]
  
  # Calculate nu (nu = alpha - d/2)
  nu <- alpha - 0.5
  
  # Calculate A
  A <- sqrt(2) * kappa^(2*nu) * gamma(nu + 0.5) / gamma(nu)
  
  # Base denominator term
  base_denom <- (1 + kappa^(-2) * w^2)^floor(alpha)
  
  # Calculate the first term
  first_term <- coeffs$k / base_denom
  
  # Calculate the sum terms based on order
  sum_terms <- 0
  for(i in 1:order) {
    r_col <- paste0("r", i)
    p_col <- paste0("p", i)
    sum_terms <- sum_terms + coeffs[[r_col]] / (base_denom * (1 + kappa^(-2) * w^2 - coeffs[[p_col]]))
  }
  
  # Combine all terms
  result <- A * sigma^2 * kappa^(-2*alpha) * (first_term + sum_terms)
  
  return(result)
}

matern_spectral <- function(kappa, sigma, alpha, w) {
  # Calculate nu (nu = alpha - d/2)
  nu <- alpha - 0.5
  
  # Calculate A
  A <- sqrt(2) * kappa^(2*nu) * gamma(nu + 0.5) / gamma(nu)
  
  # Calculate the spectral density
  result <- A * sigma^2 * (kappa^2 + w^2)^(-alpha)
  
  return(result)
}

library(ggplot2)
library(dplyr)
library(tidyr)

# Create function to get range from nu and kappa
get_kappa <- function(range, nu) {
  sqrt(8*nu)/range
}

# Create data for plotting
create_spectral_data <- function(nu, range, w_seq) {
  kappa <- get_kappa(range, nu)
  alpha <- nu + 0.5  # since we're in dimension 1
  
  # Calculate true spectral density
  true_spec <- sapply(w_seq, function(w) matern_spectral(kappa, 1, alpha, w))
  
  # Calculate SPDE approximations for orders 2 through 6
  spde_specs <- list()
  for(m in 1:6) {
    spde_specs[[paste0("m", m)]] <- sapply(w_seq, function(w) 
      spde_spectral(m, kappa, 1, alpha, w))
  }
  
  # Combine into data frame
  data.frame(
    w = w_seq,
    true = true_spec,
    m1 = spde_specs$m1,
    m2 = spde_specs$m2,
    m3 = spde_specs$m3,
    m4 = spde_specs$m4,
    m5 = spde_specs$m5,
    m6 = spde_specs$m6
  ) %>%
    pivot_longer(cols = -w,
                names_to = "type",
                values_to = "density")
}

# Set up parameters
w_seq <- seq(0, 5, length.out = 200)
ranges <- c(0.5, 1, 2)
nus <- c(0.6, 1.4)

# Create all combinations of parameters
plot_data <- expand.grid(nu = nus, range = ranges) %>%
  group_by(nu, range) %>%
  do(create_spectral_data(.$nu, .$range, w_seq))


spectral_plot <- ggplot(plot_data, aes(x = w, y = density, color = type, linetype = type)) +
  geom_line(size = 0.7) +
  facet_grid(nu ~ range, 
             labeller = label_bquote(
               rows = nu == .(nu),
               cols = rho == .(range)
             )) +
  scale_color_manual(values = c("red", "purple", "blue", "yellow", "green", "orange", "black"),
                    labels = c("m1" = "Order 1",
                            "m2" = "Order 2", 
                             "m3" = "Order 3",
                             "m4" = "Order 4",
                             "m5" = "Order 5",
                             "m6" = "Order 6",
                             "true" = "True")) +
  scale_linetype_manual(values = c("dashed", "dotted", "longdash", "twodash", "dotdash", "dotted", "solid"),
                       labels = c("m1" = "Order 1",
                                "m2" = "Order 2", 
                                "m3" = "Order 3",
                                "m4" = "Order 4",
                                "m5" = "Order 5",
                                "m6" = "Order 6",
                                "true" = "True")) +
  labs(x = "Frequency (w)",
       y = "Spectral Density") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "lightgray"),
    panel.grid.minor = element_line(color = "gray"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing = unit(1, "lines"),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.spacing.x = unit(0.5, "cm"),
    legend.margin = margin(t = 0, b = 0),
    legend.title.align = 0.5
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 0.7)),
    linetype = guide_legend(override.aes = list(size = 0.7))
  )

ggsave("spectral_comparison.png", 
       plot = spectral_plot,
       width = 10, 
       height = 7, 
       units = "in",
       dpi = 300)
