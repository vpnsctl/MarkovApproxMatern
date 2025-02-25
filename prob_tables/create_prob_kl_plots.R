# prob_regular <- readRDS("prob_tables/prob_errors_regular.RDS")
# kl_regular <- readRDS("kl_tables/post_kl_regular.RDS")

prob_regular <- readRDS("prob_tables/prob_errors_full.RDS")
kl_regular <- readRDS("kl_tables/post_kl.RDS")

prob_regular <- prob_regular %>%
  rename(order = m, n = N) %>%
  mutate(
    abs_error = abs(error),
    Error_Type = "Absolute Probability Error",
    order = as.character(order),
    n = as.numeric(n), 
    range = as.numeric(range),
    Method = case_when(
      Method == "Rational" ~ "Rational",
      Method == "NNGP" ~ "nnGP",
      Method == "FEM" ~ "FEM",
      Method == "nnGP" ~ "nnGP"
    )
  ) %>% 
  filter(order %in% 3:6)

# Create a separate dataframe for regular probability error
prob_regular_raw <- prob_regular %>%
  mutate(
    value = error,
    Error_Type = "Probability Error"
  )

# Prepare absolute probability error dataframe
prob_regular_abs <- prob_regular %>%
  mutate(
    value = abs_error,
    Error_Type = "Absolute Probability Error"
  )

# Prepare KL divergence data
kl_regular <- kl_regular %>%
  mutate(
    Error_Type = "KL Divergence",
    value = error,
    Method = case_when(
      Method == "Rational" ~ "Rational",
      Method == "NNGP" ~ "nnGP",
      Method == "FEM" ~ "FEM",
      Method == "nnGP" ~ "nnGP"
    )
  ) %>% 
  rename(order = m, n = N) %>% 
  mutate(order = as.character(order), n = as.numeric(n), range = as.numeric(range)) %>% 
  filter(order %in% 3:6)

# Combine all three datasets
combined_df <- bind_rows(
  prob_regular_abs,
  prob_regular_raw,
  kl_regular
) %>%
  mutate(Error_Type = factor(Error_Type, 
                            levels = c("Absolute Probability Error", 
                                     "Probability Error", 
                                     "KL Divergence")))

# Update the plotting function
create_comparison_plot <- function(data, sigma_e_val) {
  filtered_data <- data %>% 
    filter(abs(sigma_e - sigma_e_val) < 1e-10)
  
  ggplot(filtered_data, 
         aes(x = n, y = value, color = Method)) +
    geom_line(aes(linetype = factor(order)), linewidth = 0.7) +
    scale_color_manual(
      values = c("limegreen", "red", "black"),
      labels = c("FEM", "NNGP", "Rational")
    ) +
    scale_linetype_manual(
      name = "Order",
      values = c("solid", "dashed", "dotted", "twodash", "dotdash", "longdash"),
      labels = paste("Order", 3:6)
    ) +
    facet_grid(
      Error_Type ~ range,
      labeller = labeller(range = function(x) paste0("ρ = ", x)),
      scales = "free_y"
    ) +
    ggh4x::facetted_pos_scales(
      y = list(
        "Absolute Probability Error" = scale_y_continuous(),
        "Probability Error" = scale_y_continuous(),
        "KL Divergence" = scale_y_log10()
      )
    ) + 
    labs(
      y = "Posterior Errors",
      x = "N (number of prediction locations)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 16),
      panel.grid.major = element_line(color = "lightgray"),
      panel.grid.minor = element_line(color = "gray"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.spacing = unit(1, "lines"),
      strip.background = element_blank(),
      strip.placement = "outside",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.spacing.x = unit(0.5, "cm"),
      legend.margin = margin(t = 0, b = 0),
      legend.key.width = unit(1.5, "cm"),
      legend.key.height = unit(0.5, "cm")
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(linewidth = 1)),
      linetype = guide_legend(order = 2, nrow = 2)
    )
}

# Create and save the plots
p1 <- create_comparison_plot(combined_df, 0.1)
p2 <- create_comparison_plot(combined_df, 0.316227766016838)

ggsave("comparison_plot_sigma01.png", p1, width = 15, height = 12)
ggsave("comparison_plot_sigma316.png", p2, width = 15, height = 12)

# library(ggplot2)
# library(dplyr)


# # Prepare the data
# # Rename columns in prob_regular to match kl_regular
# prob_regular <- prob_regular %>%
#   rename(order = m, n = N) %>%
#   mutate(
#     value = abs(error),  # Take absolute value of errors
#     Error_Type = "Probability Error",
#     order = as.character(order),
#     n = as.numeric(n), range = as.numeric(range),
#     Method = case_when(
#       Method == "Rational" ~ "Rational",
#       Method == "NNGP" ~ "nnGP",
#       Method == "FEM" ~ "FEM",
#       Method == "nnGP" ~ "nnGP"
#     )
#   ) %>% filter( order %in% 3:6)

# kl_regular <- kl_regular %>%
#   mutate(
#     Error_Type = "KL Divergence",
#     value = error,
#     Method = case_when(
#       Method == "Rational" ~ "Rational",
#       Method == "NNGP" ~ "nnGP",
#       Method == "FEM" ~ "FEM",
#       Method == "nnGP" ~ "nnGP"
#     )
#   ) %>% rename(order = m, n = N) %>% mutate(order = as.character(order), n = as.numeric(n), range = as.numeric(range)) %>% filter( order %in% 3:6)

# # Combine the datasets and set factor levels for Error_Type
# combined_df <- bind_rows(prob_regular, kl_regular) %>%
#   mutate(Error_Type = factor(Error_Type, levels = c("Probability Error", "KL Divergence")))

# # Create separate plots for each sigma_e value
# # create_comparison_plot <- function(data, sigma_e_val) {
# #   filtered_data <- data %>% 
# #     filter(abs(sigma_e - sigma_e_val) < 1e-10)
  
# #   ggplot(filtered_data, 
# #          aes(x = n, y = value, color = Method)) +
# #     geom_line(aes(linetype = factor(order)), linewidth = 0.7) +
# #     scale_y_log10(limits = c(1e-10, NA)) +
# #     scale_color_manual(
# #       values = c("limegreen", "red", "black"),
# #       labels = c("FEM", "NNGP", "Rational")
# #     ) +
# #     scale_linetype_manual(
# #       name = "Order",
# #       values = c("solid", "dashed", "dotted", "twodash", "dotdash", "longdash"),
# #       labels = paste("Order", 1:6)
# #     ) +
# #     labs(
# #       y = "Posterior Error",
# #       x = "N (number of prediction locations)"
# #     ) +
# #     facet_grid(
# #       Error_Type ~ range,
# #       labeller = labeller(range = function(x) paste0("ρ = ", x))
# #     ) +
# #     theme_minimal() +
# #     theme(
# #       legend.position = "bottom",
# #       legend.text = element_text(size = 14),
# #       legend.title = element_text(size = 16),
# #       axis.title = element_text(size = 16),
# #       axis.text = element_text(size = 14),
# #       strip.text = element_text(size = 16),
# #       panel.grid.major = element_line(color = "lightgray"),
# #       panel.grid.minor = element_line(color = "gray"),
# #       panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
# #       panel.spacing = unit(1, "lines"),
# #       strip.background = element_blank(),
# #       strip.placement = "outside",
# #       legend.box = "horizontal",
# #       legend.box.just = "center",
# #       legend.spacing.x = unit(0.5, "cm"),
# #       legend.margin = margin(t = 0, b = 0),
# #       legend.key.width = unit(1.5, "cm"),
# #       legend.key.height = unit(0.5, "cm")
# #     ) +
# #     guides(
# #       color = guide_legend(order = 1, override.aes = list(linewidth = 1)),
# #       linetype = guide_legend(order = 2, nrow = 2)
# #     )
# # }

# create_comparison_plot <- function(data, sigma_e_val) {
#   filtered_data <- data %>% 
#     filter(abs(sigma_e - sigma_e_val) < 1e-10)
  
#   ggplot(filtered_data, 
#          aes(x = n, y = value, color = Method)) +
#     geom_line(aes(linetype = factor(order)), linewidth = 0.7) +
#     scale_color_manual(
#       values = c("limegreen", "red", "black"),
#       labels = c("FEM", "NNGP", "Rational")
#     ) +
#     scale_linetype_manual(
#       name = "Order",
#       values = c("solid", "dashed", "dotted", "twodash", "dotdash", "longdash"),
#       labels = paste("Order", 1:6)
#     ) +
#     # Different y-axis scales for each facet
#     facet_grid(
#       Error_Type ~ range,
#       labeller = labeller(range = function(x) paste0("ρ = ", x)),
#       scales = "free_y"  # This allows different y-axis scales
#     ) +
#     # Use scale_y_log10 with different limits based on Error_Type
#     # scale_y_log10() +
#     ggh4x::facetted_pos_scales(
#       y = list(
#         "Probability Error" = scale_y_continuous(),
#         "KL Divergence" = scale_y_log10()
#       )) + 
#     labs(
#       y = "Posterior Error",
#       x = "N (number of prediction locations)"
#     ) +
#     theme_minimal() +
#     theme(
#       legend.position = "bottom",
#       legend.text = element_text(size = 14),
#       legend.title = element_text(size = 16),
#       axis.title = element_text(size = 16),
#       axis.text = element_text(size = 14),
#       strip.text = element_text(size = 16),
#       panel.grid.major = element_line(color = "lightgray"),
#       panel.grid.minor = element_line(color = "gray"),
#       panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#       panel.spacing = unit(1, "lines"),
#       strip.background = element_blank(),
#       strip.placement = "outside",
#       legend.box = "horizontal",
#       legend.box.just = "center",
#       legend.spacing.x = unit(0.5, "cm"),
#       legend.margin = margin(t = 0, b = 0),
#       legend.key.width = unit(1.5, "cm"),
#       legend.key.height = unit(0.5, "cm")
#     ) +
#     guides(
#       color = guide_legend(order = 1, override.aes = list(linewidth = 1)),
#       linetype = guide_legend(order = 2, nrow = 2)
#     )
# }

# # Create plots for each sigma_e value
# p1 <- create_comparison_plot(combined_df, 0.1)
# p2 <- create_comparison_plot(combined_df, 0.316227766016838)

# # Save the plots
# ggsave("comparison_plot_sigma01.png", p1, width = 15, height = 8)
# ggsave("comparison_plot_sigma316.png", p2, width = 15, height = 8)