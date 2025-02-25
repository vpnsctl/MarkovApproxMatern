library(dplyr)
library(ggplot2)

full_df <- readRDS("kl_tables/KL_full_prior.RDS")
full_df <- full_df %>%
  rename(Method = type) %>%
  mutate(Method = case_when(
    Method == "NNGP" ~ "nnGP",
    TRUE ~ Method
  ))

unique_methods <- unique(full_df$Method)

plot <- ggplot(full_df, 
       aes(x = n, y = KL, 
           color = factor(order),
           linetype = Method)) +
    geom_line(linewidth = 0.7) +
    scale_y_log10() +
    scale_color_manual(
      name = "Order",
      values = c("red", "blue", "green", "purple", "orange", "black"),
      labels = paste("Order", 1:6)
    ) +
    scale_linetype_manual(
      name = "Method",
      values = setNames(c("solid", "dotted", "dashed"), unique_methods)  # One line type for each method
    ) +
    facet_grid(
      . ~ range,
      labeller = labeller(range = function(x) paste0("ρ = ", x))
    ) +
    labs(
      y = "KL Divergence",
      x = "N (number of non-support locations)"
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
      legend.margin = margin(t = 0, b = 0)
    ) +
    guides(
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    )

ggsave("kl_tables/prior_KL_error.png", plot, width = 14, height = 5)
