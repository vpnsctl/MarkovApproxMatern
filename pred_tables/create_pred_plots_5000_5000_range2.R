library(ggplot2)
library(dplyr)
library(tidyr)

color_plot_options <- c("black", "steelblue", "limegreen", "red", "purple", "orange", "brown")
color_plot_used <- color_plot_options

# Read the data
pred_df <- readRDS("pred_tables/pred_error_range2_5000_5000_complete.RDS") 

# Initial data preparation
pred_df <- pred_df |> 
  rename(Order = m, Method = method) |> 
  dplyr::filter(Order %in% 2:6, abs(sigma_e - 0.1) < 1e-10) # Change to sigma_e - sqrt(0.1) if needed

# Method mapping
method_mapping <- c(
  "rational" = "Rational",
  "pca" = "PCA",
  "fourier" = "Fourier",
  "nngp" = "nnGP",
  "statespace" = "State-Space",
  "fem" = "FEM",
  "taper" = "Taper"
)

pred_df$Method <- method_mapping[tolower(pred_df$Method)]

# Create long format data with both error types
pred_df <- pred_df |>
  tidyr::pivot_longer(
    cols = c(mu_error, sigma_error),
    names_to = "Error_Type",
    values_to = "Error"
  ) |>
  mutate(
    Error_Type = factor(Error_Type,
                       levels = c("mu_error", "sigma_error"),
                       labels = c(
                         expression(paste("Posterior mean error")),
                         expression(paste("Posterior std. dev. error"))
                       ))
  )

# Handle small error values
pred_df <- pred_df |> mutate(Error = ifelse(Error < 1e-10, 1e-10, Error))

# Filter data
df_filtered <- pred_df |>
  dplyr::filter((Method == "Rational") | (Order %in% c(3, 5)))

df_filtered <- df_filtered |>
  mutate(Facet_Cols = case_when(
    Method == "Rational" & Order == 3 ~ "Order3",
    Method == "Rational" & Order == 5 ~ "Order5",
    Method == "Rational" ~ "Rational",
    Order == 3 ~ "Order3",
    Order == 5 ~ "Order5",
    TRUE ~ NA_character_
  ))

df_filtered_tmp <- df_filtered |> 
  dplyr::filter(Method == "Rational", Order %in% c(3, 5)) |> 
  mutate(Facet_Cols = ifelse(Method == "Rational", "Rational", Method))

df_filtered <- bind_rows(df_filtered, df_filtered_tmp)

df_filtered <- df_filtered |> 
  mutate(Facet_Cols = factor(Facet_Cols, 
                            levels = c("Rational", "Order3", "Order5"),
                            labels = c("Rational", expression(Order == 3), expression(Order == 5))))

df_filtered$Order <- factor(df_filtered$Order)

df_filtered$Method <- factor(df_filtered$Method, 
                           levels = c("Rational", "PCA", "Fourier", "nnGP", 
                                    "State-Space", "FEM", "Taper"))

line_sizes <- c("Rational" = 0.7, "nnGP" = 0.7, "State-Space" = 0.7, 
                "Fourier" = 0.7, "PCA" = 0.7, "FEM" = 0.7, "Taper" = 0.7)

df_filtered$Order <- factor(df_filtered$Order, 
                          levels = 2:6, 
                          labels = paste("Order =", 2:6))

line_types_order <- c("Order = 2" = "dotted", 
                     "Order = 3" = "dashed", 
                     "Order = 4" = "twodash", 
                     "Order = 5" = "solid", 
                     "Order = 6" = "dotdash")

# Create the plot
p <- ggplot(df_filtered, aes(x = nu, y = Error, color = Method, 
                            linetype = Order, size = Method)) +
  geom_line() + 
  scale_y_log10(limits = c(1e-10, NA)) +
  scale_color_manual(values = color_plot_used, guide = guide_legend(order = 1)) +  
  scale_linetype_manual(values = line_types_order, guide = guide_legend(order = 2, nrow = 2)) +  
  scale_size_manual(values = line_sizes) +
  labs(y = "Prediction Error", x = expression(nu ~ "(smoothness parameter)")) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "lightgray"),
    panel.grid.minor = element_line(color = "gray"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  facet_grid(
    rows = vars(Error_Type),
    cols = vars(Facet_Cols),
    labeller = label_parsed
  ) +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y = element_text(size = 16, face = "bold"),
    strip.text.x = element_text(size = 16, face = "bold"),
    legend.box = "horizontal",  
    legend.box.just = "center", 
    legend.spacing.x = unit(0.5, "cm"),  
    legend.margin = margin(t = 0, b = 0), 
    legend.title.align = 0.5,
    legend.key.width = unit(1.5, "cm"),  
    legend.key.height = unit(0.5, "cm")
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 1)), 
    linetype = guide_legend(order = 2, override.aes = list(size = 1), nrow = 2), 
    size = "none"  
  )

print(p)
ggsave("pred_tables/pred_error_5000_5000_range2.png", p, width = 14, height = 8)