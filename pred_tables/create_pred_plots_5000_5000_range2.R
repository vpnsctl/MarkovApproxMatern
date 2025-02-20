library(ggplot2)
library(dplyr)

color_plot_options <- c("black", "steelblue", "limegreen", "red", "purple", "orange", "brown")
color_plot_used <- color_plot_options

pred_df <- readRDS("pred_tables/pred_error_range2_5000_5000_complete.RDS") 

# # First, let's create the function to check and replace outliers
# check_and_replace <- function(df) {
#   # Sort by nu
#   df <- df[order(df$nu), ]
  
#   # Initialize vectors to store cleaned values
#   cleaned_mu <- df$mu_error
#   cleaned_sigma <- df$sigma_error
  
#   # For each nu value
#   for(i in 1:nrow(df)) {
#     # Find closest nu values
#     nu_diffs <- abs(df$nu - df$nu[i])
#     nu_diffs[i] <- Inf  # Exclude current point
#     closest_indices <- which(rank(nu_diffs) <= 2)[1:2]  # Get two closest
    
#     # Calculate averages of closest points
#     avg_mu <- mean(df$mu_error[closest_indices])
#     avg_sigma <- mean(df$sigma_error[closest_indices])
    
#     # Replace if current value is more than 10 times the average
#     if(df$mu_error[i] > 1.2 * avg_mu) {
#       cleaned_mu[i] <- avg_mu
#     }
#     if(df$sigma_error[i] > 1.2 * avg_sigma) {
#       cleaned_sigma[i] <- avg_sigma
#     }
#   }
  
#   df$mu_error <- cleaned_mu
#   df$sigma_error <- cleaned_sigma
#   return(df)
# }

# # Split data frame into list by groups
# split_df <- split(pred_df, list(pred_df$method, pred_df$m, pred_df$sigma_e))

# # Apply function to each group
# cleaned_list <- lapply(split_df, function(group) {
#   if(!(unique(group$method) %in% c("pca", "fourier", "statespace"))) {
#     return(check_and_replace(group))
#   } else {
#     return(group)
#   }
# })

# # Combine back into a single data frame
# pred_df <- do.call(rbind, cleaned_list)
# rownames(pred_df) <- NULL


pred_df <- pred_df |> 
  rename(Order = m, Method = method) |> 
  dplyr::filter(Order %in% 2:6, abs(sigma_e - 0.1) < 1e-10)



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

pred_df <- pred_df |> rename(Error = sigma_error)

pred_df <- pred_df |> mutate(Error = ifelse(Error < 1e-10, 1e-10, Error))

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

df_filtered$Facet_Cols <- factor(df_filtered$Facet_Cols, 
  levels = c("Rational", "Order3", "Order5"), 
  labels = c("Rational", "Order = 3", "Order = 5"))

df_filtered$Order <- factor(df_filtered$Order)

df_filtered$Method <- factor(df_filtered$Method, levels = c("Rational", "PCA", "Fourier", "nnGP", "State-Space", "FEM", "Taper"))


line_sizes <- c("Rational" = 0.7, "nnGP" = 0.7, "State-Space" = 0.7, "Fourier" = 0.7, "PCA" = 0.7, "FEM" = 0.7, "Taper" = 0.7)

df_filtered$Order <- factor(df_filtered$Order, levels = 2:6, labels = paste("Order =", 2:6))

line_types_order <- c("Order = 2" = "dotted", "Order = 3" = "dashed", "Order = 4" = "twodash", "Order = 5" = "solid", "Order = 6" = "dotdash")


p <- ggplot(df_filtered, aes(x = nu, y = Error, color = Method, linetype = Order, size = Method)) +
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
  facet_grid(~Facet_Cols) +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y = element_text(size = 16),#, face = "bold"),
    strip.text.x = element_text(size = 16),# face = "bold"),
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
ggsave("pred_tables/pred_error_5000_5000_range2.png", p, width = 14, height = 5)
