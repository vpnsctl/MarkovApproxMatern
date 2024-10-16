library(ggplot2)
library(dplyr)

# Define the color palette for Methods
color_plot_options <- c("black", "steelblue", "limegreen", "red", "purple", "orange", "brown")
color_plot_used <- color_plot_options[1:length(unique(df_filtered$Method))]

# Load and preprocess data
prob_df <- readRDS("prob_tables/prob_errors.RDS") |> 
  rename(Order = m, Range = range) |> 
  rename(Error = prob_error) |> 
  mutate(Error = abs(Error))

# Filter relevant orders
df_filtered <- prob_df |> filter(Order > 2)

# Ensure factors and numerical conversion
df_filtered$Order <- factor(df_filtered$Order)
df_filtered$N <- as.numeric(df_filtered$N)

# Define the factor levels for Methods
df_filtered$Method <- factor(df_filtered$Method, levels = c("Rational", "nnGP", "FEM"))

# Define the factor levels for Range
df_filtered$Range <- factor(df_filtered$Range, 
                            levels = c("0.5", "1", "2"), 
                            labels = c(expression(rho ~ "= 0.5"), expression(rho ~ "= 1"), expression(rho ~ "= 2")))

# Define line types for Orders
line_types <- c("2" = "solid", 
                "3" = "dashed", 
                "4" = "dotted", 
                "5" = "dotdash", 
                "6" = "twodash")

# Create the plot with reversed aesthetics: color as Method and linetype as Order
p <- ggplot(df_filtered, aes(x = N, y = Error, color = Method, linetype = Order)) +
  geom_line(size = 0.7) + 
  scale_color_manual(values = color_plot_used, guide = guide_legend(order = 1)) +
  scale_linetype_manual(values = line_types, guide = guide_legend(order = 2)) +  
  labs(y = "Posterior Probability Error", x = "N (number of prediction locations)") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14),
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "lightgray"),
    panel.grid.minor = element_line(color = "gray"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  facet_grid(~Range, labeller = label_parsed) +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y = element_text(size = 14, face = "bold"),
    legend.box = "horizontal", 
    legend.box.just = "center",  
    legend.spacing.x = unit(0.5, "cm"),  
    legend.margin = margin(t = 0, b = 0), 
    legend.title.align = 0.5
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 1)), 
    linetype = guide_legend(order = 2, override.aes = list(size = 0.7)), 
    size = "none"
  )

# Print and save the plot
print(p)
ggsave("prob_tables/prob_error.png", p, width = 12, height = 5)
