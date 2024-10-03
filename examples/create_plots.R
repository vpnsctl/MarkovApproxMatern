library(ggplot2)
library(dplyr)

# Define color options
color_plot_options <- c("black", "steelblue", "limegreen", "red", "purple", "orange", "brown")
color_plot_used <- color_plot_options

# Load and preprocess the data
dist_df <- readRDS("distance_tables/full_dists.RDS") %>% 
  rename(Order = m, Range = range) %>% 
  filter(Order %in% 2:6, N == 5000, n_obs == 5000, Range == 2)

df_filtered <- dist_df %>%
  filter((Method == "Rational") | (Order %in% c(3, 5)))

# Adjust Dist variable
df_filtered <- df_filtered %>% mutate(Dist = ifelse(Dist == "L2", "L2", "Linf"))

# Create the Facet_Cols variable
df_filtered <- df_filtered %>%
  mutate(Facet_Cols = case_when(
    Method == "Rational" & Order == 3 ~ "Order3",
    Method == "Rational" & Order == 5 ~ "Order5",
    Method == "Rational" ~ "Rational",
    Order == 3 ~ "Order3",
    Order == 5 ~ "Order5",
    TRUE ~ NA_character_
  ))

# Include Rational method again for proper facet organization
df_filtered_tmp <- df_filtered %>% 
  filter(Method == "Rational", Order %in% c(3, 5)) %>% 
  mutate(Facet_Cols = ifelse(Method == "Rational", "Rational", Method))

df_filtered <- bind_rows(df_filtered, df_filtered_tmp)

# Set the order of Facet_Cols
df_filtered$Facet_Cols <- factor(df_filtered$Facet_Cols, levels = c("Rational", "Order3", "Order5"))

# Define line types for each method
line_types <- c("Rational" = "solid", 
                "nnGP" = "dotted", 
                "State-Space" = "dashed", 
                "Fourier" = "twodash", 
                "PCA" = "dotdash")

# Adjust Dist variable to use expressions
df_filtered <- df_filtered %>% mutate(Dist = factor(Dist, levels = c("L2", "Linf"), 
                              labels = c(expression(paste("Error in ", italic(L)[2], " (I %x% I)-norm")),
                                         expression(paste("Error in ", italic(L)[infinity], " (I %x% I)-norm")))))

# Plot with custom labeller for Dist
p <- ggplot(df_filtered, aes(x = nu, y = Error, color = Order, linetype = Method)) +
  geom_line() +
  scale_y_log10() +
  scale_color_manual(values = color_plot_used) +
  scale_linetype_manual(values = line_types) +
  labs(y = "Covariance Error", x = expression(nu ~ "(smoothness parameter)")) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),   # Increase font size for legend text
    legend.title = element_text(size = 14),   # Increase font size for legend title
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = 14),
    panel.background = element_blank(),         # Remove gray background
    panel.grid.major = element_line(color = "lightgray"), # Set major grid lines color to light gray
    panel.grid.minor = element_line(color = "gray"), # Set minor grid lines color to a darker gray
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a border around the plot
  ) +
  facet_grid(rows = vars(Dist), 
             cols = vars(Facet_Cols), 
             labeller = label_parsed) + 
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside")

# Print and save the plot
print(p)
ggsave("facet_plot.png", p, width = 12, height = 8)
