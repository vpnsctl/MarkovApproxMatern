library(ggplot2)
library(dplyr)

color_plot_options <- c("black", "steelblue", "limegreen", "red", "purple", "orange", "brown")
color_plot_used <- color_plot_options

dist_df <- readRDS("distance_tables/full_dists.RDS") %>% 
  rename(Order = m, Range = range) %>% 
  filter(Order %in% 2:6, N == 5000, n_obs == 5000, Range == 2)

df_filtered <- dist_df %>%
  filter((Method == "Rational") | (Order %in% c(3, 5)))

df_filtered <- df_filtered %>% mutate(Dist = ifelse(Dist == "L2", "L2", "Linf"))

df_filtered <- df_filtered %>%
  mutate(Facet_Cols = case_when(
    Method == "Rational" & Order == 3 ~ "Order3",
    Method == "Rational" & Order == 5 ~ "Order5",
    Method == "Rational" ~ "Rational",
    Order == 3 ~ "Order3",
    Order == 5 ~ "Order5",
    TRUE ~ NA_character_
  ))

df_filtered_tmp <- df_filtered %>% 
  filter(Method == "Rational", Order %in% c(3, 5)) %>% 
  mutate(Facet_Cols = ifelse(Method == "Rational", "Rational", Method))

df_filtered <- bind_rows(df_filtered, df_filtered_tmp)

df_filtered$Facet_Cols <- factor(df_filtered$Facet_Cols, levels = c("Rational", "Order3", "Order5"))

markers <- c("Rational" = 15, "nnGP" = 3, "State-Space" = 5, "Fourier" = 5, "PCA" = 15)

df_filtered <- df_filtered %>% 
  mutate(Dist = factor(Dist, levels = c("L2", "Linf"), 
                       labels = c(expression(paste("Error in ", italic(L)[2], "(I\u00D7I)-norm")),
                                   expression(paste("Error in ", italic(L)[infinity], "(I\u00D7I)-norm")))))

df_points <- df_filtered %>% arrange(Order) %>% slice(seq(1, n(), by = 30))

df_filtered$Order <- factor(df_filtered$Order)

p <- ggplot(df_filtered, aes(x = nu, y = Error, color = Order, shape = Method)) +
  geom_line(size = 1.2) +
  geom_point(data = df_points, size = 3) +  # Increase dodge width for sparse markers
  scale_y_log10(limits = c(1e-13, NA)) +  
  scale_color_manual(values = color_plot_used) +
  scale_shape_manual(values = markers) +
  labs(y = "Covariance Error", x = expression(nu ~ "(smoothness parameter)")) +
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
  facet_grid(rows = vars(Dist), 
             cols = vars(Facet_Cols), 
             labeller = label_parsed) + 
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside") +
  theme(strip.text.y = element_text(size = 14, face = "bold"))  # Add back facet column titles

print(p)
ggsave("facet_plot.png", p, width = 12, height = 8)
