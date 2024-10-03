library(ggplot2)
library(dplyr)

dist_df <- readRDS("distance_tables/full_dists.RDS") %>% 
  rename(Order = m, Range = range) %>% 
  filter(Order %in% 2:6, N == 5000, n_obs == 5000, Range == 2)

df_filtered <- dist_df %>%
  filter((Method == "Rational") | (Order %in% c(3, 5)))

df_filtered <- df_filtered %>%
  mutate(Facet_Cols = case_when(
    Method == "Rational" & Order == 3 ~ "Order 3",
    Method == "Rational" & Order == 5 ~ "Order 5",
    Method == "Rational" ~ "Rational",
    Order == 3 ~ "Order 3",
    Order == 5 ~ "Order 5",
    TRUE ~ NA_character_
  ))

df_filtered_tmp <- df_filtered %>% 
  filter(Method == "Rational", Order %in% c(3, 5)) %>% 
  mutate(Facet_Cols = ifelse(Method == "Rational", "Rational", Method))

df_filtered <- bind_rows(df_filtered, df_filtered_tmp)

df_filtered$Facet_Cols <- factor(df_filtered$Facet_Cols, levels = c("Rational", "Order 3", "Order 5"))

# Manually set line types for each method
line_types <- c("Rational" = "solid", 
                "nnGP" = "dotted", 
                "State-Space" = "dashed", 
                "Fourier" = "twodash", 
                "PCA" = "dotdash")

p <- ggplot(df_filtered, aes(x = nu, y = Error, color = Order, linetype = Method)) +
  geom_line() +
  scale_y_log10() +
  scale_color_manual(values = color_plot_used) +
  scale_linetype_manual(values = line_types) +  # Set line types manually
  labs(y = "Covariance Error", x = expression(nu ~ "(smoothness parameter)")) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14)) +
  facet_grid(rows = vars(Dist), 
             cols = vars(Facet_Cols), 
             labeller = as_labeller(c("L2" = expression('Error in L'[2] ~ "(\U1D49F\U00D7\U1D49F)-norm"),
                                       "Linf" = expression('Error in L'["∞"] ~ "(\U1D49F\U00D7\U1D49F)-norm")))) +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside")

print(p)

ggsave("facet_plot.png", p, width = 12, height = 8)
