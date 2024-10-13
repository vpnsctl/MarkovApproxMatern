library(tidyverse)

prob_df <- data.frame()

prob_tmp <- readRDS("prob_tables/full_results/rat_vs_nngp_range0.5_nu1.RDS")

prob_df <- prob_tmp$err_rat
colnames(prob_df) <- prob_tmp$n_vec
rownames(prob_df) <- prob_tmp$m_vec

prob_df <- as.data.frame(prob_df) %>%
  rownames_to_column(var = "m") %>%  
  pivot_longer(-m, names_to = "N", values_to = "prob_error")  

prob_df[["Method"]] <- "Rational"
prob_df[["range"]] <- "0.5"

prob_df_tmp <- prob_tmp$err_nngp
colnames(prob_df_tmp) <- prob_tmp$n_vec
rownames(prob_df_tmp) <- prob_tmp$m_vec

prob_df_tmp <- as.data.frame(prob_df_tmp) %>%
  rownames_to_column(var = "m") %>%  
  pivot_longer(-m, names_to = "N", values_to = "prob_error")  

prob_df_tmp[["Method"]] <- "nnGP"
prob_df_tmp[["range"]] <- "0.5"

prob_df <- bind_rows(prob_df, prob_df_tmp)

prob_tmp <- readRDS("prob_tables/full_results/rat_vs_nngp_range1_nu1.RDS")

prob_df_tmp <- prob_tmp$err_rat
colnames(prob_df_tmp) <- prob_tmp$n_vec
rownames(prob_df_tmp) <- prob_tmp$m_vec

prob_df_tmp <- as.data.frame(prob_df_tmp) %>%
  rownames_to_column(var = "m") %>%  
  pivot_longer(-m, names_to = "N", values_to = "prob_error")  

prob_df_tmp[["Method"]] <- "Rational"
prob_df_tmp[["range"]] <- "1"

prob_df <- bind_rows(prob_df, prob_df_tmp)

prob_df_tmp <- prob_tmp$err_nngp
colnames(prob_df_tmp) <- prob_tmp$n_vec
rownames(prob_df_tmp) <- prob_tmp$m_vec

prob_df_tmp <- as.data.frame(prob_df_tmp) %>%
  rownames_to_column(var = "m") %>%  
  pivot_longer(-m, names_to = "N", values_to = "prob_error")  

prob_df_tmp[["Method"]] <- "nnGP"
prob_df_tmp[["range"]] <- "1"

prob_df <- bind_rows(prob_df, prob_df_tmp)

prob_tmp <- readRDS("prob_tables/full_results/rat_vs_nngp_range2_nu1.RDS")

prob_df_tmp <- prob_tmp$err_rat
colnames(prob_df_tmp) <- prob_tmp$n_vec
rownames(prob_df_tmp) <- prob_tmp$m_vec

prob_df_tmp <- as.data.frame(prob_df_tmp) %>%
  rownames_to_column(var = "m") %>%  
  pivot_longer(-m, names_to = "N", values_to = "prob_error")  

prob_df_tmp[["Method"]] <- "Rational"
prob_df_tmp[["range"]] <- "2"

prob_df <- bind_rows(prob_df, prob_df_tmp)

prob_df_tmp <- prob_tmp$err_nngp
colnames(prob_df_tmp) <- prob_tmp$n_vec
rownames(prob_df_tmp) <- prob_tmp$m_vec

prob_df_tmp <- as.data.frame(prob_df_tmp) %>%
  rownames_to_column(var = "m") %>%  
  pivot_longer(-m, names_to = "N", values_to = "prob_error")  

prob_df_tmp[["Method"]] <- "nnGP"
prob_df_tmp[["range"]] <- "2"


prob_df <- bind_rows(prob_df, prob_df_tmp)

saveRDS(prob_df, "prob_tables/prob_errors.RDS")


### Adding fem part


prob_df <- readRDS("prob_tables//prob_errors.RDS")

prob_tmp <- readRDS("prob_tables/full_results//fem_range0.5_nu1.RDS")

prob_df_tmp <- prob_tmp$err_fem
colnames(prob_df_tmp) <- prob_tmp$n_vec
rownames(prob_df_tmp) <- prob_tmp$m_vec

prob_df_tmp <- as.data.frame(prob_df_tmp) %>%
  rownames_to_column(var = "m") %>%  
  pivot_longer(-m, names_to = "N", values_to = "prob_error")  

prob_df_tmp[["Method"]] <- "FEM"
prob_df_tmp[["range"]] <- "0.5"

prob_df <- bind_rows(prob_df, prob_df_tmp)

prob_tmp <- readRDS("prob_tables/full_results//fem_range1_nu1.RDS")

prob_df_tmp <- prob_tmp$err_fem
colnames(prob_df_tmp) <- prob_tmp$n_vec
rownames(prob_df_tmp) <- prob_tmp$m_vec

prob_df_tmp <- as.data.frame(prob_df_tmp) %>%
  rownames_to_column(var = "m") %>%  
  pivot_longer(-m, names_to = "N", values_to = "prob_error")  

prob_df_tmp[["Method"]] <- "FEM"
prob_df_tmp[["range"]] <- "1"

prob_df <- bind_rows(prob_df, prob_df_tmp)

prob_tmp <- readRDS("prob_tables/full_results//fem_range2_nu1.RDS")

prob_df_tmp <- prob_tmp$err_fem
colnames(prob_df_tmp) <- prob_tmp$n_vec
rownames(prob_df_tmp) <- prob_tmp$m_vec

prob_df_tmp <- as.data.frame(prob_df_tmp) %>%
  rownames_to_column(var = "m") %>%  
  pivot_longer(-m, names_to = "N", values_to = "prob_error")  

prob_df_tmp[["Method"]] <- "FEM"
prob_df_tmp[["range"]] <- "2"

prob_df <- bind_rows(prob_df, prob_df_tmp)

