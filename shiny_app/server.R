library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)


dist_df <- readRDS("../distance_tables/full_dists.RDS") |> rename(Error_not_true = Error)
pred_df <- readRDS("../pred_tables/pred_error.RDS") |> 
  mutate(across(where(is.numeric), ~ ifelse(. == 0, 1e-10, .)))  # Replace 0 with 1e-10
prob_df <- readRDS("../prob_tables/prob_errors.RDS") |> mutate(Error = abs(prob_error))


methods <- c("Rational", "State-Space", "nnGP", "Fourier", "PCA", "Taper", "FEM")
orders <- as.character(0:6)


color_mapping_methods <- setNames(c("black", "steelblue", "limegreen", "red", "purple", "orange", "brown"), methods)
linetype_mapping_orders <- setNames(c("solid", "dashed", "dotted", "aa", "longdash", "twodash", "dotdash"), orders)

color_mapping_orders <- setNames(c("black", "steelblue", "limegreen", "red", "purple", "orange", "brown"), orders)
linetype_mapping_methods <- setNames(c("solid", "dashed", "dotted", "aa", "longdash", "twodash", "dotdash"), methods)

set_factor_levels <- function(df) {
  df %>%
    mutate(
      Method = factor(Method, levels = methods),
      m = factor(as.character(m), levels = orders)
    )
}

apply_aesthetic_mappings <- function(plotStyle) {
  if (plotStyle == "Linetype (Method), color (Order)") {
    list(
      aes(linetype = Method, color = m),
      scale_color_manual(values = unlist(color_mapping_orders)),
      scale_linetype_manual(values = unlist(linetype_mapping_methods))
    )
  } else {
    list(
      aes(linetype = m, color = Method),
      scale_color_manual(values = unlist(color_mapping_methods)),
      scale_linetype_manual(values = unlist(linetype_mapping_orders))
    )
  }
}

server <- function(input, output, session) {

  observe({
    output$raterrors <- renderPlotly({
      range_val <- input$rangeParameter
      approxNorm <- ifelse(input$approxNorm == "Sup", "Linf", input$approxNorm)
      plotStyle <- input$plotStyle
      m_order <- as.character(input$orderRat)
      methods_approx <- c("Rational", input$methods_approx)
      numberLoc <- as.numeric(input$numberLoc)
      numberObs <- as.numeric(input$numberObs)
      nuRange <- input$nuRange

      error_data <- if (input$trueCalibrationCov) {
        dist_df |> mutate(Error = True_Error_nnGP_Taper)
      } else {
        dist_df |> mutate(Error = Error_not_true)
      }

      filtered_data <- error_data %>%
        filter(
          range == range_val, Dist == approxNorm, 
          nu >= nuRange[1], nu <= nuRange[2],
          N == numberLoc, n_obs == numberObs,
          Method %in% methods_approx, m %in% m_order
        ) %>%
        set_factor_levels()

      aesthetics <- apply_aesthetic_mappings(plotStyle)

      fig <- ggplot(filtered_data, aes(x = nu, y = Error)) +
        geom_line(aesthetics[[1]]) +
        aesthetics[[2]] + aesthetics[[3]] +
        labs(
          y = ifelse(approxNorm == "L2", "Error in L2-norm", "Error in Sup-norm"),
          x = "\u028b (smoothness parameter)"
        )

      if (input$logScaleCoverror) {
        fig <- fig + scale_y_log10()
      }

      ggplotly(fig, height = 755)
    })
  })

  observe({
    output$prederrors <- renderPlotly({
      range_val <- input$rangeParameter_pred
      plotStyle <- input$plotStyle_pred
      m_order <- as.character(input$orderRat_pred)
      methods_approx <- c("Rational", input$methods_approx_pred)
      numberLoc <- as.numeric(input$numberLoc_pred)
      numberObs <- min(as.numeric(input$numberObs_pred), as.numeric(input$numberLoc_pred))
      nuRange <- input$nuRange_pred

      error_data <- if (input$trueCalibrationPred) {
        pred_df |> mutate(Error = true_pred_error_nnGP_Taper)
      } else {
        pred_df |> mutate(Error = pred_error)
      }

      filtered_data <- error_data %>%
        filter(
          range == range_val, 
          nu >= nuRange[1], nu <= nuRange[2],
          N == numberLoc, n_obs == numberObs,
          Method %in% methods_approx, m %in% m_order
        ) %>%
        set_factor_levels()

      aesthetics <- apply_aesthetic_mappings(plotStyle)

      fig <- ggplot(filtered_data, aes(x = nu, y = Error)) +
        geom_line(aesthetics[[1]]) +
        aesthetics[[2]] + aesthetics[[3]] +
        labs(y = "Error in l2-norm", x = "\u028b (smoothness parameter)")

      if (input$logScalepredError) {
        fig <- fig + scale_y_log10()
      }

      ggplotly(fig, height = 755)
    })
  })


   observe({
  output$proberrors <- renderPlotly({
    range_val_prob <- input$rangeParameter_prob
    plotStyle_prob <- input$plotStyle_prob
    m_order_prob <- input$orderRat_prob
    selected_methods <- c("Rational", input$methods_prob) 
    
    filtered_data <- prob_df %>% mutate(N = as.numeric(N)) %>%
      filter(
        range == range_val_prob,
        Method %in% selected_methods,
        m %in% as.character(m_order_prob)
      ) %>%
      set_factor_levels()

    aesthetics <- apply_aesthetic_mappings(plotStyle_prob)

    fig <- ggplot(filtered_data, aes(x = N, y = Error)) +
      geom_line(aesthetics[[1]]) +
      aesthetics[[2]] + aesthetics[[3]] +
      labs(y = "Probability Error", x = "N")

    if (input$logScaleprobError) {
      fig <- fig + scale_y_log10()
    }

    ggplotly(fig, height = 755)
  })
   })

}
