library(shiny)
library(dplyr)
library(tidyr)

dist_df <- readRDS("../distance_tables/full_dists.RDS")
pred_df <- readRDS("../pred_tables/pred_error.RDS")
prob_df <- readRDS("../prob_tables/prob_errors.RDS")

dist_df <- dist_df |> rename(Error_not_true = Error)

prob_df <- prob_df |> rename(Error = prob_error)

prob_df <- prob_df |> mutate(Error = abs(Error))


server <- function(input, output, session) {
  
  observe({
    range_val <- input$rangeParameter
    approxNorm <- input$approxNorm
    plotStyle <- input$plotStyle
    if(approxNorm == "Sup"){
      approxNorm <- "Linf"
    }
    m_order <- input$orderRat
    methods_approx <- input$methods_approx
    methods_approx <- c("Rational", methods_approx)
    numberLoc <- as.numeric(input$numberLoc)
    numberObs <- as.numeric(input$numberObs)
    nuRange <- input$nuRange

    true_cal <- input$trueCalibrationCov

    if(true_cal){
      dist_df <- dist_df |> mutate(Error = True_Error_nnGP_Taper)
    } else{
      dist_df <- dist_df |> mutate(Error = Error_not_true)
    }
        
    color_plot_options <- c("black", "steelblue", 
                    "limegreen", "red",
                    "purple","orange","pink")
    
    error_table <- dist_df |> dplyr::filter(range == range_val)
        
    if(!is.null(m_order)){
      color_plot_used <- color_plot_options[as.numeric(m_order)+1]
      
      error_table <- error_table %>% dplyr::filter(m %in% as.character(m_order), Dist == approxNorm, nu <= nuRange[2], nu >= nuRange[1], N == numberLoc, n_obs == numberObs)

      idx_table <- (error_table[["Method"]]%in%methods_approx)

      error_table <- error_table[idx_table,] 
      
      if(plotStyle == "Linetype (Method), color (Order)"){
        fig <- ggplot(error_table, aes(x = nu, y = Error, linetype=Method, color = m)) + geom_line() +
          scale_color_manual(values = color_plot_used)
      } else{
        fig <- ggplot(error_table, aes(x = nu, y = Error, linetype=m, color = Method)) + geom_line() 
      }
            
      y_label_cov <- ifelse(approxNorm == "L2", "Error in L2-norm", "Error in Sup-norm")
      
      fig <- fig +labs(y= y_label_cov, x = "\u028b (smoothness parameter)")
      
      if(input$logScaleCoverror){
        fig <- fig + scale_y_log10()
      }
      
            output$downloadPlotCov <- downloadHandler(
          filename = function() { paste0("Cov_error_plot.", input$fileExtensionCov) },
          content = function(file) {
            ggsave(file, plot = fig, device = input$fileExtensionCov)
          }
        )
        
        output$downloadDataFrameCov <- downloadHandler(
          filename = function() { paste0("Cov_error_data_frame_range_",range_val,".RDS") },
          content = function(file) {
            saveRDS(error_table, file = file)
          }
        )
      
      fig_plotly <- ggplotly(fig,height = 755)
      
      
    } else{
      m_order = c(2)
      updateCheckboxGroupInput(session, "orderRat", selected= c("2"))
    }
    
    output$raterrors <- renderPlotly({
      
      fig_plotly
      
    })
    
    
  })
  

  observe({
    range_val_pred <- input$rangeParameter_pred
    approxNorm_pred <- input$approxNorm_pred
    plotStyle_pred <- input$plotStyle_pred
    m_order_pred <- input$orderRat_pred
    methods_approx_pred <- input$methods_approx_pred
    methods_approx_pred <- c("Rational_LDL", methods_approx_pred)
    numberLoc_pred <- input$numberLoc_pred
    nuRange_pred <- input$nuRange_pred
    numberObs_pred <- min(as.numeric(input$numberObs_pred), as.numeric(input$numberLoc_pred))
    methods_approx_pred <- c("Rational", methods_approx_pred)

    true_cal_pred <- input$trueCalibrationPred

    if(true_cal_pred){
      pred_df <- pred_df |> mutate(Error = true_pred_error_nnGP_Taper)
    } else{
      pred_df <- pred_df |> mutate(Error = pred_error)
    }    
    
    color_plot_options <- c("black", "steelblue", 
                    "limegreen", "red",
                    "purple","orange","pink")
    
    pred_table <- pred_df |> dplyr::filter(range == range_val_pred)
    
    if(!is.null(m_order_pred)){
      color_plot_used <- color_plot_options[as.numeric(m_order_pred)+1]
      
      pred_table <- pred_table %>% dplyr::filter(m %in% as.character(m_order_pred), nu <= nuRange_pred[2], nu >= nuRange_pred[1], N == numberLoc_pred, n_obs == numberObs_pred)

      idx_table <- (pred_table[["Method"]]%in%methods_approx_pred)

      pred_table <- pred_table[idx_table,] 
      
      if(plotStyle_pred == "Linetype (Method), color (Order)"){
        fig <- ggplot(pred_table, aes(x = nu, y = Error, linetype=Method, color = m)) + geom_line() +
          scale_color_manual(values = color_plot_used)
      } else{
        fig <- ggplot(pred_table, aes(x = nu, y = Error, linetype=m, color = Method)) + geom_line() 
      }
            
      y_label_cov <- "Error in l2-norm"
      
      fig <- fig +labs(y= y_label_cov, x = "\u028b (smoothness parameter)")
      
      if(input$logScaleCoverror){
        fig <- fig + scale_y_log10()
      }
      
            output$downloadPlotPred <- downloadHandler(
          filename = function() { paste0("Pred_error_plot.", input$fileExtensionCov) },
          content = function(file) {
            ggsave(file, plot = fig, device = input$fileExtensionCov)
          }
        )
        
        output$downloadDataFramePred <- downloadHandler(
          filename = function() { paste0("Pred_error_data_frame_range_",range_val,".RDS") },
          content = function(file) {
            saveRDS(error_table, file = file)
          }
        )
      
      fig_plotly <- ggplotly(fig,height = 755)
      
      
    } else{
      m_order = c(2)
      updateCheckboxGroupInput(session, "orderRat_pred", selected= c("2"))
    }
    
    output$prederrors <- renderPlotly({
      
      fig_plotly
      
    })
   })

  
   observe({
    range_val_prob <- input$rangeParameter_prob
    approxNorm_prob <- input$approxNorm_prob
    plotStyle_prob <- input$plotStyle_prob
    m_order_prob <- input$orderRat_prob
    nuRange_prob <- 1
     
    color_plot_options <- c("black", "steelblue", 
                    "limegreen", "red",
                    "purple","orange","pink")
    
    prob_table <- prob_df |> dplyr::filter(range == range_val_prob)
     
    prob_table[["N"]] <- as.numeric(prob_table[["N"]])
    
    if(!is.null(m_order_prob)){
      color_plot_used <- color_plot_options[as.numeric(m_order_prob)+1]
      
      prob_table <- prob_table %>% dplyr::filter(m %in% as.character(m_order_prob))
      
      if(plotStyle_prob == "Linetype (Method), color (Order)"){
        fig <- ggplot(prob_table, aes(x = N, y = Error, linetype=Method, color = m)) + geom_line() +
          scale_color_manual(values = color_plot_used)
      } else{
        fig <- ggplot(prob_table, aes(x = N, y = Error, linetype=m, color = Method)) + geom_line() 
      }
            
      y_label_cov <- "Probability error"
      
      fig <- fig +labs(y= y_label_cov, x = "N")
      
      if(input$logScaleprobError){
        fig <- fig + scale_y_log10()
      }
      
            output$downloadPlotprob <- downloadHandler(
          filename = function() { paste0("Prob_error_plot.", input$fileExtensionCov) },
          content = function(file) {
            ggsave(file, plot = fig, device = input$fileExtensionCov)
          }
        )
        
        output$downloadDataFrameprob <- downloadHandler(
          filename = function() { paste0("Prob_error_data_frame_range_",range_val,".RDS") },
          content = function(file) {
            saveRDS(error_table, file = file)
          }
        )
      
      fig_plotly <- ggplotly(fig,height = 755)
      
      
    } else{
      m_order = c(2)
      updateCheckboxGroupInput(session, "orderRat_prob", selected= c("2"))
    }
    
    output$proberrors <- renderPlotly({
      
      fig_plotly
      
    })
   })
  
}