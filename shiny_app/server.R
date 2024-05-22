library(shiny)
library(dplyr)
library(tidyr)

dist_range_02 <- readRDS("../distance_tables/dist_range_02.RDS")
dist_range_05 <- readRDS("../distance_tables/dist_range_05.RDS")
dist_range_1 <- readRDS("../distance_tables/dist_range_1.RDS")

dist_range_02_samp <- readRDS("../distance_tables/dist_range_02_samp.RDS")
dist_range_05_samp <- readRDS("../distance_tables/dist_range_05_samp.RDS")
dist_range_1_samp <- readRDS("../distance_tables/dist_range_1_samp.RDS")

pred_range02 <- readRDS("../distance_tables/pred_range02.RDS")


server <- function(input, output, session) {
  
  observe({
    which_calib <- input$whichCalib
    range_val <- input$rangeParameter
    approxNorm <- input$approxNorm
    plotStyle <- input$plotStyle
    if(approxNorm == "Sup"){
      approxNorm <- "Linf"
    }
    m_order <- input$orderRat
    methods_approx <- input$methods_approx
    methods_approx <- c("Rational", methods_approx)
    numberLoc <- input$numberLoc
    nuRange <- input$nuRange
    
    if(which_calib == "Sampling"){
      error_02 <- dist_range_02_samp
      error_05 <- dist_range_05_samp
      error_1 <- dist_range_1_samp
    } else{
      error_02 <- dist_range_02
      error_05 <- dist_range_05
      error_1 <- dist_range_1
    }
    
    color_plot_options <- c("black", "steelblue", 
                    "limegreen", "red",
                    "purple","orange","pink")
    
    if(range_val == "0.2"){
      error_table <- error_02
    } else if(range_val == "0.5"){
      error_table <- error_05
    } else{
      error_table <- error_1
    }
    
    if(!is.null(m_order)){
      color_plot_used <- color_plot_options[as.numeric(m_order)+1]
      
      error_table <- error_table %>% dplyr::filter(m %in% as.character(m_order), Dist == approxNorm, nu <= nuRange[2], nu >= nuRange[1], N == numberLoc)

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
    if(approxNorm_pred == "Max"){
      approxNorm_pred <- "max"
    }
    m_order_pred <- input$orderRat_pred
    methods_approx_pred <- input$methods_approx_pred
    methods_approx_pred <- c("Rational_LDL", methods_approx_pred)
    numberLoc_pred <- input$numberLoc_pred
    nuRange_pred <- input$nuRange_pred
    
    color_plot_options <- c("black", "steelblue", 
                    "limegreen", "red",
                    "purple","orange","pink")
    
    if(range_val_pred == "0.2"){
      pred_table <- pred_range02
    } else if(range_val_pred == "0.5"){
      pred_table <- pred_range05
    } else{
      pred_table <- pred_range1
    }
    
    if(!is.null(m_order_pred)){
      color_plot_used <- color_plot_options[as.numeric(m_order_pred)+1]
      
      pred_table <- pred_table %>% dplyr::filter(m %in% as.character(m_order_pred), Norm == approxNorm_pred, nu <= nuRange_pred[2], nu >= nuRange_pred[1], N == numberLoc_pred)

      idx_table <- (pred_table[["Method"]]%in%methods_approx_pred)

      pred_table <- pred_table[idx_table,] 
      
      if(plotStyle_pred == "Linetype (Method), color (Order)"){
        fig <- ggplot(pred_table, aes(x = nu, y = Error, linetype=Method, color = m)) + geom_line() +
          scale_color_manual(values = color_plot_used)
      } else{
        fig <- ggplot(pred_table, aes(x = nu, y = Error, linetype=m, color = Method)) + geom_line() 
      }
            
      y_label_cov <- ifelse(approxNorm_pred == "l2", "Error in l2-norm", "Error in Max-norm")
      
      fig <- fig +labs(y= y_label_cov, x = "\u028b (smoothness parameter)")
      
      if(input$logScaleCoverror){
        fig <- fig + scale_y_log10()
      }
      
            output$downloadPlotCov <- downloadHandler(
          filename = function() { paste0("Pred_error_plot.", input$fileExtensionCov) },
          content = function(file) {
            ggsave(file, plot = fig, device = input$fileExtensionCov)
          }
        )
        
        output$downloadDataFrameCov <- downloadHandler(
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
  
}