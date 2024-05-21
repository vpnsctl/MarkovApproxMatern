library(shiny)
library(dplyr)
library(tidyr)

dist_range_02 <- readRDS("../distance_tables/dist_range_02.RDS")
dist_range_05 <- readRDS("../distance_tables/dist_range_05.RDS")
dist_range_1 <- readRDS("../distance_tables/dist_range_1.RDS")

dist_range_02_samp <- readRDS("../distance_tables/dist_range_02_samp.RDS")
dist_range_05_samp <- readRDS("../distance_tables/dist_range_05_samp.RDS")
dist_range_1_samp <- readRDS("../distance_tables/dist_range_1_samp.RDS")


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
  
   })
  
}