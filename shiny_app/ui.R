library(shiny)
library(plotly)
library(shinythemes)
library(shinyWidgets)

nu_val <- seq(0.01, 2.49, 0.01)

ui <- navbarPage(
  "Markov Rational Approximation",
  
  tabPanel("Covariance error",
           sidebarLayout(position = "left",
                         sidebarPanel(width = 3,
                                      checkboxGroupInput(inputId = "methods_approx", label = "Methods to be included",
                                                         choices = c("State-Space", "nnGP", "Fourier", "PCA"), 
                                                         selected = "State-Space",
                                                         inline = TRUE),                                                         
                                      checkboxGroupInput(inputId = "orderRat", label = "Order of the approximation",
                                                         choices = 0:6, selected = 0:6,
                                                         inline = TRUE),
                                      radioButtons(inputId = "numberLoc", label = "Number of locations",
                                                   choices = c("5000", "10000"),
                                                   selected = "10000",
                                                   inline = TRUE),                                           
                                      radioButtons(inputId = "rangeParameter", label = "Range Parameter",
                                                   choices = c("0.5", "1", "2"),
                                                   selected = "0.5",
                                                   inline = TRUE),
                                      radioButtons(inputId = "approxNorm", label = "Which norm?",
                                                   choices = c("L2", "Sup"),
                                                   selected = "L2",
                                                   inline = TRUE),
                                      sliderTextInput(inputId = "nuRange", label = "Which values of nu?",
                                                       choices = nu_val,
                                                       selected = range(nu_val)),
                                      radioButtons(inputId = "plotStyle", label = "Which style?",
                                                   choices = c("Linetype (Method), color (Order)", "Linetype (Order), color (Method)"),
                                                   selected = "Linetype (Method), color (Order)",
                                                   inline = FALSE),                                                   
                                      checkboxInput("logScaleCoverror", "Use log scale?", value = TRUE),
                                      downloadButton('downloadPlotCov', 'Download Plot'),
                                      radioButtons(inputId = "fileExtensionCov", label = "File extension",
                                                   choices = c("png", "pdf"),
                                                   selected = "png",
                                                   inline = TRUE),
                                      downloadButton('downloadDataFrameCov', 'Download Data Frame')
                         ),  # Removed extra comma here
                         mainPanel(
                           width = 9,  # Set to 9 for clarity
                           plotlyOutput("raterrors")
                         )
           )
  ),
  
  tabPanel("Prediction error", 
           sidebarLayout(position = "left",
                         sidebarPanel(width = 3,
                                      checkboxGroupInput(inputId = "methods_approx_pred", label = "Methods to be included",
                                                         choices = c("State-Space", "nnGP", "Fourier", "PCA"), 
                                                         selected = "State-Space",
                                                         inline = TRUE),                                                         
                                      checkboxGroupInput(inputId = "orderRat_pred", label = "Order of the approximation",
                                                         choices = 0:6, selected = 0:6,
                                                         inline = TRUE),
                                      radioButtons(inputId = "numberLoc_pred", label = "Number of locations",
                                                   choices = c("5000", "10000"),
                                                   selected = "10000",
                                                   inline = TRUE),
                                      radioButtons(inputId = "numberObs_pred", label = "Number of observations",
                                                   choices = c("5000", "10000"),
                                                   selected = "10000",
                                                   inline = TRUE),                                                   
                                      radioButtons(inputId = "rangeParameter_pred", label = "Range Parameter",
                                                   choices = c("0.5", "1", "2"),
                                                   selected = "0.5",
                                                   inline = TRUE),
                                      sliderTextInput(inputId = "nuRange_pred", label = "Which values of nu?",
                                                       choices = nu_val,
                                                       selected = range(nu_val)),
                                      radioButtons(inputId = "plotStyle_pred", label = "Which style?",
                                                   choices = c("Linetype (Method), color (Order)", "Linetype (Order), color (Method)"),
                                                   selected = "Linetype (Method), color (Order)",
                                                   inline = FALSE),  
                                      checkboxInput("logScalepredError", "Use log scale?", value = TRUE),
                                      downloadButton('downloadPlotPred', 'Download Plot'),
                                      radioButtons(inputId = "fileExtensionPred", label = "File extension",
                                                   choices = c("png", "pdf"),
                                                   selected = "png",
                                                   inline = TRUE),
                                      downloadButton('downloadDataFramePred', 'Download Data Frame')
                         ),  # Removed extra comma here
                         mainPanel(
                           width = 9,  # Set to 9 for clarity
                           plotlyOutput("prederrors")
                         )
           )
  ),
  
  tabPanel("Posterior probability error", 
           sidebarLayout(position = "left",
                         sidebarPanel(width = 3,     
                                      checkboxGroupInput(inputId = "orderRat_prob", label = "Order of the approximation",
                                                         choices = 1:6, selected = 1:6,
                                                         inline = TRUE),                                                
                                      radioButtons(inputId = "rangeParameter_prob", label = "Range Parameter",
                                                   choices = c("0.5", "1", "2"),
                                                   selected = "0.5",
                                                   inline = TRUE),
                                      radioButtons(inputId = "plotStyle_prob", label = "Which style?",
                                                   choices = c("Linetype (Method), color (Order)", "Linetype (Order), color (Method)"),
                                                   selected = "Linetype (Method), color (Order)",
                                                   inline = FALSE),  
                                      checkboxInput("logScaleprobError", "Use log scale?", value = TRUE),
                                      downloadButton('downloadPlotprob', 'Download Plot'),
                                      radioButtons(inputId = "fileExtensionprob", label = "File extension",
                                                   choices = c("png", "pdf"),
                                                   selected = "png",
                                                   inline = TRUE),
                                      downloadButton('downloadDataFrameprob', 'Download Data Frame')
                         ),  # Removed extra comma here
                         mainPanel(
                           width = 9,  # Set to 9 for clarity
                           plotlyOutput("proberrors")
                         )
           )
  ),

  tags$style(type = "text/css", "body {padding-top: 70px;}"),
  theme = shinytheme("cosmo"),
  position = "fixed-top"
)
