library(shiny)
library(plotly)
library(shinythemes)
library(shinyWidgets)


dist_range_02 <- readRDS("../distance_tables/dist_range_02.RDS")

nu_val <- sort(unique(dist_range_02[["nu"]]))


ui <- navbarPage(
  "Markov Rational Approxiomation",
  tabPanel("Covariance error",
           sidebarLayout(position = "left",
                         sidebarPanel(width = 3,
                                      radioButtons(inputId = "whichCalib", label = "Which calibration?",
                                                         choices = c("Sampling", "Prediction"), 
                                                         selected = "Sampling",
                                                         inline=TRUE),
                                      checkboxGroupInput(inputId = "methods_approx", label = "Methods to be included",
                                                         choices = c("State-Space", "nnGP", "Fourier", "PCA"), 
                                                         selected = "State-Space",
                                                         inline=TRUE),                                                         
                                      checkboxGroupInput(inputId = "orderRat", label = "Order of the approximation",choices = 0:6, selected = 0:6,
                                                         inline=TRUE),
                                      radioButtons(inputId = "numberLoc", label = "Number of locations",
                                                   choices = c("500", "1000"),
                                                   selected = "1000",
                                                   inline=TRUE),
                                      radioButtons(inputId = "rangeParameter", label = "Range Parameter",
                                                   choices = c("0.2", "0.5", "1"),
                                                   selected = "0.2",
                                                   inline=TRUE),
                                      radioButtons(inputId = "approxNorm", label = "Which norm?",
                                                   choices = c("L2", "Sup"),
                                                   selected = "L2",
                                                   inline=TRUE),
                                      sliderTextInput(inputId = "nuRange", label = "Which values of nu?",
                                                   choices = nu_val,
                                                   selected = range(nu_val)),
                                      radioButtons(inputId = "plotStyle", label = "Which style?",
                                                   choices = c("Linetype (Method), color (Order)", "Linetype (Order), color (Method)"),
                                                   selected = "Linetype (Method), color (Order)",
                                                   inline=FALSE),                                                   
                                      checkboxInput("logScaleCoverror", "Use log scale?", value = TRUE),
                                      downloadButton('downloadPlotCov', 'Download Plot'),
                                      radioButtons(inputId = "fileExtensionCov", label = "File extension",
                                                   choices = c("png", "pdf"),
                                                   selected = "png",
                                                   inline=TRUE),
                                      downloadButton('downloadDataFrameCov', 'Download Data Frame'),),
                         mainPanel(
                           width = 12 - 3,
                           plotlyOutput("raterrors")
                         )
           ),
           
  ),
  
  tabPanel("Prediction error", sidebarLayout(position = "left",
                                             sidebarPanel(width = 3,
                                      checkboxGroupInput(inputId = "methods_approx_pred", label = "Methods to be included",
                                                         choices = c("State-Space", "nnGP", "Fourier", "PCA"), 
                                                         selected = "State-Space",
                                                         inline=TRUE),                                                         
                                      checkboxGroupInput(inputId = "orderRat_pred", label = "Order of the approximation",choices = 0:6, selected = 0:6,
                                                         inline=TRUE),
                                      radioButtons(inputId = "numberLoc_pred", label = "Number of locations",
                                                   choices = c("500", "1000"),
                                                   selected = "1000",
                                                   inline=TRUE),
                                      radioButtons(inputId = "rangeParameter_pred", label = "Range Parameter",
                                                   choices = c("0.2", "0.5", "1"),
                                                   selected = "0.2",
                                                   inline=TRUE),
                                      radioButtons(inputId = "approxNorm_pred", label = "Which norm?",
                                                   choices = c("l2", "Max"),
                                                   selected = "l2",
                                                   inline=TRUE),
                                      sliderTextInput(inputId = "nuRange_pred", label = "Which values of nu?",
                                                   choices = nu_val,
                                                   selected = range(nu_val)),
                                      radioButtons(inputId = "plotStyle_pred", label = "Which style?",
                                                   choices = c("Linetype (Method), color (Order)", "Linetype (Order), color (Method)"),
                                                   selected = "Linetype (Method), color (Order)",
                                                   inline=FALSE),  
                                      checkboxInput("logScalepredError", "Use log scale?", value = TRUE),
                                      downloadButton('downloadPlotPred', 'Download Plot'),
                                      radioButtons(inputId = "fileExtensionPred", label = "File extension",
                                                   choices = c("png", "pdf"),
                                                   selected = "png",
                                                   inline=TRUE),
                                      downloadButton('downloadDataFramePred', 'Download Data Frame'),),
                                             
                                             mainPanel(
                                               width = 12 - 3,
                                               plotlyOutput("prederrors")
                                               
                                             )
  ),),
  
  tags$style(type = "text/css", "body {padding-top: 70px;}"),
  theme = shinytheme("cosmo"),
  position = "fixed-top"
  
)
