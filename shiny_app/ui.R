library(shiny)
library(plotly)
library(shinythemes)
library(shinyWidgets)

nu_val <- seq(0.01, 2.49, 0.01)

ui <- navbarPage(
  "Markov Rational Approximation",

  tabPanel("Overview",
           fluidPage(
             titlePanel("Welcome to the Markov Rational Approximation App"),
             mainPanel(
               h4("Introduction"),
               p("This application provides tools to analyze the covariance error, prediction error, 
                  and posterior probability error for various approximation methods."),
               p("To ensure fair comparisons between methods, we performed detailed calibration procedures. 
                 Please check the 'Details on calibrations' tab for further information.")
             )
           )
  ),

  tabPanel("Covariance error",  
           sidebarLayout(position = "left",
                         sidebarPanel(width = 3,
                                      checkboxGroupInput(inputId = "methods_approx", label = "Methods to be included",
                                                         choices = c("State-Space", "nnGP", "Fourier", "PCA", "Taper", "FEM"), 
                                                         selected = "State-Space", inline = TRUE),                                                         
                                      checkboxGroupInput(inputId = "orderRat", label = "Order of the approximation",
                                                         choices = 0:6, selected = 0:6, inline = TRUE),
                                      radioButtons(inputId = "numberLoc", label = "Number of locations",
                                                   choices = c("5000", "10000"), selected = "5000", inline = TRUE),                      
                                      radioButtons(inputId = "numberObs", label = "Number of observations",
                                                   choices = c("5000", "10000"), selected = "5000", inline = TRUE),                                                                           
                                      radioButtons(inputId = "rangeParameter", label = "Range Parameter",
                                                   choices = c("0.5", "1", "2"), selected = "2", inline = TRUE),
                                      radioButtons(inputId = "approxNorm", label = "Which norm?",
                                                   choices = c("L2", "Sup"), selected = "L2", inline = TRUE),
                                      sliderTextInput(inputId = "nuRange", label = "Which values of nu?",
                                                       choices = nu_val, selected = range(nu_val)),
                                      radioButtons(inputId = "plotStyle", label = "Which style?",
                                                   choices = c("Linetype (Method), color (Order)", 
                                                               "Linetype (Order), color (Method)"),
                                                   selected = "Linetype (Method), color (Order)", inline = FALSE),                                                   
                                      checkboxInput("logScaleCoverror", "Use log scale?", value = TRUE),
                                      checkboxInput("trueCalibrationCov", "True calibration for nnGP and Taper", value = TRUE),
                                      downloadButton('downloadPlotCov', 'Download Plot'),
                                      radioButtons(inputId = "fileExtensionCov", label = "File extension",
                                                   choices = c("png", "pdf"), selected = "png", inline = TRUE),
                                      downloadButton('downloadDataFrameCov', 'Download Data Frame')
                         ),
                         mainPanel(width = 9, plotlyOutput("raterrors"))
           )
  ),

  tabPanel("Prediction error", 
           sidebarLayout(position = "left",
                         sidebarPanel(width = 3,
                                      checkboxGroupInput(inputId = "methods_approx_pred", 
                                                         label = "Methods to be included",
                                                         choices = c("State-Space", "nnGP", "Fourier", "PCA", "Taper", "FEM"), 
                                                         selected = "State-Space", inline = TRUE),                                                         
                                      checkboxGroupInput(inputId = "orderRat_pred", 
                                                         label = "Order of the approximation",
                                                         choices = 0:6, selected = 0:6, inline = TRUE),
                                      radioButtons(inputId = "numberLoc_pred", label = "Number of locations",
                                                   choices = c("5000", "10000"), selected = "5000", inline = TRUE),
                                      radioButtons(inputId = "numberObs_pred", label = "Number of observations",
                                                   choices = c("5000", "10000"), selected = "5000", inline = TRUE),                                                   
                                      radioButtons(inputId = "rangeParameter_pred", label = "Range Parameter",
                                                   choices = c("0.5", "1", "2"), selected = "2", inline = TRUE),
                                      sliderTextInput(inputId = "nuRange_pred", label = "Which values of nu?",
                                                       choices = nu_val, selected = range(nu_val)),
                                      radioButtons(inputId = "plotStyle_pred", label = "Which style?",
                                                   choices = c("Linetype (Method), color (Order)", 
                                                               "Linetype (Order), color (Method)"),
                                                   selected = "Linetype (Method), color (Order)", inline = FALSE),  
                                      checkboxInput("logScalepredError", "Use log scale?", value = TRUE),
                                      checkboxInput("trueCalibrationPred", "True calibration for nnGP and Taper", value = TRUE),
                                      downloadButton('downloadPlotPred', 'Download Plot'),
                                      radioButtons(inputId = "fileExtensionPred", label = "File extension",
                                                   choices = c("png", "pdf"), selected = "png", inline = TRUE),
                                      downloadButton('downloadDataFramePred', 'Download Data Frame')
                         ),
                         mainPanel(width = 9, plotlyOutput("prederrors"))
           )
  ),

  tabPanel("Posterior probability error", 
           sidebarLayout(position = "left",
                         sidebarPanel(width = 3,
                                      checkboxGroupInput(inputId = "orderRat_prob", label = "Order of the approximation",
                                                         choices = 1:6, selected = 3:6, inline = TRUE),
                                      radioButtons(inputId = "rangeParameter_prob", label = "Range Parameter",
                                                   choices = c("0.5", "1", "2"), selected = "0.5", inline = TRUE),
                                      radioButtons(inputId = "plotStyle_prob", label = "Which style?",
                                                   choices = c("Linetype (Method), color (Order)", 
                                                               "Linetype (Order), color (Method)"),
                                                   selected = "Linetype (Method), color (Order)", inline = FALSE),
                                      checkboxInput("logScaleprobError", "Use log scale?", value = FALSE),
                                      downloadButton('downloadPlotprob', 'Download Plot'),
                                      radioButtons(inputId = "fileExtensionprob", label = "File extension",
                                                   choices = c("png", "pdf"), selected = "png", inline = TRUE),
                                      downloadButton('downloadDataFrameprob', 'Download Data Frame')
                         ),
                         mainPanel(width = 9, plotlyOutput("proberrors"))
           )
  ),

tabPanel("Details on calibrations",
           fluidPage(
             tags$head(
               tags$script(src = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js")
             ),
             titlePanel("Details on Calibration Procedures"),
             mainPanel(
               width = 12,
               h4("Calibration Overview"),
               p("Calibration ensures that different methods are compared fairly by equalizing their computational costs. 
                  This involves balancing the time required to assemble matrices (construction cost) and compute the posterior mean (prediction cost)."),

               h5("nnGP and Taper Methods: True Calibration"),
               p("When the smoothness parameter \\(\\nu\\) is within the interval \\(0 < \\nu < 0.5\\), the nnGP and Taper methods 
                  are calibrated by fixing \\(m = 1\\). This is because the rational approximation with \\(m = 1\\) achieves 
                  lower computational costs compared to nnGP and Taper. If the 'True calibration' checkbox is selected, the app enforces this setting, 
                  fixing \\(m = 1\\) for all comparisons involving nnGP and Taper within this interval."),

               p("If the 'True calibration' checkbox is unchecked, users can explore how performance changes as \\(m\\) varies from 1 to 6. 
                 However, it is important to note that comparisons may become less fair to the rational approximation when \\(m > 1\\), 
                 since the computational costs for nnGP and Taper are already higher with \\(m = 1\\)."),

               h5("Calibration by Smoothness Parameter \\(\\nu\\)"),
               tags$ul(
                 tags$li("When \\(\\nu\\) is within the interval \\(0 < \\nu < 0.5\\), calibration is not feasible for nnGP and Taper, 
                          and \\(m = 1\\) is used by default for these methods. Other methods are calibrated normally."),
                 tags$li("For the interval \\(0.5 < \\nu < 1.5\\), separate calibration is performed for all methods."),
                 tags$li("In the interval \\(1.5 < \\nu < 2.5\\), true calibration was found to be unstable for the FEM method, 
                          so calibrations were adjusted to ensure stable predictions, while other methods were calibrated normally.")
               ),

               h5("Taper Method Calibration"),
               p("For the Taper method, the taper distance is selected such that each observation has, on average, \\(m\\) neighbors within the distance. 
                 The value of \\(m\\) is adjusted to ensure that the total computational cost matches that of the rational approximation."),

               h5("PCA and Fourier Methods"),
               p("The PCA and Fourier methods share identical prediction costs, as both depend on the number of basis functions. 
                 The value of \\(m\\) for Fourier is therefore set to match that of PCA. 
                 PCA calibration assumes that the eigenvectors of the covariance matrix are known a priori, 
                 providing a theoretical lower bound on the computational cost for any low-rank method."),

               h5("State-Space Method"),
               p("The state-space method provides an alternative Markov representation, leveraging the same computational strategies as the rational approximation. 
                 For this method, the value of \\(m\\) is determined as \\(m - \\lfloor \\alpha \\rfloor\\), where \\(\\alpha\\) is a parameter that controls the order of the approximation."),

               h5("FEM Method"),
               p("The FEM method extends the original domain \\([0, 50]\\) to \\([-4\\rho, 50 + 4\\rho]\\), where \\(\\rho\\) represents the practical correlation range. 
                 This extension minimizes boundary effects and ensures that the Matérn covariance is accurately approximated at the target locations.")
             )
           )
  )
,

  tags$style(type = "text/css", "body {padding-top: 70px;}"),
  theme = shinytheme("cosmo"),
  position = "fixed-top"
)