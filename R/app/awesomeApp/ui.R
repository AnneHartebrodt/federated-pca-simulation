#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Distributed Privacy Aware PCA"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
       sliderTextInput("owners",
                   "Number of Data Owners:",
                    choices = as.character(sliderSplit)),
       sliderTextInput("epsilon",
                   "Epsilon",
                   choices =as.character( sliderEpsilon)
       ),
       sliderTextInput("delta",
                       "Delta",
                       choices = as.character(sliderDelta)
       )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
       plotOutput("distPlot")
    )
  )
))
