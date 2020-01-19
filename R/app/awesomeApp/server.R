#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#



owners=1
epsilon = 0.001
delta = 0.01
# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
    print(input$owners)
    print(input$epsilon)
    print(input$delta)
    print(head(data[split==as.numeric(input$owners) & epsilon ==as.numeric(input$epsilon) & delta == as.numeric(input$delta)]))
    g4 <- ggplot(data[split==as.numeric(input$owners) & epsilon ==as.numeric(input$epsilon) & delta == as.numeric(input$delta)], aes(PC1, PC2, color = as.factor(index)))+geom_point()+
      ggtitle('PCA')+
      ylab('PC2')+xlab('PC1')+theme(legend.title = element_blank(), plot.title = element_text(size = 28))+
      #geom_mark_circle(aes(color=as.factor(index)), expand = unit( 2, "mm"))+
      theme(legend.position = 'none')
    g4
     })
  
})
