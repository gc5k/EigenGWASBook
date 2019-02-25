#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("EigenGWAS Power Calculator"),
   hr(),
   # Sidebar with a slider input for number of bins
   sidebarLayout(
      sidebarPanel(
#        numericInput("n", "Sample size",
#                     value = 500, min=100),
        selectInput("pop", "Population type",
                     choices = c("Outbred", "Inbred"), selected="Outbred"),
        numericInput("m", "Marker number",
                     value = 10000, min=1),
        sliderInput("w1",
                     "Subpop 1 proportion",
                     min = 0.05, max = 0.95, value = 0.5),
        sliderInput("p1", "Frequency at population 1",
                    min = 0.01, max=0.99, value=0.3),
        sliderInput("p2", "Frequency at population 2",
                    min = 0.01, max=0.99, value=0.5),
        selectInput("alpha", "Alpha",
                    choices = c("0.001", "0.005", "0.01", "0.05", "0.1"), 
                    selected="0.05"),
        selectInput("beta", "Beta",
                    choices = c("0.5", "0.6", "0.7", "0.8", "0.9"),
                    selected="0.8"),
        submitButton(text="Update power", icon("refresh"))
      ),

      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

   output$distPlot <- renderPlot({
     popType=input$pop
     m=input$m
     alpha=as.numeric(input$alpha)
     pcut=alpha/m
     chiT=qchisq(pcut, 1, lower.tail = F)

     n=c(100,  200,  500,   1000,  1500,  2000,
         5000, 7500, 10000, 15000, 20000, 50000)
     PW=matrix(0, 2, length(n))

     w1=input$w1
     w2=1-w1

     p1=input$p1
     h1=2*p1*(1-p1)
     p2=input$p2
     h2=2*p2*(1-p2)
     
     p=w1*p1+w2*p2
     H=w1*h1+w2*h2
     
     for(i in 1:length(n)) {
       if(popType == "Outbred") {
         ncpA=4*n[i]*w1*w2*(p1-p2)^2/(2*p*(1-p))
         ncpD=n[i]*w1*w2*(h1-h2)^2/H
       } else {
         ncpA=2*n[i]*w1*w2*(p1-p2)^2/(2*p*(1-p))
         ncpD=0
       }

       PW[1,i]=pchisq(chiT, 1, ncp=ncpA, lower.tail = F)
       PW[2,i]=pchisq(chiT, 1, ncp=ncpD, lower.tail = F)
     }
     colnames(PW)=n
     par(las=2)
     barplot(PW, beside = T, border = F, ylab="Statistical power", xlab="Sample size")
     abline(h=as.numeric(input$beta), lty=2, col="grey")
     legend("topleft", legend=c("Add", "Dom"), pch=15, col=c("black", "grey"), bty='n')
   })
}

# Run the application
shinyApp(ui = ui, server = server)

