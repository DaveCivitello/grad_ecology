library(shiny)
library(deSolve)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Lotka Volterra competition model
LV_comp =function(t, y, parameters) { 
  N1=y[1]; N2=y[2]
  with(as.list(parameters),{
    dN1dt=r1* (1-c1*N1-a1*N2)*N1
    dN2dt=r2 *(1-c2*N2-a2*N1)*N2
    
    res = c(dN1dt,dN2dt)
    return(list(res))
  }
  ) 
} 


### To do list ###
# Color axis labels according to N1, N2
# Add 1/c_i or 1/a_i labels to figure

# Builds a polygon for shading a ggplot (Thanks StackOverflow! https://stackoverflow.com/questions/6801571/ggplot2-shade-area-above-line)
buildPoly <- function(xr, yr, slope = 1, intercept = 0, above = TRUE){
  #Assumes ggplot default of expand = c(0.05,0)
  xrTru <- xr + 0.05*diff(xr)*c(-1,1)
  yrTru <- yr + 0.05*diff(yr)*c(-1,1)
  
  #Find where the line crosses the plot edges
  yCross <- (yrTru - intercept) / slope
  xCross <- (slope * xrTru) + intercept
  
  #Build polygon by cases
  if (above & (slope >= 0)){
    rs <- data.frame(x=-Inf,y=Inf)
    if (xCross[1] < yrTru[1]){
      rs <- rbind(rs,c(-Inf,-Inf),c(yCross[1],-Inf))
    }
    else{
      rs <- rbind(rs,c(-Inf,xCross[1]))
    }
    if (xCross[2] < yrTru[2]){
      rs <- rbind(rs,c(Inf,xCross[2]),c(Inf,Inf))
    }
    else{
      rs <- rbind(rs,c(yCross[2],Inf))
    }
  }
  if (!above & (slope >= 0)){
    rs <- data.frame(x= Inf,y= -Inf)
    if (xCross[1] > yrTru[1]){
      rs <- rbind(rs,c(-Inf,-Inf),c(-Inf,xCross[1]))
    }
    else{
      rs <- rbind(rs,c(yCross[1],-Inf))
    }
    if (xCross[2] > yrTru[2]){
      rs <- rbind(rs,c(yCross[2],Inf),c(Inf,Inf))
    }
    else{
      rs <- rbind(rs,c(Inf,xCross[2]))
    }
  }
  if (above & (slope < 0)){
    rs <- data.frame(x=Inf,y=Inf)
    if (xCross[1] < yrTru[2]){
      rs <- rbind(rs,c(-Inf,Inf),c(-Inf,xCross[1]))
    }
    else{
      rs <- rbind(rs,c(yCross[2],Inf))
    }
    if (xCross[2] < yrTru[1]){
      rs <- rbind(rs,c(yCross[1],-Inf),c(Inf,-Inf))
    }
    else{
      rs <- rbind(rs,c(Inf,xCross[2]))
    }
  }
  if (!above & (slope < 0)){
    rs <- data.frame(x= -Inf,y= -Inf)
    if (xCross[1] > yrTru[2]){
      rs <- rbind(rs,c(-Inf,Inf),c(yCross[2],Inf))
    }
    else{
      rs <- rbind(rs,c(-Inf,xCross[1]))
    }
    if (xCross[2] > yrTru[1]){
      rs <- rbind(rs,c(Inf,xCross[2]),c(Inf,-Inf))
    }
    else{
      rs <- rbind(rs,c(yCross[1],-Inf))
    }
  }
  
  return(rs)
}



# Define UI for miles per gallon app ----
ui <- pageWithSidebar(
  
  # App title ----
  headerPanel("Lotka Volterra competition"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    sliderInput(inputId = "r1", label= "r1", min=0, max=1, value = 0.1),
    sliderInput(inputId = "r2", label= "r2", min=0, max=1, value = 0.1),
    
    sliderInput(inputId = "c1", label= "c1", min=0.01, max=1, value = 0.1),
    sliderInput(inputId = "c2", label= "c2", min=0.01, max=1, value = 0.1),
   
    sliderInput(inputId = "a1", label= "a1", min=0, max=1, value = 0.1),
    sliderInput(inputId = "a2", label= "a2", min=0, max=1, value = 0.1),
    
    sliderInput(inputId = "t", label= "Time steps", min=100, max=1000, value = 100, step = 100)
    
    
    
  ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    plotOutput("PopPlot"),
    plotOutput("ZNGIs")
      
    )
  )


# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  sim <- reactive({data.frame(ode(y=c(N1=runif(1, min=0, max=1/input$c1), N2=runif(1, min=0, max=1/input$c2)), times=0:input$t, 
                                  parms = c(r1 = input$r1, r2 = input$r2, 
                                            c1 = input$c1, c2 = input$c2,
                                            a1 = input$a1, a2 = input$a2), func=LV_comp, method="lsoda"))})
  
  blueshade <- reactive({data.frame(N1 = c(0, 1/input$c1), N2 = c(1/input$a1, 0))})
  
  redshade <-reactive({data.frame(N1 = c(0, 1/input$a2), N2 = c(1/input$c2, 0))})
 
  #reactive({input$N0*exp(input$r*(0:input$t))})
  output$PopPlot <- renderPlot(ggplot(data = sim()) + 
                                 geom_line(aes(x=time, y=N1), colour="blue") + 
                                 geom_line(aes(x=time, y=N2), colour="red"))
  output$ZNGIs <- renderPlot(ggplot() + xlim(0, 1.2/min(input$c1, input$a2)) + ylim(0, 1.2/min(input$c2, input$a1)) +
                               geom_abline(intercept = 1/input$a1, slope = -input$c1/input$a1, colour="blue") +
                               geom_ribbon(data=blueshade(), aes(x=N1, ymax=N2, ymin=0), fill="blue", alpha=0.2) +
                               geom_abline(intercept = 1/input$c2, slope = -input$a2/input$c2, colour="red") +
                               geom_ribbon(data=redshade(), aes(x=N1, ymax=N2, ymin=0), fill="red", alpha=0.2) +
                               xlab("N1") + ylab("N2") +
                               geom_path(data=sim(), aes(x=N1, y=N2), arrow=arrow())
                             )
  
}


shinyApp(ui, server)
