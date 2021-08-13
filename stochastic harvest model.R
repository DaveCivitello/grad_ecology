library(shiny)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Making a stochastic version of your ODE model
# Relevant to most (all?) applied questions b/c individuals are in whole units
# Sometimes new dynamics can emerge when stochasticity interacts with other mechanisms in models

# Gillespie's direct method (Exact, simple, but slow with lots of individuals) (Gillespie 1977, Exact Stochastic Simulation of Coupled Chemical Reactions, Journal of Physical Chemistry)
#Algorithm:
# 1. Figure out all of the possible events, E_i, in your model
# 2. For each event, determine the rate at which it occurs, R_i
# 3. The rate at which any event occurs is R_tot = sum(R_i)
# 4. The time to the next reaction (of yet unknown type) is exponentially-distributed with mean rate R_tot
# 5. Which event occurs is essentially a multinomial with probabilities equal to R_i/R_tot


Direct.step.harvest = function(state, parameters){
  N = state[1]; H = state[2]; t = state[3]
  with(as.list(parameters),
       {
         if(N == 0){return(c(N, H, t+1))}
         # There are only two transitions (birth and death)
         R_birth = max(b*(1 - c*N)*N, 0)
         R_death = d*N
         R_harvested = h*N
         # Total rate is sum of all transition rates
         R_total = R_birth + R_death + R_harvested
         # Determine length of timestep (exponential distributed wait time with mean rate R_total)
         delta.t = -1/R_total*log(runif(1))
         # Determine which rate occurs
         p = runif(1)
         if(p <= R_birth/R_total){ # birth happens
           N.out = N + 1
           H.out = H}else if(p <= (R_birth + R_death)/R_total){ # else-if a death happens
             N.out = N - 1
             H.out = H} else { # else a harvest happens
               N.out = N - 1
               H.out = H + 1
             }
         state.out = c(N.out, H.out, t+delta.t)
         state.out
       }
  )
}

stochastic_harvest = function(b, d, c, h, timespan){
  pop.state = matrix(c("N" = max(round((b-d)/(b*c)), 0), "H" = 0, "t" = 0), nrow=1, ncol=3)
  colnames(pop.state) = c("N", "H","t")
  i = 1
  while (pop.state[i, "t"] < timespan){
    pop.state = rbind(pop.state, Direct.step.harvest(state = pop.state[i,], parameters = c("b"=b,"d"=d,"c"=c, "h"=h)))
    i = i + 1
  }
  pop.state
}

# Define UI for miles per gallon app ----
ui <- pageWithSidebar(
  
  # App title ----
  headerPanel("Stochastic logistic growth & harvest model"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    checkboxInput(inputId= "hide1", label = "hide plot 1?", value=TRUE),
    checkboxInput(inputId= "hide2", label = "hide plot 2?", value=TRUE),
    sliderInput(inputId = "b", label= "maximum birth rate", min=0, max=1, value = 0.2),
    sliderInput(inputId = "d", label= "death rate", min=0, max=1, value = 0.1),
    sliderInput(inputId = "c", label= "competition coefficient (on births)", min=0.005, max=.02, value = 0.01),
    sliderInput(inputId = "h", label= "harvesting effort", min=0, max=2, value = 0, step = 0.01),
    sliderInput(inputId = "t", label= "Time span", min=100, max=500, value = 100, step = 100)

  ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    plotOutput("PopPlot"),
    plotOutput("HarvestPlot")
      
    )
  )


# Define server logic to make plots 
server <- function(input, output) {
  popsize <- reactive({data.frame(stochastic_harvest(b = input$b, d = input$d, c = input$c, h = input$h, timespan = input$t))})
  
  #reactive({input$N0*exp(input$r*(0:input$t))})
  output$PopPlot <- renderPlot(
    if(input$hide1 == TRUE){ ggplot() +  ylab("Population size") + xlab("time")}else{
    ggplot(data = popsize(), aes(x=t, y = N)) + geom_line() + ylab("Population size") + xlab("time") +
        geom_abline(intercept = max(round((input$b-input$d)/(input$b*input$c)), 0), slope=0, lty=2, col="blue") +
      draw_label(ifelse(min(popsize()$N) == 0 & input$h > 0, "Nice job, greedy jerk!", ""), x = 0.5*max(popsize()$t), y = 0.5*max(round((input$b-input$d)/(input$b*input$c)), 0))})
  output$HarvestPlot <- renderPlot({
    plot_grid(
    if(input$hide2 == TRUE){ ggplot() +  ylab("Cumulative harvest") + xlab("time")} else
      if(input$h == 0){ ggplot() +  ylab("Cumulative harvest")  + xlab("time") + draw_label("No harvest")} else {
    ggplot(data = popsize(), aes(x=t, y = H)) + geom_line() + ylab("Cumulative harvest") + xlab("time") },
    if(input$hide2 == TRUE){ ggplot() +  ylab("Cumulative catch per unit effort") + xlab("time")} else 
      if(input$h == 0){ ggplot() +  ylab("Cumulative catch per unit effort")  + xlab("time") + draw_label("No harvest")} else {
        ggplot(data = subset(popsize(), t > 0), aes(x=t, y = H/(t*input$h))) + geom_line() + ylab("Cumulative catch per unit effort") + xlab("time")}, nrow=1, ncol=2)
  })
}

shinyApp(ui, server)
