library(shiny)
library(ggplot2)
library(gridExtra)

# Define UI ----
ui <- fluidPage(
  # Application title
  titlePanel("Bayes Triplot"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("a",
                  "Beta Shape Parameter a:",
                  min = 1,
                  max = 30,
                  value = 4,
                  step = 0.2),
      sliderInput("b",
                  "Beta Shape Parameter b:",
                  min = 1,
                  max = 30,
                  value = 4,
                  step = 0.2),
      sliderInput("s",
                  "Number of Successes:",
                  min = 0,
                  max = 50,
                  value = 0),
      sliderInput("f",
                  "Number of Failures:",
                  min = 0,
                  max = 50,
                  value = 0)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot", height="450px")
    )
  ))

# Define server logic ----
server <- function(input, output) {
  gtriplot <- function(a, b, s, f){

    p <- seq(0, 1, length.out = 200)
    prior <- data.frame(p = p,
                        Density = dbeta(p, a, b),
                        Type = "Prior")
    like <- data.frame(p = p,
                       Density = dbeta(p, s + 1, f + 1),
                       Type = "Likelihood")
    post <- data.frame(p = p,
                       Density = dbeta(p, a + s, b + f),
                       Type = "Posterior")
    dall <- rbind(prior, like, post)

    text <- paste("Beta(", a, ", ",
                  b, ") Prior, (s, f) = (",
                  s, ", ", f,")", sep = "")

    ggplot(dall, aes(p, Density, color = Type)) +
      geom_line(size = 1.5) +
      ggtitle(paste("Triplot: ", text)) +
      theme(plot.title = element_text(colour = "red",
                                      size = 24,  hjust = 0.5,
                                      vjust = 0.8, angle = 0)) +
      theme(text = element_text(size = 18)) +
      xlab("P") + ylab("Density")
  }

    output$distPlot <- renderPlot({
      gtriplot(input$a, input$b,
               input$s, input$f)
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
