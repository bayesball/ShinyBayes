library(shiny)
library(ggplot2)
library(ProbBayes)
require(metR)

# Define UI ----
ui <- fluidPage(
#  titlePanel("Visualizing Posterior of Two Proportions"),
  h2(id="big-heading", "Visualizing Posterior of Two Proportions"),
  tags$style(HTML("#big-heading{color: red;}")),
  fluidRow(
    column(4, wellPanel(
#      h4(id="model-heading", "Model:"),
 #     tags$style(HTML("#model-heading{color: red;}")),
      h5("s1 is Binom(s1 + f1, p1)"),
      h5("s2 is Binom(s2 + f2, p2)"),
      h5("theta1 = logit(p1) - logit(p2)"),
      h5("theta2 = logit(p1) + logit(p2)"),
      h5("(theta1, theta2) is Uniform"),
      h4(id="data-heading", "Enter Data:"),
      tags$style(HTML("#data-heading{color: red;}")),
      sliderInput("s1", "s1 = # Successes in 1st Sample:", 
                  min = 1, max = 20, 
                  value = 5,
                  step = 1),
      sliderInput("f1", "f1 = # Failures in 1st Sample:", 
                  min = 1, max = 20, 
                  value = 10,
                  step = 1),
      sliderInput("s2", "s2 = # Successes in 2nd Sample:", 
                  min = 1, max = 20, 
                  value = 3,
                  step = 1),
      sliderInput("f2", "f2 = # Failures in 2nd Sample:", 
                  min = 1, max = 20, 
                  value = 8,
                  step = 1),
      sliderInput("iter", "# Simulations:", 
                  min = 100, max = 5000, 
                  value = 1000,
                  step = 20)
    )),
    column(8,
           plotOutput("plot", height = "380px"),
           h4(id="data-heading", "Data:"),
           tags$style(HTML("#data-heading{color: red;}")),
    #       h4("Data:"),
           verbatimTextOutput("stats"),
           h4(id="stats-heading", "Simulation Posterior Summaries:"),
           tags$style(HTML("#stats-heading{color: red;}")),
    
    #       h4("Simulation Posterior Summaries:"),
           verbatimTextOutput("stats2")
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  gcontour2 <- function(logf, limits, ...){
    LOGF <- function(theta, ...) {
      if (is.matrix(theta) == TRUE) {
        val <- matrix(0, c(dim(theta)[1], 1))
        for (j in 1:dim(theta)[1]){ 
          val[j] <- logf(theta[j, ], ...)
        }
      }
      else val <- logf(theta, ...)
      return(val)
    }
    ng <- 50
    df <- expand.grid(
      x = seq(limits[1], limits[2], length = ng),
      y = seq(limits[3], limits[4], length = ng)
    )
    df$Z <- unlist(LOGF(as.matrix(df), ...))
    df$Z <- df$Z - max(df$Z)
#    BR <- seq(-6.9, -1.15, 1.15)
    BR <- c(-6.9, -4.6, -2.3)
    ggplot(df)  +
      geom_contour_fill(aes(x=x,
                            y=y,
                            z=Z),
                        breaks=BR,
                        size=1.5) +
      scale_fill_distiller(palette="Spectral") +
      theme(text=element_text(size=18)) 
  }
  data <- reactive({
    logpost <- function(theta){
      d <- c(input$s1, input$f1, 
             input$s2, input$f2)
      logctablepost(theta, d)
    }
    fit <- laplace(logpost, c(0, 0))
    mo <- fit$mode
    sds <- sqrt(diag(fit$var))
    mdata <- matrix(c(input$s1, input$f1, 
                      input$s2, input$f2), 
                    2, 2, byrow = TRUE)
    dimnames(mdata)[[2]] <- c("Successes", "Failures")
    dimnames(mdata)[[1]] <- c("Sample 1", "Sample 2")
    summ <- data.frame(Parameter = 
                         c("theta1", "theta2"),
                       Mean = mo,
                       Standard_Deviation = sds)
    limits <- c(mo[1] - 5 * sds[1], 
                mo[1] + 5 * sds[1],
                mo[2] - 5 * sds[2],
                mo[2] + 5 * sds[2])
    
    sim_post <- simcontour(logctablepost, limits,
        c(input$s1, input$f1, input$s2, input$f2),
        input$iter)
    new_mo <- as.character(round(c(mean(sim_post$x), 
                mean(sim_post$y)), 3))
    new_sd <- as.character(round(c(sd(sim_post$x), 
                sd(sim_post$y)), 3))
    new_cor <- as.character(round(cor(sim_post$x, 
                      sim_post$y), 3))
    prob <- as.character(round(mean(sim_post$x > 
                        0), 3))
    newsumm <- data.frame(Parameter = 
                         c("theta1", "theta2",
                           "(theta1, theta2)",
                           "P(theta1 > 0)"),
              Mean = c(new_mo, "", prob),
            Stan_Dev = c(new_sd, "", ""),
            Correlation = c("", "", new_cor, ""))
    
    sim_post <- data.frame(x = sim_post$x,
                           y = sim_post$y)
    gplot <- gcontour2(logpost, limits) +
      geom_point(data = sim_post, aes(x, y),
                 alpha = 0.5, color = "black")
    list(data = mdata,
         summary = newsumm,
         plot = gplot)
  })
  output$plot <- renderPlot({
    out <- data()
    t0 <- paste("Data: (", input$s1, ", ",
                input$f1, ", ", input$s2, ", ",
                input$f2, "), ", sep="")
    data()$plot +
      increasefont() +
      ggtitle("Posterior of (theta1, theta2)") +
      theme(plot.title = element_text(colour = "blue", 
          size = 24, 
         hjust = 0.5, vjust = 0.8, angle = 0)) +
      xlab("theta1 = Difference in Logits") +
      ylab("theta2 = Sum of Logits") +
      theme(
        axis.title = element_text(color = "blue", 
                    size = 18))
  })
  output$stats <- renderPrint({
    data()$data
  })
  output$stats2 <- renderPrint({
    data()$summary
  })
}
  
# Run the app ----
shinyApp(ui = ui, server = server)