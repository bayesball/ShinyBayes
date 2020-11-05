library(shiny)
library(ggplot2)
library(ProbBayes)
require(metR)
require(gridExtra)

# Define UI ----
ui <- fluidPage(
#  titlePanel("Visualizing Posterior of Two Proportions"),
  h2(id="big-heading", "Visualizing Posterior of Two Proportions"),
  tags$style(HTML("#big-heading{color: red;}")),
  fluidRow(
    column(4, wellPanel(
      h4(id="data-heading", "Enter Successes and Failures:"),
      tags$style(HTML("#data-heading{color: red;}")),
      h5(id="data-heading1", "Sample 1:"),
      fluidRow(
        column(2,  HTML('<h5><b>s1:</b></h5>')),
        column(4, numericInput("s1", "", min = 0, max = 100, value = 5)),
        column(2, HTML('<h5><b>f1:</b></h5>')),
        column(4, numericInput("f1", "", min = 0, max = 100, value = 5))
      ),
      h5(id="data-heading1", "Sample 2:"),
      fluidRow(
        column(2,  HTML('<h5><b>s2:</b></h5>')),
        column(4, numericInput("s2", "", min = 0, max = 100, value = 5)),
        column(2, HTML('<h5><b>f2:</b></h5>')),
        column(4, numericInput("f2", "", min = 0, max = 100, value = 5))
      ),
      sliderInput("iter", "# Simulations:",
                  min = 100, max = 5000,
                  value = 1000,
                  step = 20)
    )),
    column(8,
           tabsetPanel(type = "tabs",
                       tabPanel("Model",
                                h4("Description"),
                                p('We have two independent samples
                                  where y1 is binomial(n1, p1) and
                                  y2 is binomial(n2, p2).  We
                                  reparametrize to the difference in
                                  logits theta1 = logit(p1) - logit(p2)
                                  and the sum of logits theta2 =
                                  logit(p1) + logit(p2).  We assume that
                                  the prior on (theta1, theta2) is
                                  uniform.'),
                                h4("Using the App"),
                                p('One enters in the number of successes
                                  and the number of failures for each
                                  of the two samples.  Also one inputs
                                  the number of posterior simulations.'),
                                p("The Data tab shows the observed data,
                                  the Contour tab displays a contour
                                  plot of the joint posterior density
                                  of (theta1, theta2).  The Marginals tab
                                  displays density estimates of the
                                  marginal posterior densities of theta1
                                  and theta2.  The Summaries tab gives
                                  inferential summaries of the posterior.
                                  ")),
                       tabPanel("Data",
                                tableOutput("stats")),
                       tabPanel("Contour",
                                plotOutput("plot",
                                           height = "400px")),
                       tabPanel("Marginals",
                                plotOutput("mplot",
                                           height = "400px")),
                       tabPanel("Summaries",
                                tableOutput("stats2"))
           )
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
         plot = gplot,
         sim_post = sim_post)
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
  output$mplot <- renderPlot({
    out <- data()
    p1 <- ggplot(out$sim_post, aes(x)) +
      geom_density(size = 2, color = "red") +
      xlab("theta1") +
      ggtitle("Difference in Logits") +
      increasefont() + centertitle()
    p2 <- ggplot(out$sim_post, aes(y)) +
      geom_density(size = 2, color = "red") +
      xlab("theta2") +
      ggtitle("Sum of Logits") +
      increasefont() + centertitle()
    grid.arrange(p1, p2, ncol = 1)
  })
  output$stats <- renderTable({
    as.data.frame(data()$data) -> df
    df$Sample <- c("Sample 1", "Sample 2")
    df$'Sample Size' <- df[, 1] + df[, 2]
    df[, c(3, 1, 2, 4)]
  })
  output$stats2 <- renderTable({
    data()$summary
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
