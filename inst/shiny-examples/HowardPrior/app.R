library(shiny)
library(ggplot2)
library(LearnBayes)
library(metR)

ui <- fluidPage(
  fluidRow(
    column(4, wellPanel(
      sliderInput("sigma", "Prior parameter sigma:",
                  min = 0.001, max = 2, value = 0.5),
      sliderInput("s1", "# of Successes in Sample 1:",
              min = 0, max = 40, value = 0, step = 1),
      sliderInput("f1", "# of Failures in Sample 1:",
                  min = 0, max = 40, value = 0, step = 1),
      sliderInput("s2", "# of Successes in Sample 2:",
                  min = 0, max = 40, value = 0, step = 1),
      sliderInput("f2", "# of Failures in Sample 2:",
                  min = 0, max = 40, value = 0, step = 1)
    )),
    column(8,
      tabsetPanel(type = "tabs",
          tabPanel("Model",
                   h3("Comparing Two Proportions Using an Independence (Howard's) Prior"),
                   h4("Description"),
                   p("We have two independent samples
                                  where y1 is binomial(n1, p1) and
                                  y2 is binomial(n2, p2). We assume Howard's prior
                     where (p1, p2) has the density exp(-u ^ 2 / 2) where u =
                     (logit p1 - logit p2) / sigma.  The prior reflects the belief
                     that the proportions are equal and the parameter sigma
                     reflects the sureness of this belief.  As sigma approaches 0,
                     one is placing all of the prior on the line where p1 = p2."),
                   h4("Using the App"),
                   p('One enters in the number of successes
                                  and the number of failures for each
                                  of the two samples.  Also one inputs
                                  the value of the prior standard deviation sigma.'),
                   p("The App tab shows the contour
                                  plot of the joint posterior density
                                  of (p1, p2). The posterior probability P(p1 > p2)
                                  is displayed at the bottom of the graph.
                                  ")),
          tabPanel("Posterior",
          h3(id="big-heading",
         "     Posterior Using Howard's Prior for 2 x 2 Table"),
          tags$style(HTML("#big-heading{color: red;}")),
          plotOutput("plot", height = "500px"))
      ))
))
server <- function(input, output, session) {
  output$plot <- renderPlot({
    gcontour2 <- function(logf, limits, ...){
      require(ggplot2)
      require(metR)
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
      BR <- c(seq(-6.9, -1.15, by=1.15), 0)
      ggplot(df)  +
        geom_contour_fill(aes(x=x,
                              y=y,
                              z=Z),
                          breaks=BR,
                          size=1.5) +
        scale_fill_distiller(palette="Spectral")
    }
    plo <- .0001; phi <- .9999
    s <- simcontour(howardprior, c(plo, phi, plo, phi),
                 + c(1 + input$s1, 1 + input$f1,
                     1 + input$s2, 1 + input$f2,
                     input$sigma), 50000)
    prob <- round(sum(s$x > s$y) / 50000, 3)
    mytitle <- paste(", Data: (", input$s1, ",", input$f1,
                     "), (", input$s2, ",", input$f2,
                     ")", sep = "")
    gcontour2(howardprior,
              c(plo, phi, plo, phi),
              c(1 + input$s1, 1 + input$f1,
                1 + input$s2, 1 + input$f2 ,
                input$sigma)) +
      labs(title = "Contours of Posterior of (p1, p2)",
           subtitle = paste("sigma =",
                    input$sigma, mytitle),
           caption = paste("Prob(p1 > p2) =", prob)) +
      xlab("p1") + ylab("p2") +
      xlim(0, 1) + ylim(0, 1) +
      geom_abline() +
      coord_fixed() +
      annotate(geom = "text", x = 0.9, y = 0.05,
               label = "p1 > p2", size = 5, color = "blue") +
      theme(plot.title = element_text(colour = "blue",
                                      size = 16,
                                      hjust = 0.5),
            plot.subtitle = element_text(colour = "red",
                                      size = 12,
                                      hjust = 0.5),
            plot.caption = element_text(colour = "blue",
                                        size = 16,
                                        hjust = 0.5))
  }, res = 96)
}

shinyApp(ui = ui, server = server)
