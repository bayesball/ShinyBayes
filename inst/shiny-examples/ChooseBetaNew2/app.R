library(shiny)
library(ggplot2)
library(gridExtra)

ui <- fluidPage(
  theme = bslib::bs_theme(version = 4,
                          bootswatch = "spacelab"),
  h2(id="big-heading", "Using a Beta(a, b) Prior: Batting Average Edition"),
  tags$style(HTML("#big-heading{color: blue;}")),
  sidebarLayout(
    sidebarPanel(
      h4(id="heading2", "Beta Prior:"),
      tags$style(HTML("#heading2{color: red;}")),
      sliderInput("PriorMean",
                  "Choose Mean on True BA:",
                  min = .15,
                  max = .45,
                  value = 0.25),
      sliderInput("PriorPrecision",
                  "Choose Precision on True BA:",
                  min = 1,
                  max = 2000,
                  value = 200),
      hr(),
      sliderInput("Prob",
                  "Middle Probability Content:",
                  min = 0.50,
                  max = 0.95,
                  value = 0.9),
      hr(),
      sliderInput("n",
                  "Number of Future At-Bats:",
                  min = 50,
                  max = 500,
                  value = 100)
    ),

    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Description",
                           br(),
                           h4('Description'),
                           p('This app is a tool for constructing a beta
                             prior for a true batting probability.'),
                           h4('Using the App'),
                           p("One uses the slider to input one's opinion
                             about the location of the mean eta and
                             precision K of the true BA."),
                           p("The beta density will have shape parameters
                             equal to a = K eta and b = K (1 - eta)."),
                           p('The app displays the matching beta curve and
                             the values of the beta shape parameters.
                             The interval that covers the middle probability
                             content of the proportion is displayed where
                             one can choose the value of the probability
                             content.'),
                           p('The app also displays the corresponding
                             predictive distribution of the
                             batting average for a future number of at-bats (AB).
                             The value of AB can be selected by
                             the slider.')),
                  tabPanel("Choose Prior",
                           plotOutput("distPlot",
                             height = "550px"))
      )
    )
  )
)


server <- function(input, output) {
  discint <- function (dist, prob){
    x = dist[, 1]
    p = dist[, 2]
    n = length(x)
    sp = sort(p, index.return = TRUE)
    ps = sp$x
    i = sp$ix[seq(n, 1, -1)]
    ps = p[i]
    xs = x[i]
    cp = cumsum(ps)
    ii = 1:n
    j = ii[cp >= prob]
    j = j[1]
    eprob = cp[j]
    set = sort(xs[1:j])
    v = list(prob = eprob, set = set)
    return(v)
  }

  output$distPlot <- renderPlot({

    a <- input$PriorPrecision * input$PriorMean
    b <- input$PriorPrecision * (1 - input$PriorMean)
    p_lo <- (1 - input$Prob) / 2
    p_hi <- 1 - (1 - input$Prob) / 2
    q <- round(qbeta(c(p_lo, p_hi), a, b), 3)

    x0 <- seq(q[1], q[2], length.out = 100)
    y0 <- dbeta(x0, a, b)
    xx <- c(q[1], x0, q[2], q[1])
    yy <- c(0, y0, 0, 0)

    t1 <- paste("Prob( ", q[1], " < p < ", q[2], " ) = ",
                input$Prob,
                sep = "")

    mytitle <- paste("Beta(a, b) Prior with a = ",
                     a, ", b = ", b, sep="")
    p1 <- ggplot() +
      labs(title = mytitle,
           subtitle = t1) +
      theme(plot.title = element_text(colour = "red",
                                      size = 22,
                                      hjust = 0.5,
                                      vjust = 0.8,
                                      angle = 0),
            plot.subtitle = element_text(color = "blue",
                                      size = 18,
                                      hjust = 0.5,
                                      vjust = 0.8)) +
      geom_segment(mapping = aes(x = q,
                       y = 0,
                       xend = q,
                       yend = dbeta(q, a, b)),
                   size = 2, color = "blue") +
      geom_polygon(data = data.frame(xx, yy),
                   mapping = aes(xx, yy),
                   fill = "yellow") +
      stat_function(data = data.frame(x = c(0, 1)),
                    mapping = aes(x),
                    fun = dbeta,
                    args = list(shape1 = a, shape2 = b),
                    color = "red",
                    size = 1.5) +
      theme(text = element_text(size = 18)) +
      xlab("True Batting Average") + ylab("Density") +
      xlim(0.1, 0.5)

    y <- 0:input$n
    phat <- y / input$n  #new
    p <- a / (a + b)
    py <- dbeta(p, a, b) *
      dbinom(y, size = input$n, prob = p) /
      dbeta(p, a + y, b + input$n - y)
    pred_dist <- data.frame(pHat = phat,
                            Probability = py)  #new
    ps <- discint(pred_dist, input$Prob)

    mytitle <- paste("Predictive Distribution for BA for",
                     input$n, "Future At-Bats")
    t1 <- paste("Prob(", round(min(ps$set), 3), " ≤ Y ≤ ",
                round(max(ps$set), 3), ") = ",
                round(ps$prob, 3), sep = "")
    p2 <- ggplot() +
      labs(title = mytitle,
           subtitle = t1) +
      theme(plot.title = element_text(colour = "red",
                                      size = 22,
                                      hjust = 0.5,
                                      vjust = 0.8,
                                      angle = 0),
            plot.subtitle = element_text(color = "blue",
                                         size = 18,
                                         hjust = 0.5,
                                         vjust = 0.8)) +
      geom_segment(data = pred_dist,
                   mapping = aes(x = pHat,
                       y = Probability,
                       xend = pHat,
                       yend = 0),
                   size = 2, lineend = "butt",
                   color = "blue") +
      geom_segment(data = dplyr::filter(pred_dist,
                        phat %in% ps$set),
                   mapping = aes(x = pHat,
                                 y = Probability,
                                 xend = pHat,
                                 yend = 0),
                   size = 2, lineend = "butt",
                   color = "red") +
      theme(text = element_text(size = 18)) +
      xlab("Batting Average") +
      ylab("Probability") +
      xlim(0.1, 0.5)

    grid.arrange(p1, p2, ncol = 1)

  })
}


shinyApp(ui = ui, server = server)
