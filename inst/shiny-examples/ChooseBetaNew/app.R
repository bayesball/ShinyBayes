library(shiny)
library(ggplot2)
library(gridExtra)

ui <- fluidPage(
  h2(id="big-heading", "Constructing a Beta(a, b) Prior: Batting Average Edition"),
  tags$style(HTML("#big-heading{color: blue;}")),
  sidebarLayout(
    sidebarPanel(
      sliderInput("Prior",
                  "Choose Median and 90th Percentile of Prior on True BA:",
                  min = .2,
                  max = .4,
                  value = c(0.25, 0.30)),
      hr(),
      sliderInput("Prob",
                  "Probability Content for Middle Interval:",
                  min = 0.50,
                  max = 0.95,
                  value = 0.9),
      hr(),
      sliderInput("n",
                  "Choose Future At-Bats:",
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
                             about the location of the median and the 90th
                             percentile of the true BA."),
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
  beta.select <- function(quantile1, quantile2){
    betaprior1 <- function(K, x, p){
      m.lo <- 0
      m.hi <- 1
      flag <- 0
      while(flag==0){
        m0 <- (m.lo + m.hi) / 2
        p0 <- pbeta(x, K * m0, K * (1 - m0))
        if(p0 < p) m.hi <- m0 else m.lo <- m0
        if(abs(p0 - p) < .0001) flag <- 1
      }
      return(m0)
    }

    p1 <- quantile1$p
    x1 <- quantile1$x
    p2 <- quantile2$p
    x2 <- quantile2$x

    logK <- seq(-3, 8, length=100)
    K <- exp(logK)
    m <- sapply(K, betaprior1, x1, p1)

    prob2 <- pbeta(x2, K * m, K * (1 - m))
    ind <- ((prob2 > 0) & (prob2 < 1))
    app <- approx(prob2[ind], logK[ind], p2)
    K0 <- exp(app$y)
    m0 <- betaprior1(K0, x1, p1)

    return(round(K0 * c(m0, (1 - m0)),   2))
  }

  output$distPlot <- renderPlot({

    ab <- beta.select(list(x=input$Prior[1],
                           p=0.5),
                      list(x=input$Prior[2],
                           p=0.9))
    a <- ab[1]
    b <- ab[2]
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
                                      size = 18,
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
                     input$n, "Future AB")
    t1 <- paste("Prob(", round(min(ps$set), 3), " ≤ Y ≤ ",
                round(max(ps$set), 3), ") = ",
                round(ps$prob, 3), sep = "")
    p2 <- ggplot() +
      labs(title = mytitle,
           subtitle = t1) +
      theme(plot.title = element_text(colour = "blue",
                                      size = 18,
                                      hjust = 0.5,
                                      vjust = 0.8,
                                      angle = 0),
            plot.subtitle = element_text(color = "red",
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
