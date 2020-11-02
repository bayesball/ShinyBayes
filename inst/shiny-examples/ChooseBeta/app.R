
ui <- fluidPage(
  h1(id="big-heading", "Constructing a Beta(a, b) Prior From Two Quantiles"),
  tags$style(HTML("#big-heading{color: blue;}")),
  sidebarLayout(
    sidebarPanel(
      sliderInput("Prior",
                  "Choose Median and 90th Percentile of Prior:",
                  min = 0,
                  max = 1,
                  value = c(0.5, 0.8)),
      hr(),
      sliderInput("Prob",
                  "Choose Probability Content for Middle Interval:",
                  min = 0.50,
                  max = 0.95,
                  value = 0.9),
      hr(),
      sliderInput("n",
                  "Choose Future Sample Size N:",
                  min = 10,
                  max = 50,
                  value = 20)
    ),

    mainPanel(
      plotOutput("distPlot",
                 height = "550px")
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
      xlab("P") + ylab("Density")
    y <- 0:input$n
    p <- a / (a + b)
    py <- dbeta(p, a, b) *
      dbinom(y, size = input$n, prob = p) /
      dbeta(p, a + y, b + input$n - y)
    pred_dist <- data.frame(Y = y,
                            Probability = py)
    ps <- discint(pred_dist, input$Prob)

    mytitle <- paste("Predictive Distribution for Sample of Size",
                     input$n)
    t1 <- paste("Prob(", min(ps$set), " ≤ Y ≤ ",
                max(ps$set), ") = ",
                round(ps$prob, 3), sep = "")
    p2 <- ggplot(data = pred_dist,
                 aes(Y, Probability)) +
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
      geom_segment(aes(xend = Y, yend = 0),
                   size = 2, lineend = "butt",
                   color = "blue") +
      geom_segment(data = dplyr::filter(pred_dist,
                        Y %in% ps$set),
                   aes(xend = Y, yend = 0),
                   size = 2, lineend = "butt",
                   color = "red") +
      theme(text = element_text(size = 18)) +
      xlab("Number of Successes") + ylab("Probability")

    grid.arrange(p1, p2, ncol = 1)

  })

}


shinyApp(ui = ui, server = server)
