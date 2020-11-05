library(shiny)
library(ggplot2)
library(ProbBayes)
require(gridExtra)
require(dplyr)

# Define UI ----
ui <- fluidPage(
#  titlePanel("Visualizing Posterior of Two Proportions"),
  h2(id="big-heading", "Learning About Two Proportions Using a Discrete Prior"),
  tags$style(HTML("#big-heading{color: red;}")),
  fluidRow(
    column(4, wellPanel(
    h4(id="prior-heading", "Select Prior:"),
    tags$style(HTML("#prior-heading{color: red;}")),
    sliderInput("P",
                  "Choose Limits of Each P:",
                  min = 0,
                  max = 1,
                  value = c(0.05, 0.95)),
    sliderInput("N", "# of Values:",
                min = 5, max = 50,
                value = 10,
                step = 5),
    radioButtons("type", "Prior Type:",
                 c("Uniform" = "unif",
                   "Testing" = "test")),
      h4(id="data-heading", "Enter Data:"),
      tags$style(HTML("#data-heading{color: red;}")),
      h4(id="data-heading1", "Sample 1:"),
  #    tags$style(HTML("#data-heading1{color: red;}")),
      fluidRow(
        column(2,  HTML('<h5><b>s1:</b></h5>')),
        column(4, numericInput("s1", "", min = 0, max = 100, value = 0)),
        column(2, HTML('<h5><b>f1:</b></h5>')),
        column(4, numericInput("f1", "", min = 0, max = 100, value = 0))
      ),
      h4(id="data-heading2", "Sample 2:"),
 #     tags$style(HTML("#data-heading2{color: red;}")),
      fluidRow(
        column(2,  HTML('<h5><b>s2:</b></h5>')),
        column(4, numericInput("s2", "", min = 0, max = 100, value = 0)),
        column(2, HTML('<h5><b>f2:</b></h5>')),
        column(4, numericInput("f2", "", min = 0, max = 100, value = 0))
      )
    )),
    column(8,
           tabsetPanel(type = "tabs",
                       tabPanel("Story",
          br(),
          h4('Description'),
          p('This app considers the sampling model where s1 is binomial
          with sample size s1 + f1 and probability of success p1 and
          s2 (independent of s1) is binomial with sample size s2 + f2 and probability of
          success p2.'),
          p('Consider a prior where each of the proportions p1 and p2 are assumed to
            take on N equally spaced values between LO and HI. One can
            either have a UNIFORM prior where all pairs (p1, p2) are
            given the same probability, or a TESTING prior where
            P(p1 = p2) = 0.5 and equal probabilities are
            assigned to all off-diagonal pairs where p1 ≠ p2,
            and uniform probabilities for the pairs along the diagonal
            where p1 = p2.'),
          h4('Using the App'),
          p('One inputs the prior by choosing Limits for Each P (the LO and HI
          values), the # of Values N,
            and the Prior Type (UNIFORM or TESTING).  In addition, one
            inputs the number of successes and failures for each sample (s1, f1,
          s2, f2).'),
          p('The Data tab shows the table of observed data.'),
          p('The Prior tab displays the joint prior of (p1, p2).'),
          p('The Post tab displays the joint posterior of (p1, p2).'),
          p('The Prior/Post of Diff tab displays the prior and posterior
            probability functions of the difference in probabilities p1 - p2.'),
          p('The Post CDF of Diff tab displays the cumulative distribution
            function of the posterior of the difference p1 - p2.')

          ),
                       tabPanel("Data",
                                tableOutput("stats")),
                       tabPanel("Prior",
                                plotOutput("prior1",
                                           height = "500px")),
                       tabPanel("Post",
                                plotOutput("post1",
                                           height = "500px")),
                       tabPanel("Prior/Post of Diff",
                                plotOutput("post2",
                                           height = "500px")),
                       tabPanel("Post CDF of Diff",
                                plotOutput("post3",
                                           height = "500px"))
           )
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  data <- reactive({
    bivariate_prob <- function(lo, hi, n_values,
                               pequal = 0.5,
                               uniform = FALSE,
                               s1f1 = c(0, 0),
                               s2f2 = c(0, 0)){
      n_diagonal <- n_values
      n_off_diag <- n_values ^ 2 - n_values

      p1 <- seq(lo, hi, length = n_values)
      p2 <- p1
      grid <- expand.grid(p1, p2)
      names(grid) <- c("p1", "p2")
      grid$p1Mp2 <- round(grid$p1 - grid$p2, 3)

      if(uniform == TRUE){
        grid$prior <- 1 / (n_values ^ 2)
      } else {
        grid$prior <- (grid$p1 == grid$p2) *
          pequal / n_values +
          (grid$p1 != grid$p2) *
          (1 - pequal) / n_off_diag
      }

      ## update

      grid$posterior <- log(grid$prior) +
        dbinom(s1f1[1],
               size = sum(s1f1),
               prob = grid$p1,
               log = TRUE) +
        dbinom(s2f2[1],
               size = sum(s2f2),
               prob = grid$p2,
               log = TRUE)
      grid$posterior <- exp(grid$posterior -
                              max(grid$posterior))
      grid$posterior <- grid$posterior /
        sum(grid$posterior)

      ###########

      grid %>%
        group_by(p1Mp2) %>%
        summarize(Prior = sum(prior),
                  Posterior = sum(posterior),
                  .groups = "drop") -> S

      S$Prior <- round(S$Prior, 4)
      S$Posterior <- round(S$Posterior, 4)
      names(S)[1] <- "p1-p2"

      list(grid = grid,
           Summary = S)
    }

    flag <- ifelse(input$type == "unif", TRUE, FALSE)
    out <- bivariate_prob(input$P[1], input$P[2],
                          input$N,
                    pequal = 0.5,
                    uniform = flag,
                    s1f1 = c(input$s1, input$f1),
                    s2f2 = c(input$s2, input$f2))

    mdata <- matrix(c(input$s1, input$f1,
                      input$s2, input$f2),
                    2, 2, byrow = TRUE)
    dimnames(mdata)[[2]] <- c("Successes", "Failures")
    dimnames(mdata)[[1]] <- c("Sample 1", "Sample 2")

    list(grid = out$grid,
         S = out$Summary,
         data = mdata)
  })
  output$prior1 <- renderPlot({
    out <- data()
    ggplot(out$grid,
           aes(p1, p2, size = prior)) +
      geom_point() +
      increasefont() +
      centertitle() +
      ggtitle("Joint Prior Distribution of P1 and P2")
  })
  output$prior2 <- renderPlot({
    S <- data()$S
    names(S)[1] <- "p12"
    prob_plot(select(S, p12, Prior)) +
      increasefont() +
      ggtitle("Prior of p1 MINUS p2") +
      centertitle() +
      xlab("p1 - p2") + ylab("Probability")
  })
  output$post1 <- renderPlot({
    ggplot(data()$grid,
           aes(p1, p2, size = posterior)) +
      geom_point() +
      increasefont() +
      centertitle() +
      ggtitle("Joint Posterior Distribution of P1 and P2")
  })
  output$post2 <- renderPlot({
    S <- data()$S
    names(S)[1] <- "p12"
    M <- max(c(S$Prior, S$Posterior))
    p1 <- prob_plot(select(S, p12, Prior)) +
      increasefont() +
      ggtitle("Prior of p1 MINUS p2") +
      centertitle() +
      ylim(0, M * 1.1) +
      xlab("p1 - p2") + ylab("Probability")
    p2 <- prob_plot(select(S, p12, Posterior)) +
      increasefont() +
      ggtitle("Posterior of p1 MINUS p2") +
      centertitle() +
      ylim(0, M * 1.1) +
      xlab("p1 - p2") + ylab("Probability")
    grid.arrange(p1, p2, ncol = 1)
  })
  output$post3 <- renderPlot({
    S <- data()$S
    names(S)[1] <- "p12"
    S$cprob <- cumsum(S$Posterior)
    ggplot(S, aes(p12, cprob)) +
      geom_line(size = 1.5, color = "red") +
      increasefont() +
      ggtitle("CDF of P1 - P2") +
      centertitle() +
      xlab("p1 - p2") + ylab("Cumulative Probability")
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
    data()$S
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
