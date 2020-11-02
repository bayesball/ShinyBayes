library(shiny)
library(ProbBayes)

ui <- fluidPage(
  pageWithSidebar(
  headerPanel(
    "Discrete Bayes for a Proportion",
    ),
  sidebarPanel(
    tags$head(
      tags$style(type="text/css", "label{ display: table-cell; text-align: center; vertical-align: middle; } .form-group { display: table-row;}")
    ),
    h4('Enter Prior Weights:'),
    numericInput("p0", "p = 0: ", min = 0, max = 10, value = 1),
    numericInput("p1", "p = 0.1: ", min = 0, max = 10, value = 1),
    numericInput("p2", "p = 0.2: ", min = 0, max = 10, value = 1),
    numericInput("p3", "p = 0.3: ", min = 0, max = 10, value = 1),
    numericInput("p4", "p = 0.4: ", min = 0, max = 10, value = 1),
    numericInput("p5", "p = 0.5: ", min = 0, max = 10, value = 1),
    numericInput("p6", "p = 0.6: ", min = 0, max = 10, value = 1),
    numericInput("p7", "p = 0.7: ", min = 0, max = 10, value = 1),
    numericInput("p8", "p = 0.8: ", min = 0, max = 10, value = 1),
    numericInput("p9", "p = 0.9: ", min = 0, max = 10, value = 1),
    numericInput("p10", "p = 1: ", min = 0, max = 10, value = 1),
    br(),
    h4("Enter Data:"),
    numericInput("s", "# of Successes", min = 0, max = 30, value = 0),
    numericInput("f", "# of Failures", min = 0, max = 30, value = 0),
    br(),
    tags$head(
      tags$style(HTML('#goButton{background-color:orange}'))
    ),
    actionButton("goButton", "Update"),
  ),
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Table",
                         dataTableOutput("table")),
                tabPanel("Graph",
                         plotOutput("myplot",
                                height = "500px")),
                tabPanel("Summary",
                         dataTableOutput("summary"))
    )
  )
)
)

server <- function(input, output) {
  # builds a reactive expression that only invalidates
  # when the value of input$goButton becomes out of date
  # (i.e., when the button is pressed)
  constructprior <- eventReactive(input$goButton, {
    wt <- c(input$p0, input$p1, input$p2, input$p3,
            input$p4, input$p5, input$p6,
            input$p7, input$p8, input$p9, input$p10)
    d <- data.frame(
       P = seq(0, 1, by = 0.1),
       Weight = wt,
       Prior = wt / sum(wt)
    )
    d$Likelihood <- dbinom(input$s,
                           size = input$s + input$f,
                           prob = d$P)
    round(bayesian_crank(d), 3)
  })

  summarize_post <- eventReactive(input$goButton, {
    wt <- c(input$p0, input$p1, input$p2, input$p3,
            input$p4, input$p5, input$p6,
            input$p7, input$p8, input$p9, input$p10)
    d <- data.frame(
      P = seq(0, 1, by = 0.1),
      Weight = wt,
      Prior = wt / sum(wt)
    )
    d$Likelihood <- dbinom(input$s,
                           size = input$s + input$f,
                           prob = d$P)
    b_table <- bayesian_crank(d)
    postm <- sum(b_table$P * b_table$Posterior)
    posts <- sqrt(sum((b_table$P - postm) ^ 2 *
                        b_table$Posterior))
    bint <- discint(b_table[, c("P", "Posterior")],
                    0.9)
    values1 <- c(as.character(
                    c(round(c(postm, posts), 3))),
                     as.character(min(bint$set)))
    values2 <- c("", "",
                 as.character(max(bint$set)))
    data.frame(Summary = c("Posterior Mean",
                           "Posterior Standard Deviation",
                           "90% Probability Interval"),
               Value = values1,
               Value2 = values2)
  })

  priorpostplot <- eventReactive(input$goButton, {
    wt <- c(input$p0, input$p1, input$p2, input$p3,
            input$p4, input$p5, input$p6,
            input$p7, input$p8, input$p9, input$p10)
    d <- data.frame(
      P = seq(0, 1, by = 0.1),
      Weight = wt,
      Prior = wt / sum(wt)
    )
    d$Likelihood <- dbinom(input$s,
                           size = input$s + input$f,
                           prob = d$P)
    d <- bayesian_crank(d)
    prior_post_plot(d) + xlab("P") +
      theme(text = element_text(size = 18))
  })

  output$table <- renderDataTable({
    constructprior()
  })
  output$myplot <- renderPlot({
    priorpostplot()
  })
  output$summary <-  renderDataTable({
    summarize_post()
  })
}

shinyApp(ui, server)
