library(shiny)
library(ggplot2)
library(ProbBayes)
library(dplyr)
library(tidyr)
library(gridExtra)
library(runjags)
library(coda)

# Define UI ----
ui <- fluidPage(
  h1(id="big-heading", "Normal-Normal Multilevel Model - JAGS"),
  tags$style(HTML("#big-heading{color: red;}")),
  fluidRow(
    column(4, wellPanel(
      h4(id="model-heading", "Model:"),
      tags$style(HTML("#model-heading{color: red;}")),
      h5("y_1, ..., y_N indep, y_i is normal(mu_j[group_i], sigma)"),
      h5("mu_1, ..., mu_J iid normal(mu, tau)"),
      h5("mu is t(3, mu_0, tau_0)"),
      h5("tau is t(3, 0, tau_1) I(tau > 0)"),
      h5("sigma is t(3, 0, tau_2) I(sigma > 0)"),
      hr(),
      h4(id="data-heading", "Read in Data [y, group]:"),
      tags$style(HTML("#data-heading{color: red;}")),
      fileInput("file1", "CSV File",
                accept = ".csv"),
      checkboxInput("header", "Header", TRUE),
      h4(id = 'enter-prior', 'Prior Parameters:'),
      tags$style(HTML("#enter-prior{color: red;}")),
      tags$head(
        tags$style(type="text/css", "label{ display: table-cell; text-align: center; vertical-align: middle; } .form-group { display: table-row;}")
      ),
      fluidRow(
        column(2,  HTML('<h5><b>mu_0:</b></h5>')),
        column(4, numericInput("mu0", "", min = 0, max = 100, value = 0)),
        column(2, HTML('<h5><b>tau_0:</b></h5>')),
        column(4, numericInput("tau0", "", min = 0, max = 100, value = 1))
      ),
      fluidRow(
        column(2,  HTML('<h5><b>tau_1:</b></h5>')),
        column(4, numericInput("tau1", "", min = .001, max = 100, value = 12.6)),
        column(2, HTML('<h5><b>tau_2:</b></h5>')),
        column(4, numericInput("tau2", "", min = .001, max = 100, value = 12.6))
      ),
      h4(id = 'enter-iter', '# of Simulations:'),
      tags$style(HTML("#enter-iter{color: red;}")),
      tags$head(
        tags$style(type="text/css", "label{ display: table-cell; text-align: center; vertical-align: middle; } .form-group { display: table-row;}")
      ),
      fluidRow(
        column(1,  HTML('<h5><b>iter:</b></h5>')),
        column(5, numericInput("iter", "", min = 1000, max = 10000, value = 1000)),
      ),
      tags$head(
        tags$style(HTML('#goButton{background-color:orange}'))
      ),
      actionButton("goButton", "UPDATE")
#      tableOutput("contents")
    )),
    column(8,
      tabsetPanel(type = "tabs",
            tabPanel("Data",
                     tableOutput("contents")),
            tabPanel("JAGS Script",
                     verbatimTextOutput("jagscript")),
            tabPanel("Marginals",
                     plotOutput("mplot",
                                height = "500px")),
            tabPanel("Summaries",
                     dataTableOutput("stats2")),
            tabPanel("Shrinkage",
                     plotOutput("shrink",
                                height = "500px"))
)
    )
  )
)
# Define server logic ----
server <- function(input, output) {
  the_data <- reactive({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    read.csv(file$datapath, header = input$header)
  })
  output$contents <- renderTable({
    the_data()
  })
  data <- eventReactive(input$goButton, {
    the_data() %>%
        select(1:2) -> df
    names(df)[1] <- "y"
    df$group <- as.numeric(factor(df[, 2]))
    N <- length(df$y)
    J <- max(df$group)
    the_data <- list("y" = df$y, "group" = df$group,
                     "N" = N, "J" = J,
                     "mu0" = input$mu0,
                     "invtau_0" =  1 / input$tau0 ^ 2,
                     "invtau_1" = 1 / input$tau1 ^ 2,
                     "invtau_2" = 1 / input$tau2 ^ 2)
    runjags.options(silent.jags = TRUE,
                    silent.runjags = TRUE)
    options(warn = -1)
    ModelScript <- "model {
for (i in 1:N){
y[i] ~ dnorm(mu_j[group[i]], invsigma2)
}
for (j in 1:J){
  mu_j[j] ~ dnorm(mu, invtau2)
}
invsigma2 <- 1 / sigma / sigma
sigma ~ dt(0, invtau_2, 3) T(0, )
mu ~ dt(mu0, invtau_0, 3)
invtau2 <- 1 / tau / tau
tau ~ dt(0, invtau_1, 3) T(0, )
}"
    posterior <- run.jags(ModelScript,
                          n.chains = 1,
                          data = the_data,
                          monitor = c("mu", "tau", "sigma"),
                          adapt = 1000,
                          burnin = 1000,
                          sample = input$iter)
    S <- round(summary(posterior), 3)
    newdf <- data.frame(Parameter = row.names(S),
                        S[, 1:5])
    row.names(newdf) <- NULL

    post <- as.data.frame(as.mcmc(posterior))
    post$Iteration <- 1:input$iter
    p1a <- ggplot(post, aes(Iteration, mu)) +
      geom_line() + ggtitle("MU") +
      increasefont() + centertitle()
    p1b <- ggplot(post, aes(mu)) +
      geom_density() + ggtitle("MU") +
      increasefont() + centertitle() + xlab("")
    p2a <- ggplot(post, aes(Iteration, sigma)) +
      geom_line() + ggtitle("SIGMA") +
      increasefont() + centertitle()
    p2b <- ggplot(post, aes(sigma)) +
      geom_density() + ggtitle("SIGMA") +
      increasefont() + centertitle()+ xlab("")
    p3a <- ggplot(post, aes(Iteration, tau)) +
      geom_line() + ggtitle("TAU") +
      increasefont() + centertitle()
    p3b <- ggplot(post, aes(tau)) +
      geom_density() + ggtitle("TAU") +
      increasefont() + centertitle()+ xlab("")
    mcmcplot <- grid.arrange(p1a, p1b, p2a, p2b, p3a, p3b,
                 ncol = 2)

    df$group <- as.numeric(factor(df[, 2]))
    df %>%
      group_by(group) %>%
      summarize(N = n(),
                Individual = mean(y),
                .groups = "drop") -> S

    pm <- apply(post, 2, mean)

    S$Multilevel <- (S$Individual * S$N / pm[3] ^ 2 +
                     pm[1] / pm[2] ^ 2) /
      (S$N / pm[3] ^ 2 + 1 / pm[2] ^ 2)

    estimates <- select(S, group, Individual, Multilevel)
    names(estimates)[1] <- "Player"
    list(summary = newdf,
         mcmcplot = mcmcplot,
         estimates = estimates)
  })
  output$mplot <- renderPlot({
    data()$mcmcplot
    #####
  })
  output$stats2 <- renderDataTable({
    data()$summary
  })
  output$jagscript <- renderText({
"model {
for (i in 1:N){
y[i] ~ dnorm(mu_j[group[i]], invsigma2)
}
for (j in 1:J){
  mu_j[j] ~ dnorm(mu, invtau2)
}
invsigma2 <- 1 / sigma / sigma
sigma ~ dt(0, invtau_2, 3) T(0, )
mu ~ dt(mu0, invtau_0, 3)
invtau2 <- 1 / tau / tau
tau ~ dt(0, invtau_1, 3) T(0, )
}"
  })
  output$shrink <- renderPlot({
    shrinkage_plot <- function (S, N = 15) {
      S2 <- gather(sample_n(S, size = N),
                   Type, AVG, -Player)
      S2a <- filter(S2, Type == "Individual")
      ggplot(S2, aes(Type, AVG, group = Player)) +
        geom_line(color = "red") +
        geom_point(color = "blue") +
        ggplot2::annotate("text", x = 0.75,
                          y = S2a$AVG,
                          label = S2a$Player,
                          color = "blue") +
        ggtitle("Shrinkage Plot") +
        theme(plot.title =
                element_text(colour = "blue",
                             size = 18, hjust = 0.5))
    }
    out <- data()$estimates
    N1 <- min(dim(out)[1], 15)
    shrinkage_plot(out, N = N1) +
      increasefont() +
      xlab("Mean Estimate") +
      ylab("Estimate")
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
