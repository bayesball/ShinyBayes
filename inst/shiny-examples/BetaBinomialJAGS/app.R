library(shiny)
library(ggplot2)
library(ProbBayes)
library(metR)
library(dplyr)
library(tidyr)
library(gridExtra)
library(runjags)
library(coda)

# Define UI ----
ui <- fluidPage(
  h1(id="big-heading", "Binomial-Beta Multilevel Model - JAGS"),
  tags$style(HTML("#big-heading{color: red;}")),
  fluidRow(
    column(4, wellPanel(
      h4(id="data-heading", "Read in Data [y, n]:"),
      tags$style(HTML("#data-heading{color: red;}")),
      fileInput("file1", "CSV File",
                accept = ".csv"),
      checkboxInput("header", "Header", TRUE),
      h4(id = 'enter-prior', 'Enter Prior Parameters:'),
      tags$style(HTML("#enter-prior{color: red;}")),
      tags$head(
        tags$style(type="text/css", "label{ display: table-cell; text-align: center; vertical-align: middle; } .form-group { display: table-row;}")
      ),
      numericInput("a", "a: ", min = 0, max = 100, value = 1),
      numericInput("b", "b: ", min = 0, max = 100, value = 1),
      numericInput("logn", "logn: ", min = 0, max = 1000, value = 0),
      hr(),
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
#           plotOutput("plot", height = "380px"),
#           h4(id="stats-heading", "Simulation Posterior Summaries:"),
#           tags$style(HTML("#stats-heading{color: red;}")),
#           verbatimTextOutput("stats2")
      tabsetPanel(type = "tabs",
                  tabPanel("Story",
                           br(),
                           h4('Description'),
                           p('This app fits the following Beta/Binomial multilevel
                      model using the JAGS software.'),
                      p('We observe independent observations y1, ..., yN,
                      where yi is binomial(ni, pi).  Assume p1, ..., pN
                      are a random sample from a beta distribution with
                      shape parameters K eta and K (1 - eta).  At the
                      last stage, eta is assumed beta(a, b) and log K
                      is logistic with location logn and scale 1.'),
                      h4('Using the App'),
                      p('One inputs the data from a csv file where the first
                      column contains the yi and the second column contains
                      the ni.  One inputs values of the hyperparameters a, b,
                      and logn, and the number of simulations from the
                      posterior distribution.  The Update button will start
                      MCMC sampling from the posterior using JAGS.'),
                      p('The outputs are Data, the listing of the observed data,
                      JAGS Script, a listing of the JAGS model
                      script, Contour, a contour graph of the joint posterior of
                      logit eta and log K, Marginals, density estimates of the
                      marginal posterior densities of logit eta and log K.  The
                      tab Summaries provides summaries of the posterior density
                      and Shrinkage shows how the observed rates yi / ni are
                      shrunk to the posterior means of the pi.')
                  ),
            tabPanel("Data",
                     tableOutput("contents")),
            tabPanel("JAGS Script",
                     verbatimTextOutput("jagscript")),
            tabPanel("Contour",
                     plotOutput("plot",
                                height = "500px")),
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

   data <- eventReactive(input$goButton, {
     bblogpost <- function(theta, datapar){
        y <- datapar$df[, 1]
        n <- datapar$df[, 2]
        ab <- datapar$ab
        logn <- datapar$logn
        eta <- theta[1]
        K <- exp(theta[2])
        logf <- function(y, n, K, eta){
          lbeta(K * eta + y, K * (1 - eta) + n - y) -
            lbeta(K * eta, K * (1 - eta))
        }
        loglike <- sum(logf(y, n, K, eta))
        loglike + dbeta(eta, ab[1], ab[2], log = TRUE) +
          + dlogis(theta[2], logn, 1, log = TRUE)
    }
    the_data() %>%
        select(1:2) -> df
    datapar <- list(df = df,
                    ab = c(input$a, input$b),
                    logn = input$logn)
    fit <- laplace(bblogpost, c(0.5, 0), datapar)
    mo <- fit$mode
    sds <- sqrt(diag(fit$var))
    limits <- c(mo[1] - 5 * sds[1],
                mo[1] + 5 * sds[1],
                mo[2] - 5 * sds[2],
                mo[2] + 5 * sds[2])

#    sim_post <- data.frame(x = sim_post$x,
#                           y = sim_post$y)
    the_data <- list("y" = df[, 1],
                      "n" = df[, 2],
                      "N" = dim(df)[1],
                      "mua" = input$a,
                      "mub" = input$b,
                      "logn" = input$logn)
    runjags.options(silent.jags = TRUE,
                    silent.runjags = TRUE)
    options(warn = -1)
    posterior <- run.jags(JAGS_script("beta_binomial"),
                          n.chains = 1,
                          data = the_data,
                          monitor = c("mu", "logeta"),
                          adapt = 1000,
                          burnin = 1000,
                          sample = input$iter)
    sim_post <- as.data.frame(as.mcmc(posterior))
    names(sim_post) <- c("x", "y")

    eta <- sim_post$x
    K <- exp(sim_post$y)
    new_mo <- as.character(round(c(mean(sim_post$x),
                                   mean(sim_post$y)), 3))
    new_sd <- as.character(round(c(sd(sim_post$x),
                                   sd(sim_post$y)), 3))
    new_cor <- as.character(round(cor(sim_post$x,
                                      sim_post$y), 3))
    eta_m <- round(mean(eta), 4)
    eta_sd <- round(sd(eta), 4)
    K_m <- round(mean(K), 1)
    K_sd <- round(sd(K), 1)
    new_mo2 <- as.character(c(K_m))
    new_sd2 <- as.character(c(K_sd))
    newsumm <- data.frame(Parameter =
                            c("eta", "log K",
                              "(eta, log K)",
                              "K"),
                          Mean = c(new_mo, "", new_mo2),
                          Stan_Dev = c(new_sd, "", new_sd2),
                          Correlation = c("", "", new_cor, ""))

    gplot <- gcontour2(bblogpost, limits, datapar) +
      geom_point(data = sim_post, aes(x, y),
                 alpha = 0.5, color = "black")
    p1 <- ggplot(sim_post, aes(x)) +
      geom_density(size = 1.5, color = "red") +
      ggtitle("Marginal Density of eta") +
      increasefont() + xlab("") +
      centertitle()
    p2 <- ggplot(sim_post, aes(y)) +
      geom_density(size = 1.5, color = "red") +
      ggtitle("Marginal Density of Log K") +
      increasefont() + xlab("") +
      centertitle()
    mgplot <- grid.arrange(p1, p2, ncol = 1)
    N <- dim(df)[1]
    estimates <- data.frame(Player = 1:N,
              Individual = df[, 1] / df[, 2],
              Multilevel = (df[, 1] + K_m * eta_m) /
                (df[, 2] + K_m))
    list(summary = newsumm,
         plot = gplot,
         mplot1 = p1,
         mplot2 = p2,
         estimates = estimates)
  })
  output$plot <- renderPlot({
    data()$plot +
      increasefont() +
      ggtitle("Posterior of (eta, log K)") +
      theme(plot.title = element_text(colour = "blue",
          size = 24,
         hjust = 0.5, vjust = 0.8, angle = 0)) +
      xlab("eta") +
      ylab("log K") +
      theme(
        axis.title = element_text(color = "blue",
                    size = 18))
  })
  output$mplot <- renderPlot({
    out <- data()
    grid.arrange(out$mplot1,
                 out$mplot2)
  })
  output$shrink <- renderPlot({
    shrinkage_plot <- function (S, N = 15) {
      S2 <- gather(sample_n(S, size = N),
                   Type, AVG, -Player)
      S2a <- filter(S2, Type == "Individual")
      ggplot(S2, aes(Type, AVG, group = Player)) + geom_line() +
        geom_point() +
        ggplot2::annotate("text", x = 0.75,
                          y = S2a$AVG,
                          label = S2a$Player) +
        ggtitle("Shrinkage Plot") +
        theme(plot.title =
                element_text(colour = "blue",
                  size = 18, hjust = 0.5))
    }
    out <- data()$estimates
    N1 <- min(dim(out)[1], 15)
    shrinkage_plot(out, N = N1) +
      increasefont() +
      xlab("Proportion Estimate") +
      ylab("Estimate")
  })
  output$stats2 <- renderDataTable({
    data()$summary
  })
  output$jagscript <- renderText({
    JAGS_script("beta_binomial")
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
