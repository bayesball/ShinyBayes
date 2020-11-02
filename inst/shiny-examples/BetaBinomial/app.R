
# Define UI ----
ui <- fluidPage(
  h1(id="big-heading", "Binomial-Beta Multilevel Model"),
  tags$style(HTML("#big-heading{color: red;}")),
  fluidRow(
    column(4, wellPanel(
      h4(id="model-heading", "Model:"),
      tags$style(HTML("#model-heading{color: red;}")),
      h5("y1, ..., yn indep, yi is binomial(ni, pi)"),
      h5("p1, ..., pn iid beta(K eta, K (1 - eta))"),
      h5("eta is beta(a, b)"),
      h5("log K is logistic(logn, 1)"),
      hr(),
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
            tabPanel("Data",
                     tableOutput("contents")),
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
        eta <- exp(theta[1])/(1 + exp(theta[1]))
        K <- exp(theta[2])
        logf <- function(y, n, K, eta){
          lbeta(K * eta + y, K * (1 - eta) + n - y) -
            lbeta(K * eta, K * (1 - eta))
        }
        loglike <- sum(logf(y, n, K, eta))
        loglike + dbeta(eta, ab[1], ab[2], log = TRUE) +
          + dlogis(theta[2], logn, 1, log = TRUE) +
          log(eta * (1 - eta))
    }
    the_data() %>%
        select(1:2) -> df
    datapar <- list(df = df,
                    ab = c(input$a, input$b),
                    logn = input$logn)
    fit <- laplace(bblogpost, c(0, 0), datapar)
    mo <- fit$mode
    sds <- sqrt(diag(fit$var))
    limits <- c(mo[1] - 5 * sds[1],
                mo[1] + 5 * sds[1],
                mo[2] - 5 * sds[2],
                mo[2] + 5 * sds[2])
    sim_post <- simcontour(bblogpost,
                limits, datapar, input$iter)
    eta <- exp(sim_post$x) / (1 + exp(sim_post$x))
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
    new_mo2 <- as.character(c(eta_m, K_m))
    new_sd2 <- as.character(c(eta_sd, K_sd))
    newsumm <- data.frame(Parameter =
                         c("logit eta", "log K",
                           "(logit eta, log K)",
                           "eta", "K"),
              Mean = c(new_mo, "", new_mo2),
            Stan_Dev = c(new_sd, "", new_sd2),
            Correlation = c("", "", new_cor, "", ""))

    sim_post <- data.frame(x = sim_post$x,
                           y = sim_post$y)
    gplot <- gcontour2(bblogpost, limits, datapar) +
      geom_point(data = sim_post, aes(x, y),
                 alpha = 0.5, color = "black")
    p1 <- ggplot(sim_post, aes(x)) +
      geom_density(size = 1.5, color = "red") +
      ggtitle("Marginal Density of Logit eta") +
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
      ggtitle("Posterior of (logit eta, log K)") +
      theme(plot.title = element_text(colour = "blue",
          size = 24,
         hjust = 0.5, vjust = 0.8, angle = 0)) +
      xlab("logit eta") +
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
}

# Run the app ----
shinyApp(ui = ui, server = server)
