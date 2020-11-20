library(shiny)

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Project App"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      tags$div(
        tags$p (" Welcome to This Nonparametric Tests Shiny App!"),
        br(),
        tags$p("The app is designed for users to use Kruskal-Wallis test,Fisher’s Comparison Test,and Jonckheere-Terpstra Test"),
        br(),
        tags$p("Before you use this app, it is required that file uploaded should be in the txt format. Additionally, it should be formatted 
               in such a way that groups should be in rows and blocks should be in columns.It's also required to have equal sample size in each group."),
        br(),
        tags$p ("Have Fun !")
        ),
      
      # Input: Select a file ----
      fileInput("file1", "Choose File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".txt")),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Checkbox if file has header ----
      checkboxInput("header", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      # Input: Select quotes ----
      radioButtons("quote", "Quote",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = '"'),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Select number of rows to display ----
      radioButtons("disp", "Display",
                   choices = c(Head = "head",
                               All = "all"),
                   selected = "head"),
      
  
      sliderInput("samplesize", " Simulation Sample Size:",
                  min=500, max=5000, value=500)
      
    ),
    
   
  
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      tableOutput("contents"),
      
      # Output: Boxplots ----
      plotOutput("boxplots"),
      plotOutput("histogram_JT"),
      
      
      # Output: P-value Permutation Test ----
      textOutput("kwp_value"),
      textOutput("jtp_value"),
      textOutput("list")
      
    )
    
  )
  )


# Define server logic to read selected file ----
server <- function(input, output) {
  
  
  data <- reactive({
    
    req(input$file1)
    
    df <- read.table(input$file1$datapath,
                     header = input$header,
                     sep = input$sep,
                     quote = input$quote)[,-1]
    
    if(input$disp == "head") {
      
      return(head(df))
      
    } else {
      
      return(df)
      
    }
    
  })
  
  
  
  output$contents <- renderTable({
    
    data()
    
  })
  
  output$boxplots <- renderPlot({
    
    boxplot(data(), col = "grey80")
    
  })
  
  
  ###########################################################################################################
  output$kwp_value <- renderText({
      KW <- function(X) {
      combined <- c(as.matrix(X))
      N <- length(combined)
      ranks <- rank(combined)
      rank_mat <- matrix(ranks, dim(X)[1], dim(X)[2])
      K <- dim(X)[1]
      varranks <- rep(NA, K)
      exp_rank <- (N + 1) / 2
      
      for (i in 1:K){
        varranks[i] <- dim(X)[2]*(mean(rank_mat[i,]) - exp_rank)^2 
      }
      KWobs <- (12 / (N * (N + 1))) * sum(varranks)
      KWobs
    }
    KWobs <- KW(data())
    alpha <- 0.05
    df <- K - 1
    R <- input$samplesize
    set.seed(123)
    KWstar <- NA
    for(i in 1:R) {
      xstar <- matrix(sample(combined, N), dim(data())[1], dim(data())[2])
      KWstar[i] <- KW(xstar)   
    }
    
    support <- seq(0, max(KWstar), length.out = 1000)
    p_value <- sum(KWstar >= KWobs) / R
    paste("The p-value for the KW test is", p_value)})
  
  
  ################################################################################################
  
  ## JT Test
  output$jtp_value <- renderText({
    n <- dim(data())[2]
    K <- dim(data())[1]
    m <- K * (K - 1) /2
    # Function for permutation
    JT <- function(X) {
      JTobs <- 0
      comparisons <- expand.grid(1:K, 1:K)
      comparisons <- comparisons[comparisons[, 1] < comparisons[, 2], ]
      m <- dim(comparisons)[1]
      
      for(i in 1:m) {
        
        pairs_temp <- expand.grid(X[, as.numeric(comparisons[i,1])], X[, as.numeric(comparisons[i,2])])
        JTobs <- JTobs + sum(apply(pairs_temp, 1, function(x) x[1] < x[2]))
      }
      JTobs
    }
    # Compute JT statistic
    JTobs <- JT(data())
    # Permutation null distribution
    R <- input$samplesize
    set.seed(123)
    combined <- c(as.matrix(data()))
    N <- length(combined)
    JTstar <- NA
    set.seed(123)
    
    for(i in 1:R) {
      
      perm <- matrix(sample(combined, N), dim(data())[1], dim(data())[2])
      JTstar[i] <- JT(perm)
      
    }
    
    
    jtp_value <- sum(JTstar >= JTobs) / R
    paste("The p-value for the JT test is", jtp_value)
  })
  
  output$histogram_JT <- renderPlot({
    n <- dim(data())[2]
    K <- dim(data())[1]
    m <- K * (K - 1) /2
    
    # Function for permutation
    JT <- function(X) {
      
      JTobs <- 0
      comparisons <- expand.grid(1:K, 1:K)
      comparisons <- comparisons[comparisons[, 1] < comparisons[, 2], ]
      m <- dim(comparisons)[1]
      
      for(i in 1:m) {
        
        pairs_temp <- expand.grid(X[, as.numeric(comparisons[i,1])], X[, as.numeric(comparisons[i,2])])
        JTobs <- JTobs + sum(apply(pairs_temp, 1, function(x) x[1] < x[2]))
        
      }
      
      JTobs
      
    }
    
    JTobs <- JT(data())
    # Permutation null distribution
    R <- input$samplesize
    set.seed(123)
    combined <- c(as.matrix(data()))
    N <- length(combined)
    
    JTstar <- NA
    
    set.seed(123)
    
    for(i in 1:R) {
      
      perm <- matrix(sample(combined, N), dim(data())[1], dim(data())[2])
      JTstar[i] <- JT(perm)
      
    }
    
    
    
    hist(JTstar, col = "grey80", freq = FALSE)
    abline(v = JTobs, col = "red")})
  
 
}






# Run the app ----
shinyApp(ui, server)

