library(shiny)
ui <- fluidPage( titlePanel("F test "),
                 
                 # Sidebar layout with input and output definitions ----
                 sidebarLayout(
                   
                   # Sidebar panel for inputs ----
                   sidebarPanel(
                     
                     # Input: Select a file ----
                     fileInput("file1", "Choose CSV File For F test",
                               multiple = TRUE,
                               accept = c("text/csv",
                                          "text/comma-separated-values,text/plain",
                                          ".csv"))
                     
                     
                   ),
                   
                   # Main panel for displaying outputs ----
                   mainPanel(
                     tabsetPanel(
                       tabPanel("F test", 
                                tableOutput("contents"),
                                plotOutput("perm_test_plot"),
                                textOutput("p_value"),
                                textOutput("Tp_value")
                       ),
                       tabPanel("Chi-Square Goodness of Fit", 
                                br(),
                                actionButton("chisq","Chi-Square Goodness of Fit"),
                                
                                tableOutput("chisqcontents"),
                                textOutput("chistat"),
                                textOutput("chicrit"),
                                textOutput("chipval"),
                                textOutput("chiconclusion"),
                                plotOutput("chisqplot")),
                       tabPanel("Wilcoxon", 
                                br(),
                                actionButton("wilcoxon","Wilcoxon"),
                                
                                tableOutput("wilcoxtable"),
                                plotOutput("wilcoxplot"),
                                textOutput("wilcoxperms"),
                                textOutput("wilcoxpvalue"))
                     )
                   )
                 )
)

server <- function(input, output) {
  
  data <- reactive({
    
    req(input$file1)
    
    df <- read.csv(input$file1$datapath, header = TRUE, sep=",")
    
  })
  
  output$contents <- renderTable({
    
    data()
    
  })
  
  
  output$p_value <- renderText({
    df <- data()
    group1 <- df[,'Group1']
    group2 <- df[,'Group2']
    group3 <- df[,'Group3']
    
    K <- 3
    n <- 3
    N <- n*K # since all groups have n observations
    
    # F-Statistic
    Xbar <-  (mean(group1) + mean(group2) + mean(group3))/3
    mean1 <- mean(group1) - Xbar
    mean2 <- mean(group2) - Xbar
    mean3 <- mean(group3) - Xbar
    
    sigma1 <- var(group1)
    sigma2 <- var(group2)
    sigma3 <- var(group3)
    
    SST <- 3*(mean1^2 + mean2^2 + mean3^2) # since n is the same for all groups
    MST <- SST / (K - 1)
    
    SSE <- 2*(sigma1 + sigma2 + sigma3) # since n is the same for all groups so n - 1 = 2
    MSE <- SSE / (N - K)
    
    Fobs <- MST / MSE
    
    Fstat <- function(group1, group2, group3) {
      
      Xbar <- (mean(group1) + mean(group2) +mean(group3))/3
      mean1 <- mean(group1) - Xbar
      mean2 <- mean(group2) - Xbar
      mean3 <- mean(group3) - Xbar
      
      sigma1 <- var(group1)
      sigma2 <- var(group2)
      sigma3 <- var(group3)
      
      SST <- 3*(mean1^2 + mean2^2 + mean3^2) # since n is the same for all groups
      MST <- SST / (K - 1)
      
      SSE <- 2*(sigma1 + sigma2 + sigma3) # since n is the same for all groups so n - 1 = 2
      MSE <- SSE / (N - K)
      
      Fobs <- MST / MSE
      Fobs
      
    }
    
    # Random Permutation Sampling
    R <- 1000
    Fstar <- NA
    combined_obs <- c(group1, group2, group3)
    
    set.seed(123)
    
    for(i in 1:R) {
      
      xstar <- sample(combined_obs, N)
      Fstar[i] <- Fstat(xstar[1:3], xstar[4:6], xstar[7:9])
      
    }
    
    #support <- seq(0, max(Fstar), length.out = 1000)
    #hist(Fstar, col = "grey80", freq = FALSE)
    #abline(v = Fobs, col = "red")
    #abline(v = qf(0.95, K - 1, N - K), col = "blue")
    #lines(support, df(support, K - 1, N - K), col = "blue")
    
    p_val <- sum(Fstar >= Fobs) / R
    #p_val
    paste("The p-value for the permutation test is",p_val )
    
    # Theoretical F-test
    #p_val_norm <- pf(Fobs, K - 1, N - K, lower.tail = FALSE)
    #p_val_norm
    
    
  })
  output$Tp_value <- renderText({
    df <- data()
    group1 <- df[,'Group1']
    group2 <- df[,'Group2']
    group3 <- df[,'Group3']
    
    K <- 3
    n <- 3
    N <- n*K # since all groups have n observations
    
    # F-Statistic
    Xbar <-  (mean(group1) + mean(group2) + mean(group3))/3
    mean1 <- mean(group1) - Xbar
    mean2 <- mean(group2) - Xbar
    mean3 <- mean(group3) - Xbar
    
    sigma1 <- var(group1)
    sigma2 <- var(group2)
    sigma3 <- var(group3)
    
    SST <- 3*(mean1^2 + mean2^2 + mean3^2) # since n is the same for all groups
    MST <- SST / (K - 1)
    
    SSE <- 2*(sigma1 + sigma2 + sigma3) # since n is the same for all groups so n - 1 = 2
    MSE <- SSE / (N - K)
    
    Fobs <- MST / MSE
    
    Fstat <- function(group1, group2, group3) {
      
      Xbar <- (mean(group1) + mean(group2) +mean(group3))/3
      mean1 <- mean(group1) - Xbar
      mean2 <- mean(group2) - Xbar
      mean3 <- mean(group3) - Xbar
      
      sigma1 <- var(group1)
      sigma2 <- var(group2)
      sigma3 <- var(group3)
      
      SST <- 3*(mean1^2 + mean2^2 + mean3^2) # since n is the same for all groups
      MST <- SST / (K - 1)
      
      SSE <- 2*(sigma1 + sigma2 + sigma3) # since n is the same for all groups so n - 1 = 2
      MSE <- SSE / (N - K)
      
      Fobs <- MST / MSE
      Fobs
      
    }
    
    # Random Permutation Sampling
    R <- 1000
    Fstar <- NA
    combined_obs <- c(group1, group2, group3)
    
    set.seed(123)
    
    for(i in 1:R) {
      
      xstar <- sample(combined_obs, N)
      Fstar[i] <- Fstat(xstar[1:3], xstar[4:6], xstar[7:9])
      
    }
    
    #support <- seq(0, max(Fstar), length.out = 1000)
    #hist(Fstar, col = "grey80", freq = FALSE)
    #abline(v = Fobs, col = "red")
    #abline(v = qf(0.95, K - 1, N - K), col = "blue")
    #lines(support, df(support, K - 1, N - K), col = "blue")
    
    #p_val <- sum(Fstar >= Fobs) / R
    #p_val
    
    
    # Theoretical F-test
    p_val_norm <- pf(Fobs, K - 1, N - K, lower.tail = FALSE)
    #p_val_norm
    paste("The Theoretical p-value for the permutation test is",p_val_norm )
    
    
  })
  output$perm_test_plot <- renderPlot({
    df <- data()
    group1 <- df[,'Group1']
    group2 <- df[,'Group2']
    group3 <- df[,'Group3']
    
    K <- 3
    n <- 3
    N <- n*K # since all groups have n observations
    
    # F-Statistic
    Xbar <-  (mean(group1) + mean(group2) + mean(group3))/3
    mean1 <- mean(group1) - Xbar
    mean2 <- mean(group2) - Xbar
    mean3 <- mean(group3) - Xbar
    
    sigma1 <- var(group1)
    sigma2 <- var(group2)
    sigma3 <- var(group3)
    
    SST <- 3*(mean1^2 + mean2^2 + mean3^2) # since n is the same for all groups
    MST <- SST / (K - 1)
    
    SSE <- 2*(sigma1 + sigma2 + sigma3) # since n is the same for all groups so n - 1 = 2
    MSE <- SSE / (N - K)
    
    Fobs <- MST / MSE
    
    Fstat <- function(group1, group2, group3) {
      
      Xbar <- (mean(group1) + mean(group2) +mean(group3))/3
      mean1 <- mean(group1) - Xbar
      mean2 <- mean(group2) - Xbar
      mean3 <- mean(group3) - Xbar
      
      sigma1 <- var(group1)
      sigma2 <- var(group2)
      sigma3 <- var(group3)
      
      SST <- 3*(mean1^2 + mean2^2 + mean3^2) # since n is the same for all groups
      MST <- SST / (K - 1)
      
      SSE <- 2*(sigma1 + sigma2 + sigma3) # since n is the same for all groups so n - 1 = 2
      MSE <- SSE / (N - K)
      
      Fobs <- MST / MSE
      Fobs
      
    }
    
    # Random Permutation Sampling
    R <- 1000
    Fstar <- NA
    combined_obs <- c(group1, group2, group3)
    
    set.seed(123)
    
    for(i in 1:R) {
      
      xstar <- sample(combined_obs, N)
      Fstar[i] <- Fstat(xstar[1:3], xstar[4:6], xstar[7:9])
      
    }
    
    support <- seq(0, max(Fstar), length.out = 1000)
    hist(Fstar, col = "azure2", freq = FALSE)
    abline(v = Fobs, col = "aquamarine")
    abline(v = qf(0.95, K - 1, N - K), col = "blue1")
    lines(support, df(support, K - 1, N - K), col = "blue1")
    
    #p_val <- sum(Fstar >= Fobs) / R
    #p_val
    #paste("The p-value for the permutation test is",p_val )
    
    # Theoretical F-test
    #p_val_norm <- pf(Fobs, K - 1, N - K, lower.tail = FALSE)
    #p_val_norm
  })
  
  ### CHI SQUARE ###
  observeEvent(input$chisq, {
    
    colors <- c("Brown", "Yellow", "Orange", "Green", "Coffee")
    observed <- c(224, 119, 130, 48, 59)
    probs_null <- c(0.4, 0.2, 0.2, 0.1, 0.1)
    expected <- sum( observed ) * probs_null
    table <- as.data.frame(rbind(observed, probs_null, expected))
    
    names(table) <- colors
    table[, 6] <- rowSums(table)
    names(table)[6] <- "Sum"
    row.names(table) <- c("Observed", "Null Probabilities", "Expected")
    
    df <- length(observed) - 1
    alpha <- 0.05
    
    
    output$chisqcontents <- renderTable({
      table
    })
    
    output$chistat <- renderText({
      
      chi.stat <- sum ( (observed - expected)^2 / expected )
      paste("The test statistic for the Chi-Square test is",round(chi.stat,digits=2))
      
    })
    
    output$chicrit <- renderText ({
      crit.stat <- qchisq(alpha,df,lower.tail = FALSE)
      paste("The critical value  for the Chi-Square test is",round(crit.stat,digits=2))
    })
    
    output$chipval <- renderText ({
      chi.stat <- sum ( (observed - expected)^2 / expected )
      p_val <- pchisq(chi.stat, df, lower.tail = FALSE)
      paste("The p-value for the Chi-Square test is",round(p_val,digits=2))
    })
    
    output$chiconclusion <- renderText ({
      chi.stat <- sum ( (observed - expected)^2 / expected )
      p_val <- pchisq(chi.stat, df, lower.tail = FALSE)
      paste("With a p-value of",round(p_val,digits=2),"we fail to reject the null hypothesis and conclude that the assumed probabilities are a good fit for this dataset")
    })
    
    output$chisqplot <- renderPlot ({
      chi.stat <- sum ((observed - expected)^2 / expected)
      crit.stat <- qchisq(alpha,df,lower.tail = FALSE)
      x <- rchisq(100, 4)
      hist(x, prob=TRUE, col = "azure2")
      abline(v = crit.stat, col = "aquamarine4")
      abline(v = chi.stat, col = "blue4")
      
    })
  })
  
  observeEvent(input$wilcoxon, {
    x = c(37, 49, 55, 57)
    y = c(23, 31, 46, 40)
    m = length(x)
    n = length(y)
    size = m + n
    start = m + 1
    combined = c(x,y)
    ranks = rank(combined)
    wobs = sum(ranks[1:m])
    
    # Used to choose upper or lower test
    yobs = sum(ranks[start:size])
    nperm = choose(size,m)
    
    # Uses an upper test if X ranks are larger than Y ranks
    # Calculates an upper tail p value using permutations
    if (wobs >= yobs){
      id = "higher"
      if (nperm <= 10000){
        perms = t(combn(ranks, m))
        wperm = apply(perms, 1, sum)
        wsum = sum(wperm >= wobs)
        pval <- sum(wperm >= wobs) / length(wperm)
      }
      # Calculates p value using normal approximation
      else{
        wexp = (m*(size + 1))/2
        wsd = sqrt(m*n*(size + 1)/12)
        zstat = (wobs - wexp) / wsd
        pval <- pnorm(zstat, lower.tail = FALSE)
      }
    }
    
    # Uses a lower test if X ranks are smaller than Y ranks
    # Calculates a lower tail p value using permutations
    else{
      id = "lower"
      if (nperm <= 10000){
        perms = t(combn(ranks, m))
        wperm = apply(perms, 1, sum)
        wsum = sum(wperm <= wobs)
        pval <- sum(wperm <= wobs) / length(wperm)
      }
      
      # Calculates p value using normal approximation
      else{
        wexp = (m*(size + 1))/2
        wsd = sqrt(m*n*(size + 1)/12)
        zstat = (wobs - wexp) / wsd
        pval <- pnorm(zstat, lower.tail = TRUE)
      }
    }
    
    # Output information for wilcoxon
    wrow1 = c(x,sum(x))
    wrow2 = c(y,sum(y))
    wrow3 = c(ranks[1:m],wobs)
    wrow4 = c(ranks[start:size],yobs)
    wtable = as.data.frame(rbind(wrow1,wrow2,wrow3,wrow4))
    names(wtable) = c("Obs1","Obs2","Obs3","Obs4","Sum")
    row.names(wtable) = c("X","Y","RankX","RankY")
    
    output$wilcoxpvalue <- renderText({
      paste("The Wilcoxon p-value is ",round(pval,digits=3))
    })
    
    output$wilcoxtable <- renderTable({
      wtable
    },rownames=TRUE)
    
    output$wilcoxperms <- renderText({
      paste(wsum," of ",nperm ," wilcoxon rank permutations are ", id, "than
        
                the observed rank sum of X.")
    })
    output$wilcoxplot <- renderPlot ({
      
      x<-rnorm(100000,mean=0, sd=1)
      hist(x,freq=FALSE,col="azure2")
      curve(dnorm(x), col='chartreuse4', add=TRUE)
      abline(v = qnorm(1 - (pval/2)), col = "blue4")
      abline(v = qnorm(1 - .05/2), col = "chocolate2")
    })
    
  })
  
}

shinyApp(ui = ui, server = server)