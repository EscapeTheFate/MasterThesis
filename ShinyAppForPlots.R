
library(shiny)
library(ggplot2)
library(dplyr)
library(readr)

# Define UI for application
ui <- fluidPage(
  titlePanel("Error Visualization - Model Estimation"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("rho", "Select Rho:", min = -0.9, max = 0.9, step = 0.1, value = 0),
      sliderInput("Tfull", "Select Tfull:", min = 2, max = 4, step = 1, value = 2),
      sliderInput("n_reps", "Select Number of Repetitions:", min = 500, max = 1000, step = 500, value = 500),
      radioButtons("plot_type", "Select Plot Type:",
                   choices = c("Error Metrics" = "error", "Average SD" = "sd"),
                   selected = "error"),
      actionButton("plot_button", "Create Plot"),
      # Add a download button
      downloadButton("download_plot", "Download Plot")
    ),
    mainPanel(
      plotOutput("error_plot")
    ),
    position = "right"
  )
)

# Define server logic
server <- function(input, output) {
  
  # Load data (adjust file path if necessary)
  dataset <- reactive({
    tryCatch({
      read_csv(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
    }, error = function(e) {
      # Handle file reading errors gracefully
      showNotification(paste("Error reading CSV:", e$message), type = "error")
      NULL # Return NULL to prevent further errors
    })
  })
  
  # Reactive expression for filtering data
  filtered_data <- reactive({
    req(input$plot_button)
    req(dataset()) # Make sure data is loaded
    
    dataset() %>%
      filter(rho == as.numeric(input$rho),
             Tfull == as.numeric(input$Tfull),
             n_reps == as.numeric(input$n_reps))
  })
  
  # Generate plot
  output$error_plot <- renderPlot({
    req(filtered_data())
    
    plot_data <- filtered_data()
    
    title_text <- ifelse(input$plot_type == "error", "Error Metrics vs. Sample Size (N)", "Average SD vs. Sample Size (N)")
    legend_title <- ifelse(input$plot_type == "error", "Error Metric", "SD Metric")
    
    p <- ggplot(plot_data, aes(x = N)) +
      {
        if (input$plot_type == "error") {
          list(
            geom_line(aes(y = beta_error, color = "Beta Error")),
            geom_line(aes(y = omega_error, color = "Omega Error")),
            geom_line(aes(y = rho_error, color = "Rho Error")),
            geom_point(aes(y = beta_error, color = "Beta Error"), size = 1),
            geom_point(aes(y = omega_error, color = "Omega Error"), size = 1),
            geom_point(aes(y = rho_error, color = "Rho Error"), size = 1)
          )
        } else {
          list(
            geom_line(aes(y = beta_sd_avg, color = "Beta SD Avg")),
            geom_line(aes(y = omega_sd_avg, color = "Omega SD Avg")),
            geom_line(aes(y = rho_sd_avg, color = "Rho SD Avg")),
            geom_point(aes(y = beta_sd_avg, color = "Beta SD Avg"), size = 1),
            geom_point(aes(y = omega_sd_avg, color = "Omega SD Avg"), size = 1),
            geom_point(aes(y = rho_sd_avg, color = "Rho SD Avg"), size = 1)
          )
        }
      } +
      labs(title = title_text,
           x = "Sample Size (N)",
           y = "Value",
           color = legend_title) +
      theme_minimal() +
      theme(legend.position = c(0.85, 0.85)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
      scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50))
    
    print(p) # Explicitly print the plot
  })
  
  # Download Handler
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("error_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(filtered_data())
      png(file, width = 800, height = 600, units = "px") # Adjust dimensions as needed
      print(output$error_plot()) # Use the already rendered plot
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)


# -------------------------------------------------------------------------

library(shiny)
library(ggplot2)
library(dplyr)
library(readr)

# Define UI for application
ui <- fluidPage(
  titlePanel("Error Visualization - Model Estimation"),
  sidebarLayout(
    sidebarPanel(width = 2, # Sidebar takes 2/12 of the width
                 sliderInput("rho", "Select Rho:", min = -0.9, max = 0.9, step = 0.1, value = 0),
                 sliderInput("Tfull", "Select Tfull:", min = 2, max = 4, step = 1, value = 2),
                 sliderInput("n_reps", "Select Number of Repetitions:", min = 500, max = 1000, step = 500, value = 500),
                 radioButtons("plot_type", "Select Plot Type:",
                              choices = c("Error Metrics" = "error", "Average SD" = "sd"),
                              selected = "error"),
                 actionButton("plot_button", "Create Plot"),
                 downloadButton("download_plot", "Download Plot")
    ),
    mainPanel(width = 10,  # Main panel takes 10/12 of the width
              plotOutput("error_plot")
    ),
    position = "right"
  )
)


shinyApp(ui = ui, server = server)


# -------------------------------------------------------------------------

library(shiny)
library(ggplot2)
library(dplyr)
library(readr)

# Define UI for application
ui <- fluidPage(
  titlePanel("Error Visualization - Model Estimation"),
  sidebarLayout(
    sidebarPanel(width = 2,
                 sliderInput("rho", "Select Rho:", min = -0.9, max = 0.9, step = 0.1, value = 0),
                 sliderInput("Tfull", "Select Tfull:", min = 2, max = 4, step = 1, value = 2),
                 sliderInput("n_reps", "Select Number of Repetitions:", min = 500, max = 1000, step = 500, value = 500),
                 radioButtons("plot_type", "Select Plot Type:",
                              choices = c("Error Metrics" = "error", "Average SD" = "sd"),
                              selected = "error"),
                 checkboxInput("filter_n20", "Filter N = 20", value = FALSE), # Toggle button
                 actionButton("plot_button", "Create Plot"),
                 downloadButton("download_plot", "Download Plot")
    ),
    mainPanel(width = 10,
              plotOutput("error_plot")
    ),
    position = "right"
  )
)

# Define server logic
server <- function(input, output) {
  
  dataset <- reactive({
    tryCatch({
      read_csv(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
    }, error = function(e) {
      showNotification(paste("Error reading CSV:", e$message), type = "error")
      NULL
    })
  })
  
  filtered_data <- reactive({
    req(input$plot_button)
    req(dataset())
    
    df <- dataset()
    
    if (input$filter_n20) {
      df <- df %>% filter(N != 20)
    }
    
    df %>%
      filter(rho == as.numeric(input$rho),
             Tfull == as.numeric(input$Tfull),
             n_reps == as.numeric(input$n_reps))
  })
  
  output$error_plot <- renderPlot({
    req(filtered_data())
    
    plot_data <- filtered_data()
    
    title_text <- ifelse(input$plot_type == "error", "Error Metrics vs. Sample Size (N)", "Average SD vs. Sample Size (N)")
    legend_title <- ifelse(input$plot_type == "error", "Error Metric", "SD Metric")
    
    p <- ggplot(plot_data, aes(x = N)) +
      {
        if (input$plot_type == "error") {
          list(
            geom_line(aes(y = beta_error, color = "Beta Error")),
            geom_line(aes(y = omega_error, color = "Omega Error")),
            geom_line(aes(y = rho_error, color = "Rho Error")),
            geom_point(aes(y = beta_error, color = "Beta Error"), size = 1),
            geom_point(aes(y = omega_error, color = "Omega Error"), size = 1),
            geom_point(aes(y = rho_error, color = "Rho Error"), size = 1)
          )
        } else {
          list(
            geom_line(aes(y = beta_sd_avg, color = "Beta SD Avg")),
            geom_line(aes(y = omega_sd_avg, color = "Omega SD Avg")),
            geom_line(aes(y = rho_sd_avg, color = "Rho SD Avg")),
            geom_point(aes(y = beta_sd_avg, color = "Beta SD Avg"), size = 1),
            geom_point(aes(y = omega_sd_avg, color = "Omega SD Avg"), size = 1),
            geom_point(aes(y = rho_sd_avg, color = "Rho SD Avg"), size = 1)
          )
        }
      } +
      labs(title = title_text,
           x = "Sample Size (N)",
           y = "Value",
           color = legend_title) +
      theme_minimal() +
      theme(legend.position = c(0.85, 0.85)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
      scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50))
    
    print(p)
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("error_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(filtered_data())
      png(file, width = 800, height = 600, units = "px")
      print(output$error_plot())
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)


# Newest (TODO: make it work) ----------------------------------------------------------------()
library(shiny)
library(ggplot2)
library(dplyr)
library(readr)

# Define UI for application
ui <- fluidPage(
  titlePanel("Error Visualization - Model Estimation"),
  sidebarLayout(
    sidebarPanel(width = 2,
                 sliderInput("rho", "Select Rho:", min = -0.9, max = 0.9, step = 0.1, value = 0),
                 sliderInput("Tfull", "Select Tfull:", min = 2, max = 4, step = 1, value = 2),
                 sliderInput("n_reps", "Select Number of Repetitions:", min = 500, max = 1000, step = 500, value = 500),
                 radioButtons("plot_type", "Select Plot Type:",
                              choices = c("Error Metrics" = "error", "Average SD" = "sd", "CI Overlap" = "ci"),
                              selected = "error"),
                 checkboxInput("filter_n20", "Filter N = 20", value = FALSE),
                 actionButton("plot_button", "Create Plot"),
                 downloadButton("download_plot", "Download Plot")
    ),
    mainPanel(width = 10,
              plotOutput("error_plot")
    ),
    position = "right"
  )
)

# Define server logic
server <- function(input, output) {
  
  error_data <- reactive({
    tryCatch({
      read_csv(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
    }, error = function(e) {
      showNotification(paste("Error reading ErrorAvgSd CSV:", e$message), type = "error")
      NULL
    })
  })
  
  ci_overlap_data <- reactive({
    tryCatch({
      read_csv(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/CI_Overlap_Results.csv"))
    }, error = function(e) {
      showNotification(paste("Error reading CI_Overlap CSV:", e$message), type = "error")
      NULL
    })
  })
  
  filtered_data <- reactive({
    req(input$plot_button)
    req(error_data())
    req(ci_overlap_data())
    
    df_error <- error_data()
    df_ci <- ci_overlap_data()
    
    if (input$filter_n20) {
      df_error <- df_error %>% filter(N != 20)
      df_ci <- df_ci %>% filter(N != 20)
    }
    
    filtered_error <- df_error %>%
      filter(rho == as.numeric(input$rho),
             Tfull == as.numeric(input$Tfull),
             n_reps == as.numeric(input$n_reps))
    
    filtered_ci <- df_ci %>%
      filter(rho == as.numeric(input$rho),
             Tfull == as.numeric(input$Tfull),
             n_reps == as.numeric(input$n_reps))
    
    list(error = filtered_error, ci = filtered_ci)
  })
  
  output$error_plot <- renderPlot({
    req(filtered_data())
    plot_data <- filtered_data()
    
    if (input$plot_type == "error" || input$plot_type == "sd") {
      df <- plot_data$error
      y_var_label <- "Value"
      hline_intercept <- 0
      title_text <- ifelse(input$plot_type == "error", "Error Metrics vs. Sample Size (N)", "Average SD vs. Sample Size (N)")
      legend_title <- ifelse(input$plot_type == "error", "Error Metric", "SD Metric")
      
      p <- ggplot(df, aes(x = N)) +
        {
          if (input$plot_type == "error") {
            list(
              geom_line(aes(y = beta_error, color = "Beta Error")),
              geom_line(aes(y = omega_error, color = "Omega Error")),
              geom_line(aes(y = rho_error, color = "Rho Error")),
              geom_point(aes(y = beta_error, color = "Beta Error"), size = 1),
              geom_point(aes(y = omega_error, color = "Omega Error"), size = 1),
              geom_point(aes(y = rho_error, color = "Rho Error"), size = 1)
            )
          } else {
            list(
              geom_line(aes(y = beta_sd_avg, color = "Beta SD Avg")),
              geom_line(aes(y = omega_sd_avg, color = "Omega SD Avg")),
              geom_line(aes(y = rho_sd_avg, color = "Rho SD Avg")),
              geom_point(aes(y = beta_sd_avg, color = "Beta SD Avg"), size = 1),
              geom_point(aes(y = omega_sd_avg, color = "Omega SD Avg"), size = 1),
              geom_point(aes(y = rho_sd_avg, color = "Rho SD Avg"), size = 1)
            )
          }
        }
    } else if (input$plot_type == "ci") {
      df <- plot_data$ci
      y_var_label <- "Success Rate"
      hline_intercept <- 0.95
      title_text <- "CI Overlap Success Rate vs. Sample Size (N)"
      legend_title <- "Estimate"
      
      p <- ggplot(df, aes(x = N)) +
        geom_line(aes(y = beta, color = "Beta")) +
        geom_line(aes(y = omega, color = "Omega")) +
        geom_line(aes(y = rho, color = "Rho")) +
        geom_point(aes(y = beta, color = "Beta"), size = 1) +
        geom_point(aes(y = omega, color = "Omega"), size = 1) +
        geom_point(aes(y = rho, color = "Rho"), size = 1)
    }
    
    p <- p +
      labs(title = title_text,
           x = "Sample Size (N)",
           y = y_var_label,
           color = legend_title) +
      theme_minimal() +
      theme(legend.position = c(0.85, 0.85)) +
      geom_hline(yintercept = hline_intercept, linetype = "dotted", color = "black") +
      scale_x_continuous(breaks = seq(0, max(df$N, na.rm = TRUE), by = 50))
    
    print(p)
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("error_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(filtered_data())
      png(file, width = 800, height = 600, units = "px")
      print(output$error_plot())
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)

