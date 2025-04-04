library(shiny)
library(ggplot2)
library(dplyr)
library(readr)

# Define UI for application
ui <- fluidPage(
  titlePanel("Visualisation of Error, Average Standard Deviation"),
  sidebarLayout(
    sidebarPanel(width = 2,
                 numericInput("beta_free", "Select Beta:", min = 0, max = 1, value = 1, step = 0.1),
                 selectInput("omega_var", "Select Omega Variance:",
                             choices = c(0.1^2, 0.2^2, 1, 4^2, 5^2),
                             selected = 1),
                 sliderInput("rho", "Select Rho:", min = -0.9, max = 0.9, step = 0.3, value = 0),
                 sliderInput("Tfull", "Select Tfull:", min = 2, max = 4, step = 1, value = 2),
                 sliderInput("n_reps", "Select Number of Repetitions:", min = 500, max = 1000, step = 500, value = 500),
                 selectInput("y_scale_fix", "Fix Y-Scale:",
                             choices = c("None", "Compare Time", "Compare Omega", "Compare Rho"),
                             selected = "None"),
                 radioButtons("plot_type", "Select Plot Type:",
                              choices = c("Error Metrics" = "error", "Average SD" = "sd"),
                              selected = "error"),
                 checkboxInput("show_beta_ci", "Show Beta CI", value = TRUE),
                 checkboxInput("show_omega_ci", "Show Omega CI", value = TRUE),
                 checkboxInput("show_rho_ci", "Show Rho CI", value = TRUE),
                 actionButton("plot_button", "Create Plot"),
                 downloadButton("download_plot", "Download Plot")
    ),
    mainPanel(width = 10,
              plotOutput("error_plot")
    ),
    position = "right"
  )
)

# server <- function(input, output) {
#   
#   dataset <- reactive({
#     tryCatch({
#       read_csv(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
#     }, error = function(e) {
#       showNotification(paste("Error reading CSV:", e$message), type = "error")
#       NULL
#     })
#   })
#   
#   filtered_data <- reactive({
#     req(input$plot_button)
#     req(dataset())
#     
#     df <- dataset()
#     
#     if (input$plot_type %in% c("error", "sd")) {
#       df <- df %>% filter(N != 20) # Always filter N = 20 for these plot types
#     }
#     
#     df <- df %>%
#       filter(rho == as.numeric(input$rho),
#              Tfull == as.numeric(input$Tfull),
#              n_reps == as.numeric(input$n_reps),
#              beta_free == as.numeric(input$beta_free),
#              omega_var == as.numeric(input$omega_var))
#     
#     if (input$plot_type == "error") {
#       df <- df %>%
#         mutate(
#           beta_lower = beta_error - 1.96 * beta_estimate_sd,
#           beta_upper = beta_error + 1.96 * beta_estimate_sd,
#           omega_lower = omega_error - 1.96 * omega_estimate_sd,
#           omega_upper = omega_error + 1.96 * omega_estimate_sd,
#           rho_lower = rho_error - 1.96 * rho_estimate_sd,
#           rho_upper = rho_error + 1.96 * rho_estimate_sd
#         )
#     }
#     
#     df
#   })
#   
#   output$error_plot <- renderPlot({
#     req(filtered_data())
#     req(nrow(filtered_data()) > 0)
#     
#     plot_data <- filtered_data()
#     
#     title_text <- ifelse(input$plot_type == "error", "Error Metrics vs. Sample Size (N)", "Average SD vs. Sample Size (N)")
#     legend_title <- ifelse(input$plot_type == "error", "Error Metric", "SD Metric")
#     
#     p <- ggplot(plot_data, aes(x = N)) +
#       {
#         if (input$plot_type == "error") {
#           ci_ribbons <- list()
#           if (input$show_beta_ci) {
#             ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper, fill = "Beta CI"), alpha = 0.2)))
#           }
#           if (input$show_omega_ci) {
#             ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = "Omega CI"), alpha = 0.2)))
#           }
#           if (input$show_rho_ci) {
#             ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = rho_lower, ymax = rho_upper, fill = "Rho CI"), alpha = 0.2)))
#           }
#           
#           list(
#             ci_ribbons,
#             geom_line(aes(y = beta_error, color = "Beta Error")),
#             geom_line(aes(y = omega_error, color = "Omega Error")),
#             geom_line(aes(y = rho_error, color = "Rho Error")),
#             geom_point(aes(y = beta_error, color = "Beta Error"), size = 1),
#             geom_point(aes(y = omega_error, color = "Omega Error"), size = 1),
#             geom_point(aes(y = rho_error, color = "Rho Error"), size = 1)
#           )
#         } else {
#           list(
#             geom_line(aes(y = beta_sd_avg, color = "Beta SD Avg")),
#             geom_line(aes(y = omega_sd_avg, color = "Omega SD Avg")),
#             geom_line(aes(y = rho_sd_avg, color = "Rho SD Avg")),
#             geom_point(aes(y = beta_sd_avg, color = "Beta SD Avg"), size = 1),
#             geom_point(aes(y = omega_sd_avg, color = "Omega SD Avg"), size = 1),
#             geom_point(aes(y = rho_sd_avg, color = "Rho SD Avg"), size = 1)
#           )
#         }
#       } +
#       labs(title = title_text,
#            x = "Sample Size (N)",
#            y = "Value",
#            color = legend_title,
#            fill = "Confidence Interval") +
#       theme_minimal() +
#       theme(legend.position = c(0.85, 0.85)) +
#       geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
#       scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50))
#     
#     print(p)
#   })
#   
#   output$download_plot <- downloadHandler(
#     filename = function() {
#       paste0("error_plot_", Sys.Date(), ".png")
#     },
#     content = function(file) {
#       req(filtered_data())
#       png(file, width = 800, height = 600, units = "px")
#       print(output$error_plot())
#       dev.off()
#     }
#   )
# }

# shinyApp(ui = ui, server = server)

# With y-fix in visual ----------------------------------------------------

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
    
    if (input$plot_type %in% c("error", "sd")) {
      df <- df %>% filter(N != 20) # Always filter N = 20 for these plot types
    }
    
    df <- df %>%
      filter(rho == as.numeric(input$rho),
             Tfull == as.numeric(input$Tfull),
             n_reps == as.numeric(input$n_reps),
             beta_free == as.numeric(input$beta_free),
             omega_var == as.numeric(input$omega_var))
    
    if (input$plot_type == "error") {
      df <- df %>%
        mutate(
          beta_lower = beta_error - 1.96 * beta_estimate_sd,
          beta_upper = beta_error + 1.96 * beta_estimate_sd,
          omega_lower = omega_error - 1.96 * omega_estimate_sd,
          omega_upper = omega_error + 1.96 * omega_estimate_sd,
          rho_lower = rho_error - 1.96 * rho_estimate_sd,
          rho_upper = rho_error + 1.96 * rho_estimate_sd
        )
    }
    
    df
  })
  
  output$error_plot <- renderPlot({
    req(filtered_data())
    req(nrow(filtered_data()) > 0)
    
    plot_data <- filtered_data()
    
    title_text <- ifelse(input$plot_type == "error", "Error Metrics vs. Sample Size (N)", "Average SD vs. Sample Size (N)")
    legend_title <- ifelse(input$plot_type == "error", "Error Metric", "SD Metric")
    
    # For y-axis scale fix
    # foo <- dataset()
    # foo <- foo %>% filter(beta_free == input$beta_free,
    #                       rho == input$rho,
    #                       Tfull == input$Tfull,
    #                       n_reps == input$n_reps) %>% filter(!(N == 20))
    
    p <- ggplot(plot_data, aes(x = N)) +
      {
        if (input$plot_type == "error") {
          ci_ribbons <- list()
          if (input$show_beta_ci) {
            ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper, fill = "Beta CI"), alpha = 0.3)))
          }
          if (input$show_omega_ci) {
            ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = "Omega CI"), alpha = 0.3)))
          }
          if (input$show_rho_ci) {
            ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = rho_lower, ymax = rho_upper, fill = "Rho CI"), alpha = 0.3)))
          }
          
          list(
            ci_ribbons,
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
           color = legend_title,
           fill = "Confidence Interval") +
      theme_minimal() +
      theme(legend.position = c(0.85, 0.85)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
      scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50)) +
      
      # **Custom Fill Colors for the Ribbons**
      scale_fill_manual(values = c("Beta CI" = "#E69F00",   # Orange
                                   "Omega CI" = "#7570B3",  # Purple
                                   "Rho CI" = "#1B9E77"))   # Teal
    
    foo <- dataset()
    
    # Conditional y-axis limit adjustment
    if (input$y_scale_fix == "Compare Omega" && input$plot_type == "error") {
      foo <- foo %>% filter(beta_free == input$beta_free,
                            rho == input$rho,
                            Tfull == input$Tfull,
                            n_reps == input$n_reps) %>% filter(!(N == 20)) %>% filter(omega_error < 5)
      p <- p + coord_cartesian(ylim = c(min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE ),
                                        max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE ))
      )
    }
    
    if (input$y_scale_fix == "Compare Time" && input$plot_type == "error") {
      # For y-axis scale fix
      foo <- foo %>% filter(beta_free == input$beta_free,
                            rho == input$rho,
                            n_reps == input$n_reps) %>% filter(!(N == 20))
      p <- p + coord_cartesian(ylim = c(min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE ),
                                        max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE ))
      )
    }
    
    if (input$y_scale_fix == "Compare Rho" && input$plot_type == "error") {
      # For y-axis scale fix
      foo <- foo %>% filter(beta_free == input$beta_free,
                            Tfull == input$Tfull,
                            n_reps == input$n_reps) %>% filter(!(N == 20))
      p <- p + coord_cartesian(ylim = c(min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE ),
                                        max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE ))
      )
    }
    
    if (input$y_scale_fix == "Compare Omega" && input$plot_type == "sd") {
      foo <- foo %>% filter(beta_free == input$beta_free,
                            rho == input$rho,
                            Tfull == input$Tfull,
                            n_reps == input$n_reps) %>% filter(!(N == 20))
      p <- p + coord_cartesian(ylim = c(min( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE ),
                                        max( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE ))
      )
    }
    
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




# Extra shiny app for error plots -----------------------------------------
## Define UI for application -------
ui <- fluidPage(
  titlePanel("Visualisation of Estimate Error (Bias), Avg. Estimate SD, and %-Rate of Estimates within 95%-CI"),
  sidebarLayout(
    sidebarPanel(width = 2,
                 selectInput("n_reps", "Select # of reps:",
                             choices = c(500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000),
                             selected = 1),
                 selectInput("omega_var", "Select Omega Variance:",
                             choices = c(0.25^2, 1, 4^2),
                             selected = 1),
                 sliderInput("rho", "Select Rho:", min = -0.9, max = 0.9, step = 0.3, value = 0),
                 sliderInput("Tfull", "Select Tfull:", min = 2, max = 4, step = 1, value = 2),
                 selectInput("y_scale_fix", "Fix Y-Scale for:",
                             choices = c("None", "Compare Time", "Compare Omega", "Compare Rho", "Compare n_reps"),
                             selected = "None"),
                 checkboxInput("show_beta_ci", "Show Beta CI", value = TRUE),
                 checkboxInput("show_omega_ci", "Show Omega CI", value = TRUE),
                 checkboxInput("show_rho_ci", "Show Rho CI", value = TRUE),
                 actionButton("plot_button", "Create Plot"),
                 downloadButton("download_plot", "Download Plot")
    ),
    mainPanel(width = 10,
              plotOutput("error_plot")
    ),
    position = "right"
  )
)

server <- function(input, output) {
  
  dataset <- reactive({
    tryCatch({
      read_csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results.csv"))
    }, error = function(e) {
      showNotification(paste("Error reading CSV:", e$message), type = "error")
      NULL
    })
  })
  
  filtered_data <- reactive({
    req(input$plot_button)
    
    df <- dataset()
    
    df <- df %>%
      filter(rho == as.numeric(input$rho),
             Tfull == as.numeric(input$Tfull),
             n_reps == as.numeric(input$n_reps),
             omega_var == as.numeric(input$omega_var))
    
    df <- df %>% mutate(beta_lower = beta_error - 1.96 * beta_estimate_sd,
                        beta_upper = beta_error + 1.96 * beta_estimate_sd,
                        omega_lower = omega_error - 1.96 * omega_estimate_sd,
                        omega_upper = omega_error + 1.96 * omega_estimate_sd,
                        rho_lower = rho_error - 1.96 * rho_estimate_sd,
                        rho_upper = rho_error + 1.96 * rho_estimate_sd)
    df
  })
  
  output$error_plot <- renderPlot({
    
    # Requirements and load-in
    req(filtered_data())
    req(input$plot_button)
    
    df <- filtered_data()
    p <- ggplot(df, aes(x = N)) + list(
      geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper, fill = "Beta CI"), alpha = 0.3),
      geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = "Omega CI"), alpha = 0.3),
      geom_ribbon(aes(ymin = rho_lower, ymax = rho_upper, fill = "Rho CI"), alpha = 0.3),
      geom_line(aes(y = beta_error, color = "Beta Error")),
      geom_line(aes(y = omega_error, color = "Omega Error")),
      geom_line(aes(y = rho_error, color = "Rho Error")),
      geom_point(aes(y = beta_error, color = "Beta Error"), size = 1),
      geom_point(aes(y = omega_error, color = "Omega Error"), size = 1),
      geom_point(aes(y = rho_error, color = "Rho Error"), size = 1)) +
      labs(title = "title_text",
           x = "Sample Size (N)",
           y = "Value",
           fill = "Confidence Interval") +
      theme_minimal() +
      theme(legend.position = c(0.85, 0.85)) +
      scale_fill_manual(values = c("Beta CI" = "#E69F00",   # Orange
                                   "Omega CI" = "#009900",  # Purple
                                   "Rho CI" = "#0066ff"))   # Teal
    
    if(input$y_scale_fix == "None"){
      foo <- filtered_data()
      foo <- foo %>% filter(beta_error < 10, beta_error > -10, omega_error < 10,
                            omega_error > -10, rho_error < 10, rho_error > -10)
      min = min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
      max = max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
      p = p + coord_cartesian(ylim = c(min, max))
    } else {
      p = p + coord_cartesian(ylim = c(0, 10)) 
    }
    
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

