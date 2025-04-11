
# Testing for appropriate n_reps sample size ------------------------------
library(shiny)
library(ggplot2)
library(dplyr)
library(readr)

## Define UI for application -------
ui<- fluidPage(
  titlePanel("Visualisation of Estimate Error (Bias), Avg. Estimate SD, and %-Rate of Estimates within 95%-CI"),
  sidebarLayout(
    sidebarPanel(width = 2,
                 radioButtons("plot_type", "Select Plot Type:",
                              choices = c("Error Metrics" = "error",
                                          "SMSE Metrics" = "smse",
                                          "Average SD" = "sd",
                                          "CI-Overlap rate" = "ci_succ"),
                              selected = "error"),
                 selectInput("beta", "Select Beta:",
                             choices = c(0, 1, 2, 3, 4),
                             selected = 1),
                 # dynamic nreps input
                 uiOutput("n_reps_ui"),
                 
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
                 downloadButton("download_plot", "Download Plot"),
                 
                 # Toggle advanced options
                 checkboxInput("advanced", "Enable Advanced nreps Options", FALSE),
                 checkboxInput("advancedView", "Enable Advanced view Options", FALSE),
                 checkboxInput("advancedSE", "Use empirical SE", FALSE)
    ),
    mainPanel(width = 10,
              plotOutput("error_plot")
    ),
    position = "right"
  )
)

server <- function(input, output, session) {
  
  output$n_reps_ui <- renderUI({
    if(input$advanced){
      sliderInput("n_reps", "Select # of reps:", 
                  min = 500, max = 3000, step = 2500, value = 500)
    } else {
      selectInput("n_reps", "Select # of reps:",
                  choices = c(100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000),
                  selected = 500)
    }
  })
  
  
  dataset <- reactive({
    tryCatch({
      read_csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results_allB_v2.csv"))
    }, error = function(e) {
      showNotification(paste("Error reading CSV:", e$message), type = "error")
      NULL
    })
  })
  
  filtered_data_for_20 <- reactive({
    req(input$plot_button)
    req(dataset())
    
    df <- dataset()
    df <- df %>% filter(N != 20)
    
    df
  })
  
  filtered_data <- reactive({
    req(input$plot_button)
    req(filtered_data_for_20())
    
    if (input$plot_type %in% c("error", "sd")) {
      df <- dataset() # filtered_data_for_20()
    } else {
      df <- dataset()
    }
    
    df <- df %>%
      filter(rho == as.numeric(input$rho),
             Tfull == as.numeric(input$Tfull),
             n_reps == as.numeric(input$n_reps),
             omega_var == as.numeric(input$omega_var),
             beta_free == as.numeric(input$beta))
    
    alpha = 0.05
    z_value = qnorm(1-alpha/2)
    
    if(input$plot_type == "error"){
      df <- df %>%
        mutate(beta_lower = beta_error - z_value * beta_estimate_sd,
               beta_upper = beta_error + z_value * beta_estimate_sd,
               omega_lower = omega_error - z_value * omega_estimate_sd,
               omega_upper = omega_error + z_value * omega_estimate_sd,
               rho_lower = rho_error - z_value * rho_estimate_sd,
               rho_upper = rho_error + z_value * rho_estimate_sd)
    }
    
    if(input$plot_type == "sd"){
      df <- df %>%
        mutate(beta_lower = beta_sd_avg - z_value * sd_of_beta_sds,
               beta_upper = beta_sd_avg + z_value * sd_of_beta_sds,
               omega_lower = omega_sd_avg - z_value * sd_of_omega_sds,
               omega_upper = omega_sd_avg + z_value * sd_of_omega_sds,
               rho_lower = rho_sd_avg - z_value * sd_of_rho_sds,
               rho_upper = rho_sd_avg + z_value * sd_of_rho_sds)
    }
    
    df
  })
  
  output$error_plot <- renderPlot({
    
    # Requirements and load-in
    req(filtered_data())
    
    plot_data <- filtered_data()
    
    ## Plot title & legend --------
    title_text <- switch(input$plot_type,
                         "error" = "Error Metrics (Bias) vs. Sample Size",
                         "smse" = "Standardised Mean Squared Error vs. Sample Size",
                         "sd" = "Average SD of Estimate vs. Sample Size",
                         "Success Rate of lying within 95%-CI vs. Sample Size")
    
    legend_title <- switch(input$plot_type,
                           "error" = "Error Metric",
                           "smse" = "Error Metric",
                           "sd" = "SD Metric",
                           "Estimates within 95%-CI")
    
    ## Create entire plot body ----------
    if(input$plot_type == "error"){
      p <- ggplot(plot_data, aes(x = N))
      p <- p + {
        ci_ribbons <- list()
        if(input$show_beta_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper, fill = "Beta CI"), alpha = 0.3)))
        }
        if(input$show_omega_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = "Omega CI"), alpha = 0.3)))
        }
        if(input$show_rho_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = rho_lower, ymax = rho_upper, fill = "Rho CI"), alpha = 0.3)))
        }
        list(ci_ribbons)
      }
      p <- p +
        list(geom_line(aes(y = beta_error, color = "Beta Error")),
             geom_line(aes(y = omega_error, color = "Omega Error")),
             geom_line(aes(y = rho_error, color = "Rho Error")),
             geom_point(aes(y = beta_error, color = "Beta Error"), size = 1),
             geom_point(aes(y = omega_error, color = "Omega Error"), size = 1),
             geom_point(aes(y = rho_error, color = "Rho Error"), size = 1)) +
        labs(title = title_text,
             x = "Sample Size (N)",
             y = "Value",
             color = legend_title,
             fill = "Confidence Interval") +
        theme_minimal() +
        theme(legend.position = c(0.85, 0.85)) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
        scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50)) + 
        coord_cartesian(ylim = c(0, 50)) +
        scale_color_manual(values = c("Beta Error" = "red",   
                                      "Omega Error" = "#006400",  
                                      "Rho Error" = "blue")) +
        scale_fill_manual(values = c("Beta CI" = "#E69F00",   
                                     "Omega CI" = "#009900",  
                                     "Rho CI" = "#0066ff"))  
    }
    
    if(input$plot_type == "smse"){
      p <- ggplot(plot_data, aes(x = N))
      if(input$advancedSE){
        p <- p +
          list(geom_line(aes(y = beta_mse_with_emp_se, color = "Beta Error")),
               geom_line(aes(y = omega_mse_with_emp_se, color = "Omega Error")),
               geom_line(aes(y = rho_mse_with_emp_se, color = "Rho Error")),
               geom_point(aes(y = beta_mse_with_emp_se, color = "Beta Error"), size = 1),
               geom_point(aes(y = omega_mse_with_emp_se, color = "Omega Error"), size = 1),
               geom_point(aes(y = rho_mse_with_emp_se, color = "Rho Error"), size = 1)) +
          labs(title = title_text,
               x = "Sample Size (N)",
               y = "Value",
               color = legend_title,
               fill = "Confidence Interval") +
          theme_minimal() +
          theme(legend.position = c(0.85, 0.85)) +
          geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
          scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50)) + 
          coord_cartesian(ylim = c(0.5, 1.5)) +
          scale_color_manual(values = c("Beta Error" = "red",   
                                        "Omega Error" = "#006400",  
                                        "Rho Error" = "blue"))
      } else {
      p <- p +
        list(geom_line(aes(y = beta_mse_with_model_se, color = "Beta Error")),
             geom_line(aes(y = omega_mse_with_model_se, color = "Omega Error")),
             geom_line(aes(y = rho_mse_with_model_se, color = "Rho Error")),
             geom_point(aes(y = beta_mse_with_model_se, color = "Beta Error"), size = 1),
             geom_point(aes(y = omega_mse_with_model_se, color = "Omega Error"), size = 1),
             geom_point(aes(y = rho_mse_with_model_se, color = "Rho Error"), size = 1)) +
        labs(title = title_text,
             x = "Sample Size (N)",
             y = "Value",
             color = legend_title,
             fill = "Confidence Interval") +
        theme_minimal() +
        theme(legend.position = c(0.85, 0.85)) +
        geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
        scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50)) + 
        coord_cartesian(ylim = c(0.5, 1.5)) +
        scale_color_manual(values = c("Beta Error" = "red",   
                                      "Omega Error" = "#006400",  
                                      "Rho Error" = "blue"))
      }
    }
    
    if(input$plot_type == "sd"){
      p <- ggplot(plot_data, aes(x = N))
      p <- p + {
        ci_ribbons <- list()
        if(input$show_beta_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper, fill = "Beta CI"), alpha = 0.3)))
        }
        if(input$show_omega_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = "Omega CI"), alpha = 0.3)))
        }
        if(input$show_rho_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = rho_lower, ymax = rho_upper, fill = "Rho CI"), alpha = 0.3)))
        }
        list(ci_ribbons)
      }
      p <- p + 
        list(geom_line(aes(y = beta_sd_avg, color = "Beta SD Avg")),
             geom_line(aes(y = omega_sd_avg, color = "Omega SD Avg")),
             geom_line(aes(y = rho_sd_avg, color = "Rho SD Avg")),
             geom_point(aes(y = beta_sd_avg, color = "Beta SD Avg"), size = 1),
             geom_point(aes(y = omega_sd_avg, color = "Omega SD Avg"), size = 1),
             geom_point(aes(y = rho_sd_avg, color = "Rho SD Avg"), size = 1)) +
        labs(title = title_text,
             x = "Sample Size (N)",
             y = "Value",
             color = legend_title,
             fill = "Confidence Interval") +
        theme_minimal() +
        theme(legend.position = c(0.85, 0.85)) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
        scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50)) +
        scale_color_manual(values = c("Beta SD Avg" = "red",   
                                      "Omega SD Avg" = "#006400",  
                                      "Rho SD Avg" = "blue")) +
        scale_fill_manual(values = c("Beta CI" = "#E69F00",   # Orange
                                     "Omega CI" = "#009900",  # Purple
                                     "Rho CI" = "#0066ff"))   # Teal
    }
    
    if(input$plot_type == "ci_succ"){
      p <- ggplot(plot_data, aes(x = N))
      p <- p + {
        ci_ribbons <- list()
        if(input$show_beta_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = bino_ci_exact_lower_beta, ymax = bino_ci_exact_upper_beta, fill = "Beta CI"), alpha = 0.3)))
        }
        if(input$show_omega_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = bino_ci_exact_lower_omega, ymax =  bino_ci_exact_upper_omega, fill = "Omega CI"), alpha = 0.3)))
        }
        if(input$show_rho_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = bino_ci_exact_lower_rho, ymax = bino_ci_exact_upper_rho, fill = "Rho CI"), alpha = 0.3)))
        }
        list(ci_ribbons)
      }
      p <- p +
        list(geom_line(aes(y = p_beta, color = "Beta")),
             geom_line(aes(y = p_omega, color = "Omega")),
             geom_line(aes(y = p_rho, color = "Rho")),
             geom_point(aes(y = p_beta, color = "Beta"), size = 1),
             geom_point(aes(y = p_omega, color = "Omega"), size = 1),
             geom_point(aes(y = p_rho, color = "Rho"), size = 1)) +
        labs(title = title_text,
             x = "Sample Size (N)",
             y = "Value",
             color = legend_title,
             fill = "Confidence Interval") +
        theme_minimal() +
        theme(legend.position = c(0.85, 0.25)) +
        geom_hline(yintercept = 0.95, linetype = "dotted", color = "black") +
        scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50)) + 
        scale_color_manual(values = c("Beta" = "red",   
                                      "Omega" = "#006400",  
                                      "Rho" = "blue")) +
        scale_fill_manual(values = c("Beta CI" = "#E69F00",   # Orange
                                     "Omega CI" = "#009900",  # Purple
                                     "Rho CI" = "#0066ff"))   # Teal
    }
    
    ## The y-axis limit adjustment for better comparison of values/info ------
    ## 1.) Extract the correct info for limiting the scale (so what is it based on)
    if(input$y_scale_fix == "None"){
      foo <- filtered_data()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_error < 10, beta_error > -10, omega_error < 10,
                              omega_error > -10, rho_error < 10, rho_error > -10) %>% 
          filter(N != 20)
        min = min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        max = max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(beta_sd_avg < 10, omega_sd_avg < 10, rho_sd_avg < 10) %>% 
          filter(N != 20)
        min = min( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        max = max( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "ci_succ"){
        min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
    }
    
    if(input$y_scale_fix == "Compare Time"){
      foo <- dataset()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_error < 10, beta_error > -10, omega_error < 10,
                              omega_error > -10, rho_error < 10, rho_error > -10) %>% 
          filter(rho == input$rho,omega_var == input$omega_var,n_reps == input$n_reps) %>% 
          filter(N != 20)
        min = min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        max = max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(beta_sd_avg < 10, omega_sd_avg < 10, rho_sd_avg < 10) %>% filter(rho == input$rho,
                                                  omega_var == input$omega_var,
                                                  n_reps == input$n_reps)  %>% 
          filter(N != 20)
        min = min( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        max = max( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        
        
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "ci_succ"){
         
        foo <- foo %>% filter(rho == input$rho,
                              omega_var == input$omega_var,
                              n_reps == input$n_reps)
        min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
    }
    
    if(input$y_scale_fix == "Compare Omega"){
      foo <- dataset()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_error < 10, beta_error > -10, omega_error < 10,
                              omega_error > -10, rho_error < 10, rho_error > -10) %>% filter(rho == input$rho,
                                                  Tfull == input$Tfull,
                                                  n_reps == input$n_reps) %>% 
          filter(N != 20)
        min = min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        max = max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(beta_sd_avg < 10, omega_sd_avg < 10, rho_sd_avg < 10) %>% filter(rho == input$rho,
                                                  Tfull == input$Tfull,
                                                  n_reps == input$n_reps)  %>% 
          filter(N != 20)
        min = min( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        max = max( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      } 
      if(input$plot_type == "ci_succ"){
        
        foo <- foo %>% filter(rho == input$rho,
                              Tfull == input$Tfull,
                              n_reps == input$n_reps)
        min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
    }
    
    if(input$y_scale_fix == "Compare Rho"){
      foo <- dataset()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_error < 10, beta_error > -10, omega_error < 10,
                              omega_error > -10, rho_error < 10, rho_error > -10) %>% filter(omega_var == input$omega_var,
                                                  Tfull == input$Tfull,
                                                  n_reps == input$n_reps) %>% 
          filter(N != 20)
        min = min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        max = max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(beta_sd_avg < 10, omega_sd_avg < 10, rho_sd_avg < 10) %>% filter(omega_var == input$omega_var,
                                                  Tfull == input$Tfull,
                                                  n_reps == input$n_reps)  %>% 
          filter(N != 20)
        min = min( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        max = max( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "ci_succ"){
        
        foo <- foo %>% filter(omega_var == input$omega_var,
                              Tfull == input$Tfull,
                              n_reps == input$n_reps)
        min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
    }
    
    if(input$y_scale_fix == "Compare n_reps"){
      foo <- dataset()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_error < 10, beta_error > -10, omega_error < 10,
                              omega_error > -10, rho_error < 10, rho_error > -10) %>% filter(omega_var == input$omega_var,
                                                  Tfull == input$Tfull,
                                                  rho == input$rho) %>% 
          filter(N != 20)
        min = min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        max = max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(beta_sd_avg < 10, omega_sd_avg < 10, rho_sd_avg < 10) %>% filter(omega_var == input$omega_var,
                                                  Tfull == input$Tfull,
                                                  rho == input$rho)  %>% 
          filter(N != 20)
        min = min( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        max = max( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "ci_succ"){
        
        foo <- foo %>% filter(omega_var == input$omega_var,
                              Tfull == input$Tfull,
                              rho == input$rho)
        min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
    }
    
    if(input$advancedView){
      p <- p + coord_cartesian(ylim = c(0, 30))
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




# Comparing all inbetween beta's ----------------------------------------------------------
library(shiny)
library(ggplot2)
library(dplyr)
library(readr)

## Define UI for application -------
ui<- fluidPage(
  titlePanel("Visualisation of Estimate Error (Bias), Avg. Estimate SD, and %-Rate of Estimates within 95%-CI"),
  sidebarLayout(
    sidebarPanel(width = 2,
                 radioButtons("plot_type", "Select Plot Type:",
                              choices = c("Error Metrics" = "error",
                                          "Average SD" = "sd",
                                          "CI-Overlap rate" = "ci_succ"),
                              selected = "error"),
                 sliderInput("beta", "Select Beta:", min = 0, max = 4, step = 0.1, value = 1),
                 
                 selectInput("omega_var", "Select Omega Variance:",
                             choices = c(0.25^2, 1, 4^2),
                             selected = 1),
                 sliderInput("rho", "Select Rho:", min = -0.9, max = 0.9, step = 0.3, value = 0),
                 sliderInput("Tfull", "Select Tfull:", min = 2, max = 4, step = 1, value = 2),
                 selectInput("y_scale_fix", "Fix Y-Scale for:",
                             choices = c("None", "Compare Time", "Compare Omega", "Compare Rho"),
                             selected = "None"),
                 checkboxInput("show_beta_ci", "Show Beta CI", value = TRUE),
                 checkboxInput("show_omega_ci", "Show Omega CI", value = TRUE),
                 checkboxInput("show_rho_ci", "Show Rho CI", value = TRUE),
                 actionButton("plot_button", "Create Plot"),
                 downloadButton("download_plot", "Download Plot"),
                 
                 # Toggle advanced options
                 checkboxInput("advancedView", "Enable Advanced view Options", FALSE),
    ),
    mainPanel(width = 10,
              plotOutput("error_plot")
    ),
    position = "right"
  )
)

server <- function(input, output, session) {
  
  dataset <- reactive({
    tryCatch({
      read_csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_500reps_Results.csv"))
    }, error = function(e) {
      showNotification(paste("Error reading CSV:", e$message), type = "error")
      NULL
    })
  })
  
  filtered_data_for_20 <- reactive({
    req(input$plot_button)
    req(dataset())
    
    df <- dataset()
    df <- df %>% filter(N != 20)
    
    df
  })
  
  filtered_data <- reactive({
    req(input$plot_button)
    req(filtered_data_for_20())
    
    if (input$plot_type %in% c("error", "sd")) {
      df <- dataset() # filtered_data_for_20()
    } else {
      df <- dataset()
    }
    
    df <- df %>%
      filter(rho == as.numeric(input$rho),
             Tfull == as.numeric(input$Tfull),
             omega_var == as.numeric(input$omega_var),
             beta_free == as.numeric(input$beta))
    
    alpha = 0.05
    z_value = qnorm(1-alpha/2)
    
    if(input$plot_type == "error"){
      df <- df %>%
        mutate(beta_lower = beta_error - z_value * beta_estimate_sd,
               beta_upper = beta_error + z_value * beta_estimate_sd,
               omega_lower = omega_error - z_value * omega_estimate_sd,
               omega_upper = omega_error + z_value * omega_estimate_sd,
               rho_lower = rho_error - z_value * rho_estimate_sd,
               rho_upper = rho_error + z_value * rho_estimate_sd)
    }
    
    if(input$plot_type == "sd"){
      df <- df %>%
        mutate(beta_lower = beta_sd_avg - z_value * sd_of_beta_sds,
               beta_upper = beta_sd_avg + z_value * sd_of_beta_sds,
               omega_lower = omega_sd_avg - z_value * sd_of_omega_sds,
               omega_upper = omega_sd_avg + z_value * sd_of_omega_sds,
               rho_lower = rho_sd_avg - z_value * sd_of_rho_sds,
               rho_upper = rho_sd_avg + z_value * sd_of_rho_sds)
    }
    
    df
  })
  
  output$error_plot <- renderPlot({
    
    # Requirements and load-in
    req(filtered_data())
    
    plot_data <- filtered_data()
    
    ## Plot title & legend --------
    title_text <- switch(input$plot_type,
                         "error" = "Error Metrics (Bias) vs. Sample Size",
                         "sd" = "Average SD of Estimate vs. Sample Size",
                         "Success Rate of lying within 95%-CI vs. Sample Size")
    
    legend_title <- switch(input$plot_type,
                           "error" = "Error Metric",
                           "sd" = "SD Metric",
                           "Estimates within 95%-CI")
    
    ## Create entire plot body ----------
    if(input$plot_type == "error"){
      p <- ggplot(plot_data, aes(x = N))
      p <- p + {
        ci_ribbons <- list()
        if(input$show_beta_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper, fill = "Beta CI"), alpha = 0.3)))
        }
        if(input$show_omega_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = "Omega CI"), alpha = 0.3)))
        }
        if(input$show_rho_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = rho_lower, ymax = rho_upper, fill = "Rho CI"), alpha = 0.3)))
        }
        list(ci_ribbons)
      }
      p <- p +
        list(geom_line(aes(y = beta_error, color = "Beta Error")),
             geom_line(aes(y = omega_error, color = "Omega Error")),
             geom_line(aes(y = rho_error, color = "Rho Error")),
             geom_point(aes(y = beta_error, color = "Beta Error"), size = 1),
             geom_point(aes(y = omega_error, color = "Omega Error"), size = 1),
             geom_point(aes(y = rho_error, color = "Rho Error"), size = 1)) +
        labs(title = title_text,
             x = "Sample Size (N)",
             y = "Value",
             color = legend_title,
             fill = "Confidence Interval") +
        theme_minimal() +
        theme(legend.position = c(0.85, 0.85)) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
        scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50)) + 
        coord_cartesian(ylim = c(0, 50)) +
        scale_color_manual(values = c("Beta Error" = "red",   
                                      "Omega Error" = "#006400",  
                                      "Rho Error" = "blue")) +
        scale_fill_manual(values = c("Beta CI" = "#E69F00",   
                                     "Omega CI" = "#009900",  
                                     "Rho CI" = "#0066ff"))  
    }
    
    if(input$plot_type == "sd"){
      p <- ggplot(plot_data, aes(x = N))
      p <- p + {
        ci_ribbons <- list()
        if(input$show_beta_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper, fill = "Beta CI"), alpha = 0.3)))
        }
        if(input$show_omega_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = "Omega CI"), alpha = 0.3)))
        }
        if(input$show_rho_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = rho_lower, ymax = rho_upper, fill = "Rho CI"), alpha = 0.3)))
        }
        list(ci_ribbons)
      }
      p <- p + 
        list(geom_line(aes(y = beta_sd_avg, color = "Beta SD Avg")),
             geom_line(aes(y = omega_sd_avg, color = "Omega SD Avg")),
             geom_line(aes(y = rho_sd_avg, color = "Rho SD Avg")),
             geom_point(aes(y = beta_sd_avg, color = "Beta SD Avg"), size = 1),
             geom_point(aes(y = omega_sd_avg, color = "Omega SD Avg"), size = 1),
             geom_point(aes(y = rho_sd_avg, color = "Rho SD Avg"), size = 1)) +
        labs(title = title_text,
             x = "Sample Size (N)",
             y = "Value",
             color = legend_title,
             fill = "Confidence Interval") +
        theme_minimal() +
        theme(legend.position = c(0.85, 0.85)) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
        scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50)) +
        scale_color_manual(values = c("Beta SD Avg" = "red",   
                                      "Omega SD Avg" = "#006400",  
                                      "Rho SD Avg" = "blue")) +
        scale_fill_manual(values = c("Beta CI" = "#E69F00",   # Orange
                                     "Omega CI" = "#009900",  # Purple
                                     "Rho CI" = "#0066ff"))   # Teal
    }
    
    if(input$plot_type == "ci_succ"){
      p <- ggplot(plot_data, aes(x = N))
      p <- p + {
        ci_ribbons <- list()
        if(input$show_beta_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = bino_ci_exact_lower_beta, ymax = bino_ci_exact_upper_beta, fill = "Beta CI"), alpha = 0.3)))
        }
        if(input$show_omega_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = bino_ci_exact_lower_omega, ymax =  bino_ci_exact_upper_omega, fill = "Omega CI"), alpha = 0.3)))
        }
        if(input$show_rho_ci){
          ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = bino_ci_exact_lower_rho, ymax = bino_ci_exact_upper_rho, fill = "Rho CI"), alpha = 0.3)))
        }
        list(ci_ribbons)
      }
      p <- p +
        list(geom_line(aes(y = p_beta, color = "Beta")),
             geom_line(aes(y = p_omega, color = "Omega")),
             geom_line(aes(y = p_rho, color = "Rho")),
             geom_point(aes(y = p_beta, color = "Beta"), size = 1),
             geom_point(aes(y = p_omega, color = "Omega"), size = 1),
             geom_point(aes(y = p_rho, color = "Rho"), size = 1)) +
        labs(title = title_text,
             x = "Sample Size (N)",
             y = "Value",
             color = legend_title,
             fill = "Confidence Interval") +
        theme_minimal() +
        theme(legend.position = c(0.85, 0.25)) +
        geom_hline(yintercept = 0.95, linetype = "dotted", color = "black") +
        scale_x_continuous(breaks = seq(0, max(plot_data$N, na.rm = TRUE), by = 50)) + 
        scale_color_manual(values = c("Beta" = "red",   
                                      "Omega" = "#006400",  
                                      "Rho" = "blue")) +
        scale_fill_manual(values = c("Beta CI" = "#E69F00",   # Orange
                                     "Omega CI" = "#009900",  # Purple
                                     "Rho CI" = "#0066ff"))   # Teal
    }
    
    ## The y-axis limit adjustment for better comparison of values/info ------
    ## 1.) Extract the correct info for limiting the scale (so what is it based on)
    if(input$y_scale_fix == "None"){
      foo <- filtered_data()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_error < 10, beta_error > -10, omega_error < 10,
                              omega_error > -10, rho_error < 10, rho_error > -10) %>% 
          filter(N != 20)
        min = min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        max = max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(beta_sd_avg < 10, omega_sd_avg < 10, rho_sd_avg < 10) %>% 
          filter(N != 20)
        min = min( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        max = max( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "ci_succ"){
        min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
    }
    
    if(input$y_scale_fix == "Compare Time"){
      foo <- dataset()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_error < 10, beta_error > -10, omega_error < 10,
                              omega_error > -10, rho_error < 10, rho_error > -10) %>% 
          filter(rho == input$rho,omega_var == input$omega_var) %>% 
          filter(N != 20)
        min = min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        max = max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(beta_sd_avg < 10, omega_sd_avg < 10, rho_sd_avg < 10) %>% filter(rho == input$rho,
                                                                                               omega_var == input$omega_var)  %>% 
          filter(N != 20)
        min = min( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        max = max( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        
        
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "ci_succ"){
        
        foo <- foo %>% filter(rho == input$rho,
                              omega_var == input$omega_var)
        min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
    }
    
    if(input$y_scale_fix == "Compare Omega"){
      foo <- dataset()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_error < 10, beta_error > -10, omega_error < 10,
                              omega_error > -10, rho_error < 10, rho_error > -10) %>% filter(rho == input$rho,
                                                                                             Tfull == input$Tfull) %>% 
          filter(N != 20)
        min = min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        max = max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(beta_sd_avg < 10, omega_sd_avg < 10, rho_sd_avg < 10) %>% filter(rho == input$rho,
                                                                                               Tfull == input$Tfull)  %>% 
          filter(N != 20)
        min = min( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        max = max( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      } 
      if(input$plot_type == "ci_succ"){
        
        foo <- foo %>% filter(rho == input$rho,
                              Tfull == input$Tfull)
        min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
    }
    
    if(input$y_scale_fix == "Compare Rho"){
      foo <- dataset()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_error < 10, beta_error > -10, omega_error < 10,
                              omega_error > -10, rho_error < 10, rho_error > -10) %>% filter(omega_var == input$omega_var,
                                                                                             Tfull == input$Tfull) %>% 
          filter(N != 20)
        min = min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        max = max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(beta_sd_avg < 10, omega_sd_avg < 10, rho_sd_avg < 10) %>% filter(omega_var == input$omega_var,
                                                                                               Tfull == input$Tfull)  %>% 
          filter(N != 20)
        min = min( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        max = max( c(foo$beta_sd_avg, foo$omega_sd_avg, foo$rho_sd_avg), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "ci_succ"){
        
        foo <- foo %>% filter(omega_var == input$omega_var,
                              Tfull == input$Tfull)
        min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
    }
    
    
    if(input$advancedView){
      p <- p + coord_cartesian(ylim = c(0, 30))
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


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Standardised Mean Squared Error Plot ------------------------------------
library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Model Parameter vs. SMSE"),
  
  sidebarLayout(
    sidebarPanel(
      width = 2,
      selectInput("plot_type", "Select Plot Type:",
                  choices = c("beta", "rho", "omega", "T")),
      # Dynamic UI for the selected plot type*s corresponding variable input
      uiOutput("filter_ui"),
      checkboxInput("advancedSEopt", "Switch to Empirical SE", value = FALSE)
    ),
    mainPanel(width = 10,
      plotOutput("main_plot")
    ),
    position = "right"
  )
)


server <- function(input, output, session) {
  
  # Load dataset
  dataset <- reactive({
  req(input$plot_type)
  df <- tryCatch({
    read_csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results_allB_v2.csv"))
  }, error = function(e) {
    showNotification(paste("Error reading CSV:", e$message), type = "error")
    NULL
  })
  
  df
})
  
  # Custom allowed values for each parameter
  param_values <- list(
    beta = c(0, 1, 2, 3, 4),
    rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    omega = c(0.25^2, 1, 4^2),
    T = c(2, 3, 4)
  )
  
  # Mapping from UI param names to dataset column names
  param_mapping <- list(
    beta = "beta_free",
    omega = "omega_var",
    rho = "rho",
    T = "Tfull"
  )
  
  # Dynamic UI creation
  output$filter_ui <- renderUI({
    req(input$plot_type)
    other_params <- setdiff(names(param_values), input$plot_type)
    
    ui_list <- lapply(other_params, function(param) {
      vals <- param_values[[param]]
      
      # If the parameter is "omega", use a selectInput instead of a sliderInput
      if (param == "omega") {
        selectInput(
          inputId = paste0("input_", param),
          label = paste("Select", param),
          choices = vals,
          selected = vals[1]
        )
      } else {
        # For other parameters, use sliderInput
        step_size <- if (length(unique(diff(vals))) == 1) {
          unique(diff(vals))
        } else {
          min(diff(sort(unique(vals))))
        }
        
        sliderInput(
          inputId = paste0("input_", param),
          label = paste("Select", param),
          min = min(vals),
          max = max(vals),
          value = vals[1],
          step = step_size,
          round = TRUE
        )
      }
    })
    
    do.call(tagList, ui_list)
  })
  
  # Reactive for filtered data
  filtered_data <- reactive({
    req(input$plot_type, dataset())
    
    df <- dataset()
    plot_param <- input$plot_type
    param_names <- setdiff(names(param_values), plot_param)
    
    for (param in param_names) {
      col_name <- param_mapping[[param]] # correct var name
      value <- input[[paste0("input_", param)]]  # input value
      
      df <- df[df[[col_name]] == value, ]
    }
    
    df <- df %>% dplyr::filter(n_reps == 3000, N == 500)
    
    df
  })
  
  output$main_plot <- renderPlot({
    req(filtered_data())
    
    df <- filtered_data()
    
    # Use param_values to define x-axis order
    x_levels <- param_values[[input$plot_type]]
    
    p <- ggplot(df, aes(x = x_levels)) 
    
    if(!(input$advancedSEopt)){
      p <- p + geom_line(aes(y = beta_mse_with_emp_se, color = "Beta Error")) +
        geom_point(aes(y = beta_mse_with_emp_se, color = "Beta Error"), size = 1) +
        geom_line(aes(y = omega_mse_with_emp_se, color = "Omega Error")) +
        geom_point(aes(y = omega_mse_with_emp_se, color = "Omega Error"), size = 1) +
        geom_line(aes(y = rho_mse_with_emp_se, color = "Rho Error")) +
        geom_point(aes(y = rho_mse_with_emp_se, color = "Rho Error"), size = 1)
    } else {
      p <- p + geom_line(aes(y = beta_mse_with_model_se, color = "Beta Error")) +
        geom_point(aes(y = beta_mse_with_model_se, color = "Beta Error"), size = 1) +
        geom_line(aes(y = omega_mse_with_model_se, color = "Omega Error")) +
        geom_point(aes(y = omega_mse_with_model_se, color = "Omega Error"), size = 1) +
        geom_line(aes(y = rho_mse_with_model_se, color = "Rho Error")) +
        geom_point(aes(y = rho_mse_with_model_se, color = "Rho Error"), size = 1)
    }
     
    p <- p + labs(
        title = paste("SMSE vs.", input$plot_type),
        x = input$plot_type,
        y = "SMSE",
        color = "Error Type"
      ) +
      theme_minimal() +
      theme(
        legend.position = c(0.85, 0.85),
        plot.title = element_text(hjust = 0.5)
      ) +
      geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
      scale_x_continuous(breaks = x_levels) +
      coord_cartesian(ylim = c(0.5, 1.5)) +
      scale_color_manual(values = c(
        "Beta Error" = "red",
        "Omega Error" = "#006400",
        "Rho Error" = "blue"
      ))
    
    print(p)
  })
}


shinyApp(ui = ui, server = server)


