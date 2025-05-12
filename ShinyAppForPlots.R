
# Testing for appropriate n_reps sample size ------------------------------
library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)

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
                 checkboxInput("UseModelbasedSE", "True = Rescaled", T),
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
                 checkboxInput("advancedBiasCI", "Enable divide by sqrt(n_reps)", FALSE),
                 # Creates input fields if advanced option is true
                 conditionalPanel(
                   condition = "input.advancedView == true",
                   numericInput("advancedupper", "Advanced Parameter 1", value = 0.3),
                   numericInput("advancedlower", "Advanced Parameter 2", value = 0)
                 )
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
      read_csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results_allB_v3.csv"))
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
    # Now, an additional option for the 'correct' CIs because I was wrong intially,
    # but plotting the SD is still informative, so we keep them around for now
    if(input$advancedBiasCI == T){
      sqrtN_to_divide_by = sqrt(as.numeric(input$n_reps))
    } else {
      sqrtN_to_divide_by = 1
    }
    z_value = qnorm(1-alpha/2)
    
    
    if(input$plot_type == "error"){
      df <- df %>%
        mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
               beta_upper = beta_bias + z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
               omega_lower = omega_bias - z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
               omega_upper = omega_bias + z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
               rho_lower = rho_bias - z_value * sd_of_rho_estimates/sqrtN_to_divide_by,
               rho_upper = rho_bias + z_value * sd_of_rho_estimates/sqrtN_to_divide_by)
    }
    
    if(input$plot_type == "sd"){
      df <- df %>%
        mutate(beta_lower = avg_of_beta_se - z_value * sd_of_beta_se/sqrtN_to_divide_by,
               beta_upper = avg_of_beta_se + z_value * sd_of_beta_se/sqrtN_to_divide_by,
               omega_lower = avg_of_omega_se - z_value * sd_of_omega_se/sqrtN_to_divide_by,
               omega_upper = avg_of_omega_se + z_value * sd_of_omega_se/sqrtN_to_divide_by,
               rho_lower = avg_of_rho_se - z_value * sd_of_rho_se/sqrtN_to_divide_by,
               rho_upper = avg_of_rho_se + z_value * sd_of_rho_se/sqrtN_to_divide_by)
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
        list(geom_line(aes(y = beta_bias, color = "Beta Error")),
             geom_line(aes(y = omega_bias, color = "Omega Error")),
             geom_line(aes(y = rho_bias, color = "Rho Error")),
             geom_point(aes(y = beta_bias, color = "Beta Error"), size = 1),
             geom_point(aes(y = omega_bias, color = "Omega Error"), size = 1),
             geom_point(aes(y = rho_bias, color = "Rho Error"), size = 1)) +
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
        list(geom_line(aes(y = avg_of_beta_se, color = "Beta SD Avg")),
             geom_line(aes(y = avg_of_omega_se, color = "Omega SD Avg")),
             geom_line(aes(y = avg_of_rho_se, color = "Rho SD Avg")),
             geom_line(aes(y = sd_of_rho_estimates, color = "black")),
             geom_point(aes(y = avg_of_beta_se, color = "Beta SD Avg"), size = 1),
             geom_point(aes(y = avg_of_omega_se, color = "Omega SD Avg"), size = 1),
             geom_point(aes(y = avg_of_rho_se, color = "Rho SD Avg"), size = 1),
             geom_point(aes(y = sd_of_rho_estimates), shape = 17, color = "black", size = 1)) +
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
      
      if(input$UseModelbasedSE == T){
        
        p <- p + {
          ci_ribbons <- list()
          if(input$show_beta_ci){
            ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = bino_ci_exact_lower_beta_rescaled, ymax = bino_ci_exact_upper_beta_rescaled, fill = "Beta CI"), alpha = 0.3)))
          }
          if(input$show_omega_ci){
            ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = bino_ci_exact_lower_omega_rescaled, ymax =  bino_ci_exact_upper_omega_rescaled, fill = "Omega CI"), alpha = 0.3)))
          }
          if(input$show_rho_ci){
            ci_ribbons <- append(ci_ribbons, list(geom_ribbon(aes(ymin = bino_ci_exact_lower_rho_rescaled, ymax = bino_ci_exact_upper_rho_rescaled, fill = "Rho CI"), alpha = 0.3)))
          }
          list(ci_ribbons)
        }
        p <- p +
          list(geom_line(aes(y = p_beta_rescaled, color = "Beta")),
               geom_line(aes(y = p_omega_rescaled, color = "Omega")),
               geom_line(aes(y = p_rho_rescaled, color = "Rho")),
               geom_point(aes(y = p_beta_rescaled, color = "Beta"), size = 1),
               geom_point(aes(y = p_omega_rescaled, color = "Omega"), size = 1),
               geom_point(aes(y = p_rho_rescaled, color = "Rho"), size = 1)) +
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
        
      } else {
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
    }
    ## The y-axis limit adjustment for better comparison of values/info ------
    ## 1.) Extract the correct info for limiting the scale (so what is it based on)
    if(input$y_scale_fix == "None"){
      foo <- filtered_data()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_bias < 10, beta_bias > -10, omega_bias < 10,
                              omega_bias > -10, rho_bias < 10, rho_bias > -10) %>% 
          filter(N != 20)
        min = min( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        max = max( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(avg_of_beta_se < 10, avg_of_omega_se < 10, avg_of_rho_se < 10) %>% 
          filter(N != 20)
        min = min( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        max = max( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "ci_succ"){
        if(input$UseModelbasedSE == T){
          min = min( c(foo$p_beta_rescaled, foo$p_omega_rescaled, foo$p_rho_rescaled), na.rm = TRUE )
          max = max( c(foo$p_beta_rescaled, foo$p_omega_rescaled, foo$p_rho_rescaled), na.rm = TRUE )
          p <- p + coord_cartesian(ylim = c(min, max))
        } else {
          min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
          max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
          p <- p + coord_cartesian(ylim = c(min, max))
        }
      }
    }
    
    if(input$y_scale_fix == "Compare Time"){
      foo <- dataset()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_bias < 10, beta_bias > -10, omega_bias < 10,
                              omega_bias > -10, rho_bias < 10, rho_bias > -10) %>% 
          filter(rho == input$rho,
                 omega_var == input$omega_var,
                 n_reps == input$n_reps,
                 beta_free == input$beta) %>% 
          filter(N != 20)
        min = min( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        max = max( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(avg_of_beta_se < 10, avg_of_omega_se < 10, avg_of_rho_se < 10) %>% filter(rho == input$rho,
                                                  omega_var == input$omega_var,
                                                  n_reps == input$n_reps)  %>% 
          filter(N != 20)
        min = min( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        max = max( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        
        
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "ci_succ"){
        foo <- foo %>% filter(rho == input$rho,
                              omega_var == input$omega_var,
                              n_reps == input$n_reps)
        
         if(input$UseModelbasedSE == T){
           min = min( c(foo$p_beta_rescaled, foo$p_omega_rescaled, foo$p_rho_rescaled), na.rm = TRUE )
           max = max( c(foo$p_beta_rescaled, foo$p_omega_rescaled, foo$p_rho_rescaled), na.rm = TRUE )
           p <- p + coord_cartesian(ylim = c(min, max))
         } else {
           min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
           max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
           p <- p + coord_cartesian(ylim = c(min, max))
         }
      }
    }
    
    if(input$y_scale_fix == "Compare Omega"){
      foo <- dataset()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_bias < 10, beta_bias > -10, omega_bias < 10,
                              omega_bias > -10, rho_bias < 10, rho_bias > -10) %>% filter(rho == input$rho,
                                                  Tfull == input$Tfull,
                                                  n_reps == input$n_reps) %>% 
          filter(N != 20)
        min = min( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        max = max( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(avg_of_beta_se < 10, avg_of_omega_se < 10, avg_of_rho_se < 10) %>% filter(rho == input$rho,
                                                  Tfull == input$Tfull,
                                                  n_reps == input$n_reps)  %>% 
          filter(N != 20)
        min = min( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        max = max( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      } 
      if(input$plot_type == "ci_succ"){
        
        foo <- foo %>% filter(rho == input$rho,
                              Tfull == input$Tfull,
                              n_reps == input$n_reps)
        if(input$UseModelbasedSE){
          min = min( c(foo$p_beta_rescaled, foo$p_omega_rescaled, foo$p_rho_rescaled), na.rm = TRUE )
          max = max( c(foo$p_beta_rescaled, foo$p_omega_rescaled, foo$p_rho_rescaled), na.rm = TRUE )
          p <- p + coord_cartesian(ylim = c(min, max))
        } else {
          min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
          max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
          p <- p + coord_cartesian(ylim = c(min, max)) 
        }
      }
    }
    
    if(input$y_scale_fix == "Compare Rho"){
      foo <- dataset()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_bias < 10, beta_bias > -10, omega_bias < 10,
                              omega_bias > -10, rho_bias < 10, rho_bias > -10) %>% filter(omega_var == input$omega_var,
                                                  Tfull == input$Tfull,
                                                  n_reps == input$n_reps,
                                                  beta_free == input$beta) %>% 
          filter(!(N %in% c(20, 30)))
        min = min( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        max = max( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(avg_of_beta_se < 10, avg_of_omega_se < 10, avg_of_rho_se < 10) %>% filter(omega_var == input$omega_var,
                                                  Tfull == input$Tfull,
                                                  n_reps == input$n_reps,
                                                  beta_free == input$beta) %>% 
          filter(!(N %in% c(20, 30)))
        min = min( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        max = max( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "ci_succ"){
        
        foo <- foo %>% filter(omega_var == input$omega_var,
                              Tfull == input$Tfull,
                              n_reps == input$n_reps)
        if(input$UseModelbasedSE == T){
          min = min( c(foo$p_beta_rescaled, foo$p_omega_rescaled, foo$p_rho_rescaled), na.rm = TRUE )
          max = max( c(foo$p_beta_rescaled, foo$p_omega_rescaled, foo$p_rho_rescaled), na.rm = TRUE )
          p <- p + coord_cartesian(ylim = c(min, max))
        } else {
          min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
          max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
          p <- p + coord_cartesian(ylim = c(min, max))
        }
      }
    }
    
    if(input$y_scale_fix == "Compare n_reps"){
      foo <- dataset()
      
      if(input$plot_type == "error"){
        
        foo <- foo %>% filter(beta_bias < 10, beta_bias > -10, omega_bias < 10,
                              omega_bias > -10, rho_bias < 10, rho_bias > -10) %>% filter(omega_var == input$omega_var,
                                                  Tfull == input$Tfull,
                                                  rho == input$rho) %>% 
          filter(N != 20)
        min = min( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        max = max( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(avg_of_beta_se < 10, avg_of_omega_se < 10, avg_of_rho_se < 10) %>% filter(omega_var == input$omega_var,
                                                  Tfull == input$Tfull,
                                                  rho == input$rho)  %>% 
          filter(N != 20)
        min = min( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        max = max( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "ci_succ"){
        
        foo <- foo %>% filter(omega_var == input$omega_var,
                              Tfull == input$Tfull,
                              rho == input$rho)
        if(input$UseModelbasedSE == T){
          min = min( c(foo$p_beta_rescaled, foo$p_omega_rescaled, foo$p_rho_rescaled), na.rm = TRUE )
          max = max( c(foo$p_beta_rescaled, foo$p_omega_rescaled, foo$p_rho_rescaled), na.rm = TRUE )
          p <- p + coord_cartesian(ylim = c(min, max))
        } else {
          min = min( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
          max = max( c(foo$p_beta, foo$p_omega, foo$p_rho), na.rm = TRUE )
          p <- p + coord_cartesian(ylim = c(min, max))
        }
      }
    }
    
    if(input$advancedView){
      if(input$plot_type == "ci_succ"){
        p <- p + coord_cartesian(ylim = c(0.7, 1))
      } else if (input$plot_type == "error"){
        p <- p + coord_cartesian(ylim = c(input$advancedlower, input$advancedupper))
      } else {
        p <- p + coord_cartesian(ylim = c(input$advancedlower, input$advancedupper)) # p + coord_cartesian(ylim = c(-2, 15))
      }
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
        mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates,
               beta_upper = beta_bias + z_value * sd_of_beta_estimates,
               omega_lower = omega_bias - z_value * sd_of_omega_estimates,
               omega_upper = omega_bias + z_value * sd_of_omega_estimates,
               rho_lower = rho_bias - z_value * sd_of_rho_estimates,
               rho_upper = rho_bias + z_value * sd_of_rho_estimates)
    }
    
    if(input$plot_type == "sd"){
      df <- df %>%
        mutate(beta_lower = avg_of_beta_se - z_value * sd_of_beta_se,
               beta_upper = avg_of_beta_se + z_value * sd_of_beta_se,
               omega_lower = avg_of_omega_se - z_value * sd_of_omega_se,
               omega_upper = avg_of_omega_se + z_value * sd_of_omega_se,
               rho_lower = avg_of_rho_se - z_value * sd_of_rho_se,
               rho_upper = avg_of_rho_se + z_value * sd_of_rho_se)
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
        list(geom_line(aes(y = beta_bias, color = "Beta Error")),
             geom_line(aes(y = omega_bias, color = "Omega Error")),
             geom_line(aes(y = rho_bias, color = "Rho Error")),
             geom_point(aes(y = beta_bias, color = "Beta Error"), size = 1),
             geom_point(aes(y = omega_bias, color = "Omega Error"), size = 1),
             geom_point(aes(y = rho_bias, color = "Rho Error"), size = 1)) +
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
        list(geom_line(aes(y = avg_of_beta_se, color = "Beta SD Avg")),
             geom_line(aes(y = avg_of_omega_se, color = "Omega SD Avg")),
             geom_line(aes(y = avg_of_rho_se, color = "Rho SD Avg")),
             geom_point(aes(y = avg_of_beta_se, color = "Beta SD Avg"), size = 1),
             geom_point(aes(y = avg_of_omega_se, color = "Omega SD Avg"), size = 1),
             geom_point(aes(y = avg_of_rho_se, color = "Rho SD Avg"), size = 1)) +
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
        
        foo <- foo %>% filter(beta_bias < 10, beta_bias > -10, omega_bias < 10,
                              omega_bias > -10, rho_bias < 10, rho_bias > -10) %>% 
          filter(N != 20)
        min = min( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        max = max( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(avg_of_beta_se < 10, avg_of_omega_se < 10, avg_of_rho_se < 10) %>% 
          filter(N != 20)
        min = min( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        max = max( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
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
        
        foo <- foo %>% filter(beta_bias < 10, beta_bias > -10, omega_bias < 10,
                              omega_bias > -10, rho_bias < 10, rho_bias > -10) %>% 
          filter(rho == input$rho,omega_var == input$omega_var) %>% 
          filter(N != 20)
        min = min( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        max = max( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(avg_of_beta_se < 10, avg_of_omega_se < 10, avg_of_rho_se < 10) %>% filter(rho == input$rho,
                                                                                               omega_var == input$omega_var)  %>% 
          filter(N != 20)
        min = min( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        max = max( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        
        
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
        
        foo <- foo %>% filter(beta_bias < 10, beta_bias > -10, omega_bias < 10,
                              omega_bias > -10, rho_bias < 10, rho_bias > -10) %>% filter(rho == input$rho,
                                                                                             Tfull == input$Tfull) %>% 
          filter(N != 20)
        min = min( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        max = max( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(avg_of_beta_se < 10, avg_of_omega_se < 10, avg_of_rho_se < 10) %>% filter(rho == input$rho,
                                                                                               Tfull == input$Tfull)  %>% 
          filter(N != 20)
        min = min( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        max = max( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
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
        
        foo <- foo %>% filter(beta_bias < 10, beta_bias > -10, omega_bias < 10,
                              omega_bias > -10, rho_bias < 10, rho_bias > -10) %>% filter(omega_var == input$omega_var,
                                                                                             Tfull == input$Tfull) %>% 
          filter(N != 20)
        min = min( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        max = max( c(foo$beta_bias, foo$omega_bias, foo$rho_bias), na.rm = TRUE )
        p <- p + coord_cartesian(ylim = c(min, max))
        
      }
      if(input$plot_type == "sd"){
        
        foo <- foo %>% filter(avg_of_beta_se < 10, avg_of_omega_se < 10, avg_of_rho_se < 10) %>% filter(omega_var == input$omega_var,
                                                                                               Tfull == input$Tfull)  %>% 
          filter(N != 20)
        min = min( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
        max = max( c(foo$avg_of_beta_se, foo$avg_of_omega_se, foo$avg_of_rho_se), na.rm = TRUE )
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



# Standardised Mean Squared Error Plot ------------------------------------
library(shiny)
library(ggplot2)

#####
ui <- fluidPage(
  titlePanel("Model Parameter vs. SMSE/MSE/Relative MSE"),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      selectInput("plot_type", "Plot Type:",
                  choices = c("SMSE", "MSE"), selected = "SMSE"),
      
      selectInput("xaxis", "Select x-axis of plot:",
                  choices = c("beta", "rho", "omega", "T")),
      # Dynamic UI for the selected plot type*s corresponding variable input
      uiOutput("filter_ui"),
      conditionalPanel(
        condition = "input.plot_type == 'SMSE'",
        selectInput("switchToEmpirical", "Switch appropriate SE", 
                    choices = c("modelSE", "empiricalSE", "avgmodelSE"))
      ),
      conditionalPanel(
        condition = "input.plot_type != 'SMSE'",
        selectInput("reference", "Reference Options:",
                    choices = c("Non", "Relative to set MSE"),
                    selected = "Non")  
      ),
      conditionalPanel(
        condition = "input.plot_type == 'MSE' && input.reference == 'Relative to set MSE'",
        selectInput("whatIsRef", "Divide by:",
                    choices = c("MSE/Ref.MSE", "MSE/Ref.Emp.SE", "MSE/Ref.Avg.SE", "(MSE/Emp.SE)/Ref.NMSE",
                                "(MSE/Avg.SE)/Ref.NMSE"),
                    selected = "MSE")  
      ),
      conditionalPanel(
        condition = "input.plot_type == 'MSE' && input.reference == 'Relative to set MSE'",
        numericInput("ref_value", "Set ref:", value = 0, min = -1, max = 10)
      )
    ),
    
    
    mainPanel(width = 10,
      plotOutput("main_plot")
    ),
    position = "right"
  )
)

#####
server <- function(input, output, session) {
  
  # Load dataset
  dataset <- reactive({
  req(input$xaxis)
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
    req(input$xaxis)
    other_params <- setdiff(names(param_values), input$xaxis)
    
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
    req(input$xaxis, dataset())
    
    df <- dataset()
    plot_param <- input$xaxis
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
    x_levels <- param_values[[input$xaxis]]
    
    if(input$plot_type == "SMSE"){
      p <- ggplot(df, aes(x = x_levels)) 
      
      #SMSE plot logic
      if(input$switchToEmpirical == "empiricalSE"){
        p <- p + geom_line(aes(y = beta_smse_with_emp_se, color = "Beta Error")) +
          geom_point(aes(y = beta_smse_with_emp_se, color = "Beta Error"), size = 1) +
          geom_line(aes(y = omega_smse_with_emp_se, color = "Omega Error")) +
          geom_point(aes(y = omega_smse_with_emp_se, color = "Omega Error"), size = 1) +
          geom_line(aes(y = rho_smse_with_emp_se, color = "Rho Error")) +
          geom_point(aes(y = rho_smse_with_emp_se, color = "Rho Error"), size = 1)
      } 
      if(input$switchToEmpirical == "modelSE"){
        p <- p + geom_line(aes(y = beta_smse_with_model_se, color = "Beta Error")) +
          geom_point(aes(y = beta_smse_with_model_se, color = "Beta Error"), size = 1) +
          geom_line(aes(y = omega_smse_with_model_se, color = "Omega Error")) +
          geom_point(aes(y = omega_smse_with_model_se, color = "Omega Error"), size = 1) +
          geom_line(aes(y = rho_smse_with_model_se, color = "Rho Error")) +
          geom_point(aes(y = rho_smse_with_model_se, color = "Rho Error"), size = 1)
      }
      if(input$switchToEmpirical == "avgmodelSE"){
        p <- p + geom_line(aes(y = beta_smse_with_avg_model_se, color = "Beta Error")) +
          geom_point(aes(y = beta_smse_with_avg_model_se, color = "Beta Error"), size = 1) +
          geom_line(aes(y = omega_smse_with_avg_model_se, color = "Omega Error")) +
          geom_point(aes(y = omega_smse_with_avg_model_se, color = "Omega Error"), size = 1) +
          geom_line(aes(y = rho_smse_with_avg_model_se, color = "Rho Error")) +
          geom_point(aes(y = rho_smse_with_avg_model_se, color = "Rho Error"), size = 1)
      }
      
      p <- p + labs(
        title = paste(input$plot_type, " vs. ", input$xaxis),
        x = input$xaxis,
        y = input$plot_type,
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
      
      if(input$switchToEmpirical == "empiricalSE"){
        p <- p + coord_cartesian(ylim = c(min(c(df$beta_smse_with_emp_se, df$omega_smse_with_emp_se, df$rho_smse_with_emp_se, 1)),
                                          max(c(df$beta_smse_with_emp_se, df$omega_smse_with_emp_se, df$rho_smse_with_emp_se, 1))))
      } else if (input$switchToEmpirical == "avgmodelSE"){
        p <- p + coord_cartesian(ylim = c(min(c(df$beta_smse_with_avg_model_se, df$omega_smse_with_avg_model_se, df$rho_smse_with_avg_model_se, 1)),
                                          max(c(df$beta_smse_with_avg_model_se, df$omega_smse_with_avg_model_se, df$rho_smse_with_avg_model_se, 1))))
      }
      
    } else if (input$plot_type == "MSE"){
      if(input$reference == "Relative to set MSE"){
        ref_value = as.numeric(input$ref_value)
        ref_df = df[df[[param_mapping[[input$xaxis]]]] == ref_value,]
        if(input$whatIsRef == "MSE/Ref.MSE"){
          df <- df %>% rowwise %>% mutate(beta_mse = beta_mse/ref_df[["beta_mse"]],
                                          omega_mse = omega_mse/ref_df[["omega_mse"]],
                                          rho_mse = rho_mse/ref_df[["rho_mse"]])
        }
        if(input$whatIsRef == "MSE/Ref.Emp.SE"){
          df <- df %>% rowwise %>% mutate(beta_mse = beta_mse/(ref_df[["sd_of_beta_estimates"]])^2,
                                          omega_mse = omega_mse/(ref_df[["sd_of_omega_estimates"]])^2,
                                          rho_mse = rho_mse/(ref_df[["sd_of_rho_estimates"]])^2)
          
        }
        if(input$whatIsRef == "MSE/Ref.Avg.SE"){
          df <- df %>% rowwise %>% mutate(beta_mse = beta_mse/(ref_df[["avg_of_beta_se"]])^2,
                                          omega_mse = omega_mse/(ref_df[["avg_of_omega_se"]])^2,
                                          rho_mse = rho_mse/(ref_df[["avg_of_rho_se"]])^2)
        }
        if(input$whatIsRef == "(MSE/Emp.SE)/Ref.NMSE"){
          df <- df %>% rowwise %>% mutate(beta_mse = (beta_mse/sd_of_beta_estimates^2)/(ref_df[["beta_mse"]]/ref_df[["sd_of_beta_estimates"]]^2),
                                          omega_mse = (omega_mse/sd_of_omega_estimates^2)/(ref_df[["omega_mse"]]/ref_df[["sd_of_omega_estimates"]]^2),
                                          rho_mse = (rho_mse/sd_of_rho_estimates^2)/(ref_df[["rho_mse"]]/ref_df[["sd_of_rho_estimates"]]^2))
        }
        if(input$whatIsRef == "(MSE/Avg.SE)/Ref.NMSE"){
          df <- df %>% rowwise %>% mutate(beta_mse = (beta_mse/avg_of_beta_se^2)/(ref_df[["beta_mse"]]/ref_df[["avg_of_beta_se"]]^2),
                                          omega_mse = (omega_mse/avg_of_omega_se^2)/(ref_df[["omega_mse"]]/ref_df[["avg_of_omega_se"]]^2),
                                          rho_mse = (rho_mse/avg_of_rho_se^2)/(ref_df[["rho_mse"]]/ref_df[["avg_of_rho_se"]]^2))
        }
      } 
      
      
      p <- ggplot(df, aes(x = x_levels)) 
      
      p <- p + geom_line(aes(y = beta_mse, color = "Beta Error")) +
          geom_point(aes(y = beta_mse, color = "Beta Error"), size = 1) +
          geom_line(aes(y = omega_mse, color = "Omega Error")) +
          geom_point(aes(y = omega_mse, color = "Omega Error"), size = 1) +
          geom_line(aes(y = rho_mse, color = "Rho Error")) +
          geom_point(aes(y = rho_mse, color = "Rho Error"), size = 1)
        
      p <- p + labs(
          title = paste(input$plot_type, " vs. ", input$xaxis),
          x = input$xaxis,
          y = input$plot_type,
          color = "Error Type"
        ) +
          theme_minimal() +
          theme(
            legend.position = c(0.85, 0.85),
            plot.title = element_text(hjust = 0.5)
          ) +
          scale_x_continuous(breaks = x_levels) +
          scale_color_manual(values = c(
            "Beta Error" = "red",
            "Omega Error" = "#006400",
            "Rho Error" = "blue"
          ))
      
      if(input$whatIsRef %in% c("MSE/Ref.Emp.SE", "MSE/Ref.Avg.SE")){
        p <- p + coord_cartesian(ylim = c(0, 1.25)) 
      } else if (input$whatIsRef == "MSE/Ref.MSE"){
        p <- p + coord_cartesian(ylim = c(0, 1.25)) + geom_hline(yintercept = 1, linetype = "dotted", color = "black")
      } else if (input$whatIsRef == "(MSE/Avg.SE)/Ref.NMSE"){
        p <- p + coord_cartesian(ylim = c(0.8, 1.6)) + geom_hline(yintercept = 1, linetype = "dotted", color = "black")
      } else if (input$whatIsRef == "(MSE/Emp.SE)/Ref.NMSE"){
        p <- p + coord_cartesian(ylim = c(0.9, 1.1)) + geom_hline(yintercept = 1, linetype = "dotted", color = "black")
      }
      
      } # <- bracket is from if(plot_type == "MSE)
      print(p)
  })
}


shinyApp(ui = ui, server = server)


# Testing -----------------------------------------------------------------
ref_df = df = read.csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results_allB_v2.csv"))
# Choose a ref
ref_df = ref_df %>% filter(Tfull == 2, rho == 0, beta_free == 4, omega_var == 0.0625, N == 500, n_reps == 3000) 
  
df <- df %>% filter(Tfull == 2, N == 500, n_reps == 3000, omega_var == 0.0625) 
df <- df %>% rowwise %>% mutate(beta_mse = (beta_mse/avg_of_beta_se^2)/(ref_df[["beta_mse"]]/ref_df[["avg_of_beta_se"]]^2),
                                omega_mse = (omega_mse/avg_of_omega_se^2)/(ref_df[["omega_mse"]]/ref_df[["avg_of_omega_se"]]^2),
                                rho_mse = (rho_mse/avg_of_rho_se^2)/(ref_df[["rho_mse"]]/ref_df[["avg_of_rho_se"]]^2))

ref_text <- paste0(
  "N = ", ref_df$N, "\n",
  "n_reps = ", ref_df$n_reps, "\n",
  "Reference:\nbeta = ", ref_df$beta_free,
  ", T = ", ref_df$Tfull,
  ", rho = ", ref_df$rho,
  ", omega = ", ref_df$omega_var, "\n",
  "Active In-plot filter:\nT = ", df$Tfull,
  ", omega = ", df$omega_var
)

p <- ggplot(df, aes(x = rho, y = beta_mse, group = beta_free, color = beta_free)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for beta param vs. Rho",
       x = "Rho",
       y = "NMSE / ReferenceNMSE") +
  coord_cartesian(ylim = c(0, 1.5)) +
  scale_x_continuous(breaks = unique(df$rho)) +
  scale_color_gradient(low = "blue", high = "red", name = "Beta") + 
  annotate("text", x = min(df$rho), y = 0.3, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

# Does this make much sense in it's current form? The division of the one rho ----
# reference value for all other possible rho values doesn't make much sense when
# interpreting, i.e. I can only compare the relative NMSE ranking for other rho 
# values with each other, but I can't tell anything compared to the reference
# So, should the entire rho-realm be used as a reference, meaning that each rho
# value is divided by the NMSE with respect to the same rho value?

## In-plot: Rho and Beta, Omega varies among plots ----
plotted_omega = c(0.25^2, 0.5^2, 0.75^2, 1, 2^2, 3^2, 4^2)
plotted_beta = c(rep.int(2, length(plotted_omega))) # reference beta used
plotted_T = 2 # T used for plots
plotted_rho = 0 # reference rho used in additional plots

for (i in 1:length(plotted_omega)){
ref_df = df = read.csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results_allB_v4.csv"))
# Choose a ref
ref_df = ref_df %>% filter(Tfull == plotted_T, beta_free == plotted_beta[i], omega_var == plotted_omega[i], N == 500, n_reps == 3000)

ref_df_filtered <- ref_df %>%
  dplyr::select(rho, beta_mse, omega_mse, rho_mse, avg_of_beta_se, avg_of_omega_se, avg_of_rho_se) %>%
  rename(
    ref_beta_mse = beta_mse,
    ref_omega_mse = omega_mse,
    ref_rho_mse = rho_mse,
    ref_avg_beta_se = avg_of_beta_se,
    ref_avg_omega_se = avg_of_omega_se,
    ref_avg_rho_se = avg_of_rho_se
  )

df_filtered <- df %>%
  filter(Tfull == plotted_T, N == 500, n_reps == 3000, omega_var == plotted_omega[i])

df_final <- df_filtered %>%
  left_join(ref_df_filtered, by = "rho") %>%
  mutate(
    beta_mse = (beta_mse / avg_of_beta_se^2) / (ref_beta_mse / ref_avg_beta_se^2),
    omega_mse = (omega_mse / avg_of_omega_se^2) / (ref_omega_mse / ref_avg_omega_se^2),
    rho_mse   = (rho_mse   / avg_of_rho_se^2)   / (ref_rho_mse   / ref_avg_rho_se^2)
  )

ref_text <- paste0(
  "N = ", ref_df$N, "\n",
  "n_reps = ", ref_df$n_reps, "\n",
  "Reference:\nbeta = ", ref_df$beta_free,
  ", T = ", ref_df$Tfull,
  ", rho = (-0.9, +0.9)",
  ", omega = ", ref_df$omega_var, "\n",
  "Non-reference filter:\nT = ", df_final$Tfull,
  ", omega = ", df_final$omega_var
)

p <- ggplot(df_final, aes(x = rho, y = beta_mse, group = beta_free, color = beta_free)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for beta param vs. Rho",
       x = "Rho",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.5, 1.5)) +
  scale_x_continuous(breaks = unique(df$rho)) +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 2, name = "Beta") + 
  #scale_color_gradient(low = "yellow", high = "red", name = "Beta") + 
  annotate("text", x = min(df$rho), y = 0.55, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)
}

## In-plot: Rho and Omega, Beta varies among plots ----
plotted_omega = c(1)
plotted_beta = c(0, 1, 2, 3, 4) # reference beta used
plotted_T = 2

for (i in 1:length(plotted_beta)){
  ref_df = df = read.csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results_allB_v4.csv"))
  # Choose a ref
  ref_df = ref_df %>% filter(Tfull == plotted_T, beta_free == plotted_beta[i], omega_var == plotted_omega, N == 500, n_reps == 3000)
  
  ref_df_filtered <- ref_df %>%
    dplyr::select(rho, beta_mse, omega_mse, rho_mse, avg_of_beta_se, avg_of_omega_se, avg_of_rho_se) %>%
    rename(
      ref_beta_mse = beta_mse,
      ref_omega_mse = omega_mse,
      ref_rho_mse = rho_mse,
      ref_avg_beta_se = avg_of_beta_se,
      ref_avg_omega_se = avg_of_omega_se,
      ref_avg_rho_se = avg_of_rho_se
    )
  
  df_filtered <- df %>%
    filter(Tfull == plotted_T, N == 500, n_reps == 3000, beta_free == plotted_beta[i])
  
  df_final <- df_filtered %>%
    left_join(ref_df_filtered, by = "rho") %>%
    mutate(
      beta_mse = (beta_mse / avg_of_beta_se^2) / (ref_beta_mse / ref_avg_beta_se^2),
      omega_mse = (omega_mse / avg_of_omega_se^2) / (ref_omega_mse / ref_avg_omega_se^2),
      rho_mse   = (rho_mse   / avg_of_rho_se^2)   / (ref_rho_mse   / ref_avg_rho_se^2)
    )
  
  ref_text <- paste0(
    "N = ", ref_df$N, "\n",
    "n_reps = ", ref_df$n_reps, "\n",
    "Reference:\nbeta = ", ref_df$beta_free,
    ", T = ", ref_df$Tfull,
    ", rho = (-0.9, +0.9)",
    ", omega = ", ref_df$omega_var, "\n",
    "Non-reference filter:\nT = ", df_final$Tfull,
    ", Corresp. omega & rho, beta =", df_filtered$beta_free
  )
  
  p <- ggplot(df_final, aes(x = rho, y = beta_mse, group = factor(omega_var), color = factor(omega_var))) + 
    geom_line() +
    geom_point(size = 1) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) + 
    labs(title = "Relative NMSE Ranking for beta param vs. Rho & Omega",
         x = "Rho",
         y = "NMSE_i / NMSE_ref_i") +
    coord_cartesian(ylim = c(0.5, 1.5)) +
    scale_x_continuous(breaks = unique(df$rho)) +
    scale_color_manual(
      values = c( 
        "0.0625" = "purple",  
        "0.25"   = "#2166ac",  # deep blue
        "0.5625" = "#67a9cf",  # light blue
        "1"      = "#ffff33",  # bright yellow (midpoint)
        "4"      = "#fdae61",  # orange
        "9"      = "#f46d43",  # red-orange
        "16"     = "#d73027"   # deep red
      ),,
      name = "Omega"
    ) +
    annotate("text", x = min(df$rho), y = 0.55, 
             label = ref_text, hjust = 0, vjust = 0, size = 3.5)
  print(p)
}

## In-plot: Beta and Omega, Rho varies among plots ----
plotted_omega = c(1)
plotted_rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9) # reference rho used
plotted_T = 2

for (i in 1:length(plotted_rho)){
  ref_df = df = read.csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results_allB_v4.csv"))
  # Choose a ref
  ref_df = ref_df %>% filter(Tfull == plotted_T, rho == plotted_rho[i], omega_var == plotted_omega, N == 500, n_reps == 3000)
  
  ref_df_filtered <- ref_df %>%
    dplyr::select(beta_free, beta_mse, omega_mse, rho_mse, avg_of_beta_se, avg_of_omega_se, avg_of_rho_se) %>%
    rename(
      ref_beta_mse = beta_mse,
      ref_omega_mse = omega_mse,
      ref_rho_mse = rho_mse,
      ref_avg_beta_se = avg_of_beta_se,
      ref_avg_omega_se = avg_of_omega_se,
      ref_avg_rho_se = avg_of_rho_se
    )
  
  df_filtered <- df %>%
    filter(Tfull == plotted_T, N == 500, n_reps == 3000, rho == plotted_rho[i])
  
  df_final <- df_filtered %>%
    left_join(ref_df_filtered, by = "beta_free") %>%
    mutate(
      beta_mse = (beta_mse / avg_of_beta_se^2) / (ref_beta_mse / ref_avg_beta_se^2),
      omega_mse = (omega_mse / avg_of_omega_se^2) / (ref_omega_mse / ref_avg_omega_se^2),
      rho_mse   = (rho_mse   / avg_of_rho_se^2)   / (ref_rho_mse   / ref_avg_rho_se^2)
    )
  
  ref_text <- paste0(
    "N = ", ref_df$N, "\n",
    "n_reps = ", ref_df$n_reps, "\n",
    "Reference:\nbeta = (0, 4)",
    ", T = ", ref_df$Tfull,
    ", rho = ", ref_df$rho,
    ", omega = ", ref_df$omega_var, "\n",
    "Non-reference filter:\nT = ", df_final$Tfull,
    ", rho = ", ref_df$rho, ", Corresp. Omega"
  )
  
  p <- ggplot(df_final, aes(x = beta_free, y = beta_mse, group = factor(omega_var), color = factor(omega_var))) + 
    geom_line() +
    geom_point(size = 1) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) + 
    labs(title = "Relative NMSE Ranking for beta param vs. Beta & Omega",
         x = "Beta",
         y = "NMSE_i / NMSE_ref_i") +
    coord_cartesian(ylim = c(0.75, 1.25)) +
    scale_x_continuous(breaks = unique(df$beta_free)) +
    scale_color_manual(
      values = c( 
        "0.0625" = "purple",  
        "0.25"   = "#2166ac",  # deep blue
        "0.5625" = "#67a9cf",  # light blue
        "1"      = "#ffff33",  # bright yellow (midpoint)
        "4"      = "#fdae61",  # orange
        "9"      = "#f46d43",  # red-orange
        "16"     = "#d73027"   # deep red
      ),
      name = "Omega"
    ) +
    annotate("text", x = min(df$beta_free), y = 0.8, 
             label = ref_text, hjust = 0, vjust = 0, size = 3.5)
  print(p)
}

# Build onto previous plot: -----------------------------------------------
# This plot uses the wrong reference
df_final_non_ref = df_final # %>% filter(beta_free != unique(ref_df$beta_free)) #%>% filter(rho %in% c(-0.9, -0.6, -0.3, 0))


p <- ggplot(df_final_non_ref, aes(x = beta_free, y = beta_mse, group = rho, color = rho)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for beta param vs. Beta",
       x = "Beta",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.75, 1.5)) +
  scale_x_continuous(breaks = unique(df_final_non_ref$beta_free)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_color_gradientn(colors = c("blue", "purple", "red"), name = "Rho") + 
  annotate("text", x = min(df_final_non_ref$beta_free), y = 0.8, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

# Use another color scheme?

p <- ggplot(df_final_non_ref, aes(x = beta_free, y = beta_mse, group = rho, color = rho)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for beta param vs. Beta",
       x = "Beta",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.75, 1.5)) +
  scale_x_continuous(breaks = unique(df_final_non_ref$beta_free)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_color_gradient2(
    low = "red",      
    mid = "yellow",     
    high = "blue",     
    midpoint = 0,
    name = "Rho"
  ) + 
  annotate("text", x = min(df_final_non_ref$beta_free), y = 0.8, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

# . -----------------------------------------------------------------------
# This one uses varying beta values as reference, so it divides by the NMSE_beta=i
ref_df2 = df2 = read.csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results_allB_v2.csv"))
# Choose a ref
ref_df2 = ref_df2 %>% filter(Tfull == 2, rho == plotted_rho, omega_var == plotted_omega, N == 500, n_reps == 3000)

ref_df_filtered2 <- ref_df2 %>%
  dplyr::select(beta_free, beta_mse, omega_mse, rho_mse, avg_of_beta_se, avg_of_omega_se, avg_of_rho_se) %>%
  rename(
    ref_beta_mse = beta_mse,
    ref_omega_mse = omega_mse,
    ref_rho_mse = rho_mse,
    ref_avg_beta_se = avg_of_beta_se,
    ref_avg_omega_se = avg_of_omega_se,
    ref_avg_rho_se = avg_of_rho_se
  )

df_filtered2 <- df2 %>%
  filter(Tfull == 2, N == 500, n_reps == 3000, omega_var == plotted_omega)

df_final2 <- df_filtered2 %>%
  left_join(ref_df_filtered2, by = "beta_free") %>%
  mutate(
    beta_mse = (beta_mse / avg_of_beta_se^2) / (ref_beta_mse / ref_avg_beta_se^2),
    omega_mse = (omega_mse / avg_of_omega_se^2) / (ref_omega_mse / ref_avg_omega_se^2),
    rho_mse   = (rho_mse   / avg_of_rho_se^2)   / (ref_rho_mse   / ref_avg_rho_se^2)
  )

ref_text <- paste0(
  "N = ", ref_df2$N, "\n",
  "n_reps = ", ref_df2$n_reps, "\n",
  "Reference:\nbeta = (0, 1, 2, 3, 4)",
  ", T = ", ref_df2$Tfull,
  ", rho = 0",
  ", omega = ", ref_df2$omega_var, "\n",
  "Active In-plot filter:\nT = ", df_final2$Tfull,
  ", omega = ", df_final2$omega_var
)

p <- ggplot(df_final2, aes(x = beta_free, y = beta_mse, group = rho, color = rho)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for beta param vs. Beta",
       x = "Beta",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.75, 1.5)) +
  scale_x_continuous(breaks = unique(df_final2$beta_free)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_color_gradientn(colors = c("blue", "purple", "red"), name = "Rho") + 
  annotate("text", x = min(df_final2$beta_free), y = 0.8, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)


#. -------------------------------------------------------------------------


# What about the other parameters? ----------------------------------------

p <- ggplot(df_final, aes(x = rho, y = omega_mse, group = beta_free, color = beta_free)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for omega param vs. Rho",
       x = "Rho",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0, 1.5)) +
  scale_x_continuous(breaks = unique(df$rho)) +
  scale_color_gradient(low = "blue", high = "red", name = "Beta") + 
  annotate("text", x = min(df$rho), y = 0.3, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

p <- ggplot(df_final_non_ref, aes(x = beta_free, y = omega_mse, group = rho, color = rho)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for omega param vs. Beta",
       x = "Beta",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.75, 1.5)) +
  scale_x_continuous(breaks = unique(df_final_non_ref$beta_free)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_color_gradient2(
    low = "red",      
    mid = "blue",     
    high = "red",     
    midpoint = 0,
    name = "Rho"
  ) + 
  annotate("text", x = min(df_final_non_ref$beta_free), y = 0.8, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

# . -----------------------------------------------------------------------
# This one uses varying beta values as reference, so it divides by the NMSE_beta=i
p <- ggplot(df_final2, aes(x = beta_free, y = omega_mse, group = rho, color = rho)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for omega param vs. Beta",
       x = "Beta",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.75, 1.5)) +
  scale_x_continuous(breaks = unique(df_final2$beta_free)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_color_gradientn(colors = c("blue", "purple", "red"), name = "Rho") + 
  annotate("text", x = min(df_final2$beta_free), y = 0.8, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)


#. -------------------------------------------------------------------------


p <- ggplot(df_final, aes(x = rho, y = rho_mse, group = beta_free, color = beta_free)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for rho param vs. Rho",
       x = "Rho",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0, 1.5)) +
  scale_x_continuous(breaks = unique(df$rho)) +
  scale_color_gradient(low = "blue", high = "red", name = "Beta") + 
  annotate("text", x = min(df$rho), y = 0.3, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

p <- ggplot(df_final_non_ref, aes(x = beta_free, y = rho_mse, group = rho, color = rho)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for rho param vs. Beta",
       x = "Beta",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.75, 1.5)) +
  scale_x_continuous(breaks = unique(df_final_non_ref$beta_free)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_color_gradient2(
    low = "red",      
    mid = "blue",     
    high = "red",     
    midpoint = 0,
    name = "Rho"
  ) + 
  annotate("text", x = min(df_final_non_ref$beta_free), y = 0.8, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

# . -----------------------------------------------------------------------
# This one uses varying beta values as reference, so it divides by the NMSE_beta=i
p <- ggplot(df_final2, aes(x = beta_free, y = rho_mse, group = rho, color = rho)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for rho param vs. Beta",
       x = "Beta",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.75, 1.5)) +
  scale_x_continuous(breaks = unique(df_final2$beta_free)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_color_gradientn(colors = c("blue", "purple", "red"), name = "Rho") + 
  annotate("text", x = min(df_final2$beta_free), y = 0.8, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)


#. -------------------------------------------------------------------------

# The same thing as before, but different SE selected (HERE: SD of estimates) ---------------------
ref_df = df = read.csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results_allB_v2.csv"))
# Choose a ref
ref_df = ref_df %>% filter(Tfull == 2, beta_free == 4, omega_var == 0.0625, N == 500, n_reps == 3000)

ref_df_filtered <- ref_df %>%
  dplyr::select(rho, beta_mse, omega_mse, rho_mse, sd_of_beta_estimates, sd_of_omega_estimates, sd_of_rho_estimates) %>%
  rename(
    ref_beta_mse = beta_mse,
    ref_omega_mse = omega_mse,
    ref_rho_mse = rho_mse,
    ref_avg_beta_se = sd_of_beta_estimates,
    ref_avg_omega_se = sd_of_omega_estimates,
    ref_avg_rho_se = sd_of_rho_estimates
  )

df_filtered <- df %>%
  filter(Tfull == 2, N == 500, n_reps == 3000, omega_var == 0.0625)

df_final <- df_filtered %>%
  left_join(ref_df_filtered, by = "rho") %>%
  mutate(
    beta_mse = (beta_mse / sd_of_beta_estimates^2) / (ref_beta_mse / ref_avg_beta_se^2),
    omega_mse = (omega_mse / sd_of_omega_estimates^2) / (ref_omega_mse / ref_avg_omega_se^2),
    rho_mse   = (rho_mse   / sd_of_rho_estimates^2)   / (ref_rho_mse   / ref_avg_rho_se^2)
  )

ref_text <- paste0(
  "N = ", ref_df$N, "\n",
  "n_reps = ", ref_df$n_reps, "\n",
  "Reference:\nbeta = ", ref_df$beta_free,
  ", T = ", ref_df$Tfull,
  ", rho = (-0.9, +0.9)",
  ", omega = ", ref_df$omega_var, "\n",
  "Active In-plot filter:\nT = ", df_final$Tfull,
  ", omega = ", df_final$omega_var
)

p <- ggplot(df_final, aes(x = rho, y = beta_mse, group = beta_free, color = beta_free)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for beta param vs. Rho",
       x = "Rho",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.75, 1.25)) +
  scale_x_continuous(breaks = unique(df$rho)) +
  scale_color_gradient(low = "blue", high = "red", name = "Beta") + 
  annotate("text", x = min(df$rho), y = 0.3, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

df_final_non_ref = df_final #%>% filter(beta_free != unique(ref_df$beta_free)) #%>% filter(rho %in% c(-0.9, -0.6, -0.3, 0))

p <- ggplot(df_final_non_ref, aes(x = beta_free, y = beta_mse, group = rho, color = rho)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for beta param vs. Beta",
       x = "Beta",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.75, 1.25)) +
  scale_x_continuous(breaks = unique(df_final_non_ref$beta_free)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_color_gradientn(colors = c("blue", "purple", "red"), name = "Rho") + 
  annotate("text", x = min(df_final_non_ref$beta_free), y = 0.8, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

# Use another color scheme?

p <- ggplot(df_final_non_ref, aes(x = beta_free, y = beta_mse, group = rho, color = rho)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for beta param vs. Beta",
       x = "Beta",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.75, 1.5)) +
  scale_x_continuous(breaks = unique(df_final_non_ref$beta_free)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_color_gradient2(
    low = "red",      
    mid = "blue",     
    high = "red",     
    midpoint = 0,
    name = "Rho"
  ) + 
  annotate("text", x = min(df_final_non_ref$beta_free), y = 0.8, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

p <- ggplot(df_final, aes(x = rho, y = omega_mse, group = beta_free, color = beta_free)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for omega param vs. Rho",
       x = "Rho",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.95, 1.05)) +
  scale_x_continuous(breaks = unique(df$rho)) +
  scale_color_gradient(low = "blue", high = "red", name = "Beta") + 
  annotate("text", x = min(df$rho), y = 0.3, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

p <- ggplot(df_final_non_ref, aes(x = beta_free, y = omega_mse, group = rho, color = rho)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for omega param vs. Beta",
       x = "Beta",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.95, 1.05)) +
  scale_x_continuous(breaks = unique(df_final_non_ref$beta_free)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_color_gradient2(
    low = "red",      
    mid = "blue",     
    high = "red",     
    midpoint = 0,
    name = "Rho"
  ) + 
  annotate("text", x = min(df_final_non_ref$beta_free), y = 0.8, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

p <- ggplot(df_final, aes(x = rho, y = rho_mse, group = beta_free, color = beta_free)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for rho param vs. Rho",
       x = "Rho",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.95, 1.05)) +
  scale_x_continuous(breaks = unique(df$rho)) +
  scale_color_gradient(low = "blue", high = "red", name = "Beta") + 
  annotate("text", x = min(df$rho), y = 0.3, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

p <- ggplot(df_final_non_ref, aes(x = beta_free, y = rho_mse, group = rho, color = rho)) + 
  geom_line() +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(title = "Relative NMSE Ranking for omega param vs. Beta",
       x = "Beta",
       y = "NMSE_i / NMSE_ref_i") +
  coord_cartesian(ylim = c(0.75, 1.5)) +
  scale_x_continuous(breaks = unique(df_final_non_ref$beta_free)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
  scale_color_gradient2(
    low = "red",      
    mid = "blue",     
    high = "red",     
    midpoint = 0,
    name = "Rho"
  ) + 
  annotate("text", x = min(df_final_non_ref$beta_free), y = 0.8, 
           label = ref_text, hjust = 0, vjust = 0, size = 3.5)
print(p)

