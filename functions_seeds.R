# Define custom function to select peak and filter data
select_peak_and_filter <- function(data, population, max_heat_treatment) {
  # Print population
  cat("Population:", population, "\n")
  
  # Get unique maximum heat treatment value for the population
  unique_max_heat_treatment <- unique(max_heat_treatment)
  
  # Print the maximum heat treatment value (if it's unique)
  if (length(unique_max_heat_treatment) == 1) {
    cat("Max Treatment:", unique_max_heat_treatment, "\n")
  }
  
  # Calculate mean germination percentage for each heat treatment
  mean_germination <- data %>%
    group_by(species, population, heat_treatment) %>%
    summarise(mean_gp = mean(gp_standardized))
  
  # Plot mean germination percentage
  plot(mean_germination$heat_treatment, mean_germination$mean_gp, 
       type = "o", 
       xlab = "Heat Treatment", 
       ylab = "Mean Germination Percentage",
       main = paste("Mean Germination Percentage vs. Heat Treatment for Population:", population, "\n", "Max Treatment:", unique_max_heat_treatment))
  
  # Customize x-axis to split by 5 degrees
  axis(1, at = seq(0, max(mean_germination$heat_treatment), by = 5))
  
  # Allow user to select peak interactively
  peak <- readline(prompt = "Enter peak (y/n): ")
  
  if (tolower(peak) == "y") {
    peak_value <- as.numeric(readline(prompt = "Enter peak value: "))
    
    # Filter data based on peak heat treatment
    filtered_data <- data[data$heat_treatment <= peak_value, ]
    return(filtered_data)
  } else {
    # If peak is not selected, use the original max_heat_treatment value for filtering
    max_heat_treatment <- max(max_heat_treatment)
    filtered_data <- data[data$heat_treatment <= max_heat_treatment, ]
    return(filtered_data)
  }
}


# function for 'ed'
ed <- function(x) {
  ED(x, c(0.2, 0.3, 0.4, 0.5, 0.6), type = "absolute", interval = "delta")
}