##### RQ 1 #####

library(dplyr)       # For data manipulation
library(ggplot2)     # For data visualization
library(tidyr)       # For data tidying
library(readr)       # For reading and writing data
library(purrr)       # For functional programming
library(tibble)      # For creating and manipulating tibbles
library(stringr)     # For string manipulation
library(forcats)     # For working with factors

unlist_keep_na <- function(lst) {
  lst[sapply(lst, is.null)] <- NA
  result <- unlist(lst, use.names = FALSE)
  return(result)
}

parse_abroca_data <- function(abroca_matrix) {
  if ('matrix' %in% class(abroca_matrix)) {
    return(
      as.data.frame(t(abroca_matrix)) %>%
        mutate(abroca = unlist_keep_na(abroca), 
               auc0 = unlist_keep_na(auc0), 
               auc1 = unlist_keep_na(auc1)) %>%
        mutate(auc_delta = auc1-auc0)
    )
  }
  abroca_matrix <- abroca_matrix[!sapply(abroca_matrix, function(x) all(is.na(x)))]
  result_df <- data.frame(
    abroca = numeric(),
    auc0 = numeric(),
    auc1 = numeric(),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(abroca_matrix)) {
    element_data <- abroca_matrix[[i]]
    abroca_value <- NA
    auc0_value <- NA
    auc1_value <- NA
    if (!is.null(element_data$abroca)) {
      abroca_value <- element_data$abroca[[1]]
    }
    if (!is.null(element_data$auc0)) {
      auc0_value <- element_data$auc0[[1]]
    }
    if (!is.null(element_data$auc1)) {
      auc1_value <- element_data$auc1[[1]]
    }
    result_df <- rbind(result_df, data.frame(
      abroca = abroca_value,
      auc0 = auc0_value,
      auc1 = auc1_value
    ))
  }
  return(result_df)
}

process_auc_sim <- function(f='myoutput.rds') {
  ans <- readRDS(f) %>%
    imap_dfr(function(m, i) {
      m %>%
        parse_abroca_data() %>%
        mutate(ref=i)
    }) %>%
    separate(ref, into=c('ratio_minority', 'ratio_pos_case', 'auc_minority', 'n'), sep='-') %>%
    tibble() %>%
    mutate(auc_delta = auc1-auc0)
  return(ans)
}

d <- process_auc_sim('input.rds')

# Convert columns to numeric where applicable
data <- d %>%
  mutate(
    ratio_minority = as.numeric(ratio_minority),
    ratio_pos_case = as.numeric(ratio_pos_case),
    auc_minority = as.numeric(auc_minority),
    n = as.numeric(n)
  )

# Load necessary libraries
library(gamlss)       # For generalized gamma fitting
library(gamlss.dist)  # For generalized gamma functions
library(fitdistrplus) # For fitting distributions
library(tidyverse)

# Subset the data for n = 500
data_filtered <- data %>%
  filter(n == 5000) %>%
  dplyr::select(abroca)

abroca_values <- data_filtered$abroca

# Required libraries
library(MASS)
library(ggplot2)
library(gridExtra)
library(fitdistrplus)

normalize_data <- function(data, m, s) {
  (data - m) / s
}

# Function to generate Q-Q plot with KS test
generate_qqplot <- function(data, dist_name, q_func, p_func, fit_params, method = "mle") {
  # Conditional handling for the t-distribution
  if (dist_name == "t") {
    # Normalize the data using location (m) and scale (s)
    normalized_data <- (data - fit_params$m) / fit_params$s
    
    # Perform KS test
    ks_test <- ks.test(
      normalized_data, 
      function(x) do.call(p_func, c(list(x), list(df = fit_params$df)))
    )
    
    # Generate theoretical quantiles
    theoretical_quantiles <- do.call(q_func, c(list(ppoints(length(normalized_data))), list(df = fit_params$df)))
    
    # Use normalized data for Q-Q plot
    sample_quantiles <- sort(normalized_data)
    xlab_text <- paste("Theoretical Quantiles (", dist_name, ")", sep = "")
    ylab_text <- "Sample Quantiles (Normalized)"
  } else {
    # Perform KS test for other distributions
    ks_test <- ks.test(data, function(x) do.call(p_func, c(list(x), fit_params)))
    
    # Generate theoretical quantiles
    theoretical_quantiles <- do.call(q_func, c(list(ppoints(length(data))), fit_params))
    
    # Use original data for Q-Q plot
    sample_quantiles <- sort(data)
    xlab_text <- paste("Theoretical Quantiles (", dist_name, ")", sep = "")
    ylab_text <- "Sample Quantiles"
  }
  
  # Generate Q-Q plot
  plot <- ggplot() +
    geom_point(aes(x = theoretical_quantiles, y = sample_quantiles), color = "blue") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    ggtitle(paste("Q-Q Plot:", dist_name)) +
    xlab(xlab_text) +
    ylab(ylab_text) +
    theme_minimal() +
    annotate(
      "text", 
      x = min(theoretical_quantiles), 
      y = max(sample_quantiles), 
      label = paste("KS Statistic:", round(ks_test$statistic, 3), "\np =", format(round(ks_test$p.value, 3), nsmall=3)), 
      hjust = 0, vjust = 1
    )
  
  return(list(plot = plot, ks_test = ks_test))
}

# Normalize data for T-distribution
normalize_data <- function(data, m, s) {
  (data - m) / s
}

# Generate the T-distribution Q-Q plot
t_qq <- generate_qqplot_t(abroca_values, t_params)

# Display the plot
print(t_qq$plot)

# Weibull Distribution
fit_weibull <- fitdistr(abroca_values, "weibull") # Shift data to positive
weibull_params <- list(shape = fit_weibull$estimate["shape"], scale = fit_weibull$estimate["scale"])
weibull_qq <- generate_qqplot(abroca_values, "Weibull", qweibull, pweibull, weibull_params)

# Normal Distribution
fit_normal <- fitdistr(abroca_values, "normal")
normal_params <- list(mean = fit_normal$estimate["mean"], sd = fit_normal$estimate["sd"])
normal_qq <- generate_qqplot(abroca_values, "Normal", qnorm, pnorm, normal_params)

# T Distribution
fit_t <- fitdistr(abroca_values, "t")
t_params <- list(m = fit_t$estimate['m'], 
                 s = fit_t$estimate['s'], 
                 df = fit_t$estimate['df'])
t_qq <- generate_qqplot(abroca_values, "t", qt, pt, t_params)

# F Distribution
log_likelihood_f <- function(params, data) {
  df1 <- params[1]
  df2 <- params[2]
  if (df1 <= 0 || df2 <= 0) return(Inf)  # Penalize invalid parameters
  -sum(df(data, df1 = df1, df2 = df2, log = TRUE))
}

# Fit the F-distribution using optim
start_params <- c(df1 = 1, df2 = 1)  
fit_f <- optim(
  par = start_params,
  fn = log_likelihood_f,
  data = adjusted_values,
  method = "L-BFGS-B",
  lower = c(1e-6, 1e-6)  # Constraints for degrees of freedom
)

f_params <- list(df1 = fit_f$par[1], df2 = fit_f$par[2])

f_qq <- generate_qqplot(abroca_values, "F", qf, pf, f_params)

## Arrange all plots in a grid
# Export the plots to a PDF file
pdf("qq_plots-5000-imbalanced-null-hypothesis.pdf", width = 8, height = 6)  # Set the desired dimensions

# Arrange the plots in a grid and render them
grid.arrange(weibull_qq$plot, normal_qq$plot, t_qq$plot, f_qq$plot, ncol = 2)

dev.off()

##### RQ2+RQ3 Randomization test #####

source('functions.R')

library(tidyverse)
library(progress)

# Function to simulate predictions and compute ABROCA p-values
simulate_power <- function(n_majority, n_minority, auc_majority, auc_minority, 
                           ratio_pos_case, test_set_size, n_iter_power=100, n_iter_test = 1000) {
  # Generate predictions
  d <- generate_predictions(
    n_majority = n_majority, 
    n_minority = n_minority,
    auc_majority = auc_majority, 
    auc_minority = auc_minority,
    ratio_pos_case = ratio_pos_case, 
    test_set_size = test_set_size
  )
  
  # Perform randomization test
  abroca_randomization_test <- function(d, n_iter_test) {
    abrocas <- replicate(n_iter_test, {
      d2 <- d %>% mutate(group = sample(group, size = nrow(d), replace = FALSE))
      tryCatch(
        {
          calculate_abroca(d2)
        },
        error = function(e) {
          # Log the error
          write(paste("Error in simulate_power with n_majority =", n_majority, 
                      "n_minority =", n_minority, ": ", e$message, sep = " "), 
                file = "simulation_errors.log", append = TRUE)
          NA # Return NA to indicate skipped iteration
        }
      )
    })
    
    # Filter out NA values and calculate p-value
    valid_abrocas <- abrocas[!is.na(abrocas)]
    if (length(valid_abrocas) == 0) {
      stop("All iterations resulted in errors or NA values.")
    }
    p_value <- mean(calculate_abroca(d) > valid_abrocas, na.rm = TRUE)
    return(p_value)
  }
  
  p_values <- replicate(n_iter_power, {abroca_randomization_test(d, n_iter_test)})
  #print(p_values)
  return(mean(p_values < 0.05)) # Return TRUE if significant
}


# Simulation parameters
sample_sizes <- seq(100, 5000, by = 100) # Varying sample sizes
n_sim_power <- 2 # Number of simulations per sample size

# Track results
results <- data.frame(SampleSize = integer(), Power = numeric())

# Start progress bar
pb <- progress_bar$new(
  format = "Simulating [:bar] :percent (:elapsed/:eta)", total = length(sample_sizes), clear = FALSE
)

# Start timing
start_time <- Sys.time()

# Run simulations
for (n in sample_sizes) {

  pb$tick()
  
  power <- suppressWarnings(suppressMessages(
      simulate_power(n_majority = n / 2, n_minority = n / 2,
                     auc_majority = 0.8, auc_minority = 0.4,
                     ratio_pos_case = 0.5, test_set_size = 0.2,
                     n_iter_power=n_sim_power, n_iter_test = 100)
    ))
  
  results <- rbind(results, data.frame(SampleSize = n, Power = power))
  
  # Save plot to PDF
  pdf_filename <- paste0("Power_Analysis_SampleSize_", n, ".pdf")
  p <- ggplot(results, aes(x = SampleSize, y = Power)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    labs(
      title = "Power Analysis for Varying Sample Sizes",
      x = "Sample Size",
      y = "Power",
      caption = "Power calculated using ABROCA randomization test"
    ) +
    theme_minimal(base_size = 14) +
    scale_y_continuous(labels = scales::percent_format())
  ggsave(filename = pdf_filename, plot = p, device = "pdf", width = 8, height = 6)
}
# End timing
end_time <- Sys.time()
cat("Simulation completed in:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

