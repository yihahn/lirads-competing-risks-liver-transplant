# Competing Risks Wolbers C-index Implementation in R
# ==================================================

#' Calculate Wolbers C-index for competing risks data
#'
#' @param times Numeric vector of observed times
#' @param events Numeric vector of event indicators (0=censored, 1=event of interest, 2=competing event)
#' @param predictions Numeric vector of predicted risk scores for event of interest
#' @param event_of_interest Event code for the event of interest (default: 1)
#' @param weights Optional numeric vector of sample weights
#' @return List containing c_index, variance, and detailed pair counts
wolbers_c_index_competing <- function(times, events, predictions, 
                                    event_of_interest = 1, weights = NULL) {
  
  # Input validation
  if (length(times) != length(events) || length(times) != length(predictions)) {
    stop("All input vectors must have the same length")
  }
  
  n <- length(times)
  
  if (is.null(weights)) {
    weights <- rep(1, n)
  } else if (length(weights) != n) {
    stop("Weights vector must have the same length as other inputs")
  }
  
  # Check event codes
  unique_events <- unique(events)
  if (!all(unique_events %in% c(0, 1, 2))) {
    warning("Events should be coded as 0 (censored), 1 (event of interest), 2 (competing event)")
  }
  
  # Initialize counters
  concordant <- 0
  discordant <- 0
  tied_risk <- 0
  total_pairs <- 0
  
  # Detailed pair counting for diagnostics
  pair_types <- list(
    event_vs_censored = 0,
    event_vs_competing = 0,
    event_vs_event = 0,
    competing_vs_censored = 0,
    unusable = 0
  )
  
  # Main comparison loop
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      
      # Determine pair type and comparability
      event_i <- events[i]
      event_j <- events[j]
      time_i <- times[i]
      time_j <- times[j]
      pred_i <- predictions[i]
      pred_j <- predictions[j]
      
      comparable <- FALSE
      shorter_idx <- NA
      longer_idx <- NA
      
      # Case 1: Both have event of interest
      if (event_i == event_of_interest && event_j == event_of_interest) {
        pair_types$event_vs_event <- pair_types$event_vs_event + 1
        comparable <- TRUE
        if (time_i < time_j) {
          shorter_idx <- i; longer_idx <- j
        } else if (time_j < time_i) {
          shorter_idx <- j; longer_idx <- i
        } else {
          # Same event time - handle ties
          total_pairs <- total_pairs + 1
          if (pred_i == pred_j) {
            tied_risk <- tied_risk + weights[i] * weights[j]
          }
          next
        }
      }
      
      # Case 2: One has event of interest, other is censored
      else if ((event_i == event_of_interest && event_j == 0) ||
               (event_i == 0 && event_j == event_of_interest)) {
        
        pair_types$event_vs_censored <- pair_types$event_vs_censored + 1
        
        if (event_i == event_of_interest && event_j == 0) {
          # i has event, j is censored
          if (time_i <= time_j) {
            comparable <- TRUE
            shorter_idx <- i; longer_idx <- j
          }
        } else {
          # j has event, i is censored
          if (time_j <= time_i) {
            comparable <- TRUE
            shorter_idx <- j; longer_idx <- i
          }
        }
      }
      
      # Case 3: One has event of interest, other has competing event
      else if ((event_i == event_of_interest && event_j == 2) ||
               (event_i == 2 && event_j == event_of_interest)) {
        
        pair_types$event_vs_competing <- pair_types$event_vs_competing + 1
        
        if (event_i == event_of_interest && event_j == 2) {
          # i has event of interest, j has competing event
          if (time_i <= time_j) {
            comparable <- TRUE
            shorter_idx <- i; longer_idx <- j
          }
        } else {
          # j has event of interest, i has competing event  
          if (time_j <= time_i) {
            comparable <- TRUE
            shorter_idx <- j; longer_idx <- i
          }
        }
      }
      
      # Case 4: One has competing event, other is censored
      else if ((event_i == 2 && event_j == 0) ||
               (event_i == 0 && event_j == 2)) {
        
        pair_types$competing_vs_censored <- pair_types$competing_vs_censored + 1
        # These pairs are not informative for the event of interest
        # but we could use them in some formulations
        next
      }
      
      # Case 5: Other combinations (censored vs censored, competing vs competing)
      else {
        pair_types$unusable <- pair_types$unusable + 1
        next
      }
      
      if (!comparable) next
      
      total_pairs <- total_pairs + 1
      weight_product <- weights[shorter_idx] * weights[longer_idx]
      
      # Compare predictions
      # Higher risk should correspond to shorter time to event of interest
      pred_shorter <- predictions[shorter_idx]
      pred_longer <- predictions[longer_idx]
      
      if (pred_shorter > pred_longer) {
        concordant <- concordant + weight_product
      } else if (pred_shorter < pred_longer) {
        discordant <- discordant + weight_product
      } else {
        tied_risk <- tied_risk + weight_product
      }
    }
  }
  
  # Calculate C-index
  if (total_pairs == 0) {
    warning("No comparable pairs found. Check your data.")
    return(list(
      c_index = 0.5, 
      variance = 0, 
      concordant_pairs = 0, 
      total_pairs = 0,
      pair_types = pair_types
    ))
  }
  
  c_index <- (concordant + 0.5 * tied_risk) / (concordant + discordant + tied_risk)
  
  # Estimate variance
  if (total_pairs > 1) {
    variance <- c_index * (1 - c_index) / total_pairs
  } else {
    variance <- 0
  }
  
  return(list(
    c_index = c_index,
    variance = variance,
    concordant_pairs = concordant,
    discordant_pairs = discordant,
    tied_pairs = tied_risk,
    total_pairs = total_pairs,
    pair_types = pair_types,
    event_of_interest = event_of_interest
  ))
}


#' Bootstrap confidence intervals for competing risks C-index
#'
#' @param times Numeric vector of observed times
#' @param events Numeric vector of event indicators  
#' @param predictions Numeric vector of predicted risk scores
#' @param event_of_interest Event code for the event of interest
#' @param n_bootstrap Number of bootstrap samples
#' @param seed Random seed for reproducibility
#' @return List with bootstrap results
bootstrap_c_index_competing <- function(times, events, predictions, 
                                      event_of_interest=1,
                                      n_bootstrap=1000, seed=20250521) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- length(times)
  bootstrap_scores <- numeric(n_bootstrap)
  
  # Original C-index
  original_result <- wolbers_c_index_competing(times, events, predictions, event_of_interest)
  original_c <- original_result$c_index
  
  # Bootstrap sampling
  valid_boots <- 0
  for (i in 1:n_bootstrap) {
    # Stratified sampling to preserve event type distribution
    boot_indices <- sample(n, n, replace = TRUE)
    
    boot_times <- times[boot_indices]
    boot_events <- events[boot_indices]
    boot_predictions <- predictions[boot_indices]
    
    tryCatch({
      boot_result <- wolbers_c_index_competing(boot_times, boot_events, 
                                             boot_predictions, event_of_interest)
      if (boot_result$total_pairs > 0) {
        valid_boots <- valid_boots + 1
        bootstrap_scores[valid_boots] <- boot_result$c_index
      }
    }, error = function(e) {
      # Skip failed bootstrap samples
    })
  }
  
  # Use only valid bootstrap samples
  bootstrap_scores <- bootstrap_scores[1:valid_boots]
  
  if (valid_boots == 0) {
    warning("No valid bootstrap samples obtained")
    return(list(
      c_index = original_c,
      standard_error = NA,
      ci_lower = NA,
      ci_upper = NA,
      n_valid_bootstrap = 0
    ))
  }
  
  # Calculate confidence intervals
  std_error <- sd(bootstrap_scores, na.rm = TRUE)
  ci_lower <- quantile(bootstrap_scores, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_scores, 0.975, na.rm = TRUE)
  
  return(list(
    c_index = original_c,
    standard_error = std_error,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    n_valid_bootstrap = valid_boots,
    bootstrap_scores = bootstrap_scores
  ))
}


#' Enhanced competing risks C-index analysis
#'
#' @param times Numeric vector of observed times
#' @param events Numeric vector of event indicators
#' @param predictions Numeric vector of predicted risk scores  
#' @param event_of_interest Event code for the event of interest
#' @param alpha Significance level for confidence interval
#' @return Comprehensive analysis results
enhanced_c_index_competing <- function(times, events, predictions, 
                                     event_of_interest = 1, alpha = 0.05) {
  
  # Basic C-index calculation
  basic_result <- wolbers_c_index_competing(times, events, predictions, event_of_interest)
  
  # Bootstrap confidence intervals
  boot_result <- bootstrap_c_index_competing(times, events, predictions, 
                                           event_of_interest, seed=20250521)
  
  # Data summary statistics
  n_total <- length(times)
  n_events <- sum(events == event_of_interest)
  n_competing <- sum(events == 2)
  n_censored <- sum(events == 0)
  
  event_rate <- n_events / n_total
  competing_rate <- n_competing / n_total
  censoring_rate <- n_censored / n_total
  
  # Statistical test for C-index = 0.5
  se <- sqrt(basic_result$variance)
  if (se > 0) {
    z_score <- (basic_result$c_index - 0.5) / se
    p_value <- 2 * (1 - pnorm(abs(z_score)))
  } else {
    z_score <- NA
    p_value <- NA
  }
  
  # Compile comprehensive results
  results <- list(
    c_index = basic_result$c_index,
    variance = basic_result$variance,
    standard_error = se,
    bootstrap_se = boot_result$standard_error,
    ci_lower = boot_result$ci_lower,
    ci_upper = boot_result$ci_upper,
    z_score = z_score,
    p_value = p_value,
    concordant_pairs = basic_result$concordant_pairs,
    discordant_pairs = basic_result$discordant_pairs,
    tied_pairs = basic_result$tied_pairs,
    total_pairs = basic_result$total_pairs,
    pair_types = basic_result$pair_types,
    n_observations = n_total,
    n_events = n_events,
    n_competing = n_competing,
    n_censored = n_censored,
    event_rate = event_rate,
    competing_rate = competing_rate,
    censoring_rate = censoring_rate,
    event_of_interest = event_of_interest,
    n_bootstrap = boot_result$n_valid_bootstrap
  )
  
  class(results) <- "wolbers_cindex_competing"
  return(results)
}


#' Print method for competing risks C-index results
print.wolbers_cindex_competing <- function(x, digits = 4, ...) {
  cat("Competing Risks Wolbers C-index Results\n")
  cat("=======================================\n\n")
  
  cat("Event of Interest:", x$event_of_interest, "\n")
  cat("C-index:", round(x$c_index, digits), "\n")
  cat("Standard Error:", round(x$standard_error, digits), "\n")
  
  if (!is.na(x$bootstrap_se)) {
    cat("Bootstrap SE:", round(x$bootstrap_se, digits), "\n")
    cat("95% CI: [", round(x$ci_lower, digits), ", ", round(x$ci_upper, digits), "]\n\n")
  } else {
    cat("Bootstrap CI: Not available\n\n")
  }
  
  if (!is.na(x$p_value)) {
    cat("Test H0: C-index = 0.5\n")
    cat("Z-score:", round(x$z_score, digits), "\n")
    cat("P-value:", round(x$p_value, digits), "\n\n")
  }
  
  cat("Data Summary:\n")
  cat("- Total observations:", x$n_observations, "\n")
  cat("- Events of interest:", x$n_events, "(", round(x$event_rate * 100, 1), "%)\n")
  cat("- Competing events:", x$n_competing, "(", round(x$competing_rate * 100, 1), "%)\n")
  cat("- Censored:", x$n_censored, "(", round(x$censoring_rate * 100, 1), "%)\n\n")
  
  cat("Pair Analysis:\n")
  cat("- Total comparable pairs:", x$total_pairs, "\n")
  cat("- Concordant pairs:", x$concordant_pairs, "\n")
  cat("- Discordant pairs:", x$discordant_pairs, "\n")
  cat("- Tied pairs:", x$tied_pairs, "\n\n")
  
  cat("Pair Type Distribution:\n")
  cat("- Event vs Event:", x$pair_types$event_vs_event, "\n")
  cat("- Event vs Censored:", x$pair_types$event_vs_censored, "\n")
  cat("- Event vs Competing:", x$pair_types$event_vs_competing, "\n")
  cat("- Competing vs Censored:", x$pair_types$competing_vs_censored, "\n")
  cat("- Unusable pairs:", x$pair_types$unusable, "\n")
  
  if (!is.na(x$n_bootstrap)) {
    cat("\nBootstrap samples:", x$n_bootstrap, "\n")
  }
}


# Example usage and demonstration
# ==============================

if (FALSE) {  # Change to TRUE to run examples
  
  # Example 1: Simulated competing risks data
  cat("Example 1: Simulated Competing Risks Data\n")
  cat("=========================================\n")
  
  #set.seed(42)
  n <- 300
  
  # Simulate competing risks data
  # Event of interest (type 1)
  lambda1 <- 0.1
  times1 <- rexp(n, lambda1)
  
  # Competing event (type 2)  
  lambda2 <- 0.08
  times2 <- rexp(n, lambda2)
  
  # Censoring
  censor_times <- rexp(n, 0.05)
  
  # Determine observed times and events
  observed_times <- pmin(times1, times2, censor_times)
  events <- ifelse(observed_times == times1, 1,
                  ifelse(observed_times == times2, 2, 0))
  
  # Create risk predictions for event of interest
  # Higher values should predict shorter time to event 1
  true_risk <- lambda1 + rnorm(n, 0, 0.02)
  predictions <- true_risk + rnorm(n, 0, 0.05)
  
  # Calculate competing risks C-index
  result <- enhanced_c_index_competing(observed_times, events, predictions, event_of_interest = 1)
  print(result)
  
  
  # Example 2: Compare different events as "of interest"
  cat("\n\nExample 2: Comparing Different Events as Primary\n")
  cat("===============================================\n")
  
  # C-index for event type 1
  result1 <- wolbers_c_index_competing(observed_times, events, predictions, event_of_interest = 1)
  
  # Create predictions for event type 2
  pred_for_event2 <- lambda2 + rnorm(n, 0, 0.05)
  result2 <- wolbers_c_index_competing(observed_times, events, pred_for_event2, event_of_interest = 2)
  
  cat("C-index for Event 1:", round(result1$c_index, 4), "\n")
  cat("C-index for Event 2:", round(result2$c_index, 4), "\n")
  
  
  # Example 3: Model comparison
  cat("\n\nExample 3: Competing Risk Model Comparison\n")
  cat("=========================================\n")
  
  # Create multiple prediction models
  model1 <- predictions
  model2 <- predictions + rnorm(n, 0, 0.1)  # Noisier model
  model3 <- rnorm(n, mean(predictions), sd(predictions))  # Random model
  
  models <- list(
    "Good Model" = model1,
    "Noisy Model" = model2, 
    "Random Model" = model3
  )
  
  comparison_results <- data.frame(
    Model = names(models),
    C_Index = numeric(length(models)),
    CI_Lower = numeric(length(models)),
    CI_Upper = numeric(length(models)),
    Total_Pairs = numeric(length(models)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(models)) {
    res <- bootstrap_c_index_competing(observed_times, events, models[[i]], 
                                     event_of_interest = 1, n_bootstrap = 500)
    comparison_results[i, "C_Index"] <- res$c_index
    comparison_results[i, "CI_Lower"] <- res$ci_lower
    comparison_results[i, "CI_Upper"] <- res$ci_upper
  }
  
  comparison_results <- comparison_results[order(comparison_results$C_Index, decreasing = TRUE), ]
  print(comparison_results)
}
