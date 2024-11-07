#' Sample Size Calculation for Case-Control Studies
#' 
#' @param alpha Significance level
#' @param power Desired power (optional if n_cases is provided)
#' @param n_cases Number of cases (optional if power is provided)
#' @param p0 Proportion of exposure in controls (required if method = "rho")
#' @param OR Odds Ratio to detect
#' @param CI_upper Upper confidence interval limit of OR
#' @param CI_lower Lower confidence interval limit of OR (optional if SE is provided)
#' @param SE Standard error of coefficient (optional if CI is provided)
#' @param n_previous Sample size of previous study (required for method "se")
#' @param r Ratio of controls to cases, default is 1
#' @param method Calculation method: "rho" or "se" (standard error)
#' @param rho_values Vector of rho values for sample size calculation (if method = "rho")
#' @return A data frame with results according to chosen method
#' @export
SampleCCNoMatchedMultiENG <- function(alpha, power = NULL, n_cases = NULL,
                                    p0 = NULL, OR,
                                    CI_upper = NULL, CI_lower = NULL,
                                    SE = NULL, n_previous = NULL,
                                    r = 1, 
                                    method = "rho",
                                    rho_values = seq(0, 0.9, by = 0.1)) {
  
  # Verifications according to method
  if (method == "rho") {
    if (is.null(p0)) {
      stop("For 'rho' method, p0 is required")
    }
  } else if (method == "se") {
    if (is.null(SE) && (is.null(CI_upper) || is.null(CI_lower))) {
      stop("For 'se' method, SE or both CI limits are required")
    }
    if (is.null(n_previous)) {
      stop("For 'se' method, n_previous is required")
    }
  } else {
    stop("Method must be 'rho' or 'se'")
  }
  
  # If CI provided, calculate SE
  if (!is.null(CI_upper) && !is.null(CI_lower)) {
    SE <- (log(CI_upper) - log(OR)) / qnorm(0.975)
    cat("Calculated Standard Error:", SE, "
")
  }
  
  if (method == "rho") {
    p1 <- (OR * p0) / (1 - p0 + (OR * p0))
    
    z_alpha <- qnorm(1 - alpha/2)
    if (!is.null(power)) {
      z_beta <- qnorm(power)
      
      calculate_n <- function(rho) {
        numerator <- (z_alpha * sqrt((1 + 1/r) * p0 * (1 - p0)) + 
                       z_beta * sqrt(p1 * (1 - p1) + p0 * (1 - p0) / r))^2
        denominator <- (p1 - p0)^2 * (1 - rho^2)
        
        n_cases <- ceiling(numerator / denominator)
        n_controls <- ceiling(n_cases * r)
        
        return(c(n_cases, n_controls))
      }
      
      results <- data.frame(
        rho = rho_values,
        n_cases = sapply(rho_values, function(rho) calculate_n(rho)[1]),
        n_controls = sapply(rho_values, function(rho) calculate_n(rho)[2])
      )
      
      results$total_size <- results$n_cases + results$n_controls
      
    } else {
      # Power calculation with rho method
      calculate_power <- function(rho) {
        numerator <- (p1 - p0)^2 * (1 - rho^2)
        denominator <- sqrt((1 + 1/r) * p0 * (1 - p0))
        z_beta <- sqrt(n_cases * numerator) / denominator - z_alpha
        return(pnorm(z_beta))
      }
      
      results <- data.frame(
        rho = rho_values,
        power = sapply(rho_values, calculate_power)
      )
    }
  } else {
    # Standard error based method
    z_alpha <- qnorm(1 - alpha/2)
    if (!is.null(power)) {
      z_gamma <- qnorm(power)
      n <- ceiling((z_alpha + z_gamma)^2 * n_previous * SE^2 / (log(OR))^2)
      n_cases <- ceiling(n / (1 + r))
      n_controls <- ceiling(n_cases * r)
      results <- data.frame(
        n_cases = n_cases,
        n_controls = n_controls,
        total_size = n_cases + n_controls
      )
    } else {
      z_gamma <- sqrt(n_cases * (log(OR))^2 / (n_previous * SE^2)) - z_alpha
      power <- pnorm(z_gamma)
      results <- data.frame(
        power = power,
        n_cases = n_cases,
        n_controls = n_cases * r,
        total_size = n_cases * (1 + r)
      )
    }
  }
  
  return(results)
}

#' Logistics for Case-Control Study
#'
#' @param n_cases Number of required cases
#' @param n_controls Number of required controls
#' @param case_identification_rate Rate of case identification per day
#' @param control_selection_rate Rate of control selection per day
#' @param case_rejection_rate Expected rejection rate for cases
#' @param control_rejection_rate Expected rejection rate for controls
#' @param working_days_month Number of working days per month
#' @return A list with study logistics calculations
#' @export
case_control_study_logistics <- function(n_cases, 
                                       n_controls,
                                       case_identification_rate,
                                       control_selection_rate,
                                       case_rejection_rate,
                                       control_rejection_rate,
                                       working_days_month) {
  
  cases_to_contact <- n_cases / (1 - case_rejection_rate)
  controls_to_contact <- n_controls / (1 - control_rejection_rate)
  
  days_for_cases <- ceiling(cases_to_contact / case_identification_rate)
  days_for_controls <- ceiling(controls_to_contact / control_selection_rate)
  
  total_days <- max(days_for_cases, days_for_controls)
  duration_months <- total_days / working_days_month
  
  results <- list(
    required_cases = n_cases,
    required_controls = n_controls,
    cases_to_contact = ceiling(cases_to_contact),
    controls_to_contact = ceiling(controls_to_contact),
    estimated_days = total_days,
    estimated_months = round(duration_months, 2)
  )
  
  # Print summary
  cat("
Study logistics summary:
")
  cat("----------------------------------------
")
  cat("Required cases:", n_cases, "
")
  cat("Required controls:", n_controls, "
")
  cat("Cases to contact:", ceiling(cases_to_contact),
      "(", case_rejection_rate*100, "% rejection)
")
  cat("Controls to contact:", ceiling(controls_to_contact),
      "(", control_rejection_rate*100, "% rejection)
")
  cat("Days needed:", total_days, "
")
  cat("Months needed:", round(duration_months, 2),
      "(", working_days_month, "working days per month)
")
  cat("----------------------------------------
")
  
  return(invisible(results))
}
