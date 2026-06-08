# ============================================================================
# Main Functions for COXNet and AFTNet
# ============================================================================

#' Negative Partial Log-Likelihood for \code{COXNet} Model
#'
#' @keywords internal
#' @noRd

  nll_COX <- function(
      Y, eta, delta
  ) {
    
    risk <- risk_fun(Y = Y, eta = eta, delta = delta)
    loss <- sum(risk$delta_sorted * (log(risk$w_sorted) - log(risk$W_sorted)))
    nll  <- - loss / length(Y)
    
    return(nll)
  }

# ============================================================================

#' Gradient of the Negative Partial Log-Likelihood for \code{COXNet} Model
#' 
#' @keywords internal
#' @noRd

  gradient_COX <- function(
      X, Y, eta, delta
  ) {
    
    n <- length(Y)
    
    risk   <- risk_fun(Y = Y, eta = eta, delta = delta)
    Y_risk <- outer(Y, risk$Y_sorted, function(ti, tj) as.numeric(ti >= tj))
    w_mat  <- matrix(risk$w, nrow = n, ncol = n)
    W_vec  <- matrix(risk$W_sorted, nrow = n, ncol = n, byrow = TRUE)
    
    pi_mat <- (Y_risk * w_mat) / W_vec
    pi_mat[is.nan(pi_mat) | is.infinite(pi_mat)] <- 0
    
    res <- delta - as.numeric(pi_mat %*% delta)
    score_beta <- t(X) %*% res
    
    grad_beta <- - score_beta / n
    
    return(list(
      pi_mat    = pi_mat,
      grad_beta = grad_beta))
  }

# ============================================================================

#' Penalized Hessian Matrix and Lipschitz Constant for \code{COXNet} Model 
#' 
#' @keywords internal
#' @noRd

  hessian_COX <- function(
      X, Y, eta, delta,
      alpha, lambda,
      L = NULL
  ) {
    
    n <- length(Y)
    
    if (is.null(L)) {
      p <- ncol(X)
      L <- matrix(0, p, p)
    }
    
    grad      <- gradient_COX(X = X, Y = Y, eta = eta, delta = delta)
    delta_mat <- matrix(delta, nrow = n, ncol = n, byrow = TRUE)
    W_diag    <- rowSums(grad$pi_mat * (1 - grad$pi_mat) * delta_mat)
    
    E <- t(X) %*% diag(W_diag) %*% X
    H <- (1/n) * E + lambda * (1 - alpha) * L
    
    eig     <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max.eig <- max(Re(eig))
    
    return(list(
      H = H,
      M = max.eig))
  }

# ============================================================================

#' Risk Set Quantities for \code{COXNet} Model
#' 
#' @keywords internal
#' @noRd

  risk_fun <- function(
      Y, eta, delta
  ) {
    
    eta <- pmin(pmax(eta, -500), 500)
    
    order_Y      <- order(Y)
    Y_sorted     <- Y[order_Y]
    delta_sorted <- delta[order_Y]
    
    w <- exp(eta)
    w_sorted <- w[order_Y]
    
    w_rev <- rev(w_sorted)
    W_rev <- cumsum(w_rev)
    
    W_sorted <- rev(W_rev)
    W_sorted[W_sorted == 0] <- .Machine$double.eps
    
    return(list(
      Y_sorted     = Y_sorted,
      delta_sorted = delta_sorted,
      w            = w,
      w_sorted     = w_sorted,
      W_sorted     = W_sorted))
  }

# ============================================================================
# ============================================================================

#' Negative Log-Likelihood for \code{AFTNet} Model Under Several Distributions
#' 
#' @keywords internal
#' @noRd

  nll_AFT <- function(
      Y, eta, delta,
      sigma,
      dist = NULL
  ) {
    
    dist <- match.arg(dist, choices = c("weibull", "lognormal", "loglogistic"))
    n    <- length(Y)
    
    res <- as.numeric((Y - eta) / sigma)
    res <- pmin(pmax(res, -30), 30)
    
    loss <- switch(dist,
                   
                   "weibull" = as.numeric(crossprod(delta, (-log(sigma) + res)) - sum(exp(res))),
                   
                   "lognormal" = {
                     surv  <- pmax(1 - stats::pnorm(res), 1e-16)
                     term1 <- -0.5 * log(2 * pi * sigma^2) * sum(delta)
                     term2 <- -0.5 * as.numeric(crossprod(delta, res^2))
                     term3 <- as.numeric(crossprod(1 - delta, log(surv)))
                     term1 + term2 + term3},
                   
                   "loglogistic" = {
                     exp_res <- exp(res)
                     term1 <- as.numeric(crossprod(delta, (-log(sigma) + res)))
                     term2 <- as.numeric(crossprod(1 + delta, log1p(exp_res)))
                     term1 - term2},
                   
                   stop("Invalid distribution specified."))
    
    nll <- - loss / n
    
    return(nll)
  }

# ============================================================================

#' Gradient of the Negative Log-Likelihood for \code{AFTNet} Model Under Several Distributions
#' 
#' @keywords internal
#' @noRd

  gradient_AFT <- function(
      X, Y, eta,
      delta, sigma,
      dist = NULL
  ) {
    
    dist <- match.arg(dist, choices =  c("weibull", "lognormal", "loglogistic"))
    n    <- length(Y)
    
    res <- as.numeric((Y - eta) / sigma)
    
    exp_res <- exp(res)
    
    a <- switch(dist,
                
                "weibull" = exp_res - delta,
                
                "lognormal" = {
                  dens <- stats::dnorm(res)
                  surv <- pmax(1 - stats::pnorm(res), 1e-16)
                  (delta * res) + ((1 - delta) * dens / surv)},
                
                "loglogistic" = {
                  (1 + delta) * exp_res / (1 + exp_res) - delta},
                
                stop("Invalid distribution specified."))
    
    grad_beta <- as.vector(- (1 / (n * sigma)) * crossprod(X, a))
    
    return(grad_beta)
  }

# ============================================================================

#' Penalized Hessian Matrix and Lipschitz Constant for \code{AFTNet} Model
#' 
#' @keywords internal
#' @noRd

  hessian_AFT <- function(
      X, Y, eta,
      delta, sigma,
      dist = NULL,
      alpha, lambda,
      L = NULL
  ) {
    
    dist <- match.arg(dist, choices = c("weibull", "lognormal", "loglogistic"))
    n    <- length(Y)
    
    if (is.null(L)) {
      p <- ncol(X)
      L <- matrix(0, p, p)
    }
    
    res <- as.numeric((Y - eta) / sigma)
    res <- pmin(pmax(res, -30), 30)
    
    exp_res <- exp(res)
    
    E <- switch(dist,
                
                "weibull" = diag(exp_res),
                
                "lognormal" = {
                  dens <- stats::dnorm(res)
                  surv <- pmax(1 - stats::pnorm(res), 1e-16)
                  diag(delta + (1 - delta) * (dens * (- res * surv + dens)) / surv^2)},
                
                "loglogistic" = {
                  diag((1 + delta) * exp_res / ((1 + exp_res)^2))},
                
                stop("Invalid distribution specified."))
    
    H <- (1 / n) * crossprod(X, E %*% X) + lambda * (1 - alpha) * L
    
    eig     <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max.eig <- max(Re(eig))
    
    return(list(
      H = H,
      M = max.eig))
  }

# ============================================================================
# ============================================================================

#' LASSO Proximal Operator (Soft-Thresholding)
#' 
#' @keywords internal
#' @noRd

  prox_l1 <- function(
      u, mu
  ) {
    
    uhat  <- abs(u) - mu
    ubind <- cbind(rep(0, length(u)), uhat)
    prox  <- sign(u) * apply(ubind, 1, max)
    
    return(prox)
  }

# ============================================================================
# ============================================================================
  
#' Matrix standardization
#'
#' @keywords internal
#' @noRd
  
  standardize <- function(mat) {
    
    scaled_mat <- scale(mat)
    
    return(list(
      std   = scaled_mat,
      means = attr(scaled_mat, "scaled:center"),
      sds   = attr(scaled_mat, "scaled:scale")))
  }

# ============================================================================
# ============================================================================  
  
#' Performance Metrics for Survival Models
#'
#' @description
#' Computes a variety of performance metrics for survival model 
#' supporting both real-data evaluation and simulation studies.
#'
#' @param Y_train Numeric vector of observed training survival times
#'        (log-transformed under \code{"AFTNet"}).
#' @param delta_train Integer vector of training censoring indicators 
#'        (1 = event, 0 = censored).
#' @param X_test Numeric matrix of testing covariates standardized using the training data.
#' @param Y_test Numeric vector of observed testing survival times 
#'        (log-transformed under \code{"AFTNet"}).
#' @param delta_test Integer vector of testing censoring indicators
#'        (1 = event, 0 = censored).
#' @param beta_est Numeric vector of estimated regression coefficients obtained from the training set.
#' @param beta_true Optional numeric vector of true regression coefficients.
#'                  Required only for simulation-based metrics (FPR, FNR, PMSE).
#' @param model Character string specifying the fitted survival model
#'              (\code{"COXNet"}, or \code{"AFTNet"}).
#' @param p_active Integer scalar specifying the number of truly active covariates,
#'                 required only when \code{metrics} includes \code{"FPR"} or \code{"FNR"}
#'                 and \code{beta_true} is supplied.
#' @param times_auc Optional numeric vector of time points at which the time-dependent AUC is evaluated.
#'                  If \code{NULL} (default), empirical quantiles of \code{Y_test} are used.
#' @param metrics Character vector specifying the performance measures to compute.
#'                Allowed values:
#'                \itemize{
#'                  \item \code{"PredRisk"} - Predicted Risk or expected survival time,
#'                  \item \code{"CIndex"} - Harrell's concordance index,
#'                  \item \code{"FPR"} - False Positive Rate,
#'                  \item \code{"FNR"} - False Negative Rate,
#'                  \item \code{"NSR"} - Number of Selected variables Rate,
#'                  \item \code{"PMSE"} - Predictive Mean Square Error,
#'                  \item \code{"AUC"} - time-dependent AUC.
#'                }
#'
#' @details
#' The predicted quantity depends on the model type:
#' \itemize{
#'   \item For \code{COXNet}, \code{PredRisk} is the hazard ratio.
#'   \item For \code{AFTNet}, \code{PredRisk} is proportional to the expected survival time.
#' }
#' 
#' Harrell's concordance index is computed using \code{\link[Hmisc]{rcorr.cens}}.
#' The time-dependent AUC is computed using Uno's estimator via
#' \code{\link[survAUC]{AUC.uno}} at the specified time points.
#' 
#' The metrics \code{FPR}, \code{FNR}, and \code{PMSE} are defined only in
#' simulation settings because they require knowledge of the true regression
#' coefficients. When \code{beta_true} is not provided, these metrics are
#' returned as \code{NA} if requested. 
#' All other metrics can be computed for both simulated and real datasets.
#' 
#' @return A named list containing the requested performance metrics.
#' 
#' @note Scalar metrics are returned as numeric values, \code{PredRisk} as
#' a numeric vector of predicted risk scores, and time-dependent AUC values
#' as separate list elements with names of the form \code{"AUC_t_<time>"}.
#'
#' @name Metrics
#'
#' @export
  
  Metrics <- function(
      Y_train     = NULL,
      delta_train = NULL,
      X_test      = NULL,
      Y_test      = NULL,
      delta_test  = NULL,
      beta_est,
      beta_true   = NULL,
      model       = NULL,
      p_active    = NULL,
      times_auc   = NULL,
      metrics     = NULL
  ) {
    
  # --- Input checks ---
    
    model   <- match.arg(model, choices = c("COXNet", "AFTNet"))
    
    allowed_metrics <- c("PredRisk", "CIndex", "FPR", "FNR", "NSR", "PMSE", "AUC")
    
    metrics <- unique(metrics)
    
    if (!all(metrics %in% allowed_metrics)) {
      stop("Invalid values in 'metrics'.")
    }
    
    p <- length(beta_est)
    
    results <- list()
    
  # -----  
    
    lp_est <- as.vector(X_test %*% beta_est)
    
    if (model == "COXNet") {
      Surv_train <- survival::Surv(Y_train, delta_train)
      Surv_test  <- survival::Surv(Y_test, delta_test)
      time_var   <- Y_test
    } else {
      Surv_train <- survival::Surv(exp(Y_train), delta_train)
      Surv_test  <- survival::Surv(exp(Y_test), delta_test)
      time_var   <- exp(Y_test)
    }
    
  # --- Predicted Risk ---
    
    if ("PredRisk" %in% metrics) {
      if (model == "COXNet") {
        results$PredRisk <- exp(lp_est)    # hazard ratio
      } else if (model == "AFTNet") {
        results$PredRisk <- exp(-lp_est)   # expected time
      }}
    
  # --- C-index ---
    
    if ("CIndex" %in% metrics) {
      lp_for_cindex <- if(model == "COXNet") -lp_est else lp_est
      results$CIndex <- as.numeric(Hmisc::rcorr.cens(lp_for_cindex, Surv_test)["C Index"])
    }
  
  # --- Simulation-only metrics ---
    
    if (!is.null(beta_true)) {
      
      lp_true <- as.vector(X_test %*% beta_true)
      
    # - FNR -
      
      if ("FNR" %in% metrics) {
        if (!is.null(p_active)) {
          results$FNR <- sum(beta_est[beta_true != 0] == 0) / p_active
        } else {
          results$FNR <- NA
          stop("'p_active' must be provided when computing FNR.")
        }}
      
    # - FPR -
      
      if ("FPR" %in% metrics) {
        if (!is.null(p_active)) {
          denom <- p - p_active
          results$FPR <- if (denom > 0) sum(beta_est[beta_true == 0] != 0) / denom else NA
        } else {
          results$FPR <- NA
          stop("'p_active' must be provided when computing FPR.")
        }}
      
    # - PMSE -
      
      if ("PMSE" %in% metrics) {
        results$PMSE <- (sum((lp_true - lp_est)^2))/n
      }
      
    } else {
      for (m in c("FNR", "FPR", "PMSE")) {
        if (m %in% metrics) results[[m]] <- NA}
    }
    
  # --- NSR ---
    
    if ("NSR" %in% metrics) {
      results$NSR <- sum(beta_est != 0)/p
    }
    
  # --- AUC(t) ---
    
    if ("AUC" %in% metrics) {
      
      if (is.null(times_auc)) {
        times_auc <- stats::quantile(time_var, probs = c(0.25, 0.5, 0.75))
      }
      
      lp_for_auc  <- if (model == "COXNet") lp_est else -lp_est
      
      AUC_t <- sapply(times_auc, function(t) {
        tryCatch({
          survAUC::AUC.uno(Surv_train, Surv_test, lp_for_auc, times = t)$auc
        }, error = function(e) NA)})
      names(AUC_t) <- paste0("AUC_t_", round(times_auc, 2))
      results      <- c(results, as.list(AUC_t))
    }
    
    return(results)
  }
