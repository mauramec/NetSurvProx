#' Proximal Gradient Descent for \code{COXNet} and \code{AFTNet}
#'
#' @description
#' Estimate the regression coefficients in \code{COXNet} and \code{AFTNet} models
#' using a proximal gradient descent algorithm. The objective function combines
#' the normalized negative (partial) log-likelihood with an \eqn{\ell_1} penalty, 
#' and a Laplacian regularization term.
#'
#' @param X Numeric matrix of standardized covariates.
#' @param Y Numeric vector of observed survival times (log-transformed under \code{AFTNet}).
#' @param delta Integer vector of censoring indicators (1 = event, 0 = censored).
#' @param L Optional positive semi-definite, symmetric, and diagonally dominant
#'          Laplacian matrix encoding prior network information
#'          (see \code{\link{CreateNetwork}} for details). If \code{NULL},
#'          no network-based penalization is applied.
#' @param beta0 Numeric vector of initial regression coefficients.
#' @param alpha Numeric parameter controlling the convex combination of the two
#'              penalty terms (value in \code{[0,1]}).
#' @param lambda Non-negative regularization parameter.
#' @param model Character string specifying the fitted survival model
#'              (\code{"COXNet"}, or \code{"AFTNet"}).
#' @param dist Character string specifying the error distribution in 
#'             \code{AFTNet} model. Must be one of \code{"weibull"},
#'             \code{"lognormal"}, or \code{"loglogistic"}.
#' @param sigma Positive numeric scalar representing the scale parameter of the
#'              error distribution in \code{AFTNet} model.
#' @param value Numeric scalar greater than 1 specifying the multiplicative
#'              factor used to increase the step-size constant during
#'              backtracking line search (\code{default: 2}).
#' @param niter Maximum number of iterations (\code{default: 1000}).
#' @param conv Convergence tolerance (\code{default: 1e-3}).
#'
#' @details
#' The algorithm minimizes the objective function:
#' \deqn{\mathcal{L}(\beta) = - \frac{1}{n} \ell(\beta) + \lambda\alpha \|\beta\|_1 +
#'  \lambda(1-\alpha)\beta^\top \mathbf{L} \beta}
#' where \eqn{\ell(\beta)} is the log-likelihood (partial for \code{COXNet}),
#' \eqn{\|\beta\|_1} is the LASSO penalty, \eqn{\beta^\top \mathbf{L} \beta} is
#' the Laplacian constraint.
#'
#' At each iteration the method performs the backtracking line search to enforce
#' a sufficient decrease condition, the gradient step size adaptation (initialized
#' as Lipschitz constant), and an early stopping based on relative change in objective function.
#'
#' Convergence is reached when either the maximum number of iterations is attained,
#' or the relative change in the objective function between consecutive iterations
#' falls below the specific tolerance \code{conv}.
#'
#' @return A list with the following components
#' \itemize{
#'   \item \code{beta}: numeric vector of estimated regression coefficients.
#'   \item \code{objective}: numeric scalar, the final value of the objective function.
#'   \item \code{iterations}: number of iterations performed until convergence
#'   (or until the maximum number of iterations \code{niter} is reached).
#' }
#'
#' @export

  ProxGDNet <- function(
      X, Y, delta,
      L     = NULL,
      beta0, alpha, lambda,
      model = NULL,
      dist  = NULL,
      sigma = NULL,
      value = 2,
      niter = 1000,
      conv  = 1e-3
  ) {
    
  # --- Input checks ---  
    
    if (!is.matrix(X)) stop("'X' must be a numeric matrix.")
    
    if (length(Y) != nrow(X)) stop("'Y' must have length equal to nrow(X).")
    
    if (length(delta) != nrow(X)) stop("'delta' must have length equal to nrow(X).")
    if (!all(delta %in% c(0,1))) stop("delta must be 0/1.")
    
    if (length(beta0) != ncol(X)) stop("'beta0' must have length equal to ncol(X).")
    
    if (alpha < 0 || alpha > 1) stop("'alpha' must be in [0, 1].")
    if (lambda < 0) stop("'lambda' must be non-negative.")
    
    model <- match.arg(model, choices = c("COXNet", "AFTNet"))
    
    if (model == "AFTNet" && (is.null(dist) || is.null(sigma)))
      stop("'dist' and 'sigma' must be provided for AFTNet model.")
    
    if (model == "AFTNet" && sigma <= 0) stop("'sigma' must be positive for AFTNet.")
    
    dist  <- if (model == "AFTNet") match.arg(dist, choices = c("weibull", "lognormal", "loglogistic")) else NULL
    
    if (!is.numeric(value) || length(value) != 1 || value <= 1) stop("'value' must be a numeric scalar greater than 1.")
    
  # ----- 
    
    if (is.null(L)) {
      p <- ncol(X)
      L_mat <- matrix(0, p, p)
    } else {
      L_mat <- L
    }
    
    lambda_alpha   <- lambda * alpha
    lambda_1_alpha <- lambda * (1 - alpha)
    
    beta_current <- beta0
    eta <- as.vector(X %*% beta_current)
    
    get_nll <- function(eta_val, b_val) {
      
      ll <- if(model == "COXNet") {
        nll_COX(Y = Y, eta = eta_val, delta = delta)
      } else {
        nll_AFT(Y = Y, eta = eta_val, delta = delta, sigma = sigma, dist = dist)
      }
      
      quad <- drop(t(b_val) %*% L_mat %*% b_val)
      
      return(list(
        ll   = ll,
        quad = quad))
    }
    
    init_vals   <- get_nll(eta_val = eta, b_val = beta_current)
    obj_current <- init_vals$ll + lambda_alpha * sum(abs(beta_current)) + lambda_1_alpha * init_vals$quad
    fun_old     <- init_vals$ll + lambda_1_alpha * init_vals$quad
    
    M <- if(model == "COXNet") {
      hessian_COX(X = X, Y = Y, eta = eta, delta = delta, alpha = alpha,
                  lambda = lambda, L = L_mat)$M
    } else {
      hessian_AFT(X = X, Y = Y, eta = eta, delta = delta, sigma = sigma,
                  dist = dist, alpha = alpha, lambda = lambda, L = L_mat)$M
    }
    
  # - Main optimization loop -
    
    for (k in 1:niter) {
      
      grad_ll <- if(model == "COXNet") {
        gradient_COX(X = X, Y = Y, eta = eta, delta = delta)$grad_beta
      } else {
        gradient_AFT(X = X, Y = Y, eta = eta, delta = delta, sigma = sigma, dist = dist)
      }
      
      L_beta    <- L_mat %*% beta_current
      full_grad <- grad_ll + 2 * lambda_1_alpha * L_beta
      
    # - Line Search (Backtracking) -
      
      repeat {
        
        u_beta    <- beta_current - (1/M) * full_grad
        beta_next <- prox_l1(u_beta, lambda_alpha / M)
        
        eta_next <- as.vector(X %*% beta_next)
        new_vals <- get_nll(eta_val = eta_next, b_val = beta_next)
        fun_new  <- new_vals$ll + lambda_1_alpha * new_vals$quad
        
        beta_diff <- beta_next - beta_current
        
        if (fun_new <= fun_old + drop(t(full_grad) %*% beta_diff) + (M/2) * sum(beta_diff^2)) {
          break
        }
        
        M <- value * M
        
        if (M > 1e15) break
      }
      
      obj_new    <- fun_new + lambda_alpha * sum(abs(beta_next))
      rel_change <- abs(obj_new - obj_current) / (abs(obj_current) + 1e-8)
      
      beta_current <- beta_next
      eta          <- eta_next
      obj_current  <- obj_new
      fun_old      <- fun_new
      
      if (rel_change <= conv) break
    }
    
    return(list(
      beta       = beta_current,
      objective  = obj_current,
      iterations = k))
  }
  