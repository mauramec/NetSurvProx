#' NetSurvProx Complete Routine
#'
#' @description
#' Fits network-constrained penalized survival models (\code{COXNet} and \code{AFTNet})
#' to identify prognostic signature genes and build a Prognostic Index (PI).
#' The model is trained on a training dataset by incorporating both Laplacian
#' constraints and LASSO regularization, with optional feature standardization.
#' The tuning parameters are jointly selected through cross-validation.
#' An optimal cutoff for the PI is estimated from the training data to enable
#' prognostic stratification. Predictive performance is subsequently evaluated
#' on an independent testing dataset. Model assessment includes survival curve
#' analyses and visualization. Predictive accuracy is quantified using selected metrics.
#'
#' @param X_train Numeric matrix of training covariates standardized
#'                (possibly screened using \code{screen_vars},
#'                see \code{\link{VariableScreening}}).
#' @param Y_train Numeric vector of observed training survival times (log-transformed under \code{AFTNet}).
#' @param delta_train Integer vector of training censoring indicators (1 = event, 0 = censored).
#' @param X_test Numeric matrix of testing covariates.
#' @param Y_test Numeric vector of observed testing survival times (log-transformed under \code{AFTNet}).
#' @param delta_test Integer vector of testing censoring indicators (1 = event, 0 = censored).
#' @param L Optional positive semi-definite, symmetric, and diagonally dominant
#'          Laplacian matrix encoding prior network information (see \code{\link{CreateNetwork}}).
#'          If \code{NULL}, no network-based penalization is applied.
#' @param standardize_train Logical value indicating whether to standardize the training matrix:
#'                          if \code{TRUE} (default), each column is centered to have
#'                          mean 0 and scaled to have unit variance,
#'                          if \code{FALSE}, the matrix is assumed pre-standardized by the user.
#' @param standardize_test Logical value indicating whether to standardize \code{X_test}
#'                    with respect to \code{X_train} (\code{default: TRUE}).
#' @param model Character string specifying the fitted survival model
#'              (\code{"COXNet"}, or \code{"AFTNet"}).
#' @param dist Character string specifying the \code{AFTNet} distribution.
#'             Must be one of \code{"weibull"}, \code{"lognormal"}, or \code{"loglogistic"}.
#' @param select_lambda Logical value, if \code{TRUE} (default) uses \code{lambda.min},
#'                      otherwise \code{lambda.1se}.
#' @param alpha_grid Numeric vector specifying the candidate values for \eqn{\alpha}
#'                   in \code{[0,1]} (\code{default: c(0.3, 0.5, 0.7)}).
#' @param nlambda Numeric value specifying the number of candidate values for
#'                \eqn{\lambda} in the grid (\code{default: 50}).
#' @param lambda_ratio Numeric value giving the ratio of minimum to maximum
#'                     \eqn{\lambda} in the grid (\code{default: 0.01}).
#' @param nfolds Numeric value of folds performed for tuning optimal parameters (\code{default: 5}).
#' @param method Character string specifying the cutoff selection method
#'               (\code{"median"} or \code{"minpvalue"}, see \code{\link{OptimalPICutoff}}).
#' @param probs Vector of probabilities used when \code{method = "minpvalue"} to generate
#'              candidate cutoffs based on quantiles of the PI (\code{default: probs = seq(0.25, 0.80, by = 0.05)}).
#' @param cutoffplot Logical value indicating whether survival curves should be produced
#'                   (\code{default: FALSE}).
#' @param seed Random seed for reproducibility (\code{default: 2026}).
#' @param value Numeric scalar greater than 1 specifying the multiplicative
#'              factor used to increase the step-size constant during
#'              backtracking line search in \code{\link{ProxGDNet}} (\code{default: 2}).
#' @param niter Maximum number of iterations for \code{\link{ProxGDNet}} (\code{default: 1000}).
#' @param conv Convergence tolerance for ProxGDNet (\code{default: 1e-3}).
#' @param parallel_cv  Logical value whether to use parallel processing for
#'                     \code{\link{CvNet}} (\code{default: TRUE}).
#' @param plotCV Logical value indicating whether CV curves should be shown
#'               (\code{default: FALSE}).
#' @param colors_pcv Optional named list of colors for CV plot (see \code{\link{CvNet}}).
#' @param errorbar Logical value, if \code{TRUE} the CV plot includes vertical error
#'                 bars representing 1se of the CV error (\code{default: FALSE}).
#' @param ncore_max Maximum number of cores for parallel processing over CV (\code{default: 5}).
#' @param p_active Numeric value indicating the number of truly active covariates
#'                 (required for FPR/FNR computation in simulation settings).
#' @param times_auc Numeric vector of time points for time-dependent AUC.
#'                  If \code{NULL} (default), quantiles of \code{Y_test} are used.
#' @param beta_true Numeric vector of true coefficients (used only for simulated data).
#' @param metrics Character vector specifying performance \code{\link{Metrics}} to compute.
#'                For real datasets: \code{"CIndex"}, \code{"NSR"}, \code{"AUC"}.
#'                For simulated datasets (in addition): \code{"FPR"}, \code{"FNR"}, \code{"PMSE"}.
#' @param verbose Logical value, if \code{TRUE} progress messages are
#'                printed (\code{default: FALSE}).
#' @param palette Optional character vector of length 2 specifying colors used
#'                for the survival curves. For \code{"COXNet"}, colors correspond to
#'                high- and low-risk groups. For \code{"AFTNet"}, colors correspond to
#'                short- and long-survival groups. If \code{NULL}, default colors are used.
#' @param plot_test Logical value, if \code{TRUE} returns the combined survival
#'                  plot with validation results (\code{default: FALSE}).
#' 
#' @return An object of class \code{NetSurvProx} containing:
#' \itemize{
#'    \item \code{fit_training}: training results (see \code{\link{NetSurvProx_Training}}).
#'    \item \code{fit_testing}: testing results (see \code{\link{NetSurvProx_Testing}}).
#' }
#'
#' @examples
#'   \donttest{
#'   
#'     # - Simulate 40 TFs, each regulating 10 targets with a independent structure -
#'   
#'     targets <- 10
#'     
#'     n <- 165
#'   
#'     simul_data <- Simulations(
#'         n = n, r = 40, targets = targets, p_active = 40,
#'         rho = 0.70, rate = 0.50, b_true = c(0.8, 1.2, -1.2, -0.8),
#'         nsimul = 1, model = "AFTNet", baseline = "lognormal",
#'         sigma_true = 1, shared_scheme = NULL, choice = 1,
#'         save = FALSE, save_path = NULL, seed = 2026, verbose = TRUE)
#'
#'     X     <- simul_data$X_list[[1]]
#'     Y     <- simul_data$time_list[[1]]   # generated in log-scale
#'     delta <- simul_data$delta_list[[1]]
#'     L     <- simul_data$L_list[[1]]
#'
#'     beta_true <- as.vector(unlist(simul_data$beta))
#'
#'   #  - Split the dataset (training/testing sets) -
#'   
#'     set.seed(2026)
#'     
#'     train_idx <- sample(seq_len(n), size = floor(0.7 * n))
#'
#'     X_train     <- X[train_idx,]
#'     Y_train     <- Y[train_idx]
#'     delta_train <- delta[train_idx]
#'
#'     X_test     <- X[-train_idx,]
#'     Y_test     <- Y[-train_idx]
#'     delta_test <- delta[-train_idx]
#'
#'   # - Fitting LogNormal AFTNet -
#'
#'     out <- NetSurvProx(
#'                 X_train, Y_train, delta_train, X_test, Y_test, delta_test,
#'                 L = L, standardize_train = TRUE, standardize_test = TRUE,
#'                 model = "AFTNet", dist = "lognormal", select_lambda = TRUE,
#'                 alpha_grid = 0.5, nlambda = 50, lambda_ratio = 0.1,
#'                 nfolds = 5, method = "minpvalue", probs = seq(0.25, 0.80, by = 0.05),
#'                 cutoffplot = FALSE, seed = 2026, value = 2, niter = 1000, conv = 1e-3,
#'                 parallel_cv = FALSE, plotCV = FALSE, colors_pcv = NULL, errorbar = FALSE, 
#'                 ncore_max = 1, p_active = 40, times_auc = NULL, beta_true = beta_true,
#'                 metrics = "CIndex", verbose = FALSE, palette = NULL, plot_test = FALSE)
#'   
#'   # - Results -
#'   
#'     data.frame(out$fit_testing$performance)
#'   }
#'   
#' @name NetSurvProx    
#' 
#' @export

  NetSurvProx <- function(
      X_train, Y_train, delta_train,
      X_test, Y_test, delta_test,
      L                 = NULL,
      standardize_train = TRUE,
      standardize_test  = TRUE,
      model             = NULL,
      dist              = NULL,
      select_lambda     = TRUE,
      alpha_grid        = c(0.3, 0.5, 0.7),
      nlambda           = 50,
      lambda_ratio      = 0.01,
      nfolds            = 5,
      method            = NULL,
      probs             = seq(0.25, 0.80, by = 0.05),
      cutoffplot        = FALSE,
      seed              = 2026,
      value             = 2,
      niter             = 1000,
      conv              = 1e-3,
      parallel_cv       = TRUE,
      plotCV            = FALSE,
      colors_pcv        = NULL,
      errorbar          = FALSE,
      ncore_max         = 5,
      p_active          = NULL,
      times_auc         = NULL,
      beta_true         = NULL,
      metrics           = NULL,
      verbose           = FALSE,
      palette           = NULL,
      plot_test         = FALSE
  ) {
    
    model  <- match.arg(model, choices = c("COXNet", "AFTNet"))
    method <- match.arg(method, choices = c("median", "minpvalue"))
    
  # --- Training phase ---
    
    fit_training <- NetSurvProx_Training(
      X_train = X_train, Y_train = Y_train, delta_train = delta_train,
      L = L, model = model, dist = dist, select_lambda = select_lambda,
      alpha_grid = alpha_grid, nlambda = nlambda, lambda_ratio = lambda_ratio,
      nfolds = nfolds, method = method, probs = probs, cutoffplot = cutoffplot,
      seed = seed, value = value, niter = niter, conv = conv,
      parallel = parallel_cv, plotCV = plotCV, colors_pcv = colors_pcv,
      errorbar = errorbar, ncore_max = ncore_max, standardize = standardize_train,
      verbose = verbose, palette = palette)
    
    beta <- fit_training$beta
    
    opt_cutoff <- fit_training$cutoff.opt
    
  # --- Testing phase ---
    
    fit_testing <- NetSurvProx_Testing(
      X_train = X_train, standardize = standardize_test,
      Y_train = Y_train, delta_train = delta_train,
      X_test = X_test, Y_test = Y_test, delta_test = delta_test,
      model = model, dist = dist, beta = beta, beta_true = beta_true,
      opt_cutoff = opt_cutoff, p_active = p_active,
      times_auc = times_auc, metrics = metrics, verbose = verbose,
      plot = plot_test, palette = palette)
    
    obj <- c(list(
      fit_training = fit_training,
      fit_testing = fit_testing))
    
    class(obj) = "NetSurvProx"
    
    obj
  }

# ============================================================================

#' NetSurvProx Training Routine
#'
#' @description
#' Trains penalized regression methods (\code{COXNet} or \code{AFTNet}) to incorporate
#' gene regulatory relationships and select signature genes using the training set.
#' Regularization parameters are selected via cross-validation, and an optimal
#' Prognostic Index (PI) cutoff is determined for risk stratification (\code{COXNet}) or
#' for survival time stratification (\code{AFTNet}). The procedure includes optional
#' feature standardization and simultaneous selection of the regularization
#' parameters for the Laplacian constraint and the Lasso penalty.
#'
#' @param X_train Numeric matrix of training covariates standardized
#'                (possibly screened using \code{screen_vars}).
#' @param Y_train Numeric vector of observed training survival times (log-transformed under \code{AFTNet}).
#' @param delta_train Integer vector of training censoring indicators (1 = event, 0 = censored).
#' @param L Optional positive semi-definite, symmetric, and diagonally dominant
#'          Laplacian matrix encoding prior network information. If \code{NULL},
#'          no network-based penalization is applied.
#' @param model Character string specifying the fitted survival model
#'              (\code{"COXNet"}, or \code{"AFTNet"}).
#' @param dist Character string specifying the \code{AFTNet} distribution.
#'             Must be one of \code{"weibull"}, \code{"lognormal"}, or \code{"loglogistic"}.
#' @param select_lambda Logical value, if \code{TRUE} (default) uses \code{lambda.min},
#'                      otherwise \code{lambda.1se}.
#' @param alpha_grid Numeric vector specifying the candidate values for \eqn{\alpha}
#'                   in \code{[0,1]} (\code{default: c(0.3, 0.5, 0.7)}).
#' @param nlambda Numeric value specifying the number of candidate values for
#'                \eqn{\lambda} in the grid (\code{default: 50}).
#' @param lambda_ratio Numeric value giving the ratio of minimum to maximum
#'                     \eqn{\lambda} in the grid (\code{default: 0.01}).
#' @param nfolds Number of cross-validation folds (\code{default: 5}).
#' @param method Character string specifying the cutoff selection method
#'               (\code{"median"} or \code{"minpvalue"}).
#' @param probs Vector of probabilities used when \code{method = "minpvalue"} to
#'              generate candidate cutoffs based on quantiles of the PI
#'              (\code{default: probs = seq(0.25, 0.80, by = 0.05)}).
#' @param cutoffplot Logical value indicating whether survival curves should be produced
#'                   (\code{default: FALSE}).
#' @param seed Random seed for reproducibility (\code{default: 2026}).
#' @param value Numeric scalar greater than 1 specifying the multiplicative
#'              factor used to increase the step-size constant during
#'              backtracking line search (\code{default: 2}).
#' @param niter Maximum number of iterations for ProxGDNet (\code{default: 1000}).
#' @param conv Convergence tolerance for ProxGDNet (\code{default: 1e-3}).
#' @param parallel Logical value whether to use parallel processing for CvNet (\code{default: TRUE}).
#' @param plotCV Logical value indicating whether cross-validation curves should be shown
#'               (\code{default: FALSE}).
#' @param colors_pcv Optional named list of colors:
#'                \itemize{
#'                    \item \code{line}: colorof the cross-validation error curve.
#'                    \item \code{points}: color of observed CV error evaluations.
#'                    \item \code{min}: color of the vertical line indicating \code{lambda.min}.
#'                    \item \code{one_se}: color of the vertical line indicating \code{lambda.1se}.
#'                }
#'                If \code{NULL}, a default color palette is used.
#' @param errorbar Logical value, if \code{TRUE} the CV plot includes vertical error
#'                 bars representing 1SE of the CV error (\code{default: FALSE}).
#' @param ncore_max Maximum number of cores for parallel processing over CV (\code{default: 5}).
#' @param standardize Logical value indicating whether to standardize the input matrix:
#'                    if \code{TRUE} (default), each column is centered to have mean 0 
#'                    and scaled to have unit variance,
#'                    if \code{FALSE}, the matrix is assumed pre-standardized by the user.
#' @param verbose Logical value, if \code{TRUE} progress messages are printed (\code{default: FALSE}).
#' @param palette Optional character vector of length 2 specifying colors used
#'                for the survival curves. For \code{"COXNet"}, colors correspond to
#'                high- and low-risk groups. For \code{"AFTNet"}, colors correspond to
#'                short- and long-survival groups. If \code{NULL}, default colors are used.
#' 
#' @details
#' The function performs joint tuning for regularization parameters:
#' a grid of \eqn{\alpha} values in (0, 1) is constructed, and for each candidate
#' computes corresponding \eqn{\lambda} grids via cross-validation using the negative
#' (partial for \code{COXNet}) log-likelihood's gradient.
#' 
#' Parallel computation is supported to improve efficiency.
#'
#' @return A list containing:
#' \itemize{
#'  \item \code{alpha.opt}: numeric value of optimal alpha.
#'  \item \code{lambda.opt}: numeric value of optimal lambda.
#'  \item \code{beta}: estimated regression coefficients.
#'  \item \code{index.nonzerobeta}: index of non-zero beta.
#'  \item \code{lambda.min}: value of \eqn{\lambda} minimizing the CV error.
#'  \item \code{lambda.1se}: largest \eqn{\lambda} within one standard error of the minimum.
#'  \item \code{cutoff.opt}: numeric value of optimal prognostic index cutoff.
#'  \item \code{lambda.grid}: grid of regularization parameters values.
#'  \item \code{cv.err.linPred}: cross-validated error for each value of \eqn{\lambda}.
#'  \item \code{cv.err.obj}: estimated standard error associated with each value of CV error.
#'  \item \code{full_summary}: data.frame as summary of CV results for all tested \code{alpha} values.
#' }
#'
#' @importFrom foreach %dopar%
#'
#' @seealso
#' \itemize{
#'    \item \code{\link{CreateNetwork}}: for \code{L} matrix computation.
#'    \item \code{\link{CvNet}}: for CV and parallel processing details.
#'    \item \code{\link{PlotCvNet}}: for cross-validation plot.
#'    \item \code{\link{OptimalPICutoff}}: for the optimal cutoff value to stratify observations.
#'    \item \code{\link{ProxGDNet}}: for proximal network-penalized gradient descent algorithm details.
#'    \item \code{\link{VariableScreening}}: for the \code{screen_vars} list.
#' }
#' 
#' @name NetSurvProx_Training
#'
#' @export

  NetSurvProx_Training <- function(
      X_train, Y_train, delta_train,
      L             = NULL,
      model         = NULL,
      dist          = NULL,
      select_lambda = TRUE,
      alpha_grid    = c(0.3, 0.5, 0.7),
      nlambda       = 50,
      lambda_ratio  = 0.01,
      nfolds        = 5,
      method        = NULL,
      probs         = seq(0.25, 0.80, by = 0.05),
      cutoffplot    = FALSE,
      seed          = 2026,
      value         = 2,
      niter         = 1000,
      conv          = 1e-3,
      parallel      = TRUE,
      plotCV        = FALSE,
      colors_pcv    = NULL,
      errorbar      = FALSE,
      ncore_max     = 5,
      standardize   = TRUE,
      verbose       = FALSE,
      palette       = NULL
  ) {
    
    # --- Input checks ---
      
      model <- match.arg(model, choices = c("COXNet", "AFTNet"))
      method <- match.arg(method, choices = c("median", "minpvalue"))
      
      if (model == "AFTNet" && is.null(dist)) stop("Argument 'dist' must be specified for 'AFTNet'.")
      if (!is.null(dist)) dist <- match.arg(dist, choices = c("weibull", "lognormal", "loglogistic"))
      
    # - Sigma estimation for AFTNet -
      
      if (model == "AFTNet") {
        
        fit_survreg <- survival::survreg(survival::Surv(exp(Y_train), delta_train) ~ 1, dist = dist)
        sigma <- fit_survreg$scale
        
      } else sigma = NULL
      
      if (isTRUE(verbose)) {
        message("-----------------------------------------------------")
        message(format("TRAINING PHASE", width = 53, justify = "centre"))
        message("-----------------------------------------------------")
        message(
          sprintf("Model: %s\n", model),
          if (model == "AFTNet") sprintf("Dist: %s\n", dist),
          if (model == "AFTNet") sprintf("Sigma: %.3f\n", sigma)
        )}
      
    # - Tuning parameters -
      
      results <- lapply(alpha_grid, function(a) {
        opt_data <- get_opt_params(
          X = X_train, Y = Y_train, delta = delta_train, L = L, model = model,
          dist = dist, sigma = sigma, alpha = a, nlambda = nlambda,
          lambda_ratio = lambda_ratio, select_lambda = select_lambda,
          nfolds = nfolds, seed = seed, value = value, niter = niter, conv = conv,
          parallel = parallel, ncore_max = ncore_max,
          standardize = standardize, plotCV = plotCV, colors_pcv = colors_pcv,
          errorbar = errorbar, verbose = verbose)
        
        if (isTRUE(verbose)) message(sprintf("\n Alpha: %.1f\n", a))
        
        list(
          summary = data.frame(
            alpha           = a,
            lambda.min      = opt_data$lambda.min,
            lambda.1se      = opt_data$lambda.1se,
            lambda.selected = opt_data$lambda.selected,
            index.min       = opt_data$index.min,
            index.1se       = opt_data$index.1se,
            min_cv_err      = min(opt_data$cv.err.obj, na.rm = TRUE)),
          opt_data = opt_data)})
      
    # - Best parameters -
      
      all_summaries  <- do.call(rbind, lapply(results, function(x) x$summary))
      best_alpha_idx <- which.min(all_summaries$min_cv_err)
      optimal_row    <- all_summaries[best_alpha_idx, ]
      
      if (isTRUE(verbose)) {
        message("-----------------------------------------------------")
        message(format("OPTIMAL PARAMETERS FOUND", width = 53, justify = "centre"))
        message("-----------------------------------------------------")
        message(
          sprintf("\nAlpha: %.1f | Lambda: %.6f", optimal_row$alpha, optimal_row$lambda.selected),
          sprintf("\nCV-LP Error Min: %.4f", optimal_row$min_cv_err),
          "\n"
        )}
      
      best_alpha_data <- results[[best_alpha_idx]]$opt_data
      
    # - Fit optimal model -
      
      if(standardize == TRUE){
        X_std <- standardize(X_train)$std
      } else X_std <- X_train
      
      beta0 <- rep(0, ncol(X_train))
      
      fit_Train <- ProxGDNet(
        X = X_std, Y = Y_train, delta = delta_train, L = L, beta0 = beta0,
        alpha = optimal_row$alpha, lambda = optimal_row$lambda.selected,
        model = model, dist = dist, sigma = sigma, value = value, niter = niter, conv = conv)
      
      if(standardize == TRUE){
        beta <- fit_Train$beta / standardize(X_train)$sds
      } else beta <- fit_Train$beta
      
      selected_indices <- sum(beta != 0)
      
      if (isTRUE(verbose)) message("Selected genes: ", selected_indices)
      
    # - Compute the optimal cutoff -
      
      ind_nonzerobeta <- which(beta != 0)
      
      if (length(ind_nonzerobeta) > 0) {
        
        cutoff <- OptimalPICutoff(
          X = X_train[, ind_nonzerobeta, drop = FALSE],
          Y = Y_train, delta = delta_train,
          beta = beta[ind_nonzerobeta, , drop = FALSE],
          method = method, model = model, dist = dist,
          probs = probs, plot = cutoffplot, palette = palette)
        
        cutoff_val <- cutoff$cutoff
        
      } else {
        
        warning("No variables selected.")
        cutoff_val <- NA
      }
      
      return(list(
        alpha.opt         = optimal_row$alpha,
        lambda.opt        = optimal_row$lambda.selected,
        beta              = beta,
        index.nonzerobeta = ind_nonzerobeta,
        lambda.min        = optimal_row$lambda.min,
        lambda.1se        = optimal_row$lambda.1se,
        cutoff.opt        = cutoff_val,
        lambda.grid       = best_alpha_data$lambda.grid,
        cv.err.linPred    = best_alpha_data$cv.err.linPred,
        cv.err.obj        = best_alpha_data$cv.err.obj,
        full_summary      = all_summaries))
    }

# ============================================================================

#' Joint Selection of Regularization Parameters via Cross-Validation
#' 
#' @keywords internal
#' @noRd

  get_opt_params <- function(
      X, Y, delta, L,
      model, dist, sigma,
      alpha, nlambda,
      lambda_ratio,
      select_lambda, nfolds,
      seed, value, niter, conv,
      parallel, ncore_max,
      standardize,
      plotCV, colors_pcv,
      errorbar,
      verbose
  ) {
    
    p <- ncol(X)
    
    if(standardize == TRUE){
      X_std <- standardize(X)$std
    } else X_std  <- X
    
    beta0 <- rep(0, p)
    eta0  <- as.matrix(X_std) %*% beta0
    
    zmax <- switch(model,
                   
                   "COXNet" = max(abs(gradient_COX(X = X_std, Y = Y, eta = eta0,
                                                   delta = delta)$grad_beta)),
                   
                   "AFTNet" = max(abs(gradient_AFT(X = X_std, Y = Y, eta = eta0, delta = delta,
                                                   sigma = sigma, dist = dist))))
    
    if (alpha > 0) {
      lambda_max <- (zmax / alpha) + 1e-4
    } else {
      lambda_max <- zmax + 1e-4
    }
    
    lambda_min  <- lambda_ratio * lambda_max
    lambda_grid <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
    
    cv.out <- CvNet(
      X = X_std, Y = Y, delta = delta, L = L,
      lambda = lambda_grid, alpha = alpha,
      model = model, dist = dist, sigma = sigma,
      nfolds = nfolds, seed = seed + round(100 * alpha),
      value = value, niter = niter, conv = conv, parallel = parallel,
      ncore_max = ncore_max, verbose = verbose)
    
  # - Plot CV -
    
    if(plotCV){
      
      p_cv <- PlotCvNet(cv.out = cv.out, alpha = alpha,
                        errorbar = errorbar, colors = colors_pcv)
      print(p_cv)
      
    }
    
  # - Avoid to take lambda equal to lambda_min or lambda_max -
    
    if (select_lambda) {
      
      if (cv.out$lambda.min > lambda_min ){
        
        selected_lambda <- cv.out$lambda.min
        
      } else selected_lambda <- cv.out$lambda.1se
      
    } else {
      
      if (cv.out$lambda.1se < lambda_max){
        
        selected_lambda <- cv.out$lambda.1se
        
      } else selected_lambda <- cv.out$lambda.min
    }
    
    return(list(
      alpha             = alpha,
      lambda.grid       = lambda_grid,
      lambda.min        = cv.out$lambda.min,
      lambda.1se        = cv.out$lambda.1se,
      lambda.selected   = selected_lambda,
      index.min         = cv.out$ind.lambda.min,
      index.1se         = cv.out$ind.lambda.1se,
      cv.err.linPred    = cv.out$cv.err.linPred,
      cv.err.obj        = cv.out$cv.err.obj,
      cvup              = cv.out$cvup,
      cvlo              = cv.out$cvlo))
  }

# ============================================================================

#' NetSurvProx Testing Routine
#'
#' @description
#' Evaluates predictive performance of a fitted \code{COXNet} or \code{AFTNet} model
#' on an independent testing set. The function computes the Prognostic Index (PI)
#' using the selected signature genes and the optimal cutoff obtained from the
#' training phase, generates survival curves, PI distribution plots, and calculates
#' specified performance metrics.
#' 
#' @param X_train Numeric matrix of training covariates (used only to scale
#'                \code{X_test} when \code{standardize = TRUE}).
#' @param standardize Logical value indicating whether to standardize \code{X_test}
#'                    with respect to \code{X_train} (\code{default: TRUE}).
#' @param Y_train Numeric vector of observed training survival times (log-transformed under \code{AFTNet}).
#'                Required only for time-dependent AUC computation.
#' @param delta_train Integer vector of training censoring indicators (1 = event, 0 = censored).
#'                    Required only for time-dependent AUC computation.
#' @param X_test Numeric matrix of testing covariates.
#' @param Y_test Numeric vector of observed testing survival times (log-transformed under \code{AFTNet}).
#' @param delta_test Integer vector of testing censoring indicators (1 = event, 0 = censored).
#' @param model Character string specifying the fitted survival model
#'              (\code{"COXNet"}, or \code{"AFTNet"}).
#' @param dist Character string specifying the \code{AFTNet} distribution.
#'             Must be one of \code{"weibull"}, \code{"lognormal"}, or \code{"loglogistic"}.
#' @param beta Numeric vector of regression coefficients estimated on the training set.
#' @param beta_true Numeric vector of true coefficients (used only for simulated data).
#' @param opt_cutoff Numeric value used to split the PI into two prognostic groups.
#' @param p_active Numeric value indicating the number of truly active covariates
#'                 (required for FPR/FNR computation in simulation settings).
#' @param times_auc Numeric vector of time points for time-dependent AUC.
#'                  If \code{NULL} (default), quantiles of \code{Y_test} are used.
#' @param metrics Character vector specifying performance metrics to compute.
#'                For real datasets: \code{"CIndex"}, \code{"NSR"}, \code{"AUC"}.
#'                For simulated datasets (in addition): \code{"FPR"}, \code{"FNR"}, \code{"PMSE"}.
#' @param verbose Logical value, if \code{TRUE} progress messages are printed (\code{default: FALSE}).
#' @param plot Logical value, if \code{TRUE} returns the combined survival plot (\code{default: FALSE}).
#' @param palette Optional character vector of length 2 specifying colors used
#'                for the survival curves. For \code{"COXNet"}, colors correspond to
#'                high- and low-risk groups. For \code{"AFTNet"}, colors correspond to
#'                short- and long-survival groups. If \code{NULL}, default colors are used.
#' 
#' @details
#' The testing set must be independent from the training set used in \code{NetSurvProx_Training}.
#' When \code{standardize = TRUE}, \code{X_test} is standardized using the mean and standard deviation
#' of \code{X_train}. Only covariates with non-zero coefficients in \code{beta} are retained for PI computation.
#'
#' Prognostic stratification is performed using \code{\link{ValidationPI}}, producing:
#' \itemize{
#'   \item Kaplan–Meier curves and log-rank test for \code{COXNet}.
#'   \item Parametric survival curves and likelihood ratio test for \code{AFTNet}.
#'   \item PI distribution plots by prognostic group.
#' }
#'
#' @return A list containing:
#' \itemize{
#'    \item \code{df}: data frame with \code{PI} (computed for each subject),
#'                     \code{Y}, \code{delta}, and \code{groupRisk}
#'                     (prognostic group assigned based on \code{opt_cutoff}).
#'    \item \code{p_value}: from the log-rank test (\code{COXNet}) or likelihood ratio test (\code{AFTNet}).
#'    \item \code{performance}: named list with the requested performance metrics.
#' }
#'
#' @seealso
#' \itemize{
#'    \item \code{\link{Metrics}} for available performance \code{metrics} options.
#'    \item \code{\link{NetSurvProx_Training}} for training routine.
#'    \item \code{\link{OptimalPICutoff}} for \code{opt_cutoff} estimation.
#'    \item \code{\link{ValidationPI}} for PI validation and optional plot.
#' }
#' 
#' @name NetSurvProx_Testing
#'
#' @export

  NetSurvProx_Testing <- function(
      X_train     = NULL,
      standardize = TRUE,
      Y_train     = NULL,
      delta_train = NULL,
      X_test, Y_test, delta_test,
      model       = NULL,
      dist        = NULL,
      beta,
      beta_true   = NULL,
      opt_cutoff,
      p_active    = NULL,
      times_auc   = NULL,
      metrics     = NULL,
      verbose     = FALSE,
      plot        = FALSE,
      palette     = NULL
  ) {
    
    model <- match.arg(model, choices = c("COXNet", "AFTNet"))
    
    if (isTRUE(verbose)) {
      message("-----------------------------------------------------")
      message(format("TESTING PHASE", width = 53, justify = "centre"))
      message("-----------------------------------------------------")
    }
    
    if (standardize) {
      
      mu <- colMeans(X_train)
      n  <- nrow(X_test)
      Xc <- X_test - tcrossprod(rep(1,n), mu)
      
    } else Xc <- X_test
    
    ind_nonzerobeta <- which(beta != 0)
    
    testing <- ValidationPI(X = X_test[ , ind_nonzerobeta], Y = Y_test,
                            delta = delta_test, beta = beta[ind_nonzerobeta],
                            opt_cutoff = opt_cutoff, model = model, dist = dist,
                            plot = plot, palette = palette)
    
    df   <- testing$df
    pval <- testing$p_value
    
    performance <- Metrics(Y_train = Y_train, delta_train = delta_train,
                           X_test = Xc, Y_test = Y_test, delta_test = delta_test,
                           beta_est = beta, beta_true = beta_true,
                           model = model, p_active = p_active,
                           times_auc = times_auc, metrics = metrics)
    
    if (isTRUE(verbose)) {
      message("Performance metrics:")
      print(unlist(performance))
    }
    
    return(list(
      df          = df,
      p_value     = pval,
      performance = performance))
    
  }
