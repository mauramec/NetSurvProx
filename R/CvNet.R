#' Cross-validated Linear Predictors Approach for \code{COXNet} and \code{AFTNet}
#'
#' @description
#' Performs K-fold cross-validation to select the optimal regularization
#' parameter \eqn{\lambda} for penalized survival models (\code{COXNet}, \code{AFTNet})
#' estimated via \code{\link{ProxGDNet}}. The criterion is based on cross-validated
#' linear predictors and negative (partial) log-likelihood.
#'
#' @param X Numeric matrix of standardized covariates.
#' @param Y Numeric vector of observed survival times (log-transformed under \code{AFTNet}).
#' @param delta Integer vector of censoring indicators (1 = event, 0 = censored).
#' @param L Optional positive semi-definite, symmetric, and diagonally dominant
#'          Laplacian matrix encoding prior network information
#'          (see \code{\link{CreateNetwork}} for details). If \code{NULL},
#'          no network-based penalization is applied.
#' @param lambda Numeric vector of candidate tuning parameters (in descending order).
#' @param alpha Numeric parameter controlling the convex combination of the two
#'              penalty terms (value in \code{[0,1]}).
#' @param model Character string specifying the fitted survival model
#'              (\code{"COXNet"}, or \code{"AFTNet"}).
#' @param dist Character string specifying the error distribution
#'             of \code{AFTNet} model. Must be one of \code{"weibull"}, 
#'             \code{"lognormal"}, or \code{"loglogistic"}.
#' @param sigma Positive numeric scalar representing the scale parameter of the
#'              error distribution in \code{AFTNet} model.
#' @param nfolds Number of cross-validation folds (\code{default: 5}).
#' @param seed Random seed for reproducibility (\code{default: 2026}).
#' @param value Numeric scalar greater than 1 specifying the multiplicative
#'              factor used to increase the step-size constant during
#'              backtracking line search (\code{default: 2}).
#' @param niter Maximum number of proximal gradient iterations (\code{default: 1000}).
#' @param conv Convergence tolerance for proximal gradient (\code{default: 1e-3}).
#' @param parallel Logical value, whether to use parallel processing (\code{default: TRUE}).
#' @param ncore_max Maximum number of cores for parallel processing over cross validation (\code{default: 5}).
#' @param verbose Logical value, if \code{TRUE} progress messages are printed (\code{default: FALSE}).
#'                
#' @details
#' The dataset is split into K folds. For each fold, the model is trained on K-1 folds,
#' and evaluated on the held-out fold. The cross-validated linear predictor is computed as
#' \deqn{ \hat{\eta}^{CV}_i = \boldsymbol{x}_i^\top \boldsymbol{\hat{\beta}}_\lambda^{(-k)}}
#' for \code{COXNet}, or the cross-validated standardized residual as
#' \deqn{ \hat{e}^{CV}_i = \frac{y_i - \boldsymbol{x}_i^\top \boldsymbol{\hat{\beta}}_\lambda^{(-k)}}{\hat{\sigma}}}
#' for \code{AFTNet}, and used to evaluate the cross-validation criterion over a grid of \eqn{\lambda} values.
#'
#' The optimal parameter is selected according to:
#' \itemize{
#'    \item the minimum CV error (\code{lambda.min}).
#'    \item the largest \eqn{\lambda} within one standard error
#'          of the minimum (\code{lambda.1se}).
#' }
#'
#' @note
#' Computation can be performed sequentially (\code{parallel: FALSE}), or
#' in parallel (\code{parallel: TRUE}) using \code{parLapply}.
#' The number of cores is automatically determined based on system availability,
#' number of folds and user-specified maximum \code{ncore_max}.
#'
#' @return An object of class \code{"cv.out"} containing:
#' \itemize{
#'   \item \code{cv.err.linPred}: CV error for each value of \eqn{\lambda}.
#'   \item \code{cv.err.obj}: estimated standard error associated with each value of CV error per fold.
#'   \item \code{lambda.grid}: grid of regularization parameters values.
#'   \item \code{lambda.min}: value of \eqn{\lambda} minimizing the CV error.
#'   \item \code{ind.lambda.min}: indices of \code{lambda.min}.
#'   \item \code{lambda.1se}: largest \eqn{\lambda} within one standard error of the minimum.
#'   \item \code{ind.lambda.1se}: indices of \code{lambda.1se}.
#'   \item \code{cvup}: upper error curve.
#'   \item \code{cvlo}: lower error curve.
#' }
#'
#' @importFrom foreach %dopar%
#'
#' @seealso
#' \itemize{
#'  \item \code{\link{PlotCvNet}} for visualization of the obtained cross-validation curve.
#'  \item \code{\link{ProxGDNet}} for proximal network-penalized gradient descent algorithm details.
#' }
#'
#' @name CvNet
#' 
#' @export

  CvNet <- function(
      X, Y, delta,
      L         = NULL,
      lambda, alpha,
      model     = NULL,
      dist      = NULL,
      sigma     = NULL,
      nfolds    = 5,
      seed      = 2026,
      value     = 2,
      niter     = 1000,
      conv      = 1e-3,
      parallel  = TRUE,
      ncore_max = 5,
      verbose   = FALSE
  ) {
    
  # --- Input checks --- 
    
    if (!is.matrix(X)) stop("'X' must be a numeric matrix.")
    
    if (length(Y) != nrow(X)) stop("'Y' must have length equal to nrow(X).")
    
    if (length(delta) != nrow(X)) stop("'delta' must have length equal to nrow(X).")
    if (!all(delta %in% c(0,1))) stop("delta must be 0/1.")
    
    if (alpha < 0 || alpha > 1) stop("alpha must be in [0,1].")
    
    model <- match.arg(model, choices = c("COXNet", "AFTNet"))
    
    if (model == "AFTNet" && (is.null(dist) || is.null(sigma)))
      stop("'dist' and 'sigma' must be provided for AFTNet model.")
    
    if (model == "AFTNet" && sigma <= 0) stop("'sigma' must be positive for AFTNet.")
    
    dist  <- if (model == "AFTNet") match.arg(dist, choices = c("weibull", "lognormal", "loglogistic")) else NULL
    
  # ----- 
    
    if (is.null(L)) {
      p <- ncol(X)
      L <- matrix(0, p, p)
    }
    
    set.seed(seed)
    
    n <- nrow(X)
    
    folds   <- cvTools::cvFolds(n, K = nfolds, type = "random")
    nlambda <- length(lambda)
    
    preds      <- matrix(NA, nrow = n, ncol = nlambda)
    cv.err.obj <- matrix(NA, nrow = nfolds, ncol = nlambda)
    
    ncores       <- min(parallel::detectCores() - 1, nfolds, ncore_max)
    use_parallel <- parallel && ncores > 1
    
    if (use_parallel) {
      
      if (isTRUE(verbose)) message("Running cross-validation in parallel...")
      
      cl <- parallel::makeCluster(ncores)
      
      parallel::clusterExport(
        cl,
        varlist = c("cvnet_fold","ProxGDNet","nll_COX","nll_AFT",
                    "hessian_COX", "hessian_AFT", "gradient_COX",
                    "gradient_AFT","prox_l1", "risk_fun", "L",
                    "X", "Y", "delta", "sigma", "lambda", "alpha",
                    "model", "dist", "niter", "value", "conv", "folds"), 
        envir = environment())
      
      out.fold <- parallel::parLapply(cl, 1:nfolds, function(k) {
        cvnet_fold(folds = folds, k = k, X = X, Y = Y, delta = delta, L = L,
                   sigma = sigma, lambda = lambda, alpha = alpha,
                   model = model, dist = dist, value = value,
                   niter = niter, conv = conv)})
      
      parallel::stopCluster(cl)
      
    } else {
      
      if (isTRUE(verbose)) message("Running cross-validation sequentially...")
      
      out.fold <- lapply(1:nfolds, function(k) {
        cvnet_fold(folds = folds, k = k, X = X, Y = Y, delta = delta, L = L,
                   sigma = sigma, lambda = lambda, alpha = alpha,
                   model = model, dist = dist, value = value, 
                   niter = niter, conv = conv)})
    }
    
  # - Combine results from parallel/sequential execution -
    
    for (k in seq_len(nfolds)) {
      
      index.cv   <- folds$subsets[folds$which != k]
      index.pred <- setdiff(1:n ,index.cv)
      
      preds[index.pred, ] <- out.fold[[k]]$preds
      cv.err.obj[k, ]     <- out.fold[[k]]$cv.err.obj}
    
    cv.err.linPred <- sapply(seq_len(nlambda), function(ll) {
      
      switch(model,
             
             "COXNet" = nll_COX(Y = Y, eta = preds[, ll], delta = delta),
             
             "AFTNet" = nll_AFT(Y = Y, eta = preds[, ll], delta = delta,
                                sigma = sigma, dist = dist))})
    
  # - lambda.min -
    
    ind_min    <- which.min(cv.err.linPred)
    lambda.min <- lambda[ind_min]
    
  # - lambda.1se -
    
    cv.err.linPred.se <- vapply(seq_len(nlambda), function(ll) {
      stats::sd(cv.err.obj[, ll]) / sqrt(nfolds)}, numeric(1))
    
    lambda.1se <- max(lambda[cv.err.linPred < (cv.err.linPred[ind_min] + cv.err.linPred.se[ind_min])])
    ind_1se    <- which(lambda == lambda.1se)
    
    cvup <- cv.err.linPred + cv.err.linPred.se
    cvlo <- cv.err.linPred - cv.err.linPred.se
    
    out <- list(
      cv.err.linPred = cv.err.linPred,
      cv.err.obj     = cv.err.obj,
      lambda.grid    = lambda,
      lambda.min     = lambda.min,
      ind.lambda.min = ind_min,
      lambda.1se     = lambda.1se,
      ind.lambda.1se = ind_1se,
      cvup           = cvup,
      cvlo           = cvlo)
    
    class(out) <- "cv.out"
    if (isTRUE(verbose)) message("COMPLETED")
    return(out)
  }

# ============================================================================

#' Cross-validation Computations for a Single Fold
#' 
#' @keywords internal
#' @noRd

  cvnet_fold <- function(
      folds, k, X, Y, delta,
      L, sigma, lambda,
      alpha, model,
      dist, value,
      niter, conv
  ) {
    
    index.cv <- folds$subsets[folds$which != k]
    
    if (any(index.cv <= 0)) stop("Invalid fold indices.")
    if (any(index.cv > nrow(X))) stop("Index out of range.")
    
    index.test <- setdiff(seq_len(nrow(X)), index.cv)
    
    ntrain  <- length(index.cv)
    ntest   <- length(index.test)
    nlambda <- length(lambda)
    
    preds <- matrix(NA, nrow = ntest, ncol = nlambda)
    
    standardize_with_mean_std <- function(data, mean_vals, sd_vals) {
      (data - matrix(mean_vals, nrow = nrow(data), ncol = length(mean_vals), byrow = TRUE)) /
        matrix(sd_vals, nrow = nrow(data), ncol = length(sd_vals), byrow = TRUE)}
    
  # - Standardize training and testing sets for X -
    
    mX_train  <- colMeans(X[index.cv, ])
    sdX_train <- apply(X[index.cv, ], 2, stats::sd)
    sdX_train[sdX_train == 0] <- 1
    
    X_train   <-  standardize_with_mean_std(X[index.cv, ], mX_train, sdX_train)
    X_test    <-  standardize_with_mean_std(X[index.test, ], mX_train, sdX_train)
    
    Y_train     <- Y[index.cv]
    delta_train <- delta[index.cv]
    Y_test      <- Y[index.test]
    delta_test  <- delta[index.test]
    
    beta0      <- rep(0, ncol(X_train))
    cv.err.obj <- vector()
    
  # - CV predictions loop -
    
    for (ll in seq_len(nlambda)) {
      cv.getPath <- switch(model,
                           
                           "COXNet" = ProxGDNet(X = X_train, Y = Y_train, delta = delta_train, L = L,
                                                beta0 = beta0, alpha = alpha, lambda = lambda[ll],
                                                model = model, dist = dist, sigma = sigma,
                                                value = value, niter = niter, conv = conv),
                           
                           "AFTNet" = ProxGDNet(X = X_train, Y = Y_train, delta = delta_train, L = L,
                                                beta0 = beta0, alpha = alpha, lambda = lambda[ll],
                                                model = model, dist = dist, sigma = sigma,
                                                value = value, niter = niter, conv = conv))
      beta0       <- cv.getPath$beta
      preds[, ll] <- X_test %*% beta0
      
      cv.err.obj[ll] <- if(model == "COXNet") {
        
        nll_COX(Y = Y_test, eta = preds[, ll], delta = delta_test)
        
      } else nll_AFT(Y = Y_test, eta = preds[, ll], delta = delta_test, sigma = sigma, dist = dist)
    }
    
    return(list(
      preds      = preds,
      cv.err.obj = cv.err.obj))
  }

# ============================================================================
# ============================================================================

#' Plot CV-LP Curve for \code{COXNet} and \code{AFTNet}
#' 
#' @description
#' Produces a \pkg{ggplot2} visualization of the cross-validation curve obtained
#' from \code{\link{CvNet}}. The plot displays the CV error as a function of
#' \eqn{\log(\lambda)} with optional error bars, and reference lines for
#' \code{lambda.min} and \code{lambda.1se}.
#'
#' @param cv.out Object of class \code{"cv.out"} (returned by \code{\link{CvNet}}),
#' containing at least:
#'   \itemize{
#'     \item \code{cv.err.linPred}: mean CV errors for linear predictor.
#'     \item \code{lambda.grid}: grid of \eqn{\lambda} values used as regularization path.
#'     \item \code{lambda.min}: value of \eqn{\lambda} minimizing the CV error.
#'     \item \code{lambda.1se}: largest \eqn{\lambda} within one standard error.
#'     \item \code{cvup}: upper error curve.
#'     \item \code{cvlo}: lower error curve.    
#'    }
#' @param alpha Numeric parameter controlling the convex combination of the two
#'              penalty terms (value in \code{[0,1]}), used only for plot annotation
#'              (\code{default: NULL}).
#' @param errorbar Logical value, if \code{TRUE} the plot includes vertical error
#'                 bars representing 1se of the cross-validation error at each fold
#'                 (\code{default: FALSE}).
#' @param colors Optional named list of colors:
#'                \itemize{
#'                    \item \code{line}: color of the cross-validation error curve.
#'                    \item \code{points}: color of observed CV error evaluations.
#'                    \item \code{min}: color of the vertical line indicating \code{lambda.min}.
#'                    \item \code{one_se}: color of the vertical line indicating \code{lambda.1se}.
#'                }
#'                If \code{NULL}, a default color palette is used.                
#' 
#' @return A \pkg{ggplot2} object showing the CV-LP curve.
#' 
#' @name PlotCvNet
#' 
#' @export

  PlotCvNet <- function(
      cv.out,
      alpha    = NULL,
      errorbar = FALSE,
      colors   = NULL
  ) {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.", call. = FALSE)
    
    default_colors <- list(
      line   = "#5C7997",
      points = "#F5C59F",
      min    = "#8ABFE7",
      one_se = "#7DAE9B")
    
    if (is.null(colors)) {
      colors <- default_colors
    } else {
     
      if (is.list(colors) && !is.null(names(colors)) && all(names(colors) != "")) {
        
        if (any(names(colors) %in% names(default_colors))) {
          colors <- utils::modifyList(default_colors, colors)
        } else {
          warning("None of the names in 'colors' match expected parameters, using defaults.")
          colors <- default_colors
        }
        
      } else {
        warning("'colors' must be a named list, usign defaults.")
        colors <- default_colors
      }
    }
    
    if (!is.null(alpha)) {
      plot_title <- bquote("CV-LP Curve (" * alpha == .(alpha) * ")")
    } else {
      plot_title <- "CV-LP Curve"
    }
    
    df_plot <- data.frame(
      lambda    = cv.out$lambda.grid,
      LPCVE     = cv.out$cv.err.linPred,
      LPCVEhi   = cv.out$cvup,
      LPCVElow  = cv.out$cvlo)
    
    p_cv <- ggplot2::ggplot(df_plot, ggplot2::aes(x = log(lambda), y = LPCVE)) +
      ggplot2::geom_line(col  = colors$line, alpha = 0.5) +
      ggplot2::geom_point(col = colors$points, size = 1.5) +
      ggplot2::geom_vline(xintercept = log(cv.out$lambda.min),
                          linetype = "dashed", col = colors$min, linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = log(cv.out$lambda.1se),
                          linetype = "dotted", col = colors$one_se, linewidth = 0.9) +
      ggplot2::annotate("text", x = log(cv.out$lambda.min), y = max(cv.out$cvup),
                        label = "min", col = colors$min, angle = 90, vjust = -0.5, size = 3.5) +
      ggplot2::annotate("text", x = log(cv.out$lambda.1se), y = max(cv.out$cvup),
                        label = "1se", col = colors$one_se, angle = 90, vjust = -0.5, size = 3.5) +
      ggplot2::labs(title = plot_title,
                    subtitle = sprintf("lambda.min = %.3f | lambda.1se = %.3f",
                                       cv.out$lambda.min, cv.out$lambda.1se),
                    x = expression(log(lambda)), y = "CV-LP Error") +
      ggplot2::theme_minimal() + 
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)))
    
    if(isTRUE(errorbar)){
      
      p_cv <- p_cv + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = LPCVElow, ymax = LPCVEhi),
        col = colors$line, width = 0.05)
      
    } else {
      
      p_cv
      
    }
    
    return(p_cv)
  }  
