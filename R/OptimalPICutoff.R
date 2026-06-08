#' Optimal Cutoff for Prognostic Index on Training Set
#'
#' @description
#' Identifies the optimal cutoff value of a Prognostic Index (PI)
#' to stratify subjects into prognostic groups. It supports \code{COXNet} and \code{AFTNet}
#' models with several distributions.
#'
#' @param X Numeric matrix of covariates.
#' @param Y Numeric vector of observed survival times (log-transformed under \code{AFTNet}).
#' @param delta Integer vector of censoring indicators (1 = event, 0 = censored).
#' @param beta Numeric vector of estimated regression coefficients obtained from the training set.
#' @param method Character string specifying the cutoff selection method
#'               (\code{"median"} or \code{"minpvalue"}).
#' @param model Character string specifying the fitted survival model
#'              (\code{"COXNet"}, or \code{"AFTNet"}).
#' @param dist Character string specifying the \code{AFTNet} distribution.
#'             Must be one of \code{"weibull"}, \code{"lognormal"}, or \code{"loglogistic"}.
#' @param probs Vector of probabilities used when \code{method = "minpvalue"}
#'              to generate candidate cutoffs based on quantiles of the PI 
#'              (\code{default: probs = seq(0.25, 0.80, by = 0.05)}).
#' @param plot Logical value indicating whether survival curves should be produced
#'            (\code{default: FALSE}).
#' @param palette Optional character vector of length 2 specifying colors used
#'                for the survival curves. For \code{"COXNet"}, colors correspond to
#'                high- and low-risk groups. For \code{"AFTNet"}, colors correspond to
#'                short- and long-survival groups. If \code{NULL}, default colors are used.
#' @details
#' The Prognostic Index (PI) is computed as a linear predictor.
#' Two alternative strategies are available to define the cutoff.
#' 
#' \itemize{
#' \item **Median-based cutoff**: Subjects are dichotomized as follows:
#'        \itemize{
#'            \item \code{COXNet}: PI \eqn{\geq} median is \emph{High Risk}, otherwise \emph{Low Risk}.
#'            \item \code{AFTNet}: - PI \eqn{\geq} median is \emph{Short Survival}, otherwise \emph{Long Survival}.
#'        }
#' 
#' \item **Minimum p-value approach**: A grid of candidate cutoffs is generated from the quantiles of the PI.
#'        For each candidate:
#'        \itemize{
#'            \item The cohort is dichotomized according to the model-specific direction.
#'            \item Two models are fitted (full model including the group indicator, and null model without the group indicator).
#'            \item A likelihood ratio (LR) test is performed between the two models.
#'        }
#' }    
#' 
#' Model fitting is performed using \code{survival::coxph()} for \code{COXNet}, or \code{survival::survreg()} for \code{AFTNet}.
#' 
#' The raw p-values are adjusted for multiple testing using the Benjamini–Hochberg procedure.
#' The optimal cutoff corresponds to the smallest adjusted p-value.
#' 
#' If \code{plot = TRUE}, survival curves are generated (Kaplan–Meier curves for \code{COXNet},
#' parametric survival curves based on the selected distribution for \code{AFTNet}).
#'
#' @return 
#' For \code{method = "median"}, a list with
#' \itemize{
#'    \item \code{cutoff}: numeric cutoff value.
#'    \item \code{PI.data}: data frame containing the PI, survival time, status,
#'                          and group labels.
#' }
#' 
#' For \code{method = "minpvalue"}, the list additionally contains:
#' \itemize{
#'    \item \code{summary}: table of p-values across candidate quantiles.
#'    \item \code{optimal}: optimal cutoff information (quantile, cutoff value, 
#'                          raw and adjusted p-values).
#' }
#'
#' @export

  OptimalPICutoff <- function(
      X, Y, delta, beta,
      method  = NULL,
      model   = NULL,
      dist    = NULL,
      probs   = seq(0.25, 0.80, by = 0.05),
      plot    = FALSE,
      palette = NULL
  ) {
    
  # --- Input checks ---
    
    method <- match.arg(method, choices = c("median", "minpvalue"))
    model  <- match.arg(model, choices = c("COXNet", "AFTNet"))
    
    if (model == "AFTNet") {
      if (is.null(dist)) stop("'dist' must be provided when model = 'AFTNet'.")
      dist <- match.arg(dist, c("weibull", "lognormal", "loglogistic"))
    }
    
    if (any(probs <= 0 | probs >= 1)) {
      stop("'probs' must contain values strictly between 0 and 1.")
    }
    
    if (is.null(palette)) {
      palette <- switch(
        model,
        "COXNet" = c("#CC899D", "#7DAE9B"),
        "AFTNet" = c("#CC899D", "#7A78B7")
      )}
    
  # ----- 
      
    if (plot) {
      
      required_pkgs <- c("survminer", "ggplot2")
      for (pkg in required_pkgs) {
        if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")}
    }
    
    if (model == "AFTNet" && is.null(dist)) stop("Argument 'dist' must be specified with 'AFTNet'.")
    if (!is.null(dist)) dist <- match.arg(dist, choices = c("weibull", "lognormal", "loglogistic"))
    
    if (length(Y) != nrow(as.data.frame(X)) || length(delta) != nrow(as.data.frame(X))) {
      stop("The length of Y and delta must match the number of rows of X.")}
    
    ok    <- stats::complete.cases(X) & !is.na(Y) & !is.na(delta)
    x     <- as.data.frame(X[ok, , drop = FALSE])
    y     <- Y[ok]
    delta <- delta[ok]
    
    if (length(beta) != ncol(x)) stop("Length of 'beta' must match number of columns in X.")
    
    PI <- as.numeric((switch(model,
                             "COXNet" = as.matrix(x) %*% beta,
                             "AFTNet" = -as.matrix(x) %*% beta)))
    
    time_model <- if (model == "AFTNet") exp(y) else y
    
    plot_AFT <- function(PI.data, dist, time_model, delta) {
      group <- PI.data$group
      fit   <- survival::survreg(survival::Surv(time_model, delta) ~ group, data = PI.data, dist = dist)
      
      distinct_groups <- levels(factor(group))
      lp    <- stats::predict(fit, newdata = data.frame(group = distinct_groups), type = "lp")
      scale <- fit$scale
      
      max_t <- stats::quantile(time_model[delta == 1], 0.99, na.rm = TRUE)
      if(is.na(max_t)) max_t <- max(time_model)
      tt    <- seq(0, max_t, length.out = 100)
      
      Sfun <- switch(dist,
                     "weibull"     = function(t, l, s) exp(-(t / exp(l))^(1/s)),
                     "lognormal"   = function(t, l, s) stats::plnorm(t, meanlog = l, sdlog = s, lower.tail = FALSE),
                     "loglogistic" = function(t, l, s) 1 / (1 + (t / exp(l))^(1/s)))
      
      plot_df <- do.call(rbind, lapply(seq_along(distinct_groups), function(i) {
        data.frame(time = tt, surv = Sfun(tt, lp[i], scale), group = distinct_groups[i])}))
      
      ggplot2::ggplot(plot_df, ggplot2::aes(x = time, y = surv, color = group)) +
        ggplot2::geom_line(linewidth = 1.2) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(y = "Survival Probability", x = "Time", color = "Group:") +
        ggplot2::scale_color_manual(values = c("Long Survival"= palette[2], "Short Survival"= palette[1])) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "top", 
                       axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10))) 
    }
    
  # - Median cutoff -
    
    if (method == "median") {
      
      cutoff <- stats::median(PI, na.rm = TRUE)
      group  <- if (model == "COXNet") {
        ifelse(PI >= cutoff, "High Risk", "Low Risk")
      } else {
        ifelse(PI >= cutoff, "Short Survival", "Long Survival")
      }
      
      PI.data <- data.frame(PI = PI, time = time_model, status = delta, group = group)
      
      if (plot) {
        
        if (model == "COXNet") {
          
          fit  <- survival::survfit(survival::Surv(time, status) ~ group, data = PI.data)
          labs <- gsub("group=", "", names(fit$strata))
          p    <- survminer::ggsurvplot(fit, data = PI.data, pval = FALSE, 
                                      palette = palette, legend = "top", 
                                      legend.title = "Group:", legend.labs = labs,
                                      ggtheme = ggplot2::theme_minimal())
          p$plot <- p$plot + ggplot2::theme(
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)))
          print(p)
          
        } else {
          print(plot_AFT(PI.data, dist, time_model, delta))
        }
      }
      return(list(
        cutoff = cutoff,
        PI.data = PI.data))
    }
    
  # - Minpvalue cutoff -
    
    if (method == "minpvalue") {
      
      res <- quantile_cutoff(PI = PI, Y = time_model, delta = delta,
                             dist = dist, probs = probs, model = model)
      cutoff        <- res$optimal$cutoff
      best_quantile <- res$optimal$quantile
      PI.data       <- res$data[res$data$quantile == best_quantile, ]
      
      if (plot) {
        
        if (model == "COXNet") {
          fit  <- survival::survfit(survival::Surv(time, status) ~ group, data = PI.data)
          labs <- gsub("group=", "", names(fit$strata))
          p    <- survminer::ggsurvplot(fit, data = PI.data, pval = FALSE, 
                                      palette = palette, legend = "top", 
                                      legend.title = "Group:", legend.labs = labs,
                                      ggtheme = ggplot2::theme_minimal(),
                                      grid.alpha = 0.2)
          p$plot <- p$plot + ggplot2::theme(
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)))
          print(p)
          
        } else {
          print(plot_AFT(PI.data, dist, time_model, delta))}
      }
      
      return(list(
        cutoff = cutoff,
        PI.data = PI.data,
        summary = res$summary,
        optimal = res$optimal))}
  }

# ============================================================================
# Supporting function
# ============================================================================

#' Optimal Cutoff Performing a Likelihood Ratio Test
#' 
#' @keywords internal
#' @noRd

  quantile_cutoff <- function(
      PI, Y, delta,
      dist, probs, model
  ) {
    
    cutoffs   <- stats::quantile(PI, probs, na.rm = TRUE)
    pvals_raw <- numeric(length(cutoffs))
    all_data  <- vector("list", length(cutoffs))
    
    for (i in seq_along(cutoffs)) {
      
      group <- if (model == "COXNet") {
        ifelse(PI >= cutoffs[i], "High Risk", "Low Risk")
      } else {
        ifelse(PI <= cutoffs[i], "Long Survival", "Short Survival")}
      
      if (model == "COXNet") {
        fit_full <- survival::coxph(survival::Surv(Y, delta) ~ group)
        fit_null <- survival::coxph(survival::Surv(Y, delta) ~ 1)
      } else {
        fit_full <- survival::survreg(survival::Surv(Y, delta) ~ group, dist = dist)
        fit_null <- survival::survreg(survival::Surv(Y, delta) ~ 1, dist = dist)
      }
      
      lr_stat <- 2 * (stats::logLik(fit_full) - stats::logLik(fit_null))
      pvals_raw[i] <- 1 - stats::pchisq(lr_stat, df = 1)
      
      all_data[[i]] <- data.frame(
        sample = seq_along(PI),
        PI = PI,
        time = Y,
        status = delta,
        group = group,
        quantile = rep(probs[i], length(PI)),
        cutoff = rep(cutoffs[i], length(PI)),
        p.value.raw = rep(pvals_raw[i], length(PI)))
    }
    
    all_data <- do.call(rbind, all_data)
    pvals_adj <- stats::p.adjust(pvals_raw, method = "BH")
    all_data$p.value.adj <- rep(pvals_adj, each = length(PI))
    
    summary_list <- do.call(rbind, lapply(seq_along(cutoffs), function(i) {
      dat_sub  <- all_data[all_data$quantile == probs[i], ]
      groups   <- unique(dat_sub$group)
      n_groups <- sapply(groups, function(g) sum(dat_sub$group == g))
      
      data.frame(
        quantile    = probs[i],
        cutoff      = cutoffs[i],
        group1      = groups[1],
        n_group1    = n_groups[1],
        group2      = groups[2],
        n_group2    = n_groups[2],
        p.value.raw = pvals_raw[i],
        p.value.adj = pvals_adj[i])}))
    
    best <- which.min(pvals_adj)
    
    optimal <- list(
      best_index  = best,
      quantile    = probs[best],
      cutoff      = cutoffs[best],
      p.value.raw = pvals_raw[best],
      p.value.adj = pvals_adj[best])
    
    return(list(
      data    = all_data,
      summary = summary_list,
      optimal = optimal))
  }
    