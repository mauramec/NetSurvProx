#' Prognostic Index Validation on Testing Set
#'
#' @description
#' Validates a Prognostic Index (PI) obtained from a fitted survival model
#' (\code{COXNet} or \code{AFTNet}) on an independent testing set.
#' Given the estimated regression coefficients, it computes the PI for each subject,
#' assigns prognostic groups using a pre-specified optimal cutoff, and evaluates
#' survival separation and statistical significance.
#'
#' @param X Numeric matrix of testing covariates scaled using the training data.
#' @param Y Numeric vector of observed testing survival times (log-transformed under \code{AFTNet}).
#' @param delta Integer vector of testing censoring indicators (1 = event, 0 = censored).
#' @param beta Numeric vector of estimated regression coefficients obtained from the training set.
#' @param opt_cutoff Numeric cutoff value used to split the PI into two prognostic groups.
#' @param model Character string specifying the fitted survival model
#'              (\code{"COXNet"}, or \code{"AFTNet"}).
#' @param dist Character string specifying the \code{AFTNet} distribution.
#'             Must be one of \code{"weibull"}, \code{"lognormal"}, or \code{"loglogistic"}.
#' @param plot Logical value, if \code{TRUE} returns the combined survival plot (\code{default: FALSE}).
#' @param palette Optional character vector of length 2 specifying colors used
#'                for the survival curves. For \code{"COXNet"}, colors correspond to
#'                high- and low-risk groups. For \code{"AFTNet"}, colors correspond to
#'                long- and short-survival groups. If \code{NULL}, default colors are used.
#'
#' @details
#' For \code{COXNet}, Kaplan-Meier survival curves are computed, a log-rank test
#' is performed, and the \eqn{PI = X \beta} is compared to \code{opt_cutoff}
#' to define \emph{High Risk} and \emph{Low Risk} groups.
#'
#' For \code{AFTNet}, parametric survival curves are computed using the specified
#' distribution, a likelihood ratio test is performed, and the \eqn{PI = - X \beta}
#' is compared to \code{opt_cutoff} to define \emph{Short Survival} and
#' \emph{Long Survival} groups.
#'
#' The function also produces:
#' \itemize{
#'   \item Survival curves with group-specific colors,
#'   \item Risk tables (number-at-risk) aligned with survival curves,
#'   \item Distribution plots of the PI across groups.
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{df}: data frame with columns \code{PI} (prognostic index for each subject), \code{Y}, \code{delta},
#'                    \code{groupRisk} (assigned prognostic group based on \code{opt_cutoff}),
#'   \item \code{p_value}: from the log-rank test (\code{COXNet}) or likelihood ratio test (\code{AFTNet}),
#'                         measuring survival separation between groups.
#' }
#'
#' @seealso
#' \code{\link{OptimalPICutoff}} for \code{opt_cutoff} value selection.
#'
#' @export

  ValidationPI <- function(
      X, Y, delta,
      beta,
      opt_cutoff,
      model   = NULL,
      dist    = NULL,
      plot    = FALSE,
      palette = NULL
  ) {
    
  # --- Input checks ---
    
    model <- match.arg(model, choices = c("COXNet", "AFTNet"))
    
    if (model == "AFTNet") {
      if (is.null(dist)) stop("'dist' must be provided when model = 'AFTNet'.")
      dist <- match.arg(dist, c("weibull", "lognormal", "loglogistic"))
    }
    
    for (pkg in c("survminer", "ggplot2", "ggpubr", "grid", "survival")) {
      if (!requireNamespace(pkg, quietly = TRUE)) stop(paste0("Package '", pkg, "' is required."))
    }
    
    if (is.null(palette)) palette <- if(model == "COXNet") c("#CC899D", "#7DAE9B") else c("#7A78B7", "#CC899D")
    
  # ----- 
    
    PI <- as.numeric((switch(model,
                             "COXNet" = as.matrix(X) %*% beta,
                             "AFTNet" = -as.matrix(X) %*% beta)))
    
    PI_data <- data.frame(PI = PI, Y = Y, delta = delta)
    
    if (model == "COXNet") {
      PI_data$groupRisk <- factor(ifelse(PI_data$PI >= opt_cutoff, "High Risk", "Low Risk"),
                                  levels = c("High Risk", "Low Risk"))
    } else {
      PI_data$groupRisk <- factor(ifelse(PI_data$PI >= opt_cutoff, "Short Survival", "Long Survival"),
                                  levels = c("Short Survival", "Long Survival"))
    }
    
    palette_colors <- as.character(palette)
    
    p_value <- 1
    
    if (length(unique(PI_data$groupRisk)) > 1) {
      
      if (model == "COXNet") {
        logranktest <- survival::survdiff(survival::Surv(Y, delta) ~ groupRisk, data = PI_data)
        p_value     <- signif(1 - stats::pchisq(logranktest$chisq, 1), 4)
        fit_obj     <- survival::survfit(survival::Surv(Y, delta) ~ groupRisk, data = PI_data)
      }
      
      if (model == "AFTNet") {
        fit_aft  <- survival::survreg(survival::Surv(exp(Y), delta) ~ groupRisk, data = PI_data, dist = dist)
        null_fit <- survival::survreg(survival::Surv(exp(Y), delta) ~ 1, data = PI_data, dist = dist)
        p_value  <- signif(1 - stats::pchisq(2 * (as.numeric(stats::logLik(fit_aft)) - as.numeric(stats::logLik(null_fit))), 1), 4)
        
        fit_obj  <- survival::survfit(survival::Surv(Y, delta) ~ groupRisk, data = PI_data)
        
        if (plot) {
          y_range <- range(PI_data$Y, na.rm = TRUE)
          tt <- seq(y_range[1], y_range[2], length.out = 100)
          
          lp <- stats::predict(fit_aft, newdata = data.frame(groupRisk = levels(PI_data$groupRisk)), type = "lp")
          sc <- fit_aft$scale
          
          Sfun <- switch(dist,
                         "weibull"     = function(t, l, s) exp(-(exp(t) / exp(l))^(1/s)),
                         "lognormal"   = function(t, l, s) stats::plnorm(exp(t), meanlog = l, sdlog = s, lower.tail = FALSE),
                         "loglogistic" = function(t, l, s) 1 / (1 + (exp(t) / exp(l))^(1/s)))
          
          plot_df <- do.call(rbind, lapply(seq_along(lp), function(i) {
            data.frame(time = tt, surv = Sfun(tt, lp[i], sc), groupRisk = levels(PI_data$groupRisk)[i])
          }))
          
          aft_p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time, y = surv, color = groupRisk)) +
            ggplot2::geom_line(linewidth = 1.2) +
            ggplot2::scale_color_manual(values = palette_colors) +
            ggplot2::labs(y = "Survival probability", x = "Time (log-scale)", 
                          title = paste0(toupper(dist), " - AFTNET MODEL"), 
                          subtitle = paste0("p-value = ", p_value), color = "Group:") +
            ggplot2::theme_minimal() + 
            ggplot2::theme(legend.position = "top", 
                           plot.title = ggplot2::element_text(face = "bold", size = 10),
                           axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                           axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)))
        }
      }
      
      if (plot) {
        y_limits_real <- range(PI_data$Y, na.rm = TRUE)
        
        survp <- survminer::ggsurvplot(
          fit_obj, data = PI_data, palette = palette_colors,
          risk.table = TRUE, risk.table.col = "strata",
          xlim = y_limits_real, 
          tables.theme = survminer::theme_cleantable(),
          ggtheme = ggplot2::theme_minimal(),
          legend.title = "Group:",
          legend.labs = levels(PI_data$groupRisk)
        )
        
        if (model == "AFTNet") {
          survp$plot <- aft_p
        } else {
          survp$plot <- survp$plot + 
            ggplot2::labs(title = "COXNET MODEL", subtitle = paste0("p-value = ", p_value), 
                          color = "Group:") +
            ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                           axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)))
        }
        
        table_min <- floor(min(PI_data$Y, na.rm = TRUE))
        table_max <- ceiling(max(PI_data$Y, na.rm = TRUE) * 2) / 2
        if(table_max - max(PI_data$Y) < 0.1) table_max <- table_max + 0.5
        table_breaks <- seq(table_min, table_max, by = 1)
        
        survp$table$scales$scales[[1]]$limits <- c(table_min, table_max)
        survp$table$scales$scales[[1]]$breaks <- table_breaks
        survp$table$coordinates <- ggplot2::coord_cartesian(xlim = c(table_min, table_max), expand = TRUE)
        
        x_label <- if (model == "AFTNet") "Time (log-scale)" else "Time"
        
        survp$table <- survp$table + 
          ggplot2::theme(
            axis.text.y  = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.x  = ggplot2::element_text(), 
            axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
            legend.position = "none",
            plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0)
          ) +
          ggplot2::labs(y = NULL, x = x_label, title = "SAMPLES AT RISK")
        
        PI_data_sorted <- PI_data[order(PI_data$PI), ]
        PI_data_sorted$id <- 1:nrow(PI_data_sorted)
        
        d_plot <- ggplot2::ggplot(PI_data_sorted, ggplot2::aes(x = id, y = PI, color = groupRisk)) +
          ggplot2::geom_point(alpha = 0.8) +
          ggplot2::scale_color_manual(values = palette_colors) +
          ggplot2::labs(x = "Sample Index", y = "Prognostic Index value", title = "Prognostic Groups Distribution") +
          ggplot2::theme_minimal() + 
          ggplot2::theme(legend.position = "none",
                         axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                         axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)))
        
        plot_a_arranged <- survminer::arrange_ggsurvplots(list(survp), print = FALSE, ncol = 1, nrow = 1)
        plot_a_grob <- grid::grid.grabExpr(print(plot_a_arranged))
        
        final_comp <- ggpubr::ggarrange(
          plot_a_grob, d_plot, 
          labels = c("A", "B"), 
          ncol = 2, widths = c(1.2, 0.8)
        )
        
        grid::grid.newpage()
        grid::grid.draw(final_comp)}
      
    } else {
      warning("Only one group identified: statistics and plotting skipped.")
    }
    
    return(list(
      df      = PI_data,
      p_value = p_value))
  }