#' Variables Screening Methods Based on Prior Knowledge and Marginal Utility
#'
#' @description
#' Reduces the high-dimensional feature space to a more manageable
#' subset of variables by applying one of three screening strategies:
#' \itemize{
#'   \item **BMD (Biomedical-driven)**: selects covariates based on prior
#'         biomedical knowledge about their relevance to the disease under investigation,
#'   \item **DAD (Data-driven)**: selects features using component-wise
#'         estimators obtained from the chosen penalized model,
#'   \item **BMD+DAD**: combines both biomedical knowledge and data-driven insights.
#' }
#'
#' @param X Numeric matrix of covariates.
#' @param Y Numeric vector of observed survival times (log-transformed under \code{AFTNet}).
#' @param delta Integer vector of censoring indicators (1 = event, 0 = censored).
#' @param disease_genes Character vector containing the names of genes known
#'                      to be associated with diseases.
#' @param screening Character string specifying the screening method
#'                  (\code{"BMD"}, \code{"DAD"}, or \code{"BMD+DAD"}).
#' @param model Character string specifying the fitted survival model
#'              (\code{"COXNet"}, or \code{"AFTNet"}) required for DAD-based screening.
#' @param dist Character string specifying the AFTNet distribution.
#'             Must be one of \code{"weibull"}, \code{"lognormal"}, or \code{"loglogistic"}.
#' @param rank_method Character string specifying the ranking criterion for DAD-based screening:
#'                    \code{"absmg"} (absolute marginal coefficients),
#'                    \code{"mg"} (marginal function), or \code{"mgpadj"}
#'                    (adjusted p-value from the marginal function).
#' @param d Numeric value representing the threshold for top-ranked features to select
#'          in DAD-based screening (\code{default: NULL}).
#' @param standardize Logical value indicating whether to standardize the input matrix
#'                    in DAD-based screening:
#'                    \itemize{
#'                        \item if \code{TRUE} (default) each column is centered to have
#'                              mean 0 and scaled to have unit variance.
#'                        \item if \code{FALSE} the function assumes that the matrix has
#'                              already been standardized by the user.
#'                    }
#' @param verbose Logical value, if \code{TRUE} progress messages are printed (\code{default: FALSE}).
#'
#' @details
#' The function uses marginal ranking approaches to select features based on
#' their association with survival outcomes.
#' \itemize{
#'    \item In the **BMD** approach, prior knowledge comes from literature or
#'          external biological databases such as \href{https://hb.flatironinstitute.org/}{HumanBase}.
#'    \item The **DAD** screening computes marginal regression coefficients to
#'          rank features according to their estimated importance
#'          under the selected model:
#'          \itemize{
#'              \item \code{absmg}: top \code{d} covariates by largest
#'                                  absolute marginal coefficients.
#'              \item \code{mg}: top \code{d} covariates by largest marginal
#'                               coefficients, preserving the direction.
#'              \item \code{mgpadj}: top \code{d} covariates passing significance
#'                                   thresholds based on adjusted p-values.}
#'    \item The **BMD+DAD** combines prior biological knowledge and data-driven
#'          selection for comprehensive feature screening.
#' }
#' 
#' @return A list containing selected variable names \code{screen_vars}.
#'
#' @seealso
#' \code{\link{CreateNetwork}} or \code{\link{RepositoryDisease}} for the \code{disease_genes} names.
#'
#' @export

  VariableScreening <- function(
      X, Y, delta,
      disease_genes,
      screening   = NULL,
      model       = NULL,
      dist        = NULL,
      rank_method = NULL, 
      d           = NULL,
      standardize = TRUE,
      verbose     = FALSE
  ) {
    
    screening   <- match.arg(screening, choices = c("BMD", "DAD", "BMD+DAD"))
    
    if (is.null(colnames(X))) stop("'X' must have column names.")
    
    if (!is.null(model)) {
      
      model <- match.arg(model, choices = c("COXNet", "AFTNet"))
     
      if (model == "AFTNet" && is.null(dist)) stop("Argument 'dist' must be specified for AFTNet.")
      
    }
    
    if (!is.null(dist)) dist <- match.arg(dist, choices = c("weibull", "lognormal", "loglogistic"))
    
    if (screening != "BMD") {
      
      if (is.null(model)) stop("'model' is required for DAD-based screening.")
      
      rank_method <- match.arg(rank_method, c("absmg", "mg", "mgpadj"))
      
      if (is.null(d) || !is.numeric(d) || length(d) != 1 || is.na(d) || d < 1) {
        stop("Argument 'd' must be a positive integer.")
      }
      
      d <- as.integer(d)
    }
    
    screen_vars <- NULL
    
  # --- BMD screening ---
    
    if (screening == "BMD") {
      
      screen_vars <- intersect(disease_genes, colnames(X))
      
      if (isTRUE(verbose)) message("Number of BMD-screened genes: ", length(screen_vars))
      
    }
    
  # --- DAD screening ---
    
    if (screening == "DAD") {
      
      ranking <- marginal_ranking(X = X, Y = Y, delta = delta, model = model, dist = dist,
                                  rank_method = rank_method, standardize = standardize, verbose = verbose)
    
      screen_vars <- ranking$symbol[seq_len(d)]
      
      if (isTRUE(verbose)) message("Number of DAD-screened genes: ", length(screen_vars))
      
    }
    
  # --- BMD + DAD screening ---
    
    if (screening == "BMD+DAD") {
      
      ranking <- marginal_ranking(X = X, Y = Y, delta = delta, model = model, dist = dist,
                                  rank_method = rank_method, standardize = standardize, verbose = verbose)
      
      dad_vars <- ranking$symbol[seq_len(d)]
      
      screen_vars <- intersect(union(disease_genes, dad_vars), colnames(X))
      
      if (isTRUE(verbose)) message("Number of BMD+DAD-screened genes: ", length(screen_vars))}
    
    return(list(
      screen_vars = screen_vars))
  }

# ============================================================================
# Supporting functions
# ============================================================================

#' Marginal Ranking of Variables Using Survival Models
#' 
#' @keywords internal
#' @noRd

  marginal_ranking <- function(
      X, Y, delta,
      model, dist,
      rank_method, standardize,
      verbose
  ) {
    
    if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required.", call. = FALSE)
    
    vars <- apply(X, 2, stats::var, na.rm = TRUE)
    X <- X[, vars > 1e-10, drop = FALSE]
    
    S <- survival::Surv(Y, delta)
    
    X_std <- if (isTRUE(standardize)) standardize(X)$std else X
    
    margcoef <- marg_coef(X = X_std, S, model = model,
                          dist = dist, rank_method = rank_method, verbose = verbose)
    
    if (rank_method == "mgpadj") {

      margcoef <- stats::p.adjust(margcoef, method = "fdr")
      rankcoef <- sort(margcoef, decreasing = FALSE, index.return = TRUE)
      
    } else {

      rankcoef <- sort(margcoef, decreasing = TRUE, index.return = TRUE)
      
    }
    
    ranktable <- data.frame(
      ranking = seq_along(rankcoef$ix),
      symbol = colnames(X)[rankcoef$ix],
      score = margcoef[rankcoef$ix],
      stringsAsFactors = FALSE)
    
    return(ranktable)
  }

# ============================================================================

#' Marginal Utility Calculation
#'  
#' @keywords internal
#' @noRd

  marg_coef <- function(
      X, S, model,
      dist, rank_method, verbose
  ) {
    
    p <- ncol(X)
    
    if (isTRUE(verbose)) message("Marginal coefficient computation using method: ", rank_method)
    
    FUN <- function(i) {
      mg(i, X, S, model, dist)
    }
    
    if (rank_method == "absmg") {
      
      v <- vapply(seq_len(p), FUN, numeric(1))
      return(abs(v))
      
    } else if (rank_method == "mg") {
      return(vapply(seq_len(p), FUN, numeric(1)))
      
    } else if (rank_method == "mgpadj") {
      return(vapply(seq_len(p), function(i)
        pval(i, X, S, model, dist),
        numeric(1)
      ))
    }
  }

# ============================================================================

#' Marginal Model Coefficient
#' 
#' @keywords internal
#' @noRd

  mg <- function(
      index, X, S,
      model, dist = NULL
  ) {
    
    x <- X[, index]
    
    if (model == "COXNet") {
      
      fit <- tryCatch({
        survival::coxph(S ~ x, model = FALSE)
      }, warning = function(w) {
        return(NULL)
      }, error = function(e) return(NULL))
      
      if (is.null(fit)) return(0) 
      
      res <- stats::coef(fit)

      return(if (is.finite(res) && abs(res) < 15) res else 0)
      
    } else {
      
      fit <- tryCatch({
        survival::survreg(S ~ x, dist = dist, 
                          control = survival::survreg.control(iter.max = 100))
      }, warning = function(w) {
        return(NULL)
      }, error = function(e) {
        return(NULL)
      })
      
      if (is.null(fit)) return(0)
      
      res <- stats::coef(fit)[2]
      return(if (is.finite(res)) res else 0)
    }
  }

# ============================================================================

#' Marginal Model P-Value
#' 
#' @keywords internal
#' @noRd

  pval <- function(
      index, X, S,
      model, dist = NULL
  ) {
    
    x <- X[, index]
    
    if (model == "COXNet") {
      fit <- tryCatch({
        survival::coxph(S ~ x)
      }, warning = function(w) return(NULL), 
      error = function(e) return(NULL))
      
      if (is.null(fit)) return(1)
      
      s <- summary(fit)
      p_val <- s$coefficients[1, "Pr(>|z|)"]
      
      return(if (is.finite(p_val)) p_val else 1)
      
    } else {
      
      fit <- tryCatch({
        survival::survreg(S ~ x, dist = dist, 
                          control = survival::survreg.control(iter.max = 100))
      }, warning = function(w) return(NULL), 
      error = function(e) return(NULL))
      
      if (is.null(fit)) return(1)
      
      s <- summary(fit)
      
      p_val <- s$table[2, ncol(s$table)]
      
      return(if (is.finite(p_val)) p_val else 1)
    }
  }
