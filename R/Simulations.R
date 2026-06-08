#' Simulate Transcription Factor (TF) Target Gene Networks with Survival Outcomes
#'
#' @description
#' Generates structured gene expression data based on TFs and their regulated
#' target genes, together with survival outcomes simulated from \code{COXNet}
#' or \code{AFTNet} models. The function supports both **independent** and **interconnected**
#' TF modules with user-defined shared targets via \code{shared_scheme}.
#'
#' @param n Numeric value of observations.
#' @param r Numeric value of TFs (for interconnected modules, at least 4 TFs are recommended).
#' @param targets Numeric value of target genes regulated by each TF.
#' @param p_active Numeric value of truly active predictors (non-zero coefficients).
#' @param rho Numeric value of correlation between each TF and its target (\code{default: 0.70}).
#' @param rate Numeric value of desired censoring proportion (\code{default: 0.50}).
#' @param b_true Numeric vector of length 4 \code{(pos_min, pos_max, neg_min, neg_max)}
#'               used to generate positive and negative non-zero coefficients. 
#' @param nsimul Numeric value of simulated datasets (\code{default: 10}).
#' @param model Character string specifying the survival model used for simulation
#'              (\code{"COXNet"}, or \code{"AFTNet"}).
#' @param baseline Character string specifying baseline hazard distribution.
#'                 \itemize{
#'                     \item For \code{COXNet}: exponential (\code{"exp"}), Weibull
#'                           (\code{"weibull"}), or piecewise-constant (\code{"piecewise"}).
#'                     \item For \code{AFTNet}: Weibull (\code{"weibull"}), Log-Normal
#'                           (\code{"lognormal"}), or Log-Logistic (\code{"loglogistic"}).
#'                 }
#' @param phi Numeric value of frailty parameter for \code{COXNet}'s baselines
#'            (required for \code{"exp"} and \code{"weibull"}).
#' @param sigma_true Positive numeric scalar representing the scale parameter of the
#'                   error distribution in \code{AFTNet} model (\code{default: 1}).
#' @param breaks Numeric vector of time breakpoints for piecewise exponential hazards
#'               (required if \code{baseline = "piecewise"}, \code{default: c(0, 6, 36, 60)}).
#' @param hazards Numeric vector of hazard rates corresponding to each interval in \code{breaks}
#'                (\code{default: c(0.15, 0.005, 0.1)}).
#' @param shared_scheme List defining interconnected TF modules. If \code{NULL} (default),
#'                      TFs regulate disjoint target sets (independent structure).
#'                      Otherwise, it must be a list of modules, each containing
#'                      \itemize{
#'                          \item \code{tfs}: integer vector of TF indices in the module,
#'                          \item \code{shared}: number of genes shared among those TFs,
#'                          \item \code{unique}: integer vector giving the number of TF-specific targets.
#'                      }
#' @param choice Value specifying the choice for the signs of the adjacency matrix
#'               \itemize{
#'                  \item 1 (default): for correlation-based signs,
#'                  \item 2: for ridge-based signs (see \code{\link{CreateNetwork}} for details).
#'               }
#' @param save Logical value, if \code{TRUE} each simulated dataset is saved as an
#'             \code{.rds} file in the directory specified by \code{save_path}
#'             (\code{default: FALSE}).
#' @param save_path Character string specifying an existing directory used only when
#'                  \code{save = TRUE}. No files are written by default.       
#' @param seed Random seed for reproducibility (\code{default: 2026}).
#' @param verbose Logical value, if \code{TRUE} progress and summary messages are
#'                printed during simulation (\code{default: FALSE}).
#'
#' @details
#' The total number of predictors is given by \eqn{p = r \times (targets + 1)},
#' where each TF contributes one regulatory variable in addition to its associated
#' target genes.
#' 
#' The function supports two alternative network topologies
#' \itemize{
#'      \item **Independent structure**: each TF regulates its own targets independently.
#'      \item **Interconnected structure**: TFs specified in the same \code{shared_scheme}
#'      share \code{shared} genes and additionally have their own unique genes
#'      as specified in \code{unique}.
#' }
#' 
#' These regulatory relationships are encoded in the adjacency matrix, which exhibits a
#' block-diagonal structure under independence, and introduces cross-connections between TFs
#' and shared targets when modules are specified.
#' 
#' Survival times are generated according to the chosen baseline distribution and
#' linear predictors derived from the simulated gene expression data.
#' Optional frailty effects and censoring are included, with the censoring mechanism
#' calibrated to achieve the desired censoring proportion specified by \code{rate}.
#' 
#' The function also returns the true regression coefficients, allowing the user to evaluate
#' variable selection performance using measures such as false positive and false negative rates.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{X_list}: list of simulated design matrices.
#'   \item \code{beta_list}: list of true regression coefficient vectors.
#'   \item \code{time_list}: list of observed survival times (log-transformed under \code{AFTNet}).
#'   \item \code{delta_list}: list of censoring indicators (1 = event, 0 = censored).
#'   \item \code{L_list}: list of Laplacian matrices representing the TF–gene regulatory network.
#' }
#'
#' @examples
#'   
#'   # - Simulate interconnected structure under Weibull-COXNet model -
#'   
#'     targets <- 10
#'     s1 <- 5
#'     s2 <- 3
#'     
#'     shared_scheme <- list( 
#'     list(tfs = c(1, 3), shared = s1, unique = c(targets - s1, targets - s1)),  
#'     list(tfs = c(2, 4), shared = s2, unique = c(targets - s2, targets - s2)))
#'   
#'     simul_data <- Simulations(
#'     n = 165, r = 40, targets = targets, p_active = 40, 
#'     b_true = c(0.8,1.2,-1.2,-0.8),
#'     rate = 0.3, nsimul = 1,
#'     model = "COXNet", baseline = "weibull",
#'     shared_scheme = shared_scheme,
#'     seed = 2026, verbose = FALSE)
#'         
#'   # Extract the Laplacian matrix
#'   
#'     L <- simul_data$L[[1]]
#'   
#'   # This matrix uncovers the topological overlap between TFs:
#'   # TF1 and TF3 co-regulate 5 genes, while TF2 and TF4 share 3 target genes.
#' 
#' @export

  Simulations <- function(
      n, r, targets, p_active,
      rho        = 0.70,
      rate       = 0.50,
      b_true     = c(0.8, 1.2, -1.2, -0.8),
      nsimul     = 10,
      model      = NULL,
      baseline   = NULL,
      phi        = 0.1,
      sigma_true = 1,
      breaks     = c(0, 6, 36, 60),
      hazards    = c(0.15, 0.005, 0.1),
      shared_scheme = NULL,
      choice     = 1,
      save       = FALSE,
      save_path  = NULL,
      seed       = 2026,
      verbose    = FALSE
  ) {
    
  # --- Input checks ---
    
    if (!is.null(model)) {
      model <- match.arg(model, choices = c("COXNet", "AFTNet"))
    } else {
      stop("'model' must be specified.")
    }
    
    if (!is.numeric(choice) || length(choice) != 1L || is.na(choice) || !choice %in% c(1, 2)) {
      stop("'choice' must be 1 or 2.")
    }
    
    if (!is.null(shared_scheme)) {
      for (i in seq_along(shared_scheme)) {
        mod <- shared_scheme[[i]]
        if (!all(c("tfs", "shared", "unique") %in% names(mod))) {
          stop("Each element of 'shared_scheme' must contain 'tfs', 'shared', and 'unique'.")}}
   
      if (any(mod$unique + mod$shared != targets)) {
        stop("For each TF in a shared module, unique + shared must equal 'targets'.")
      }
    }
    
    if (is.null(baseline)) stop("'baseline' must be specified.")
    
    if (!is.character(baseline) || length(baseline) != 1L) {
      stop("'baseline' must be a single character string.")
    }
    
    if (model == "COXNet") {
      baseline <- match.arg(baseline, c("exp", "weibull", "piecewise"))
    } else {
      baseline <- match.arg(baseline, c("weibull", "lognormal", "loglogistic"))
    }
    
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n <= 0)
      stop("'n' must be a positive integer.")
    
    if (!is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0)
      stop("'r' must be a positive integer.")
    
    if (!is.numeric(targets) || length(targets) != 1L || is.na(targets) || targets <= 0)
      stop("'targets' must be a positive integer.")
    
    if (!is.numeric(p_active) || length(p_active) != 1L || is.na(p_active) || p_active <= 0)
      stop("'p_active' must be a positive integer.")
    
    if (!is.numeric(rate) || length(rate) != 1L || is.na(rate) || rate <= 0 || rate >= 1) 
      stop("'rate' must be a single numeric value between 0 and 1.")
    
    if (!is.numeric(b_true) || length(b_true) != 4L || anyNA(b_true))
      stop("'b_true' must be a numeric vector of length 4.")
    
    if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || rho <= -1 || rho >= 1)
      stop("'rho' must be a single numeric value between -1 and 1.")
    
    if (!is.numeric(nsimul) || length(nsimul) != 1L || is.na(nsimul) || nsimul < 1 || nsimul != as.integer(nsimul)) 
      stop("'nsimul' must be a positive integer.")
    
    if (!is.logical(save) || length(save) != 1L || is.na(save)) stop("'save' must be TRUE or FALSE.")
    
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("'verbose' must be TRUE or FALSE.")
    
  # -----
    
    set.seed(seed)
    
    p_module <- targets + 1
    
    p <- r * p_module
    
    # - TF grouping -
    
    if (is.null(shared_scheme) || length(shared_scheme) == 0) {
      shared_tfs <- integer(0)
      structure_type <- "independent"
    } else {
      shared_tfs <- unique(unlist(lapply(shared_scheme, `[[`, "tfs")))
      structure_type <- "interconnected"
    }
    
    remaining_tfs <- setdiff(seq_len(r), shared_tfs)
    
    # - Generate regression coefficients -
    
    beta <- numeric(p)
    m <- floor(p_active / 2)
    
    beta[1:m] <- stats::runif(m, b_true[1], b_true[2])
    beta[(m + 1):p_active] <- stats::runif(p_active - m, b_true[3], b_true[4])
    
    tf_blocks <- vector("list", r)
    for (i in seq_len(r)) {
      start <- (i - 1) * p_module + 1
      tf_blocks[[i]] <- list(tf = start, genes = (start + 1):(start + targets))
    }
    
    if (!is.null(shared_scheme)) {
      for (m_idx in seq_along(shared_scheme)) {
        mod <- shared_scheme[[m_idx]]
        ref_tf <- mod$tfs[1]
        
        shared_indices <- tail(tf_blocks[[ref_tf]]$genes, mod$shared)
        for (k in 1:length(mod$tfs)) {
          target_tf <- mod$tfs[k]
          u_count <- mod$unique[k]
          
          tf_blocks[[target_tf]]$genes <- c(tf_blocks[[target_tf]]$genes[1:u_count], shared_indices)
        }
      }
    }
    
  # - Adjacency matrix -
    
    adj <- matrix(0, nrow = p, ncol = p)
    
    # Regulatory directions: +1 for activation, -1 for repression

    network_signs <- lapply(seq_len(r), function(i) sample(c(1, -1), size = targets, replace = TRUE))
    
    for (i in seq_len(r)) {
      tf_idx    <- tf_blocks[[i]]$tf
      gene_idxs <- tf_blocks[[i]]$genes
      signs     <- network_signs[[i]]

      for (g in seq_along(gene_idxs)) {
        adj[tf_idx, gene_idxs[g]] <- adj[gene_idxs[g], tf_idx] <- signs[g]
      }
    }
    
  # -- Data Generation --
    
    X_list <- beta_list <- time_list <- delta_list <- L_list <- vector("list", nsimul)

    for (k in seq_len(nsimul)) {
      
      X <- matrix(0, nrow = n, ncol = p)
      
      for (j in seq_len(n)) {
        TF_vals <- stats::rnorm(r, 0, 1)
        for (i in seq_len(r)) { X[j, tf_blocks[[i]]$tf] <- TF_vals[i] }
        
        if (length(shared_tfs) > 0) {
          for (m in seq_along(shared_scheme)) {
            
            mod <- shared_scheme[[m]]
            phi_C <- (1 / sqrt(length(mod$tfs))) * sum(TF_vals[mod$tfs])
            
            shared_signs <- tail(network_signs[[mod$tfs[1]]], mod$shared)
            
            genes_shared <- stats::rnorm(mod$shared, mean = shared_signs * rho * phi_C, sd = sqrt(1 - rho^2))
            
            for (tf_id in mod$tfs) {
              
              u_count <- mod$unique[which(mod$tfs == tf_id)]
              u_signs <- network_signs[[tf_id]][1:u_count]
              genes_unique <- stats::rnorm(u_count, (u_signs * rho) * TF_vals[tf_id], sqrt(1 - rho^2))
              
              X[j, tf_blocks[[tf_id]]$genes[1:u_count]] <- genes_unique
              X[j, tail(tf_blocks[[tf_id]]$genes, mod$shared)] <- genes_shared}}
        }
        
        for (i in remaining_tfs) {
          signs <- network_signs[[i]]
          genes <- stats::rnorm(targets, mean = (signs * rho) * TF_vals[i], sd = sqrt(1 - rho^2))
          
          X[j, tf_blocks[[i]]$genes] <- genes}
      }
      
      X <- apply(X, 2, function(col) {
        if (stats::sd(col) < 1e-7) col <- col + stats::rnorm(n, 0, 1e-5)
        return(col)
      })
      
      X <- scale(X, center = TRUE, scale = apply(X, 2, function(x) max(abs(x)) + 1e-6))
      
      surv <- generate_surv_data(X, beta, model, baseline, rate, phi, sigma_true, breaks, hazards, verbose)
      L <- signed_laplacian(choice, X, surv$times_obs, surv$delta, adj, model, baseline)
      
      X_list[[k]] <- X; beta_list[[k]] <- beta; time_list[[k]] <- surv$times_obs
      delta_list[[k]] <- surv$delta; L_list[[k]] <- L
      
      if (save) {
        
        if (is.null(save_path)) {
          stop("If save=TRUE, 'save_path' must be provided.")}
        if (!is.character(save_path) || length(save_path) != 1L || !dir.exists(save_path)) {
          stop("'save_path' must be an existing directory.")}
        
        saveRDS(
          list(
            X     = X,
            time  = time,
            delta = delta,
            L     = L_list[[k]]),
          file    = file.path(
            save_path, paste0("simul_", k, ".rds")))
      }
      if (isTRUE(verbose)) message("Simulation ", k, " of ", nsimul, " completed.\n")
    }
    
    structure_type <- if (length(shared_scheme) == 0) "independent" else "interconnected"
    
    msg <- "All data simulated."
    
    if (isTRUE(verbose)) {
      msg <- paste0(msg, "\n",
                    "    Model: ", model, " - ", baseline, "\n",
                    "    Structure: ", structure_type, "\n",
            sprintf("    Censoring rate: %.0f%%", rate * 100))
    }
    
    message(msg)
    
    return(list(
      X_list     = X_list,
      beta_list  = beta_list,
      time_list  = time_list,
      delta_list = delta_list,
      L_list     = L_list))
 }

# ============================================================================

#' Simulate Survival Outcomes under \code{COXNet} or \code{AFTNet} Models
#' 
#' @keywords internal
#' @noRd

  generate_surv_data <- function(
      X, beta,
      model, baseline,
      rate, phi, sigma_true,
      breaks, hazards,
      verbose
  ) {
    
    n       <- nrow(X)
    eta     <- as.vector(X %*% beta)
    exp_eta <- exp(eta)
    
    shape   <- 1 / sigma_true
    U       <- stats::runif(n)
    
    if (model == "COXNet") {
      
      baseline <- match.arg(baseline, c("exp", "weibull", "piecewise"))
      
      if (baseline == "exp") {
        
        T_true <- -log(U) / (phi * exp_eta)
        
      } else if (baseline == "weibull") {
        
        T_true <- phi * ((-log(U)) / exp_eta)^(1 / shape)
        
      } else if (baseline == "piecewise") {
        
        if (is.null(breaks) || is.null(hazards))
          stop("For piecewise baseline, provide both 'breaks' and 'hazards'.")
        
        K <- length(hazards)
        
        if (length(breaks) != K + 1)
          stop("Length of 'breaks' must be one more than 'hazards'.")
        
        T_true <- numeric(n)
        
        for (i in 1:n) {
          
          z <- -log(U[i]) / exp_eta[i]
          Hcum <- 0
          
          for (j in 1:K) {
            
            Hnext <- Hcum + hazards[j] * (breaks[j+1] - breaks[j])
            
            if (z <= Hnext) {
              T_true[i] <- breaks[j] + (z - Hcum) / hazards[j]
              break
            }
            Hcum <- Hnext
          }
          
          if (T_true[i] == 0)
            T_true[i] <- breaks[K + 1] + (z - Hcum) / hazards[K]}}
    }
    
    else if (model == "AFTNet") {
      
      baseline <- match.arg(baseline, c("weibull", "lognormal", "loglogistic"))
      
      if (baseline == "weibull") {
        
        T_true <- stats::rweibull(n, shape = shape, scale = exp(eta))
        
      } else if (baseline == "lognormal") {
        
        T_true <- stats::rlnorm(n, meanlog = eta, sdlog = sigma_true)
        
      } else if (baseline == "loglogistic") {
        
        if (!requireNamespace("flexsurv", quietly = TRUE)) stop("Package 'flexsurv' is required.")
        T_true <- flexsurv::rllogis(n, shape = shape, scale = exp_eta)}
    }
    
  # - Censoring -
    
    theta     <- 1 / stats::quantile(T_true, 1 - rate)
    censor    <- stats::rexp(n, rate = theta)
    delta     <- as.numeric(T_true <= censor)
    times_obs <- pmin(T_true, censor)
    
    if (model == "AFTNet") times_obs <- log(times_obs)
    
    if (isTRUE(verbose)) {
      message(sprintf("(Achieved censoring rate: %.1f%%)\n",100 * (1 - mean(delta))))
    }
    
    return(list(
      times_obs = times_obs,
      delta     = delta))
  }
