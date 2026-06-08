#' Laplacian Matrix for Prior Biological Knowledge in Network Constraint
#'
#' @description
#' Builds a Laplacian network penalty based on a prior weighted graph.
#' It encourages coefficients corresponding to connected covariates to behave similarly:
#' if two covariates are strongly connected in the network, their estimated coefficients
#' tend to be either both close to zero or both nonzero. In this way, the penalty promotes
#' smoothness and structural coherence across related variables.
#'
#' @param X Numeric matrix of standardized covariates.
#' @param Y Numeric vector of observed survival times (log-transformed under \code{AFTNet}),
#'          required for \code{choice = 2}.
#' @param delta Integer vector of censoring indicators (1 = event, 0 = censored),
#'              required for \code{choice = 2}.
#' @param doid Character string specifying Disease Ontology ID (\code{"DOID:XXXX"}),
#'             used only if \code{disease_file} is not provided.
#' @param tissue Character string specifying tissue name, used to retrieve the
#'               tissue-specific network from HumanBase, used only if
#'               \code{tissue_file} is not provided.
#' @param disease_file Character string specifying optional path to a tab-delimited
#'                     file containing disease-associated genes (columns: \code{entrez_id},
#'                     \code{standard_name}, and \code{score}).
#' @param tissue_file Character string specifying optional path to a tab-delimited
#'                    file with tissue-specific gene interactions (columns:
#'                    \code{gene1}, \code{gene2}, and \code{score}).
#' @param cache Logical value; if \code{TRUE}, downloaded HumanBase files are cached
#'              for reuse in \code{cache_dir}. If \code{FALSE} (default), files are downloaded
#'              for the current session only.
#' @param cache_dir Character string specifying a directory used to cache
#'                  downloaded HumanBase files (when \code{cache = TRUE}).
#' @param choice Value specifying the choice for the signs of the adjacency matrix
#'               \itemize{
#'                  \item \code{1} (default): for correlation-based signs.
#'                  \item \code{2}: for ridge-based signs.
#'               }
#' @param model Character string specifying the fitted survival model
#'              (\code{"COXNet"}, or \code{"AFTNet"}) required only for \code{choice = 2}.
#' @param dist Character string specifying the AFTNet distribution.
#'             Must be one of \code{"weibull"}, \code{"lognormal"}, or \code{"loglogistic"}
#'             (required only for \code{choice = 2}).
#' @param verbose Logical value, if \code{TRUE} progress messages are printed.
#'
#' @details
#' This prior network is represented by a weighted graph where each vertex
#' corresponds to a covariate and the edges describe relationships between covariates.
#' The edge weights are stored in an adjacency matrix \eqn{A}, which has zeros on its diagonal.
#' The degree matrix \eqn{D} contains on its diagonal the sum of the absolute
#' edge weights connected to each vertex. The Laplacian matrix is defined as \eqn{L = D - W},
#' where \eqn{W} is the weighted matrix estimated from \eqn{A}. 
#' Two strategies can be used.
#' \itemize{
#'    \item Correlation-based signs (\code{choice = 1}): the sign of an edge is
#'          set according to the Pearson correlation between the two corresponding covariates.
#'    \item Ridge-based signs (\code{choice = 2}): the sign of an edge is determined
#'          by the signs of ridge regression coefficients obtained from a penalized
#'          survival model. This ridge estimator provides stable coefficient estimates
#'          in high-dimensional settings. For the Cox model the ridge fit is obtained
#'          via \code{glmnet::glmnet()}), while for the AFT model via \code{survival::survreg()}).
#' }
#'
#' The framework is used to construct a disease-specific gene interaction network,
#' where edges represent biological relationships between genes relevant to a
#' given cancer and tissue type.
#'
#' Internally, the function relies on helper routines (see \code{\link{RepositoryDisease}} and \code{\link{RepositoryTissue}})
#' to retrieve biological prior information from the \href{https://hb.flatironinstitute.org/}{HumanBase} database.
#' These datasets are combined to construct a disease- and tissue-specific adjacency matrix
#' that defines the structure of the Laplacian penalty. User-provided files with
#' the same format can be supplied to bypass the download step.
#'
#' @note If tissue-specific or disease-specific files are not provided, the function
#' downloads the relevant data from HumanBase. In this case, an active internet
#' connection is required. Moreover, not all DOIDs and tissues are present in the
#' HumanBase repository. f the requested is not available, the function may return
#' an empty list.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{disease_genes}: data frame of disease genes used in the network.
#'   \item \code{L}: final Laplacian matrix.
#' }
#' 
#' @examples
#'   \donttest{
#'   
#'     data(LUADdataset)
#'   
#'     net <- CreateNetwork(
#'               LUADdataset$X_train,
#'               doid    = "DOID:1324",
#'               tissue  = "lung",
#'               choice  = 1,
#'               verbose = TRUE)
#'               
#'     L   <- net$L                          # final laplacian matrix
#'   
#'     disease_genes <- net$disease_genes    # disease genes and scores
#'   
#'   }
#'   
#' @name CreateNetwork   
#'
#' @export

  CreateNetwork <- function(
      X,
      Y            = NULL,
      delta        = NULL,
      doid         = NULL,
      tissue       = NULL,
      disease_file = NULL,
      tissue_file  = NULL,
      cache        = FALSE,
      cache_dir    = NULL,
      choice       = 1,
      model        = NULL,
      dist         = NULL,
      verbose      = FALSE
  ) {
    
  # --- Input checks ---  
    
    required_pkgs <- c("igraph", "magic", "httr","curl", "dplyr", "survival", "glmnet")
    
    for (pkg in required_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.", call. = FALSE)
    }
    
    if (!choice %in% c(1, 2)) stop("'choice' must be 1 or 2.")
    
    if (choice == 2) {
      if  (is.null(model)) stop("Parameter 'model' must be provided when 'choice = 2'.")
      if (model == "AFTNet" && is.null(dist)) stop("Argument 'dist' must be specified for AFTNet model.")
      if (!is.null(dist)) dist <- match.arg(dist, choices = c("weibull", "lognormal", "loglogistic"))
    }
  
  # -----  
  
    cache <- as.logical(cache)
    
    if (cache && !dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (isTRUE(verbose)) {message("Downloading the list of disease-associated genes...")}
    
    DiseaseGenes <- if (!is.null(disease_file)) {
      utils::read.table(disease_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    } else RepositoryDisease(doid = doid, cache = cache, cache_dir = cache_dir, verbose = verbose)
    
    tissueSpecificEdge <- if (!is.null(tissue_file)) {
      utils::read.table(tissue_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    } else RepositoryTissue(tissue = tissue, cache = cache, cache_dir = cache_dir, verbose = verbose)
    
    DiseaseGenesEntrezID <- as.character(DiseaseGenes[ , 1])
    
    idx <- which(tissueSpecificEdge[ , 1] %in% DiseaseGenesEntrezID &
                   tissueSpecificEdge[ , 2] %in% DiseaseGenesEntrezID)
    
    networkDataSubset <- tissueSpecificEdge[idx, ]
    
    if (nrow(networkDataSubset) == 0) {
      stop("No overlapping genes found between disease genes and tissue-specific network.")
    }
    
    if (isTRUE(verbose)) {
      message(sprintf("Your list contains %d genes.", length(DiseaseGenesEntrezID)))
      message(sprintf("%s of them are included in the tissue-specific network.",
                    length(unique(union(networkDataSubset[ , 1], networkDataSubset[ , 2])))))
    }
    
    g <- igraph::graph_from_data_frame(networkDataSubset, directed = FALSE)
    g <- igraph::set_edge_attr(g, "weight", value = networkDataSubset[ , 3])
    
    if (isTRUE(verbose)) message("Building the network...")
    
    NetworkMatrix <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
    
    symbol <- DiseaseGenes[match(rownames(NetworkMatrix), DiseaseGenesEntrezID), 2]
    colnames(NetworkMatrix) <- symbol
    rownames(NetworkMatrix) <- symbol
    
    NetworkMatrix <- NetworkMatrix[
      intersect(rownames(NetworkMatrix), colnames(X)),
      intersect(colnames(NetworkMatrix), colnames(X))]
    
    noncommongenes <- setdiff(colnames(X), colnames(NetworkMatrix))
    zero_matrix <- matrix(0, length(noncommongenes), length(noncommongenes))
    rownames(zero_matrix) <- noncommongenes
    colnames(zero_matrix) <- noncommongenes
    
    A <- magic::adiag(NetworkMatrix, zero_matrix)
    
    L <- signed_laplacian(
      choice = choice, X = X, Y = Y, delta = delta,
      A = A, model = model, dist = dist)
    
    if (isTRUE(verbose)) message("Process completed.")
    
    return(list(
      disease_genes = DiseaseGenes,
      L             = L))
  }

# ============================================================================
# Supporting functions
# ============================================================================

#' Disease-Specific Gene Repository from HumanBase
#'
#' @description
#' Download disease-associated gene predictions from the
#' \href{https://hb.flatironinstitute.org/}{HumanBase} resource.
#' The function retrieves gene-level association scores for a given Disease Ontology ID (DOID)
#' and returns a tidy data frame containing gene identifiers and scores.
#'
#' @param doid Character string specifying Disease Ontology ID (\code{"DOID:XXXX"}).
#' @param cache Logical value; if \code{TRUE}, downloaded HumanBase files are cached
#'              for reuse in \code{cache_dir}. If \code{FALSE} (default), files are downloaded
#'              for the current session only.
#' @param cache_dir Character string specifying a directory used to cache
#'                  downloaded HumanBase files (when \code{cache = TRUE}).
#' @param verbose Logical value, if \code{TRUE} progress messages are printed.
#' 
#' @note An active internet connection is required.
#'
#' @return A data frame with three columns:
#' \itemize{
#'   \item \code{entrez_id}: Entrez gene identifier.
#'   \item \code{standard_name}: Gene symbol.
#'   \item \code{score}: Association score from HumanBase.
#' }
#'
#' @examples
#' \donttest{
#' 
#'    # - Download disease-specific gene repository for Lung Adenocarcinoma -
#'
#'       disease_genes <- RepositoryDisease(
#'        doid      = "DOID:1324",
#'        cache     = FALSE,
#'        cache_dir = NULL,
#'        verbose   = FALSE
#'       )$standard_name
#'
#'       head(disease_genes)
#' }
#'
#' @name RepositoryDisease
#' 
#' @export

  RepositoryDisease <- function(
      doid      = NULL,
      cache     = FALSE,
      cache_dir = NULL,
      verbose   = FALSE
  ) {
    
    if (is.null(doid)) stop("'doid' must be provided when 'disease_file = NULL'.")
    
    doid_numeric <- sub("DOID:", "", doid)
    
    cache_file <- file.path(
      cache_dir,
      paste0("repositoryDisease_DOID", doid_numeric, ".tsv")
    )
    
    if (cache && file.exists(cache_file)) {
      return(utils::read.table(
        cache_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      ))
    }
    
    url <- paste0("https://humanbase.net/api/terms/disease-ontology-doid", doid_numeric, "/predictions/")
    
    response <- httr::GET(url, httr::timeout(60))
    
    httr::stop_for_status(response)
    
    predictions <- httr::content(response, as = "parsed", encoding = "UTF-8")
    if (length(predictions) == 0) stop("No predictions returned for ", doid)
    
    df <- dplyr::bind_rows(lapply(predictions, function(item) {
      data.frame(
        entrez_id        = item$gene$entrez,
        standard_name    = item$gene$standard_name,
        score            = item$score,
        stringsAsFactors = FALSE)}))
    
    # - Save the table [entrezID, gene name, score] -
    
    if (cache) {
      utils::write.table(
        df,
        file = cache_file,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
      )
      if (isTRUE(verbose)) message("Repository file saved to: ", cache_file)
    }
    
    return(df)
  }

# ============================================================================

#' Tissue-Specific Top Edge Network from HumanBase
#' 
#' @description
#' Downloads the top edge gene interaction network for a specific human tissue
#' from the \href{https://hb.flatironinstitute.org/download/tissue-networks/giant}{HumanBase} resource.
#'
#'@param tissue Character string specifying the name of the tissue to download.
#'              Spaces will automatically be converted to underscores.
#' @param cache Logical value; if \code{TRUE}, downloaded HumanBase files are cached
#'              for reuse in \code{cache_dir}. If \code{FALSE} (default), files are downloaded
#'              for the current session only.
#' @param cache_dir Character string specifying a directory used to cache
#'                  downloaded HumanBase files (when \code{cache = TRUE}).
#' @param verbose Logical value, if \code{TRUE} progress messages are printed.
#' 
#' @note An active internet connection is required.
#'
#' @return A data.frame with tissue-specific gene interactions (columns:
#'                    \code{gene1}, \code{gene2}, and \code{score}).
#'                    
#' @examples
#' \donttest{
#' 
#'    # - Download tissue-specific repository for Lung Adenocarcinoma -
#'
#'       tissue <- RepositoryTissue(
#'        tissue    = "lung",
#'        cache     = FALSE,
#'        cache_dir = NULL,
#'        verbose   = FALSE
#'       )
#'
#'       head(tissue)
#' }
#'
#' @name RepositoryTissue
#' 
#' @export

  RepositoryTissue <- function(
      tissue    = NULL, 
      cache     = FALSE, 
      cache_dir = NULL, 
      verbose   = FALSE
  ){
    
    if (is.null(tissue)) stop("'tissue' must be provided when 'tissue_file = NULL'.")
    
    tissue1 <- gsub(" ", "_", tissue)
    
    cache_file <- file.path(cache_dir, paste0("tissueSpecificEdge_", tissue1, ".tsv"))
    url <- paste0("https://s3-us-west-2.amazonaws.com/humanbase/networks/", tissue1, "_top.gz")
    
    if (cache && file.exists(cache_file)) {
      return(utils::read.table(
        cache_file,
        header = FALSE,
        sep = "\t",
        stringsAsFactors = FALSE
      ))
    }
    
    if (isTRUE(verbose)){
      message("Downloading ", tissue, "-specific gene interaction network...")
      message("(this may take a few minutes)")
    }
    
    temp_file <- tempfile(fileext = ".gz")
    
    h <- curl::new_handle()
    curl::handle_setopt(h, timeout = 1800)
    
    tryCatch({
      curl::curl_download(url, destfile = temp_file, quiet = !isTRUE(verbose), handle = h)
    }, error = function(e) {
      stop("Error downloading file: ", conditionMessage(e))
    })
    
    tissueSpecificEdge <- utils::read.table(gzfile(temp_file), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    if (cache) {
      utils::write.table(
        tissueSpecificEdge,
        file = cache_file,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE
      )
      
      if (isTRUE(verbose)) message("Tissue-specific network saved to: ", cache_file)
    }
    
    return(tissueSpecificEdge)
  }

# ============================================================================

#' Compute Signed Laplacian Matrix for a Gene Regulatory Network
#'
#' @keywords internal
#' @noRd

  signed_laplacian <- function(
      choice, X, Y,
      delta, A,
      model, dist
  ) {
    
    if (choice == 1){
  
      S <- suppressWarnings(sign(stats::cor(X, use = "pairwise.complete.obs")))
      S[is.na(S)] <- 0
      W <- S * A
      
    } else if (choice == 2) {
      
      if (is.null(Y) || is.null(delta)) stop("For 'choice = 2', Y and delta must be provided.")
      if (is.null(model)) stop("Parameter 'model' must be provided when 'choice = 2'.")
      
      fit_ridge  <- switch(model,
                           
                           "COXNet" = glmnet::cv.glmnet(X, survival::Surv(Y, delta), family = "cox", alpha = 0),
                           
                           "AFTNet" = survival::survreg(survival::Surv(exp(Y), delta) ~ survival::ridge(X, theta = 1),
                                                        dist = dist, scale = 0))
      beta_ridge <- switch(model,
                           
                           "COXNet" = as.numeric(stats::coef(fit_ridge, s = fit_ridge$lambda.min)),
                           
                           "AFTNet" = fit_ridge$coefficients[-1])
      
      beta_ridge[is.na(beta_ridge)] <- 0
      s_vec <- sign(beta_ridge)

      W <- diag(s_vec) %*% A %*% diag(s_vec)
    }
    
    # - Laplacian matrix: L = D - W
    
    D <- diag(colSums(abs(W)))
    L <- D - W
    
    return(L)
  }