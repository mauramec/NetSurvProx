# ============================================================================
# Main Functions for Pathway Analysis and Enrichment
# ============================================================================

#' Interactive Pathway Analysis Dashboard
#'
#' @description
#' Constructs interactive pathway analysis networks and generates an HTML
#' dashboard from a list of genes. Pathways can be retrieved via 
#' \href{https://www.genome.jp/kegg/}{KEGG} database or provided through a custom file.
#' 
#' @param genes_list Character vector of gene symbols, a file path to a tab-delimited file,
#'                   or a data frame where the first column contains gene symbols.
#' @param header Logical value indicating whether the input file has a header (\code{default: TRUE}).
#' @param useKeggAPI Logical value indicating whether to use the KEGG REST API
#'                   to retrieve pathways (\code{default: TRUE}).
#' @param pathway_file Optional data frame or file path containing custom pathway data.
#'                    Required if \code{useKeggAPI = FALSE}.
#'                    Must have columns: \code{pathway}, \code{gene}, optional \code{name}.
#' @param nodesCols Character vector of length 2 defining node colors.
#'                  First color for regular nodes, second for highlighted nodes (when \code{diseaseNodes = TRUE}).
#' @param diseaseNodes Logical value indicating whether to highlight
#'                     disease-associated nodes (\code{default: TRUE}).
#' @param disease_file Optional file path or data frame containing disease-associated gene scores.
#'                     Must have at least two columns: gene and score.
#' @param top_percent Numeric value indicating the percentage of top genes
#'                    to highlight based on \code{disease_file} (used with \code{diseaseNodes}, \code{default: 20}).
#' @param batch_size Numeric value indicating the batch size for KEGG API queries (\code{default: 10}).
#' @param background_genes Optional vector of background genes for enrichment analysis.
#' @param min_genes Numeric value indicating minimum number of genes in a pathway
#'                  to be considered (\code{default: 2}).
#' @param top_n Numeric value indicating the number of top pathways to display
#'              in the dashboard (\code{default: 10}).
#' @param db_name Character string specifying the Bioconductor Annotation DB name for gene
#'                mapping (\code{default: "org.Hs.eg.db"}).
#' @param organism Character string specifying KEGG organism code (\code{default: "hsa"}).
#' @param out_dir Character string specifying output directory for results.
#' @param open_browser Logical value; if \code{TRUE} and interactive session,
#'                      opens dashboard in browser (\code{default: TRUE}).
#' @param verbose Logical value, if \code{TRUE} progress messages are printed.
#' 
#' @details
#' Workflow implemented by the function:
#' \enumerate{
#'    \item Converts gene symbols to Entrez IDs for KEGG queries and maps back to gene symbols after pathway retrieval.
#'    \item Retrieves pathways using KEGG API if \code{useKeggAPI = TRUE}, otherwise uses \code{pathway_file}.
#'    \item Constructs a gene-pathway binary incidence matrix (genes as rows, pathways as columns).
#'    \item Builds an \code{igraph} network where genes are nodes and edges link genes in the same pathways.
#'    \item Assigns node colors based on connectivity and optional disease association.
#'    \item Highlights top genes by connectivity or disease association using \code{nodesCols} and \code{top_percent}.
#'    \item Saves network information in \code{network_data.rds} and optionally renders an interactive HTML dashboard
#'          (\code{Dashboard.html}).
#' }
#' 
#' The \code{network_data.rds} object contains:
#' \itemize{
#'    \item \code{g}: igraph object representing the network.
#'    \item \code{edge_info}: data frame with edges, colors, and pathway labels.
#'    \item \code{legend_info}: legend codes, colors, and counts for pathways.
#'    \item \code{all_genes}, \code{conn_genes}: all input genes and connected genes.
#'    \item \code{node_colours}: node colors and borders for plotting.
#'    \item \code{pathway_df}: data frame of pathways and genes.
#'    \item \code{background}, \code{min_genes}, \code{top_n}: parameters.
#' }
#' 
#' @note If \code{useKeggAPI = TRUE}, the function queries the KEGG REST API to
#' retrieve pathway information. An active internet connection is required in this case.
#' Moreover, gene names conversion relies on local Bioconductor Annotation DBs (e.g., org.Hs.eg.db).
#' The function returns paths to generated files but does not print to console
#' or open files unless explicitly requested.
#'
#' @return Saves:
#' \itemize{
#'    \item \code{network_data.rds}: serialized network object for later use.
#'    \item \code{Dashboard.html}: interactive dashboard showing network and enrichment panels.
#' }
#'
#' @seealso
#' \code{\link{Enrichment}} for pathway enrichment results.
#'
#' @name PathwayDashboard
#'
#' @export

  PathwayDashboard <- function(
      genes_list,
      header           = TRUE,
      useKeggAPI       = TRUE,
      pathway_file     = NULL,
      nodesCols        = c("#5C7997", "#F5C59F"),
      diseaseNodes     = FALSE,
      disease_file     = NULL,
      top_percent      = 20,
      batch_size       = 10,
      background_genes = NULL,
      min_genes        = 2,
      top_n            = 10,
      db_name          = "org.Hs.eg.db",
      organism         = "hsa",
      out_dir          = NULL,
      open_browser     = TRUE,
      verbose          = FALSE
  ) {
    
    # --- Input checks ---
    
    required_pkgs <- c("igraph", "rmarkdown", "AnnotationDbi", db_name)
    
    for (pkg in unique(required_pkgs)) {
      if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.", call. = FALSE)
    }
    
    if (is.null(out_dir)) {
      out_dir <- tempdir()
      message("Output directory not specified. Files will be saved in: ", out_dir)
    } else if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    if (is.data.frame(genes_list)) {
      genes_df <- genes_list
    } else if (is.character(genes_list)) {
      if (length(genes_list) == 1 && file.exists(genes_list)) {
        genes_df <- utils::read.table(genes_list, header = header, sep = "\t", stringsAsFactors = FALSE)
      } else {
        genes_df <- data.frame(Gene = genes_list, stringsAsFactors = FALSE)
      }
    } else {
      stop("The 'genes_list' must be a file path, a data.frame, or a character vector.")
    }
    
    gene_symbols <- unique(genes_df[[1]])
    
    if (length(gene_symbols) == 0) stop("No valid genes found.")
    
    if (diseaseNodes && !is.null(disease_file)) {
      if (is.data.frame(disease_file)) {
        disease_df <- disease_file
      } else if (is.character(disease_file)) {
        if (length(disease_file) == 1 && file.exists(disease_file)) {
          disease_df <- utils::read.table(disease_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        } else {
          disease_df <- data.frame(Gene = disease_file, stringsAsFactors = FALSE)
        }
      } else {
        stop("The 'disease_file' must be a file path, a data.frame, or a character vector.")
      }
      disease_df <- unique(as.character(disease_df[[1]]))
    }
    
    # -----  
    
    if (useKeggAPI) {
      
      if (isTRUE(verbose)) message("Loading pathways from KEGG API...")
      
      if (!requireNamespace("httr", quietly = TRUE)) {
        stop("Package 'httr' required for KEGG mode.", call. = FALSE)
      }
      
      if (!requireNamespace("curl", quietly = TRUE) || !curl::has_internet()) {
        stop("Internet connection required for KEGG mode.", call. = FALSE)
      }
      
      entrez_ids <- symbol_to_entrez(gene_symbols, db_name = db_name)
      entrez_ids <- entrez_ids[!is.na(entrez_ids)]
      
      if (length(entrez_ids) == 0) stop("No valid Entrez IDs found after local mapping.", call. = FALSE)
      
      pathway_df <- get_KEGGpathways(entrez_ids, organism = organism, batch_size = batch_size)
      
      if (nrow(pathway_df) == 0) stop("No KEGG pathways retrieved.", call. = FALSE)
      
      symbol_map      <- entrez_to_symbol(unique(pathway_df$entrez_id), db_name = db_name)
      pathway_df$gene <- symbol_map[as.character(pathway_df$entrez_id)]
      pathway_df      <- pathway_df[!is.na(pathway_df$gene), ]
      
      if (nrow(pathway_df) == 0) stop("No valid gene symbols after mapping back from Entrez IDs.", call. = FALSE)
      
      r     <- httr::GET(paste0("https://rest.kegg.jp/list/pathway/", organism))
      txt   <- httr::content(r, "text")
      lines <- strsplit(txt, "\n")[[1]]
      
      df <- do.call(rbind, lapply(lines, function(l) {
        parts <- strsplit(l, "\t")[[1]]
        if (length(parts) < 2) return(NULL)
        data.frame(
          pathway = sub(paste0("^path:", organism), "", parts[1]),
          name    = parts[2],
          stringsAsFactors = FALSE)}))
      
      pathway_df <- merge(pathway_df, df, by = "pathway", all.x = TRUE)
      
    } else {
      
      if (isTRUE(verbose)) message("Loading custom pathway file...")
      
      if (is.null(pathway_file)) stop("'pathway_file' required when 'useKeggAPI' = FALSE.")
      
      if (is.character(pathway_file) && file.exists(pathway_file)) {
        pathway_df <- utils::read.table(pathway_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      } else if (is.data.frame(pathway_file)) {
        pathway_df <- pathway_file
      } else {
        stop("The 'pathway_file' must be either a path to a valid file or a data.frame object.")}
    }
    
    gp_matrix          <- build_gp_matrix(gene_symbols, pathway_df)
    gp_matrix_filtered <- gp_matrix[ , colSums(gp_matrix) >= 2, drop = FALSE]
    conn_genes         <- rownames(gp_matrix_filtered)[rowSums(gp_matrix_filtered) > 0]
    
    if (isTRUE(verbose)) message("Connected genes: ", length(conn_genes))
    
    net_obj <- create_network(conn_genes, gp_matrix)
    g       <- net_obj$network
    
    if (!is.null(net_obj$edge_info) && igraph::ecount(g) > 0) igraph::E(g)$color <- net_obj$edge_info$color
    
    missing_genes <- setdiff(gene_symbols, igraph::V(g)$name)
    
    if(length(missing_genes) > 0) {
      g <- igraph::add_vertices(g, nv = length(missing_genes), name = missing_genes)
    }
    
    node_col_info <- define_node_colour(
      all_genes    = gene_symbols,
      conn_genes   = conn_genes,
      diseaseNodes = diseaseNodes,
      disease_file = disease_file,
      genesCol     = nodesCols,
      top_percent  = top_percent)
    
    igraph::V(g)$color       <- node_col_info$n_colours[match(igraph::V(g)$name, gene_symbols)]
    igraph::V(g)$frame.color <- node_col_info$n_borders[match(igraph::V(g)$name, gene_symbols)]
    
    if (length(igraph::E(g)) > 0) igraph::E(g)$label <- seq_along(igraph::E(g))
    
    pathway_labels <- unique(pathway_df[ , c("pathway", "name")])
    colnames(pathway_labels) <- c("code", "label")
    
    matched_labels <- merge(
      data.frame(code = net_obj$legend_info$legend_code,
                 stringsAsFactors = FALSE),
      pathway_labels, by = "code", all.x = TRUE)
    matched_labels <- matched_labels[!is.na(matched_labels$label), ]
    matched_labels <- matched_labels[match(net_obj$legend_info$legend_code, matched_labels$code), ]
    n_paths <- nrow(matched_labels)
    
    raw_cols <- RColorBrewer::brewer.pal(8, "Accent")
    
    palette_colors        <- grDevices::colorRampPalette(raw_cols)(n_paths)
    names(palette_colors) <- matched_labels$code
    
    net_obj$legend_info$legend_code   <- matched_labels$code
    net_obj$legend_info$legend_label  <- matched_labels$label
    net_obj$legend_info$legend_colour <- palette_colors[matched_labels$code]
    net_obj$legend_info$legend_count  <- seq_len(nrow(matched_labels))
    
    if (!is.null(net_obj$edge_info)) net_obj$edge_info$color <- palette_colors[net_obj$edge_info$path]
    
    rds_path <- file.path(out_dir, "network_data.rds")
    
    saveRDS(list(
      g            = g,
      edge_info    = net_obj$edge_info,
      legend_info  = net_obj$legend_info,
      all_genes    = gene_symbols,
      conn_genes   = conn_genes,
      node_colours = node_col_info,
      pathway_df   = pathway_df,
      background   = background_genes,
      min_genes    = min_genes,
      top_n        = top_n),
      file         = rds_path)
    
    if (!is.null(out_dir)) if (isTRUE(verbose)) message("Saving network informations as 'network_data.rds'...")
    
    dashboard_path <- system.file("Dashboard.Rmd", package = "NetSurvProx")
    dashboard_output <- file.path(out_dir, "Dashboard.html")
    
    rmarkdown::render(
      input = dashboard_path,
      output_file = dashboard_output,
      quiet = TRUE,
      params = list(
        inputFile    = genes_list,
        header       = header,
        useKeggAPI   = useKeggAPI,
        nodesCols    = nodesCols,
        diseaseNodes = diseaseNodes,
        disease_file = disease_file,
        top_percent  = top_percent,
        batch_size   = batch_size,
        dataset      = db_name,
        organism     = organism,
        background   = background_genes,
        min_genes    = min_genes,
        top_n        = top_n,
        rds_path     = rds_path
      ),
      envir = new.env()
    )
    
    if (open_browser && interactive()) utils::browseURL(dashboard_output)
    
    if (!is.null(out_dir)) message("Interactive Pathway Analysis saved as 'Dashboard.html.'")
  }

# ============================================================================

#' Pathway Enrichment (Over-representation Analysis)
#'
#' @description
#' Performs pathway enrichment analysis to evaluate whether a set of
#' genes is over-represented in one or more pathways compared to a background set of genes.
#' For each pathway, it calculates the number of observed genes, the Fisher's exact test
#' p-value, and FDR-adjusted p-values. Significant pathways (\code{padj < 0.05})
#' are marked with \code{Yes} in the \code{highlight} column.
#'
#' @param genes Character vector specifying the list of selected gene symbols.
#' @param pathway_df Data frame with at least the following columns:
#'   \itemize{
#'     \item \code{pathway}: pathway identifier.
#'     \item \code{gene}: gene symbol belonging to the pathway.
#'     \item \code{name}: optional descriptive name for the pathway.
#'   }
#' @param background_genes Character vector specifying background gene set.
#'                         If \code{NULL} (default), the union of \code{genes}
#'                         and all genes in \code{pathway_df} is used.
#' @param min_genes Numeric value specifying the minimum number of background
#'                  genes that a pathway must have to be considered (\code{default: 2}).
#' @param top_n Numeric value specifying the number of top pathways sorted by
#'              adjusted p-value to return (\code{default: 10}).
#' @param out_file Character string specifying the path to save the enrichment
#'                 results as an Excel file (.xlsx). If \code{NULL} (default),
#'                 the results are not written to disk.
#' @details
#' The function implements an over-representation analysis (ORA) workflow:
#' \enumerate{
#'   \item Intersects the input gene list with a background set (user-provided or derived from all pathway genes).
#'   \item Filters pathways to retain only those with at least \code{min_genes} present in the background.
#'   \item Performs Fisher's exact test for each pathway to assess over-representation.
#'   \item Adjusts p-values using the false discovery rate (FDR) method.
#'   \item Identifies significantly enriched pathways (\code{padj < 0.05}) and marks them in the \code{highlight} column.
#'   \item Selects the top \code{top_n} pathways for visualization in dashboards or plots.
#' }
#'
#' The results are automatically saved as an Excel file \code{Enrichment_results.xlsx} and are used by
#' \code{\link{PathwayDashboard}} to display enrichment results interactively
#' in the dedicated panel.
#'
#' @return A list containing:
#' \itemize{
#'    \item \code{results}: Full enrichment table with p-values and FDR correction,
#'                          including \code{pathway}, \code{nGenes}
#'                          (number of genes for pathway), \code{pval},
#'                          \code{padj}, \code{highlight}
#'                          (\code{Yes/No} if the pathway is enriched), \code{name}.
#'     \item \code{bar_data}: Top \code{top_n} enriched pathways.
#' }
#'
#' @seealso
#' \code{\link{PathwayDashboard}} for interactive visualization of enrichment results.
#'
#' @name Enrichment
#'
#' @export

  Enrichment <- function(
      genes,
      pathway_df,
      background_genes = NULL,
      min_genes        = 2,
      top_n            = 10,
      out_file         = NULL
  ) {
    
    if (!is.null(out_file)) {
      if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Package 'openxlsx' is required to save results to disk.", call. = FALSE)
    }
    
    if(is.null(background_genes)) background_genes <- unique(c(genes, pathway_df$gene))
    
    gene_list <- toupper(trimws(genes))
    
    pathway_list <- split(pathway_df$gene, pathway_df$pathway)
    pathway_list <- lapply(pathway_list, function(x) intersect(x, background_genes))
    pathway_list <- pathway_list[sapply(pathway_list, length) >= min_genes]
    
    if(length(pathway_list)==0) stop("No pathways meet the minimum background gene requirement.")
    
    n_gene_list <- length(gene_list)
    bg_set      <- background_genes
    
    res <- data.frame(
      pathway = names(pathway_list),
      nGenes  = vapply(pathway_list, function(pw) sum(pw %in% gene_list), FUN.VALUE = integer(1)),
      pval    = NA_real_)
    
    res$pval <- vapply(seq_along(pathway_list), function(i){
      pw_genes <- pathway_list[[i]]
      a <- res$nGenes[i]
      b <- n_gene_list - a
      c <- length(pw_genes) - a
      d <- max(0, length(bg_set) - (a + b + c))
      stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
    }, FUN.VALUE = numeric(1))
    
    res$padj      <- stats::p.adjust(res$pval, method = "fdr")
    res$highlight <- ifelse(res$padj < 0.05, "Yes", "No")
    res           <- res[order(res$padj),]
    
    pw_labels <- unique(pathway_df[, c("pathway", "name")])
    res$name  <- pw_labels$name[match(res$pathway, pw_labels$pathway)]
    
    if (!is.null(out_file)) {
      
      openxlsx::write.xlsx(res, file = out_file, rowNames = FALSE)
      message("Output saved as 'Enrichment_results.xlsx'.")
    }  
    
    return(list(
      results  = res,
      bar_data = utils::head(res, top_n)))
  }

# ============================================================================
# Supporting functions
# ============================================================================

#' Convert Gene Symbols to Entrez IDs using local AnnotationDbi
#'
#' @keywords internal
#' @noRd

  symbol_to_entrez <- function(
      symbols,
      db_name
  ) {
    
    db_obj <- get(db_name, envir = asNamespace(db_name))
    
    symbols <- as.character(symbols[!is.na(symbols)])
    
    res <- tryCatch({
      suppressMessages({
        AnnotationDbi::select(
          db_obj, 
          keys = symbols, 
          columns = "ENTREZID", 
          keytype = "SYMBOL"
        )
      })
    }, error = function(e) {
      return(NULL)
    })
    
    result <- stats::setNames(rep(NA_integer_, length(symbols)), symbols)
    
    if (!is.null(res) && nrow(res) > 0) {
      res <- res[!is.na(res$ENTREZID), ]
      res <- res[!duplicated(res$SYMBOL), ]
      
      common <- intersect(res$SYMBOL, symbols)
      result[common] <- as.integer(res$ENTREZID[match(common, res$SYMBOL)])
    }
    
    return(result)
  }

# ============================================================================

#' Convert Entrez IDs to Gene Symbols using local AnnotationDbi
#' 
#' @keywords internal
#' @noRd

  entrez_to_symbol <- function(
      entrez_ids,
      db_name
  ) {
    
    db_obj <- get(db_name, envir = asNamespace(db_name))
    
    entrez_ids <- as.character(entrez_ids[!is.na(entrez_ids)])
    
    res <- tryCatch({
      suppressMessages({
        AnnotationDbi::select(
          db_obj, 
          keys = entrez_ids, 
          columns = "SYMBOL", 
          keytype = "ENTREZID"
        )
      })
    }, error = function(e) {
      return(NULL)
    })
    
    result <- stats::setNames(rep(NA_character_, length(entrez_ids)), entrez_ids)
    
    if (!is.null(res) && nrow(res) > 0) {
      res <- res[!is.na(res$SYMBOL), ]
      res <- res[!duplicated(res$ENTREZID), ]
      
      common <- intersect(res$ENTREZID, entrez_ids)
      result[common] <- res$SYMBOL[match(common, res$ENTREZID)]
    }
    
    return(result)
  }

# ============================================================================    

#' Queries KEGG REST API to Get Pathways Associated with a Set of Entrez IDs
#'
#' @keywords internal
#' @noRd
  
  get_KEGGpathways <- function(
    entrez_ids,
    organism,
    batch_size
  ) {
    
    query_kegg <- function(ids_batch) {
      query <- paste0(organism, ":", ids_batch, collapse = "+")
      url   <- paste0("https://rest.kegg.jp/link/pathway/", query)
      
      res <- tryCatch(httr::GET(url, httr::timeout(10)), error = function(e) NULL)
      if (is.null(res) || httr::status_code(res) != 200) return(NULL)
      
      txt <- httr::content(res, "text")
      if (txt == "") return(NULL)
      
      lines <- strsplit(txt, "\n")[[1]]
      lines <- lines[lines != ""]
      
      df_list <- lapply(lines, function(l) {
        parts <- strsplit(l, "\t")[[1]]
        if (length(parts) == 2) {
          data.frame(
            entrez_id = sub(paste0("^", organism, ":"), "", parts[1]),
            pathway   = sub("^path:", "", parts[2]),
            stringsAsFactors = FALSE)
        } else NULL})
      
      df_list <- Filter(Negate(is.null), df_list)
      if (length(df_list) == 0) return(NULL)
      
      do.call(rbind, df_list)
    }
    
    batches <- split(entrez_ids, ceiling(seq_along(entrez_ids) / batch_size))
    
    results <- lapply(batches, query_kegg)
    results <- Filter(Negate(is.null), results)
    
    if (length(results) == 0) {
      return(data.frame(entrez_id = character(0), pathway = character(0),
                        stringsAsFactors = FALSE))}
    do.call(rbind, results)
  }

# ============================================================================

#' Create a Binary Incidence Matrix Indicating which Genes Belong to which Pathways
#'
#' @keywords internal
#' @noRd

  build_gp_matrix <- function(
      genes,
      pathway_df
  ) {
    
    df <- pathway_df[pathway_df$gene %in% genes, ]
    
    pathways <- unique(df$pathway)
    
    M <- matrix(0, nrow = length(genes), ncol = length(pathways),
                dimnames = list(genes, pathways))
    
    if(nrow(df) > 0) {
      row_idx <- match(df$gene, genes)
      col_idx <- match(df$pathway, pathways)
      M[cbind(row_idx, col_idx)] <- 1}
    
    return(M)
  }

# ============================================================================

#' Create Gene Network Based on Shared Pathways
#'
#' @keywords internal
#' @noRd

  create_network <- function(
      conn_genes,
      gp_matrix
  ) {
    
    if (length(conn_genes) == 0) {
      return(list(
        network     = igraph::make_empty_graph(),
        edge_info   = NULL,
        legend_info = list(legend_code   = character(0),
                           legend_colour = character(0))))}
    
    conn_genes <- intersect(conn_genes, rownames(gp_matrix))
    
    if (length(conn_genes) == 0) {
      return(list(
        network     = igraph::make_empty_graph(),
        edge_info   = NULL,
        legend_info = list(legend_code   = character(0),
                           legend_colour = character(0))))}
    
    path_filt <- gp_matrix[conn_genes, , drop = FALSE]
    path_filt <- path_filt[, colSums(path_filt) > 1, drop = FALSE]
    
    if (ncol(path_filt) == 0) {
      return(list(
        network     = igraph::make_empty_graph(),
        edge_info   = NULL,
        legend_info = list(legend_code   = character(0),
                           legend_colour = character(0))))}
    
    edges <- data.frame(from    = character(),
                        to      = character(),
                        pathway = character(),
                        stringsAsFactors = FALSE)
    
    for (p in colnames(path_filt)) {
      genes_in_p <- rownames(path_filt)[path_filt[ , p] == 1]
      if (length(genes_in_p) > 1) {
        comb  <- t(utils::combn(genes_in_p, 2))
        edges <- rbind(edges, data.frame(
          from = comb[ , 1],
          to   = comb[ , 2],
          pathway = p,
          stringsAsFactors = FALSE))}}
    
    if (nrow(edges) == 0) {
      return(list(
        network     = igraph::make_empty_graph(),
        edge_info   = NULL,
        legend_info = list(legend_code   = character(0),
                           legend_colour = character(0))))}
    
    all_nodes <- unique(c(edges$from, edges$to, conn_genes))
    g <- igraph::graph_from_data_frame(edges[, 1:2], directed = FALSE, vertices = all_nodes)
    unique_paths <- unique(edges$pathway)
    
    colors <- grDevices::rainbow(length(unique_paths))
    names(colors) <- unique_paths
    
    edge_info <- data.frame(
      from  = edges$from,
      to    = edges$to,
      label = seq_len(nrow(edges)),
      color = colors[edges$pathway],
      path  = edges$pathway,
      stringsAsFactors = FALSE)
    
    legend_info <- list(
      legend_code   = unique_paths,
      legend_colour = colors[unique_paths])
    
    list(
      network = g,
      edge_info = edge_info,
      legend_info = legend_info)
  }

# ============================================================================

#' Define Node Colors for Network Visualization (with Optional Disease Scores)
#'
#' @keywords internal
#' @noRd

  define_node_colour <- function(
      all_genes, conn_genes,
      diseaseNodes, disease_file,
      genesCol, top_percent
  ) {
    
    n <- length(all_genes)
    n_colours <- rep(NA, n)
    n_borders <- rep(NA, n)
    mapped_idx <- match(conn_genes, all_genes)
    
    if (length(mapped_idx) > 0) {
      
      valid_idx <- which(!is.na(mapped_idx))
      
      if (diseaseNodes && !is.null(disease_file)) {
        
        if (is.data.frame(disease_file)) {
          gname <- as.character(disease_file[[1]])
        } else {
          gname <- as.character(disease_file)
        }
        
        idx <- match(conn_genes[valid_idx], gname)
        
        top_limit <- ceiling(length(gname) * top_percent / 100)
        
        top_idx    <- which(!is.na(idx) & idx <= top_limit)
        bottom_idx <- which(!is.na(idx) & idx > top_limit)
        
        if(length(top_idx) > 0) {
          n_colours[mapped_idx[valid_idx[top_idx]]] <- genesCol[2]
          n_borders[mapped_idx[valid_idx[top_idx]]] <- genesCol[2]
        }
        
        if(length(bottom_idx) > 0) {
          n_colours[mapped_idx[valid_idx[bottom_idx]]] <- genesCol[1]
          n_borders[mapped_idx[valid_idx[bottom_idx]]] <- genesCol[1]
        }
        
        missing_in_disease <- valid_idx[is.na(idx)]
        if (length(missing_in_disease) > 0) {
          n_colours[mapped_idx[missing_in_disease]] <- genesCol[1]
          n_borders[mapped_idx[missing_in_disease]] <- genesCol[1]
        }
      } else {
        n_colours[mapped_idx[valid_idx]] <- genesCol[1]
        n_borders[mapped_idx[valid_idx]] <- genesCol[1]
      }
    }
    
    n_colours[is.na(n_colours)] <- "#FFFFFF"
    n_borders[is.na(n_borders)] <- "#e6e4df"
    
    return(list(
      n_colours = n_colours,
      n_borders = n_borders
    ))
  }