#' Example Dataset for Network-Based Survival Analysis
#'
#' A pre-processed dataset containing clinical survival information and
#' gene expression covariates for Lung Adenocarcinoma (TCGA-LUAD).
#' This dataset allows users to bypass the computationally intensive
#' download and preprocessing pipeline, providing immediate access to
#' the covariate matrix, survival outcomes, and censoring indicators.
#'
#' Gene expression data (RNA-seq) were obtained from the LinkedOmics portal
#' and processed to construct:
#' \itemize{
#'   \item screened gene expression matrix \code{X} (samples × genes),
#'   \item observed survival times \code{Y} (real scale),
#'   \item censoring indicators \code{delta} (1 = event, 0 = censored).
#' }
#' 
#' The screening was performed using the BMD method (see \code{\link{VariableScreening}}) 
#' focusing on disease-associated genes retrieved for \code{doid = "DOID:1324"} 
#' via \code{\link{RepositoryDisease}}.
#' 
#' The dataset is pre-partitioned into an 70% training set for model 
#' estimation and a 30% testing set for validation.
#'
#' @format A list with the following components.
#' \itemize{
#'   \item \code{X_train} : numeric matrix of training covariates.
#'   \item \code{X_test} : numeric matrix of testing covariates.
#'   \item \code{Y_train} : numeric vector of observed training survival times.
#'   \item \code{Y_test} : numeric vector of observed testing survival times. 
#'   \item \code{delta_train} : integer vector of training censoring indicators.
#'   \item \code{delta_test} : integer vector of testing censoring indicators.
#' }
#'
#' @source \url{https://linkedomics.org/data_download/TCGA-LUAD/}
#'
#' @usage data(LUADdataset)
"LUADdataset"
