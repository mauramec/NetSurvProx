Introduction to NetSurvProx Package
================

This package is developed by **Maura Mecchi** and **Antonella Iuliano** (University of Basilicata, Italy).

-----------------------------------------------------------------

`NetSurvProx` introduces a novel network-constrained survival analysis
environment designed to address the **challenges of high-dimensional
biomolecular data**, such as:

- limited sample size compared to the large number of variables,

- high multicollinearity among predictors,

- the need to incorporate a priori biological knowledge.

By integrating prior information from curated interaction networks like
KEGG into a double-penalty framework, the method facilitates variable
selection and estimation for both Cox Proportional Hazard (PH) and
Accelerated Failure Time (AFT) models.

The package utilizes proximal optimization algorithms and
cross-validation to ensure reliable parameter tuning and the
identification of predictive biomarkers from censored data. Ultimately,
the framework enhances interpretability through dedicated utility
functions, providing biologically coherent insights that **support
targeted therapies and the advancement of personalized medicine**.

<div style="margin-top: 1.5em;">

</div>

# <span style="color: #444444;">Package Structure</span>

In predictive modeling, the integrity of the evaluation process depends
on the rigorous separation between model development and performance
assessment. Data are typically partitioned into a **training set** and a
**testing set**, an essential split to prevent overfitting and ensure
that predictive accuracy serves as a reliable proxy for true
generalizability. To maintain this independence, all pre-processing
operations are calibrated using the training set and subsequently
applied to the testing set, ensuring the latter remains entirely
external to the model estimation process. Aligning with these
principles, the package organizes the analytical workflow into **two
distinct functional stages**:

- the <span style="color: #f5c59f;">training phase</span> to fit the
  model,

- the <span style="color: #8abfe7;">testing phase</span> to evaluate
  performance on unseen samples.

Briefly, the process begins by constructing *gene networks* from
training data to perform *variable screening*. These selected genes are
integrated into *network-penalized survival models* to derive a
prognostic gene signature, a prognostic index, and an *optimal cut-off*.

The resulting signature is validated on the testing set using *survival
analysis*, statistical significance testing, and performance metrics. To
ensure biological relevance, the workflow concludes with *pathway and
functional enrichment analyses*, translating the statistical signature
into interpretable biological insights.

## <span style="color: #7dae9b;">Regularized Estimator</span>

The package applies a double penalty to the survival model. Given a sample of size $n$, $(y_i, \delta_i, \mathbf{x}_i)$ for $i=1,\dots,n$, the package estimates the unknown parameter vector $\boldsymbol{\xi}$, where:
* $\boldsymbol{\xi}=\boldsymbol{\beta}=(\beta_1,\dots,\beta_p)^\top$ for the Cox PH model,
* $\boldsymbol{\xi}=(\boldsymbol{\beta}^\top,\sigma)^\top$ for the AFT model, where $\sigma>0$ is a scale parameter,

using the following regularized definition:

$$\hat{\boldsymbol{\xi}} = \underset{\boldsymbol{\xi}}{\text{argmin}} ~ -\frac{1}{n} ~\ell(\boldsymbol{\xi}) + \lambda \alpha | \boldsymbol{\xi}|_1 + (1-\alpha) ~\Omega(\boldsymbol{\xi})$$

Here $\lambda>0$ is the regularization parameter, $\alpha \in [0,1]$ is the parameter that balances the contribution of the two penalties, and $\ell(\boldsymbol{\xi})$ is the log-likelihood (or partial log-likelihood) of the selected model. Both parameters are selected optimally through a data-driven procedure (CV-LP).

The first term enforces variable selection through a LASSO penalty, while the second preserves gene–gene correlations by incorporating a Laplacian-based penalty $\Omega(\boldsymbol{\xi})$, ensuring that biologically relevant network structures are maintained:

$$\Omega(\boldsymbol{\xi}) = \boldsymbol{\xi}^\top \mathbf{L} \boldsymbol{\xi} = \boldsymbol{\xi}^\top (\mathbf{D} - \mathbf{A}) \boldsymbol{\xi} = \sum_{1 \leq j < k \leq p} \lvert a_{jk} \rvert \left( \xi_j - \text{sign}(a_{jk}) \xi_k \right)^2$$

with the adjacency matrix $\mathbf{A}$, the degree matrix $\mathbf{D}$, and $\text{sign}(a_{jk})$ estimated with a Pearson correlation or ridge-penalized approach.

The `NetSurvProx` algorithms use a proximal gradient discent approach
efficiently, and the proposed survival network-penalized methods are
called **COXNet** and **AFTNet**.

<div style="margin-top: 1.5em;">

</div>

# <span style="color: #444444;">Quick Start</span>

The goal of this section is to provide users with a **comprehensive
overview of the main functions and expected outputs**. We demonstrate
the workflow using both simulated data to highlight specific features in
a controlled environment, and real-world datasets, showing how the
package handles practical scenarios.

For the sake of brevity and to showcase the behavior of both models, we
will apply the `COXNet` model to the simulated data and the `AFTNet`
model to the real one. Please note that these choices are for
illustrative purposes only; both models are highly performant and can be
effectively used in either context.

By the end of this section, users will have a better idea of how to
apply these functions to their own research.

Once installed, the package can be easly loaded into R session using the
command:

``` r
library(NetSurvProx)
```

## <span style="color: #444444;">Synthetic Data</span>

To assess the package’s performance in high-dimensional genomic
contexts, we simulate scenarios where the number of predictors
significantly exceeds the training samples. These simulations
*incorporate prior biological knowledge* via a gene regulatory network
adjacency matrix, where entries represent interactions between
transcription factors and target genes.

While we examine both independent and interconnected topological
configurations, this vignette focuses on the *interconnected regulatory
modules* involving 40 active genes. This configuration provides a more
biologically representative framework, reflecting the inherently
intertwined nature of molecular pathways and regulatory circuits in
complex diseases.

### <span style="color: #444444;">Example</span>

First, we simulate genomic features $\mathbf{X}$, survival times
$\mathbf{Y}$, and the censoring indicator $\boldsymbol{\delta}$ by
employing a **shared regulatory scheme**. Using the `COXNet` model with
a Weibull baseline, we replicate a high-dimensional scenario by
involving 165 observations and a 50% censoring rate. This setup emulates
a context where dimensionality exerts a moderate influence.

``` r
 # --- SHARED REGULATORY SCHEME FOR INTERCONNECTED MODULES ---
    
    s1 <- 5
    s2 <- 3
    
    targets <- 10
    
    shared_scheme <- list( 
      list(tfs = c(1, 3), shared = s1, unique = c(targets - s1, targets - s1)),  
      list(tfs = c(2, 4), shared = s2, unique = c(targets - s2, targets - s2)))

 # --- SIMULATE DATA USING Simulations() ---
    
    simul_data <- Simulations(
        n             = 165,
        r             = 40,
        targets       = targets,
        p_active      = 40, 
        rho           = 0.7,
        rate          = 0.5,
        b_true        = c(0.8, 1.2, -1.2, -0.8),
        nsimul        = 1,
        model         = "COXNet",
        baseline      = "weibull",
        phi           = 0.1, 
        sigma_true    = NULL,
        breaks        = NULL,
        hazards       = NULL,
        shared_scheme = shared_scheme,
        choice        = 1,
        save          = FALSE,
        save_path     = NULL,
        seed          = 2026,
        verbose       = FALSE)
   
  # --- EXTRACT COMPONENTS ---

   X     <- simul_data$X_list[[1]]
   Y     <- simul_data$time_list[[1]]
   delta <- simul_data$delta_list[[1]]
   
   L <- simul_data$L_list[[1]]
   
   beta_true <- as.vector(simul_data$beta_list[[1]])
```

To evaluate the model’s predictive performance, we split the dataset
into training (50%) and testing (50%) sets.

``` r
 # --- SEED FOR REPRODUCIBILITY ---

  set.seed(2026)

 # --- SPLIT (50/50) ---

   train_frac <- 0.5
   train_idx  <- sample(seq_len(nrow(X)), size = floor(train_frac * nrow(X)), replace = FALSE)
   test_idx   <- setdiff(seq_len(nrow(X)), train_idx)
  
   X_train <- X[train_idx, , drop = FALSE]
   Y_train <- Y[train_idx]
   delta_train <- delta[train_idx]
  
   X_test <- X[test_idx, , drop = FALSE]
   Y_test <- Y[test_idx]
   delta_test <- delta[test_idx]
```

We then fit the model using the complete `NetSurvProx` routine,
executing parallelized cross-validation to determine the optimal
regularization parameters. Following estimation, the model’s performance
is rigorously assessed on the independent testing set through a suite of
selection metrics.

``` r
 # --- FIT THE ENTIRE PROCEDURE USING NetSurvProx() ---

   fit_simulation <- NetSurvProx(
          X_train, Y_train, delta_train,
          X_test, Y_test, delta_test,
          L                 = L,
          standardize_train = TRUE,
          standardize_test  = TRUE,
          model             = "COXNet",
          dist              = NULL,
          select_lambda     = TRUE,
          alpha_grid        = c(0.3, 0.5, 0.7),
          nlambda           = 50,
          lambda_ratio      = 0.01,
          nfolds            = 5,
          method            = "median",
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
          p_active          = 40,
          times_auc         = NULL,
          beta_true         = beta_true,
          metrics           = c("CIndex", "FNR", "FPR", "NSR"),
          verbose           = FALSE,
          plot_test         = TRUE)
```

The model performance yields significant predictive accuracy:

``` r
 # --- PERFORMANCE RESULTS ---

  data.frame(fit_simulation$fit_testing$performance)
```

## <span style="color: #444444;">Real Dataset</span>

In this section, we use the package to study a real genomic dataset
using the `LogLogistic - AFTNet` model. Unlike `COXNet`, this framework
allows for quantify how biomarkers directly accelerate or decelerate
survival time.

We utilize **TCGA-LUAD** (Lung Adenocarcinoma) RNA-seq data sourced from
the LinkedOmics portal (<https://www.linkedomics.org/>). The data are
processed into observed survival time $\mathbf{Y}$, censoring indicator
$\delta$. For each of the $n$ observations, $\mathbf{X}$ represents the
$p$-dimensional vector of genomic covariates.

To optimize efficiency, we bypass the time-consuming download process,
the network creation and the screening procedure by utilizing a
**pre-loaded example dataset**:

``` r
 # --- LOAD DATASET ---

  data(LUADdataset)
```

To ensure rigorous validation, the dataset is already partitioned into
an 70% training set for model fitting and a 30% testing set for
performance evaluation.

### <span style="color: #5c7997;">Gene Network Construction</span>

Prior biological knowledge, including gene–gene interactions and
co-expression relationships, is integrated into the model using
network-constrained graphs. This structural information is formalized
through the graph Laplacian matrix $\mathbf{L}$, which is constructed
using the following dedicated function on the training set:

``` r
 # --- NETWORK BUIDING USING CreateNetwork() ---

  net_result <- CreateNetwork(
      X            = LUADdataset$X_train,
      doid         = "DOID:1324",
      tissue       = "lung",
      disease_file = NULL,
      tissue_file  = NULL,
      choice       = 1,
      model        = NULL,
      dist         = NULL)
  
 # --- FINAL LAPLACIAN MATRIX ---
  
  L <- net_result$L   
  
 # --- DISEASE GENES AND SCORES ---
  
  disease_genes <- net_result$disease_genes$standard_name  
```

### <span style="color: #f5c59f;">Variable Screening</span>

To reduce dimensionality, we applied screening procedure to the training
set, shrinking the feature space from $p$ to a more tractable $d < p$.
We employed the BMD + DAD strategy, which integrates prior biological
knowledge with statistical evidence by selecting the union of covariates
identified by both methods.

``` r
 # --- SCREENING USING VariableScreening() ---

    screening_res <- VariableScreening(
        LUADdataset$X_train,
        disease_genes = disease_genes,
        screening     = "BMD")

 # --- LIST OF d SCREENED VARIABLES ---

    selected_genes <- screening_res$screen_vars   
    
    index_screen   <- match(selected_genes, colnames(LUADdataset$X_train))
    index_screen   <- index_screen[!is.na(index_screen)]
    
    X_train        <- as.matrix(LUADdataset$X_train[, index_screen])
    X_test         <- as.matrix(LUADdataset$X_test[, index_screen])
    L              <- L[index_screen, index_screen]   
```

The data in LUADdataset are screened using “BMD” approach.

### <span style="color: #f5c59f;">Training Phase</span>

We now execute the training phase to identify the prognostic gene
signature and determine the optimal cut-off value for patient
stratification into prognostic groups.

``` r
 # --- FIT THE MODEL USING NetSurvProx_Training() ---

  fit_training <- NetSurvProx_Training(
          X_train,
          LUADdataset$Y_train,
          LUADdataset$delta_train,
          L             = L,
          model         = "AFTNet",
          dist          = "loglogistic",
          select_lambda = TRUE,
          alpha_grid    = c(0.3, 0.5, 0.7),
          nlambda       = 50,
          lambda_ratio  = 0.01,
          nfolds        = 5,
          method        = "median",
          probs         = seq(0.25, 0.80, by = 0.05),
          cutoffplot    = TRUE,
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
          palette       = NULL)

 # --- ESTIMATED REGRESSION COEFFICIENTS ---

    beta <- fit_training$beta
    
 # --- OPTIMAL CUT-OFF FOR STRATIFICATION ---

    opt_cutoff <- fit_training$cutoff.opt 
```

### <span style="color: #5c7997;">Pathway Analysis</span>

We map the selected gene signature to KEGG pathways to **identify
functional relationships**. These pathways are then used to construct
subnetworks, providing a biological interpretation of the identified
prognostic genes.

``` r
 # --- PATHWAY ANALYSIS USIGN PathwayDashboard() ---

    genes_list <- rownames(beta)[beta != 0]
    
    PathwayDashboard(
        genes_list,
        useKeggAPI       = TRUE,
        pathwayFile      = NULL,
        nodesCols        = c("#5C7997", "#F5C59F"),
        diseaseNodes     = TRUE,
        disease_file     = disease_genes,
        top_percent      = 20,
        batch_size       = 10,
        background_genes = colnames(X_train),
        min_genes        = 2,
        top_n            = 10,
        db_name          = "org.Hs.eg.db",
        organism         = "hsa",
        out_dir          = NULL,
        open_browser     = TRUE,
        verbose          = FALSE)
```

The partial network visualizes only functionally connected genes, while
the full network displays the complete gene signature, including
isolated nodes. To enhance clinical relevance, the dashboard highlights
the top 20% of genes most strongly associated with the disease. This
structural view is complemented by an enrichment panel showing the top
10 pathways identified via over-representation analysis.

### <span style="color: #8abfe7;">Testing Phase</span>

This final phase evaluates the model’s predictive performance and
generalizability. Using the regression coefficients and the optimal
cut-off value from the training phase, we calculate the prognostic index
for the test set and stratify patients into prognostic groups. Model
accuracy is then assessed via specific performance metrics.

``` r
  # --- FIT THE MODEL USING NetSurvProx_Testing() ---

  fit_testing <- NetSurvProx_Testing(
                    X_train,
                    standardize = TRUE,
                    LUADdataset$Y_train,
                    LUADdataset$delta_train,
                    X_test, 
                    LUADdataset$Y_test,
                    LUADdataset$delta_test,
                    model = "AFTNet",
                    dist  = "loglogistic",
                    beta, opt_cutoff,
                    times_auc = NULL,
                    metrics   = c("CIndex", "NSR", "AUC"),
                    verbose   = FALSE,
                    plot      = FALSE,
                    palette   = NULL)
```

The model performance can be showed by the command:

``` r
 # --- PERFORMANCE RESULTS ---

  data.frame(fit_testing$performance)
```
