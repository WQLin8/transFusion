#############################################
## Data Preprocessing and Utility Functions
#############################################

#' Process scRNA-seq and stRNA-seq Data
#'
#' This function cleans and processes the scRNA-seq and stRNA-seq datasets by filtering out low-quality cells/spots,
#' selecting common genes, and applying minimal UMI count thresholds.
#'
#' @param sc_exp A matrix of scRNA-seq counts (genes x cells, raw counts).
#' @param sc_label A named vector of cell type labels (names should match sc_exp column names).
#' @param spot_exp A matrix of stRNA-seq counts (genes x spots, raw counts).
#' @param spot_loc A data.frame or matrix containing spot spatial coordinates. Row names must match spot_exp column names and include at least columns "x" and "y".
#' @param gene_det_in_min_cells_per Minimum fraction of cells in which a gene must be detected (default: 0.001).
#' @param expression_threshold Threshold to consider a gene expressed (default: 1).
#' @param nUMI Minimum UMI count for a cell or spot to be kept (default: 0).
#' @param verbose Logical; if TRUE, print processing details.
#' @param plot Logical; if TRUE, plot diagnostic histograms.
#'
#' @return A list containing:
#' \describe{
#'   \item{sc_exp}{Processed scRNA-seq expression matrix.}
#'   \item{sc_label}{Filtered cell type labels.}
#'   \item{spot_exp}{Processed stRNA-seq expression matrix.}
#'   \item{spot_loc}{Filtered spatial coordinates.}
#' }
data_process <- function(sc_exp, sc_label, spot_exp, spot_loc,
                         gene_det_in_min_cells_per = 0.001, expression_threshold = 1,
                         nUMI = 0, verbose = FALSE, plot = FALSE) {
  # Check that the number of scRNA-seq labels matches the number of cells
  if(ncol(sc_exp) != length(sc_label))
    stop("The number of cell labels does not match the number of scRNA-seq cells!")
  
  # Check that the number of spots matches the number of spatial locations
  if(ncol(spot_exp) != nrow(spot_loc))
    stop("The number of spots does not match the number of coordinate rows!")
  
  # Process scRNA-seq data (here no additional cleaning is applied; modify if needed)
  sc_matrix <- sc_exp
  sc_label <- sc_label[colnames(sc_exp)]
  
  # Process stRNA-seq data using cleanCounts
  st_matrix <- t(cleanCounts(t(as.matrix(spot_exp)),
                             gene_det_in_min_cells_per = gene_det_in_min_cells_per,
                             expression_threshold = expression_threshold,
                             nUMI = nUMI,
                             verbose = verbose, plot = plot))
  st_matrix <- as.matrix(st_matrix)
  ind_sp <- match(colnames(st_matrix), colnames(spot_exp))
  spot_loc <- spot_loc[ind_sp, ]
  
  # Keep only common genes between scRNA-seq and stRNA-seq datasets
  common_genes <- intersect(rownames(sc_matrix), rownames(st_matrix))
  sc_exp <- sc_matrix[common_genes, ]
  st_exp <- st_matrix[common_genes, ]
  
  # Filter cells and spots by total UMI counts
  idx_sc <- colSums(sc_exp) >= nUMI
  sc_exp_filter <- sc_exp[, idx_sc]
  sc_label_filter <- sc_label[idx_sc]
  
  idx_st <- colSums(st_exp) >= nUMI
  st_exp_filter <- st_exp[, idx_st]
  spot_loc_filter <- spot_loc[idx_st, ]
  
  list(sc_exp = sc_exp_filter, sc_label = sc_label_filter,
       spot_exp = st_exp_filter, spot_loc = spot_loc_filter)
}


#' Clean Counts Matrix
#'
#' This function filters the input count matrix by retaining only those genes that are expressed
#' in a sufficient number of cells and cells that have a minimal UMI count.
#'
#' @param counts Input count matrix.
#' @param gene_det_in_min_cells_per Minimum fraction of cells required to detect a gene.
#' @param expression_threshold Expression threshold to consider a gene as expressed.
#' @param nUMI Minimum total UMI count required for a cell.
#' @param verbose Logical; if TRUE, print the dimensions of the filtered matrix.
#' @param plot Logical; if TRUE, plot histograms of the gene/cell sums.
#'
#' @return A filtered count matrix.
cleanCounts <- function(counts, gene_det_in_min_cells_per = 0.001,
                        expression_threshold = 1, nUMI = 0,
                        verbose = FALSE, plot = FALSE) {
  n <- nrow(counts)
  # Keep genes detected in at least (gene_det_in_min_cells_per * n) cells
  keep_genes <- Matrix::colSums(counts >= expression_threshold) >= gene_det_in_min_cells_per * n
  # Keep cells that have a UMI sum >= nUMI (based on the filtered genes)
  keep_cells <- Matrix::rowSums(counts[, keep_genes] >= expression_threshold) >= nUMI
  counts <- counts[keep_cells, keep_genes]
  
  if (verbose) {
    message("Filtered matrix has ", nrow(counts), " cells and ", ncol(counts), " genes")
  }
  if (plot) {
    par(mfrow = c(1, 2), mar = rep(5, 4))
    hist(log10(Matrix::rowSums(counts) + 1), breaks = 20, main = 'Genes per Cell')
    hist(log10(Matrix::colSums(counts) + 1), breaks = 20, main = 'Cells per Gene')
  }
  counts
}


#' Calculate Average Expression by Cell Type
#'
#' This function groups cells by their cell type and calculates the average gene expression
#' for each group.
#'
#' @param sc_exp A scRNA-seq expression matrix.
#' @param sc_label A vector of cell type labels corresponding to the columns of sc_exp.
#'
#' @return A matrix of average gene expression per cell type (genes x cell types).
create_group_exp <- function(sc_exp, sc_label) {
  cell_types <- sort(unique(sc_label))
  group_indices <- lapply(cell_types, function(ct) which(sc_label == ct))
  group_exp <- sapply(group_indices, function(idx) Matrix::rowMeans(sc_exp[, idx]))
  group_exp <- as.matrix(group_exp)
  colnames(group_exp) <- cell_types
  group_exp
}


#############################################
## Interfaces for Individual Deconvolution Methods
#############################################

#' Run SCDC Deconvolution
#'
#' This function applies the SCDC deconvolution method.
#'
#' @param database A list containing processed scRNA-seq expression, cell labels, stRNA-seq expression, and spot locations.
#' @param iter.max Maximum number of iterations for SCDC (default: 1000).
#'
#' @return A matrix of estimated cell type proportions (spots x cell types).
SCDC_run <- function(database, iter.max = 1000) {
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_types <- sort(unique(sc_label))
  
  # Build scRNA-seq ExpressionSet
  sc_pData <- data.frame(cluster = as.character(sc_label),
                         sample = as.character(colnames(sc_exp)))
  rownames(sc_pData) <- colnames(sc_exp)
  sc_phenoData <- methods::new("AnnotatedDataFrame", data = sc_pData)
  sc_es <- Biobase::ExpressionSet(as.matrix(sc_exp), phenoData = sc_phenoData, byrow = FALSE)
  
  # Build stRNA-seq ExpressionSet
  st_pData <- data.frame(sample = as.character(colnames(spot_exp)))
  rownames(st_pData) <- colnames(spot_exp)
  st_phenoData <- methods::new("AnnotatedDataFrame", data = st_pData)
  st_es <- Biobase::ExpressionSet(as.matrix(spot_exp), phenoData = st_phenoData, byrow = FALSE)
  
  result <- suppressMessages(SCDC::SCDC_prop(bulk.eset = st_es,
                                               sc.eset = sc_es,
                                               ct.varname = "cluster",
                                               ct.sub = cell_types,
                                               iter.max = iter.max))
  result$prop.est.mvw[, cell_types]
}


#' Run RCTD Deconvolution
#'
#' This function applies the RCTD deconvolution method.
#'
#' @param database A list containing processed scRNA-seq data, cell labels, stRNA-seq data, and spot locations.
#' @param CELL_MIN_INSTANCE Minimum number of cells per cell type in the reference (default: 20).
#'
#' @return A matrix of estimated cell type proportions (spots x cell types).
RCTD_run <- function(database, CELL_MIN_INSTANCE = 20) {
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_types <- sort(unique(sc_label))
  spot_loc <- database$spot_loc
  
  # Convert to sparse matrices
  sparse_sc_exp <- as(sc_exp, "sparseMatrix")
  sparse_spot_exp <- as(spot_exp, "sparseMatrix")
  
  cell_names <- colnames(sc_exp)
  cell_types_factor <- as.factor(sc_label)
  names(cell_types_factor) <- cell_names
  sc_nUMI <- as.numeric(colSums(sc_exp))
  names(sc_nUMI) <- cell_names
  
  reference <- spacexr::Reference(sparse_sc_exp, cell_types_factor, nUMI = sc_nUMI)
  
  coords <- as.data.frame(spot_loc)
  rownames(coords) <- colnames(spot_exp)
  puck <- spacexr::SpatialRNA(coords, counts = sparse_spot_exp, nUMI = colSums(spot_exp))
  
  rctd_obj <- suppressMessages(spacexr::create.RCTD(puck, reference, max_cores = 1, CELL_MIN_INSTANCE = CELL_MIN_INSTANCE))
  rctd_obj <- suppressMessages(spacexr::run.RCTD(rctd_obj, doublet_mode = 'full'))
  
  weights <- as.matrix(rctd_obj@results$weights)
  norm_weights <- sweep(weights, 1, rowSums(weights), '/')
  norm_weights[, cell_types]
}


#' Run MuSiC Deconvolution
#'
#' This function applies the MuSiC deconvolution method (both weighted and all-gene versions).
#'
#' @param database A list containing processed scRNA-seq data, cell labels, and stRNA-seq data.
#' @param iter.max Maximum iterations for training (default: 1000).
#' @param nu Regularization parameter (default: 1e-04).
#' @param eps Convergence threshold (default: 0.01).
#'
#' @return A list containing:
#' \describe{
#'   \item{Music_weight}{Estimated proportions using the weighted approach.}
#'   \item{Music_allgene}{Estimated proportions using all genes.}
#' }
MuSiC_run <- function(database, iter.max = 1000, nu = 1e-04, eps = 0.01) {
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_types <- sort(unique(sc_label))
  
  # Convert cell type labels to numeric encoding
  sc_label_num <- as.numeric(factor(sc_label, levels = sort(unique(sc_label))))
  
  # Build scRNA-seq metadata
  sc_pData <- data.frame(sampleID = 1:ncol(sc_exp),
                         SubjectName = factor(colnames(sc_exp)),
                         cellTypeID = sc_label_num,
                         cellType = factor(sc_label))
  rownames(sc_pData) <- colnames(sc_exp)
  sc_phenoData <- methods::new("AnnotatedDataFrame", data = sc_pData)
  sc_es <- Biobase::ExpressionSet(as.matrix(sc_exp), phenoData = sc_phenoData, byrow = FALSE)
  
  # Build stRNA-seq ExpressionSet
  st_pData <- data.frame(sampleID = 1:ncol(spot_exp), SubjectName = factor(colnames(spot_exp)))
  rownames(st_pData) <- colnames(spot_exp)
  st_phenoData <- methods::new("AnnotatedDataFrame", data = st_pData)
  st_es <- Biobase::ExpressionSet(as.matrix(spot_exp), phenoData = st_phenoData, byrow = FALSE)
  
  result <- MuSiC::music_prop(bulk.eset = st_es, sc.eset = sc_es,
                              clusters = 'cellType',
                              samples = 'sampleID', select.ct = cell_types,
                              verbose = FALSE, iter.max = iter.max, nu = nu, eps = eps)
  Music_weight <- as.matrix(result$Est.prop.weighted)[, cell_types]
  Music_allgene <- as.matrix(result$Est.prop.allgene)[, cell_types]
  list(Music_weight = Music_weight, Music_allgene = Music_allgene)
}


#' Run DeconRNASeq Deconvolution
#'
#' This function applies the DeconRNASeq method for deconvolution.
#'
#' @param database A list containing processed scRNA-seq data, cell labels, and stRNA-seq data.
#' @param perc Fraction used to filter genes (default: 0.05).
#'
#' @return A matrix of estimated cell type proportions (spots x cell types).
DeconRNASeq_run <- function(database, perc = 0.05) {
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_types <- sort(unique(sc_label))
  
  # Filter genes based on total expression
  keep_genes <- rowSums(sc_exp) >= perc * ncol(sc_exp)
  sc_exp <- sc_exp[keep_genes, ]
  
  # Normalize and log-transform scRNA-seq and stRNA-seq data
  sc_exp_norm <- t(median(colSums(sc_exp)) * (t(sc_exp) / colSums(sc_exp)))
  sc_exp_norm <- log(sc_exp_norm + 1)
  spot_exp_norm <- t(median(colSums(spot_exp)) * (t(spot_exp) / colSums(spot_exp)))
  spot_exp_norm <- log(spot_exp_norm + 1)
  
  type_exp <- sapply(cell_types, function(ct) rowMeans(sc_exp_norm[, sc_label == ct]))
  type_exp <- as.matrix(type_exp)
  colnames(type_exp) <- cell_types
  
  decon_out <- DeconRNASeq::DeconRNASeq(as.data.frame(spot_exp_norm),
                                        as.data.frame(type_exp),
                                        checksig = FALSE, use.scale = TRUE)
  res <- decon_out$out.all
  rownames(res) <- colnames(spot_exp)
  res[, cell_types]
}


#' Run DWLS and SVR Deconvolution
#'
#' This function applies the Dampened Weighted Least Squares (DWLS) method and the SVR method
#' for deconvolution. It optionally runs in parallel.
#'
#' @param database A list containing processed scRNA-seq data, cell labels, and stRNA-seq data.
#' @param parallel Logical; if TRUE, use parallel processing (default: TRUE).
#' @param is_select_DEGs Logical; if TRUE, perform gene selection based on differential expression (default: TRUE).
#' @param python_env Path to the Python environment (required for some functions).
#'
#' @return A list containing:
#' \describe{
#'   \item{DampenedWLS}{Deconvolution proportions from DWLS.}
#'   \item{SVR}{Deconvolution proportions from SVR.}
#'   \item{time.DWLS}{Execution time for DWLS.}
#'   \item{time.SVR}{Execution time for SVR.}
#' }
DWLS_run <- function(database, parallel = TRUE, is_select_DEGs = TRUE, python_env) {
  # Scale the data for DWLS (by 1e-4)
  sc_exp <- database$sc_exp * 1e-4
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp * 1e-4
  cell_types <- sort(unique(sc_label))
  
  if(is_select_DEGs) {
    instrs <- Giotto::createGiottoInstructions(python_path = python_env)
    sc_obj <- Giotto::createGiottoObject(raw_exprs = sc_exp, instructions = instrs)
    sc_obj <- Giotto::normalizeGiotto(sc_obj)
    sc_obj@cell_metadata$leiden_clus <- as.character(sc_label)
    markers <- Giotto::findMarkers_one_vs_all(gobject = sc_obj,
                                               method = 'gini',
                                               expression_values = 'normalized',
                                               cluster_column = 'leiden_clus',
                                               verbose = FALSE)
    top_genes <- markers[, head(.SD, 100), by = 'cluster']
    sc_exp <- sc_exp[top_genes$genes, ]
    spot_exp <- spot_exp[top_genes$genes, ]
  }
  
  cell_type_exp <- create_group_exp(sc_exp, sc_label)
  
  if(parallel) {
    num.cores <- if(nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_", "")) &&
                      Sys.getenv("_R_CHECK_LIMIT_CORES_") == "TRUE") 2L else parallel::detectCores() - 1
    cl <- parallel::makeCluster(num.cores)
    doParallel::registerDoParallel(cl)
    
    start_dwls <- Sys.time()
    DampenedWLS <- foreach::foreach(i = 1:ncol(spot_exp), .combine = 'rbind', .inorder = TRUE) %dopar% {
      DWLS::solveDampenedWLS(cell_type_exp, spot_exp[, i])
    }
    time_dwls <- difftime(Sys.time(), start_dwls, units = "mins")
    
    start_svr <- Sys.time()
    SVR <- foreach::foreach(i = 1:ncol(spot_exp), .combine = 'rbind', .inorder = TRUE) %dopar% {
      DWLS::solveSVR(cell_type_exp, spot_exp[, i])
    }
    time_svr <- difftime(Sys.time(), start_svr, units = "mins")
    
    parallel::stopCluster(cl)
  } else {
    DampenedWLS <- matrix(NA, ncol(spot_exp), length(cell_types))
    SVR <- matrix(NA, ncol(spot_exp), length(cell_types))
    
    start_dwls <- Sys.time()
    for(i in seq_len(ncol(spot_exp))) {
      DampenedWLS[i, ] <- DWLS::solveDampenedWLS(cell_type_exp, spot_exp[, i])
    }
    time_dwls <- difftime(Sys.time(), start_dwls, units = "mins")
    
    start_svr <- Sys.time()
    for(i in seq_len(ncol(spot_exp))) {
      SVR[i, ] <- DWLS::solveSVR(cell_type_exp, spot_exp[, i])
    }
    time_svr <- difftime(Sys.time(), start_svr, units = "mins")
  }
  
  rownames(DampenedWLS) <- colnames(spot_exp)
  rownames(SVR) <- colnames(spot_exp)
  colnames(DampenedWLS) <- cell_types
  colnames(SVR) <- cell_types
  
  list(DampenedWLS = DampenedWLS, SVR = SVR, time.DWLS = time_dwls, time.SVR = time_svr)
}


#' Run SPOTlight Deconvolution
#'
#' This function applies the SPOTlight deconvolution method based on non‐negative matrix factorization (NMF).
#'
#' @param database A list containing processed scRNA-seq data, cell labels, and stRNA-seq data.
#' @param cl_n Number of cells to sample per cluster (default: 10000).
#' @param hvg Number of highly variable genes to use (default: 3000).
#' @param min_cont Minimum contribution threshold (default: 0.001).
#'
#' @return A matrix of estimated cell type proportions (spots x cell types).
SPOTlight_run <- function(database, cl_n = 10000, hvg = 3000, min_cont = 0.001) {
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_types <- sort(unique(sc_label))
  
  # Create Seurat object and perform SCTransform normalization
  sc_ref <- Seurat::CreateSeuratObject(counts = sc_exp)
  sc_ref@meta.data$subclass <- as.factor(sc_label)
  sc_ref <- Seurat::SCTransform(sc_ref, verbose = FALSE)
  Seurat::Idents(sc_ref) <- sc_ref@meta.data$subclass
  
  markers_all <- Seurat::FindAllMarkers(object = sc_ref, assay = "SCT",
                                          slot = "data", verbose = FALSE, only.pos = TRUE)
  
  se_sc_down <- SPOTlight::downsample_se_obj(se_obj = sc_ref, clust_vr = "subclass",
                                             cluster_markers = markers_all, cl_n = cl_n, hvg = hvg)
  nmf_mod_ls <- SPOTlight::train_nmf(cluster_markers = markers_all, se_sc = sc_ref,
                                     mtrx_spatial = spot_exp, clust_vr = "subclass",
                                     ntop = NULL, hvg = hvg, transf = "uv", method = "nsNMF")
  nmf_mod <- nmf_mod_ls[[1]]
  w <- NMF::basis(nmf_mod)  # Basis matrix
  h <- NMF::coef(nmf_mod)   # Coefficient matrix
  
  valid_genes <- !is.na(pmatch(rownames(spot_exp), rownames(w)))
  spot_exp_train <- spot_exp[valid_genes, ]
  
  ct_topic_profiles <- SPOTlight::topic_profile_per_cluster_nmf(h = h, train_cell_cl = nmf_mod_ls[[2]])
  decon_mtrx <- SPOTlight::mixture_deconvolution_nmf(nmf_mod = nmf_mod,
                                                     mixture_transcriptome = spot_exp_train,
                                                     transf = "uv",
                                                     reference_profiles = ct_topic_profiles,
                                                     min_cont = min_cont)
  results <- as.matrix(decon_mtrx[, -ncol(decon_mtrx)])
  rownames(results) <- colnames(spot_exp)
  colnames(results) <- cell_types
  results
}


#' Run SpatialDWLS Deconvolution
#'
#' This function applies the SpatialDWLS method using the Giotto package. It optionally selects
#' differentially expressed genes before deconvolution.
#'
#' @param database A list containing processed scRNA-seq data, cell labels, and stRNA-seq data.
#' @param my_python_path Path to the Python environment.
#' @param is_select_DEGs Logical; if TRUE, select differentially expressed genes (default: FALSE).
#' @param findmarker_method Method used for marker detection (default: "gini").
#' @param ncp_spa Number of principal components for spatial data (default: 100).
#' @param dimensions_to_use Number of dimensions used to build the KNN network (default: 10).
#' @param k Number of nearest neighbors (default: 10).
#' @param resolution Resolution parameter for clustering (default: 0.4).
#' @param n_iterations Number of iterations for clustering (default: 1000).
#' @param n_cell Number of cells per spot (default: 50).
#'
#' @return A matrix of estimated cell type proportions (spots x cell types).
spatialDWLS_run <- function(database, my_python_path, is_select_DEGs = FALSE,
                            findmarker_method = "gini", ncp_spa = 100,
                            dimensions_to_use = 10, k = 10, resolution = 0.4,
                            n_iterations = 1000, n_cell = 50) {
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_types <- sort(unique(sc_label))
  
  instrs <- Giotto::createGiottoInstructions(python_path = my_python_path)
  sc_obj <- Giotto::createGiottoObject(raw_exprs = sc_exp, instructions = instrs)
  sc_obj <- Giotto::normalizeGiotto(sc_obj)
  sc_obj@cell_metadata$leiden_clus <- as.character(sc_label)
  
  if(is_select_DEGs) {
    markers <- Giotto::findMarkers_one_vs_all(gobject = sc_obj,
                                               method = findmarker_method,
                                               expression_values = 'normalized',
                                               cluster_column = 'leiden_clus',
                                               verbose = FALSE)
    top_genes <- markers[, head(.SD, 100), by = 'cluster']
    sc_norm_exp <- 2^(sc_obj@norm_expr) - 1
    ExprSubset <- sc_norm_exp[as.character(top_genes$genes), ]
  } else {
    sc_norm_exp <- 2^(sc_obj@norm_expr) - 1
    ExprSubset <- sc_norm_exp
  }
  
  # Build signature matrix: average expression per cell type
  Sig <- sapply(unique(sc_label), function(ct) rowMeans(ExprSubset[, sc_label == ct]))
  colnames(Sig) <- unique(sc_label)
  
  grid_obj <- Giotto::createGiottoObject(raw_exprs = spot_exp, instructions = instrs)
  grid_obj <- Giotto::normalizeGiotto(grid_obj)
  grid_obj <- Giotto::calculateHVG(grid_obj, show_plot = FALSE, return_plot = FALSE)
  featgenes <- Giotto::fDataDT(grid_obj)[hvg == 'yes']$gene_ID
  grid_obj <- Giotto::runPCA(grid_obj, genes_to_use = featgenes, scale_unit = FALSE, ncp = ncp_spa)
  grid_obj <- Giotto::createNearestNetwork(grid_obj, dimensions_to_use = 1:dimensions_to_use, k = k)
  grid_obj <- Giotto::doLeidenCluster(grid_obj, resolution = resolution, n_iterations = n_iterations)
  grid_obj <- Giotto::runDWLSDeconv(grid_obj, sign_matrix = Sig, n_cell = n_cell)
  
  results <- as.matrix(grid_obj@spatial_enrichment$DWLS[, -1])
  results <- results[, cell_types]
  rownames(results) <- colnames(spot_exp)
  results
}


#' Run Stereoscope Deconvolution
#'
#' This function applies the Stereoscope deconvolution method via reticulate, using scvi-tools.
#'
#' @param database A list containing processed scRNA-seq data, cell labels, and stRNA-seq data.
#' @param python_env Path to the Python environment. If NULL, a default environment is used.
#' @param use_gpu Logical; whether to use GPU acceleration (default: FALSE).
#' @param select_HVG Logical; if TRUE, select highly variable genes (default: TRUE).
#' @param HVG_num Number of highly variable genes to select (default: 3000).
#' @param sc_training_plot Logical; if TRUE, plot training loss for scRNA-seq model.
#' @param sc_training_save_trained_model Logical; if TRUE, save the trained scRNA-seq model.
#' @param sc_max_epochs Maximum epochs for training the scRNA-seq model (default: 10000).
#' @param sc_lr Learning rate for scRNA-seq model training (default: 0.01).
#' @param st_training_plot Logical; if TRUE, plot training loss for stRNA-seq model.
#' @param st_training_save_trained_model Logical; if TRUE, save the trained stRNA-seq model.
#' @param st_max_epochs Maximum epochs for training the stRNA-seq model (default: 10000).
#' @param st_lr Learning rate for stRNA-seq model training (default: 0.01).
#'
#' @return A matrix of estimated cell type proportions (spots x cell types).
Stereoscope_run <- function(database, python_env = NULL, use_gpu = FALSE,
                            select_HVG = TRUE, HVG_num = 3000,
                            sc_training_plot = FALSE, sc_training_save_trained_model = FALSE,
                            sc_max_epochs = 10000, sc_lr = 0.01,
                            st_training_plot = FALSE, st_training_save_trained_model = FALSE,
                            st_max_epochs = 10000, st_lr = 0.01) {
  if(is.null(python_env)) {
    cat("Please set the path to your Python environment. Default Miniconda environment will be used otherwise.\n")
  } else {
    reticulate::use_python(python_env, require = TRUE)
    reticulate::py_config()
  }
  
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_types <- sort(unique(sc_label))
  
  if(nrow(sc_exp) != nrow(spot_exp))
    stop("The number of genes in scRNA-seq and stRNA-seq must be equal!")
  
  scvi <- reticulate::import("scvi")
  np <- reticulate::import("numpy")
  sc <- reticulate::import("scanpy")
  anndata <- reticulate::import("anndata")
  
  if(nrow(sc_exp) >= HVG_num && select_HVG) {
    seurat_sc <- Seurat::CreateSeuratObject(counts = sc_exp)
    seurat_sc <- Seurat::NormalizeData(seurat_sc, verbose = FALSE)
    seurat_sc <- Seurat::FindVariableFeatures(seurat_sc, nfeatures = HVG_num, verbose = FALSE)
    var.features <- seurat_sc@assays$RNA@var.features
    sc_exp_hvg <- sc_exp[var.features, ]
    spot_exp_hvg <- spot_exp[var.features, ]
    adata_sc <- anndata$AnnData(X = t(sc_exp_hvg))
    st_expr <- spot_exp_hvg
  } else {
    adata_sc <- anndata$AnnData(X = t(sc_exp))
    st_expr <- spot_exp
  }
  
  adata_sc$obs["cell_label"] <- sc_label
  scvi$external$stereoscope$RNAStereoscope$setup_anndata(adata_sc, labels_key = "cell_label")
  sc_model <- scvi$external$stereoscope$RNAStereoscope(adata_sc)
  scvi$external$stereoscope$RNAStereoscope$train(sc_model, use_gpu = use_gpu,
                                                 early_stopping_monitor = "elbo_validation",
                                                 max_epochs = as.integer(sc_max_epochs), lr = sc_lr)
  if(sc_training_plot) {
    loss_r <- np$asarray(sc_model$history["elbo_train"]$elbo_train)
    loss_df <- data.frame(epoch = 40:nrow(loss_r), loss = unlist(loss_r[40:nrow(loss_r)]))
    plot(loss_df, type = "l", main = "scRNA-seq Training Loss")
  }
  if(sc_training_save_trained_model) {
    scvi$external$stereoscope$RNAStereoscope$save(sc_model, "scmodel")
  }
  
  adata_st <- anndata$AnnData(X = t(st_expr))
  scvi$external$stereoscope$SpatialStereoscope$setup_anndata(adata_st)
  st_model <- scvi$external$stereoscope$SpatialStereoscope$from_rna_model(adata_st, sc_model, prior_weight = "minibatch")
  scvi$external$stereoscope$SpatialStereoscope$train(st_model,
                                                     early_stopping_monitor = "elbo_train", use_gpu = use_gpu,
                                                     max_epochs = as.integer(st_max_epochs), lr = st_lr)
  if(st_training_plot) {
    loss_r <- np$asarray(st_model$history["elbo_train"]$elbo_train)
    loss_df <- data.frame(epoch = 40:nrow(loss_r), loss = unlist(loss_r[40:nrow(loss_r)]))
    plot(loss_df, type = "l", main = "stRNA-seq Training Loss")
  }
  weights <- as.matrix(scvi$external$stereoscope$SpatialStereoscope$get_proportions(st_model))
  rownames(weights) <- colnames(st_expr)
  weights[, cell_types]
}


#' Run cell2location Deconvolution
#'
#' This function applies the cell2location deconvolution method via reticulate and Python.
#'
#' @param database A list containing processed scRNA-seq data, cell labels, and stRNA-seq data.
#' @param python_env Path to the Python environment.
#' @param sc_max_epoches Maximum epochs for training scRNA-seq model (default: 3000).
#' @param sc_lr Learning rate for scRNA-seq model (default: 0.002).
#' @param sc_use_gpu Logical; whether to use GPU for scRNA-seq training (default: TRUE).
#' @param st_N_cells_per_location Number of cells per spatial location (default: 0).
#' @param st_detection_alpha Regularization parameter (default: 20.00).
#' @param st_max_epoches Maximum epochs for training stRNA-seq model (default: 3000).
#'
#' @return A normalized matrix of estimated cell type proportions (spots x cell types).
cell2location_run <- function(database, python_env = NULL,
                              sc_max_epoches = 3000, sc_lr = 0.002, sc_use_gpu = TRUE,
                              st_N_cells_per_location = 0, st_detection_alpha = 20.00,
                              st_max_epoches = 3000) {
  if(is.null(python_env)) {
    cat("Please set the path to your Python environment. Default packages will be installed otherwise.\n")
  } else {
    reticulate::use_python(python_env, require = TRUE)
    reticulate::py_config()
  }
  
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_types <- sort(unique(sc_label))
  
  if(nrow(sc_exp) != nrow(spot_exp))
    stop("The number of genes in scRNA-seq and stRNA-seq must be equal!")
  
  np <- reticulate::import("numpy")
  anndata <- reticulate::import("anndata")
  cell2loc <- reticulate::import("cell2location")
  pd <- reticulate::import("pandas")
  
  # Train the scRNA-seq model
  adata_sc <- input_for_py(expr = sc_exp, cell_type = sc_label)
  cell2loc$models$RegressionModel$setup_anndata(adata_sc, labels_key = "cell_type")
  mod <- cell2loc$models$RegressionModel(adata_sc)
  batch_size <- if(ncol(sc_exp) > 10000) as.integer(2500) else as.integer(floor(ncol(sc_exp)/10))
  cell2loc$models$RegressionModel$train(mod, max_epochs = as.integer(sc_max_epoches),
                                          batch_size = batch_size, lr = sc_lr, use_gpu = sc_use_gpu)
  adata_sc <- cell2loc$models$RegressionModel$export_posterior(mod, adata_sc,
                                                                 sample_kwargs = list(batch_size = batch_size, use_gpu = sc_use_gpu))
  inf_aver <- data_frame_extracted(adata_sc)
  
  # Train the stRNA-seq model
  adata_st <- input_for_py(expr = spot_exp)
  cell2loc$models$Cell2location$setup_anndata(adata_st)
  mod_st <- cell2loc$models$Cell2location(adata_st, cell_state_df = inf_aver,
                                            N_cells_per_location = as.integer(st_N_cells_per_location),
                                            detection_alpha = st_detection_alpha)
  cell2loc$models$Cell2location$train(mod_st, max_epochs = as.integer(st_max_epoches),
                                        train_size = 1.00, use_gpu = sc_use_gpu)
  adata_st <- cell2loc$models$Cell2location$export_posterior(mod_st, adata_st,
                                                               sample_kwargs = list(batch_size = ncol(spot_exp), use_gpu = sc_use_gpu))
  
  adata_temp <- reticulate::py_to_r(adata_st)
  data_total <- adata_temp$obsm["q05_cell_abundance_w_sf"]
  factor_names <- adata_temp$uns['mod'][[1]]$factor_names
  slice_idx <- paste0("q05cell_abundance_w_sf_", factor_names)
  data_extracted <- data_total[slice_idx]
  names(data_extracted) <- factor_names
  
  weight_matrix <- matrix(unlist(data_extracted), ncol = length(data_extracted))
  colnames(weight_matrix) <- factor_names
  rownames(weight_matrix) <- rownames(data_total)
  weight_matrix <- weight_matrix[, cell_types]
  norm_matrix <- sweep(weight_matrix, 1, rowSums(weight_matrix), "/")
  norm_matrix
}


#' Run STdeconvolve Method
#'
#' This function applies the STdeconvolve method to deconvolve spatial transcriptomics data.
#'
#' @param database A list containing processed scRNA-seq data, cell labels, and stRNA-seq data.
#' @param min.lib.size Minimum library size per cell (default: 100).
#' @param min.reads Minimum reads per gene (default: 1).
#' @param nTopOD Number of top over-dispersed genes to retain (default: 1000).
#' @param betaScale Scaling factor for gene expression (default: 1000).
#'
#' @return A matrix of estimated cell type proportions (spots x cell types).
STdeconvolve_run <- function(database, min.lib.size = 100, min.reads = 1, nTopOD = 1000, betaScale = 1000) {
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_types <- sort(unique(sc_label))
  
  counts <- STdeconvolve::cleanCounts(counts = spot_exp, min.lib.size = min.lib.size, min.reads = min.reads)
  corpus <- STdeconvolve::restrictCorpus(counts = counts, nTopOD = nTopOD)
  ldas <- STdeconvolve::fitLDA(counts = t(as.matrix(corpus)), Ks = c(length(cell_types)),
                              plot = FALSE, verbose = FALSE)
  optLDA <- STdeconvolve::optimalModel(models = ldas, opt = c(length(cell_types)))
  results <- STdeconvolve::getBetaTheta(lda = optLDA, corpus = t(as.matrix(corpus)),
                                        betaScale = betaScale, verbose = FALSE)
  deconProp <- results$theta
  
  # Annotate topics by comparing to scRNA-seq reference
  deconGexp <- results$beta * 1000
  mobProxyTheta2 <- model.matrix(~0 + sc_label)
  rownames(mobProxyTheta2) <- colnames(sc_exp)
  mobProxyTheta2 <- sc_exp %*% mobProxyTheta2
  corMtx_beta <- STdeconvolve::getCorrMtx(m1 = as.matrix(deconGexp),
                                          m2 = t(as.matrix(mobProxyTheta2)),
                                          type = "b", verbose = FALSE)
  rownames(corMtx_beta) <- paste0("decon_", seq(nrow(corMtx_beta)))
  pairs_used <- STdeconvolve::lsatPairs(t(corMtx_beta))
  m <- t(corMtx_beta)[pairs_used$rowix, pairs_used$colsix]
  
  ct_annotation <- sub("sc_label", "", rownames(m))
  topic_annotation <- sub("decon_", "", colnames(m))
  deconProp_final <- deconProp[, topic_annotation]
  colnames(deconProp_final) <- ct_annotation
  deconProp_final <- deconProp_final[, cell_types]
  deconProp_final <- sweep(deconProp_final, 1, rowSums(deconProp_final), "/")
  deconProp_final
}


#' Run DestVI Deconvolution
#'
#' This function applies the DestVI deconvolution method via reticulate and scvi-tools.
#'
#' @param database A list containing processed scRNA-seq data, cell labels, and stRNA-seq data.
#' @param python_env Path to the Python environment.
#' @param n_top_genes Number of top highly variable genes to select (default: 2000).
#' @param use_gpu Logical; whether to use GPU (default: FALSE).
#' @param max_iter_sc Maximum epochs for scRNA-seq training (default: 400).
#' @param max_iter_st Maximum epochs for stRNA-seq training (default: 3000).
#'
#' @return A normalized matrix of estimated cell type proportions (spots x cell types).
DestVI_run <- function(database, python_env = NULL,
                       n_top_genes = 2000, use_gpu = FALSE,
                       max_iter_sc = 400, max_iter_st = 3000) {
  if(is.null(python_env)) {
    cat("Please set the path to your Python environment. Default packages will be installed otherwise.\n")
  } else {
    reticulate::use_python(python_env, require = TRUE)
    reticulate::py_config()
  }
  
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_types <- sort(unique(sc_label))
  
  if(nrow(sc_exp) != nrow(spot_exp))
    stop("The number of genes in scRNA-seq and stRNA-seq must be equal!")
  
  scvi <- reticulate::import("scvi")
  anndata <- reticulate::import("anndata")
  pd <- reticulate::import("pandas")
  sc <- reticulate::import("scanpy")
  
  # Process scRNA-seq data via Seurat
  meta.data <- data.frame(cell_type = sc_label)
  rownames(meta.data) <- colnames(sc_exp)
  seurat_sc <- Seurat::CreateSeuratObject(counts = sc_exp, meta.data = meta.data)
  seurat_sc <- Seurat::NormalizeData(seurat_sc, verbose = FALSE)
  seurat_sc <- Seurat::FindVariableFeatures(seurat_sc, nfeatures = n_top_genes, verbose = FALSE)
  sc_exp_norm <- as.matrix(seurat_sc@assays$RNA@data)[seurat_sc@assays$RNA@var.features, ]
  sc_label_filter <- seurat_sc$cell_type
  
  # Process stRNA-seq data via Seurat
  seurat_st <- Seurat::CreateSeuratObject(counts = spot_exp)
  seurat_st <- Seurat::NormalizeData(seurat_st, verbose = FALSE)
  spot_exp_norm <- as.matrix(seurat_st@assays$RNA@data)[rownames(sc_exp_norm), ]
  
  adata_sc <- input_for_py(expr = sc_exp_norm, cell_type = sc_label_filter)
  scvi$model$CondSCVI$setup_anndata(adata_sc, labels_key = "cell_type")
  sc_model <- scvi$model$CondSCVI(adata_sc, weight_obs = FALSE)
  scvi$model$CondSCVI$train(sc_model, use_gpu = use_gpu, max_epochs = as.integer(max_iter_sc))
  
  adata_st <- input_for_py(expr = spot_exp_norm)
  scvi$model$DestVI$setup_anndata(adata_st)
  st_model <- scvi$model$DestVI$from_rna_model(adata_st, sc_model)
  scvi$model$DestVI$train(st_model, use_gpu = use_gpu, max_epochs = as.integer(max_iter_st))
  
  weight_matrix <- as.matrix(scvi$model$DestVI$get_proportions(st_model))[, cell_types]
  sweep(weight_matrix, 1, rowSums(weight_matrix), "/")
}


#############################################
## Utility Functions for Python Interoperability
#############################################

#' Convert Expression Matrix to an AnnData Object
#'
#' This function writes an expression matrix to a CSV file, reads it using scanpy, and returns an AnnData object.
#'
#' @param expr Expression matrix.
#' @param cell_type Optional vector of cell type labels.
#'
#' @return An AnnData object.
input_for_py <- function(expr, cell_type = NULL) {
  sc <- reticulate::import("scanpy")
  anndata <- reticulate::import("anndata")
  utils::write.csv(expr, "count.csv")
  adata <- sc$read_csv("count.csv", first_column_names = TRUE, dtype = 'float32')$T
  if(!is.null(cell_type))
    adata$obs["cell_type"] <- cell_type
  unlink("count.csv", recursive = TRUE)
  adata
}


#' Extract Data Frame from AnnData Object
#'
#' This function extracts a data frame from an AnnData object (used in cell2location).
#'
#' @param adata An AnnData object.
#'
#' @return A Python object (converted from R) containing the extracted data.
data_frame_extracted <- function(adata) {
  adata <- reticulate::py_to_r(adata)
  data_total <- adata$varm['q05_per_cluster_mu_fg']
  factor_names <- adata$uns['mod'][[1]]$factor_names
  slice_idx <- paste0("q05_per_cluster_mu_fg_", factor_names)
  data_extracted <- data_total[, slice_idx]
  names(data_extracted) <- factor_names
  reticulate::r_to_py(data_extracted)
}


#' Run CARD Deconvolution
#'
#' This function applies the CARD deconvolution method.
#'
#' @param database A list containing processed scRNA-seq data, cell labels, stRNA-seq data, and spatial coordinates.
#' @param minCountGene Minimum count for each gene (default: 100).
#' @param minCountSpot Minimum count for each spot (default: 5).
#'
#' @return A matrix of estimated cell type proportions (spots x cell types).
CARD_run <- function(database, minCountGene = 100, minCountSpot = 5) {
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  spot_loc <- database$spot_loc
  cell_types <- sort(unique(sc_label))
  
  original_spot_names <- colnames(spot_exp)
  new_spot_names <- paste0("spot", seq_len(nrow(spot_loc)))
  colnames(spot_exp) <- new_spot_names
  rownames(spot_loc) <- new_spot_names
  colnames(spot_loc) <- c("x", "y")
  
  sampleInfo <- rep("sample1", ncol(sc_exp))
  names(sampleInfo) <- colnames(sc_exp)
  meta.data <- data.frame(cellID = colnames(sc_exp), cellType = sc_label, sampleInfo = sampleInfo)
  rownames(meta.data) <- colnames(sc_exp)
  
  CARD_obj <- CARD::createCARDObject(sc_count = sc_exp,
                                     sc_meta = meta.data,
                                     spatial_count = spot_exp,
                                     spatial_location = spot_loc,
                                     ct.varname = "cellType",
                                     ct.select = cell_types,
                                     sample.varname = "sampleInfo",
                                     minCountGene = minCountGene,
                                     minCountSpot = minCountSpot)
  CARD_obj <- CARD::CARD_deconvolution(CARD_object = CARD_obj)
  CARD_results <- as.matrix(CARD_obj@Proportion_CARD)
  original_names <- original_spot_names[new_spot_names %in% rownames(CARD_results)]
  CARD_results <- CARD_results[new_spot_names %in% rownames(CARD_results), cell_types]
  rownames(CARD_results) <- original_names
  CARD_results
}


#############################################
## Ensemble Functions
#############################################

#' Run Individual Deconvolution Methods
#'
#' This function runs multiple individual deconvolution methods and returns their results along with execution times.
#'
#' @param sc_exp scRNA-seq expression matrix (genes x cells).
#' @param sc_label Vector of cell type labels.
#' @param spot_exp stRNA-seq expression matrix (genes x spots).
#' @param spot_loc Spatial coordinates for spots.
#' @param gene_det_in_min_cells_per Minimum fraction of cells a gene must be detected in.
#' @param expression_threshold Expression threshold for gene detection.
#' @param nUMI Minimum UMI count for cells/spots.
#' @param verbose Logical; whether to print processing details.
#' @param plot Logical; whether to plot diagnostic plots.
#' @param python_env Path to the Python environment.
#' @param use_gpu Logical; whether to use GPU for methods that support it.
#' @param saving_results Logical; if TRUE, save the deconvolution results to a file.
#' @param SCDC Logical; if TRUE, run SCDC.
#' @param RCTD Logical; if TRUE, run RCTD.
#' @param MuSiC Logical; if TRUE, run MuSiC.
#' @param DeconRNASeq Logical; if TRUE, run DeconRNASeq.
#' @param DestVI Logical; if TRUE, run DestVI.
#' @param DWLS Logical; if TRUE, run DWLS (and SVR).
#' @param SPOTlight Logical; if TRUE, run SPOTlight.
#' @param SpatialDWLS Logical; if TRUE, run SpatialDWLS.
#' @param Stereoscope Logical; if TRUE, run Stereoscope.
#' @param cell2location Logical; if TRUE, run cell2location.
#' @param CARD Logical; if TRUE, run CARD.
#' @param STdeconvolve Logical; if TRUE, run STdeconvolve.
#' @param ... Additional parameters for individual methods.
#'
#' @return A list containing:
#' \describe{
#'   \item{Results.Deconv}{A list of deconvolution results from individual methods.}
#'   \item{times_methods}{A vector of execution times for each method.}
#' }
Decon_individual_methods <- function(sc_exp, sc_label, spot_exp, spot_loc,
                                       gene_det_in_min_cells_per = 0.001, expression_threshold = 1,
                                       nUMI = 10, verbose = FALSE, plot = FALSE, python_env = NULL,
                                       use_gpu = FALSE, saving_results = FALSE,
                                       SCDC = FALSE, RCTD = TRUE, MuSiC = FALSE, DeconRNASeq = FALSE,
                                       DestVI = TRUE, DWLS = TRUE, SPOTlight = TRUE, SpatialDWLS = TRUE,
                                       Stereoscope = TRUE, cell2location = TRUE, CARD = TRUE, STdeconvolve = TRUE,
                                       SCDC.iter.max = 1000, RCTD.CELL_MIN_INSTANCE = 10, MuSiC.iter.max = 1000,
                                       MuSiC.nu = 1e-04, MuSiC.eps = 0.01, DeconRNASeq.perc = 0.05,
                                       DWLS.parallel = TRUE, DWLS.is_select_DEGs = TRUE,
                                       SPOTlight.cl_n = 100, SPOTlight.hvg = 3000, SPOTlight.min_cont = 0.001,
                                       SpatialDWLS.findmarker_method = "gini", SpatialDWLS.ncp_spa = 100,
                                       SpatialDWLS.dimensions_to_use = 10, SpatialDWLS.k = 10,
                                       SpatialDWLS.resolution = 0.4, SpatialDWLS.n_iterations = 1000,
                                       SpatialDWLS.n_cell = 50, SpatialDWLS.is_select_DEGs = TRUE,
                                       Stereoscope.sc_training_plot = FALSE, Stereoscope.sc_training_save_trained_model = FALSE,
                                       Stereoscope.sc_max_epochs = 10000, Stereoscope.sc_lr = 0.01,
                                       Stereoscope.select_HVG = TRUE, Stereoscope.HVG_num = 5000,
                                       Stereoscope.st_training_plot = FALSE, Stereoscope.st_training_save_trained_model = FALSE,
                                       Stereoscope.st_max_epochs = 10000, Stereoscope.st_lr = 0.01,
                                       cell2location.sc_max_epoches = 1000, cell2location.sc_lr = 0.002,
                                       cell2location.st_N_cells_per_location = 10, cell2location.st_detection_alpha = 20.00,
                                       cell2location.st_max_epoches = 3000,
                                       CARD.minCountGene = 100, CARD.minCountSpot = 5,
                                       STdeconvolve.min.lib.size = 100, STdeconvolve.min.reads = 1,
                                       STdeconvolve.nTopOD = 1000, STdeconvolve.betaScale = 1000,
                                       DestVI.n_top_genes = 2000, DestVI.max_iter_sc = 400, DestVI.max_iter_st = 3000) {
  if(is.null(python_env)){
    cat("Please set the path to your Python environment; otherwise, a default will be used.\n")
  } else {
    reticulate::use_python(python_env, require = TRUE)
    reticulate::py_config()
  }
  
  # Preprocess the data
  database <- data_process(sc_exp, sc_label, spot_exp, spot_loc,
                           gene_det_in_min_cells_per, expression_threshold, nUMI, verbose, plot)
  
  Results.Deconv <- list()
  times_methods <- c()
  
  if(CARD) {
    cat("Running CARD...\n")
    t0 <- Sys.time()
    Results.Deconv$CARD <- CARD_run(database, minCountGene = CARD.minCountGene, minCountSpot = CARD.minCountSpot)
    times_methods["CARD"] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat("CARD runtime: ", times_methods["CARD"], " minutes\n")
  }
  if(cell2location) {
    cat("Running cell2location...\n")
    t0 <- Sys.time()
    Results.Deconv$cell2location <- cell2location_run(database, python_env = python_env,
                                                      sc_max_epoches = cell2location.sc_max_epoches,
                                                      sc_lr = cell2location.sc_lr, sc_use_gpu = use_gpu,
                                                      st_N_cells_per_location = cell2location.st_N_cells_per_location,
                                                      st_detection_alpha = cell2location.st_detection_alpha,
                                                      st_max_epoches = cell2location.st_max_epoches)
    times_methods["cell2location"] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat("cell2location runtime: ", times_methods["cell2location"], " minutes\n")
  }
  if(DeconRNASeq) {
    cat("Running DeconRNASeq...\n")
    t0 <- Sys.time()
    Results.Deconv$DeconRNASeq <- DeconRNASeq_run(database, perc = DeconRNASeq.perc)
    times_methods["DeconRNASeq"] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat("DeconRNASeq runtime: ", times_methods["DeconRNASeq"], " minutes\n")
  }
  if(DestVI) {
    cat("Running DestVI...\n")
    t0 <- Sys.time()
    Results.Deconv$DestVI <- DestVI_run(database, python_env = python_env,
                                        n_top_genes = DestVI.n_top_genes, use_gpu = use_gpu,
                                        max_iter_sc = DestVI.max_iter_sc, max_iter_st = DestVI.max_iter_st)
    times_methods["DestVI"] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat("DestVI runtime: ", times_methods["DestVI"], " minutes\n")
  }
  if(DWLS) {
    cat("Running DWLS and SVR...\n")
    t0 <- Sys.time()
    temp <- DWLS_run(database, parallel = DWLS.parallel, is_select_DEGs = DWLS.is_select_DEGs, python_env = python_env)
    Results.Deconv$DWLS <- temp$DampenedWLS
    Results.Deconv$SVR <- temp$SVR
    times_methods["DWLS"] <- as.numeric(temp$time.DWLS)
    times_methods["SVR"] <- as.numeric(temp$time.SVR)
    cat("DWLS runtime: ", times_methods["DWLS"], " minutes; SVR runtime: ", times_methods["SVR"], " minutes\n")
  }
  if(MuSiC) {
    cat("Running MuSiC...\n")
    t0 <- Sys.time()
    temp <- MuSiC_run(database, iter.max = MuSiC.iter.max, nu = MuSiC.nu, eps = MuSiC.eps)
    Results.Deconv$MuSiC_weight <- temp$Music_weight
    Results.Deconv$MuSiC_allgene <- temp$Music_allgene
    times_methods["MuSiC"] <- as.numeric(difftime(Sys.time(), t0, units = "mins")) / 2
    cat("MuSiC runtime: ", times_methods["MuSiC"], " minutes\n")
  }
  if(RCTD) {
    cat("Running RCTD...\n")
    t0 <- Sys.time()
    Results.Deconv$RCTD <- RCTD_run(database, CELL_MIN_INSTANCE = RCTD.CELL_MIN_INSTANCE)
    times_methods["RCTD"] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat("RCTD runtime: ", times_methods["RCTD"], " minutes\n")
  }
  if(SCDC) {
    cat("Running SCDC...\n")
    t0 <- Sys.time()
    Results.Deconv$SCDC <- suppressMessages(SCDC_run(database, iter.max = SCDC.iter.max))
    times_methods["SCDC"] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat("SCDC runtime: ", times_methods["SCDC"], " minutes\n")
  }
  if(SpatialDWLS) {
    cat("Running SpatialDWLS...\n")
    t0 <- Sys.time()
    Results.Deconv$SpatialDWLS <- suppressWarnings(spatialDWLS_run(database, my_python_path = python_env,
                                                                   findmarker_method = SpatialDWLS.findmarker_method,
                                                                   is_select_DEGs = SpatialDWLS.is_select_DEGs,
                                                                   ncp_spa = SpatialDWLS.ncp_spa,
                                                                   dimensions_to_use = SpatialDWLS.dimensions_to_use,
                                                                   k = SpatialDWLS.k, resolution = SpatialDWLS.resolution,
                                                                   n_iterations = SpatialDWLS.n_iterations,
                                                                   n_cell = SpatialDWLS.n_cell))
    times_methods["SpatialDWLS"] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat("SpatialDWLS runtime: ", times_methods["SpatialDWLS"], " minutes\n")
  }
  if(SPOTlight) {
    cat("Running SPOTlight...\n")
    t0 <- Sys.time()
    Results.Deconv$SPOTlight <- suppressWarnings(SPOTlight_run(database, cl_n = SPOTlight.cl_n, hvg = SPOTlight.hvg, min_cont = SPOTlight.min_cont))
    times_methods["SPOTlight"] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat("SPOTlight runtime: ", times_methods["SPOTlight"], " minutes\n")
  }
  if(Stereoscope) {
    cat("Running Stereoscope...\n")
    t0 <- Sys.time()
    Results.Deconv$Stereoscope <- Stereoscope_run(database, python_env = python_env, use_gpu = use_gpu,
                                                  sc_training_plot = Stereoscope.sc_training_plot,
                                                  sc_training_save_trained_model = Stereoscope.sc_training_save_trained_model,
                                                  sc_max_epochs = Stereoscope.sc_max_epochs, sc_lr = Stereoscope.sc_lr,
                                                  st_training_plot = Stereoscope.st_training_plot,
                                                  st_training_save_trained_model = Stereoscope.st_training_save_trained_model,
                                                  st_max_epochs = Stereoscope.st_max_epochs, st_lr = Stereoscope.st_lr,
                                                  select_HVG = Stereoscope.select_HVG, HVG_num = Stereoscope.HVG_num)
    times_methods["Stereoscope"] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat("Stereoscope runtime: ", times_methods["Stereoscope"], " minutes\n")
  }
  if(STdeconvolve) {
    cat("Running STdeconvolve...\n")
    t0 <- Sys.time()
    Results.Deconv$STdeconvolve <- STdeconvolve_run(database, min.lib.size = STdeconvolve.min.lib.size,
                                                    min.reads = STdeconvolve.min.reads, nTopOD = STdeconvolve.nTopOD,
                                                    betaScale = STdeconvolve.betaScale)
    times_methods["STdeconvolve"] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat("STdeconvolve runtime: ", times_methods["STdeconvolve"], " minutes\n")
  }
  names(times_methods) <- names(Results.Deconv)
  if(saving_results) {
    cat("Saving individual deconvolution results...\n")
    save(Results.Deconv, file = "Results_Deconv.RData")
  }
  list(Results.Deconv = Results.Deconv, times_methods = times_methods)
}


#' Ensemble Individual Deconvolution Results
#'
#' This function integrates the results from individual deconvolution methods using a weighted median approach.
#'
#' @param Results.Deconv A list of deconvolution result matrices (spots x cell types).
#' @param lambda Hyper-parameter to constrain method weights. If NULL, it is set automatically.
#' @param prob.quantile Quantile for setting lambda (default: 0.5).
#' @param niter Maximum number of iterations (default: 100).
#' @param epsilon Convergence threshold (default: 1e-5).
#'
#' @return A list containing:
#' \describe{
#'   \item{H_norm}{Ensembled (normalized) deconvolution result (spots x cell types).}
#'   \item{w}{Weights assigned to each individual method.}
#' }
solve_ensemble <- function(Results.Deconv, lambda = NULL, prob.quantile = 0.5,
                           niter = 100, epsilon = 1e-5) {
  num_methods <- length(Results.Deconv)
  num_spots <- nrow(Results.Deconv[[1]])
  num_cell_types <- ncol(Results.Deconv[[1]])
  
  # Initialize weights and ensemble result H (weighted mean)
  w <- rep(1/num_methods, num_methods)
  H <- Reduce("+", Map("*", Results.Deconv, w))
  
  if(is.null(lambda)) {
    cat("Automatically setting lambda...\n")
  }
  
  k <- 1
  while(k <= niter) {
    loss_prev <- if(k == 1) 0 else loss_all
    temp <- sapply(Results.Deconv, L1_norm, Y = H)
    if(k == 1) {
      lambda <- quantile(temp, probs = prob.quantile)
    }
    # Update weights
    w <- exp(-temp / lambda)
    w <- w / sum(w)
    # Compute weighted median for each spot and cell type
    array_results <- abind::abind(Results.Deconv, along = 3)
    H <- apply(array_results, c(1, 2), median_weighted, w = w)
    
    loss_main <- sum(sapply(Results.Deconv, L1_norm, Y = H) * w)
    loss_entropy <- sum(w * log(w))
    loss_all <- loss_main + lambda * loss_entropy
    
    cat("Iteration: ", k, " loss_main: ", loss_main, " loss_entropy: ", loss_entropy,
        " total loss: ", loss_all, " lambda: ", lambda, "\n")
    if(abs(loss_all - loss_prev) < epsilon)
      break
    k <- k + 1
  }
  
  colnames(H) <- colnames(Results.Deconv[[1]])
  H_norm <- sweep(H, 1, rowSums(H), "/")
  list(H_norm = H_norm, w = w)
}


#' Calculate L1 Norm Between Two Matrices
#'
#' @param X First matrix.
#' @param Y Second matrix.
#'
#' @return The L1 norm (sum of absolute differences) between X and Y.
L1_norm <- function(X, Y) {
  sum(abs(X - Y))
}


#' Compute Weighted Median
#'
#' This function computes the weighted median of a vector using the provided weights.
#'
#' @param x Numeric vector.
#' @param w Weights corresponding to x.
#'
#' @return The weighted median of x.
median_weighted <- function(x, w) {
  spatstat.geom::weighted.median(x, w)
}
