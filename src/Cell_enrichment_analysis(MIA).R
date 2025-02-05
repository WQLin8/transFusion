#' Run Marker Intersection Analysis (MIA)
#'
#' This function performs Marker Intersection Analysis (MIA) by comparing marker genes
#' identified from a scRNA-seq dataset and a spatial transcriptomics (ST) dataset. It computes
#' an enrichment score for each pair of clusters (from SC and ST) based on the overlap of marker genes.
#'
#' @param SC A Seurat object representing the scRNA-seq data.
#' @param ST A Seurat object representing the spatial transcriptomics data.
#' @param SC_assay A character string specifying the assay name to use in the SC object.
#' @param ST_assay A character string specifying the assay name to use in the ST object.
#' @param sc_min_pct Minimum percentage of cells in which a gene is expressed for SC marker detection (default: 0.25).
#' @param sc_logfc_threshold Log fold-change threshold for SC marker detection (default: 1).
#' @param st_min_pct Minimum percentage of cells in which a gene is expressed for ST marker detection (default: 0.1).
#' @param st_logfc_threshold Log fold-change threshold for ST marker detection (default: 0.1).
#'
#' @return A matrix of MIA enrichment scores with rows corresponding to SC clusters and columns corresponding to ST clusters.
#'         Positive values indicate enrichment, and negative values indicate depletion.
#'
#' @examples
#' \dontrun{
#'   mia_matrix <- run_MIA(SC = scRNA_object, ST = stRNA_object,
#'                         SC_assay = "RNA", ST_assay = "Spatial",
#'                         sc_min_pct = 0.25, sc_logfc_threshold = 1,
#'                         st_min_pct = 0.1, st_logfc_threshold = 0.1)
#' }
#'
#' @export
run_MIA <- function(SC, ST, SC_assay, ST_assay,
                    sc_min_pct = 0.25, sc_logfc_threshold = 1,
                    st_min_pct = 0.1, st_logfc_threshold = 0.1) {
  
  # Find markers in the scRNA-seq data using the specified assay and thresholds
  sc.markers <- FindAllMarkers(SC, assay = SC_assay, only.pos = TRUE,
                               test.use = "wilcox", min.pct = sc_min_pct, 
                               logfc.threshold = sc_logfc_threshold)
  
  # Find markers in the spatial transcriptomics data using the specified assay and thresholds
  st.markers <- FindAllMarkers(ST, assay = ST_assay, only.pos = TRUE,
                               test.use = "wilcox", min.pct = st_min_pct,
                               logfc.threshold = st_logfc_threshold, verbose = FALSE)
  
  # Extract cluster identities for both datasets
  st.clusts <- levels(Idents(ST))
  sc.clusts <- levels(Idents(SC))
  
  # Create marker lists for each ST cluster
  st.marker.list <- lapply(st.clusts, function(clust) {
    st.markers[st.markers$cluster == clust, "gene"]
  })
  names(st.marker.list) <- st.clusts
  
  # Create marker lists for each SC cluster
  sc.marker.list <- lapply(sc.clusts, function(clust) {
    sc.markers[sc.markers$cluster == clust, "gene"]
  })
  names(sc.marker.list) <- sc.clusts
  
  # Initialize the MIA results matrix (rows: SC clusters; columns: ST clusters)
  M <- length(sc.clusts)
  N <- length(st.clusts)
  MIA.results <- matrix(0, nrow = M, ncol = N)
  rownames(MIA.results) <- sc.clusts
  colnames(MIA.results) <- st.clusts
  
  # Define the gene universe as the number of common genes in the specified assays
  gene.universe <- length(intersect(rownames(ST[[ST_assay]]), rownames(SC[[SC_assay]])))
  
  # Loop over each pair of ST and SC clusters and compute the enrichment score
  for (i in seq_len(N)) {
    for (j in seq_len(M)) {
      genes_st <- st.marker.list[[st.clusts[i]]]
      genes_sc <- sc.marker.list[[sc.clusts[j]]]
      
      # A: overlap; B: total markers in ST cluster; C: total markers in SC cluster
      A <- length(intersect(genes_st, genes_sc))
      B <- length(genes_st)
      C <- length(genes_sc)
      
      # Calculate enrichment and depletion scores using the hypergeometric test
      enr <- -log10(phyper(A, B, gene.universe - B, C, lower.tail = FALSE))
      dep <- -log10(1 - phyper(A, B, gene.universe - B, C, lower.tail = FALSE))
      
      # Assign the score based on which value is larger; negative for depletion
      if (enr < dep) {
        MIA.results[j, i] <- -dep
      } else {
        MIA.results[j, i] <- enr
      }
    }
  }
  
  # Replace any infinite values with 0
  MIA.results[is.infinite(MIA.results)] <- 0
  
  return(MIA.results)
}
