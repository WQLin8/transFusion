#' Plot Ligand-Receptor Interaction Map (plotLR)
#'
#' This function visualizes ligand–receptor interactions on spatial transcriptomics data.
#' It computes neighborhood expression for the ligand and receptor, calculates an interaction score,
#' and returns a combined plot including a spatial group plot and a feature plot.
#'
#' @param ST A Seurat object representing the spatial transcriptomics data.
#' @param normalize Logical; if TRUE, the expression data will be normalized (default: TRUE).
#' @param knn Number of nearest neighbors to consider (default: 8).
#' @param LRpair A character vector of length 2 containing the ligand and receptor gene names.
#' @param pt.size Numeric; point size for plotting (default: 2).
#' @param alpha.min Minimum alpha value for scaling the interaction score (default: 0.1).
#' @param max.cut Quantile cutoff (default: 0.95) used to cap the interaction score.
#'
#' @return A combined ggplot object showing the spatial dimension plot and feature plot.
#' @export
plotLR <- function(ST, normalize = TRUE, knn = 8, LRpair, pt.size = 2, alpha.min = 0.1, max.cut = 0.95) {
  # Retrieve expression matrix and spatial coordinates
  expr <- GetAssayData(object = ST)
  location <- GetTissueCoordinates(ST, cols = c("row", "col"), scale = NULL)
  location <- location[colnames(expr), ]
  
  topn <- floor(0.2 * dim(location)[1])
  if (sum(rownames(expr) %in% LRpair) != 2) { 
    stop("Ligand or receptor are not expressed")
  }
  
  # Compute the nearest neighbor indices using RANN
  nnmatrix <- RANN::nn2(location, k = knn)$nn.idx
  countsum <- Matrix::colSums(expr)
  ncell <- dim(expr)[2]
  
  # Normalize expression data if required
  if (normalize == TRUE) {
    expr <- Matrix::t(log(Matrix::t(expr) / countsum * median(countsum) + 1))
  }
  
  # Extract ligand and receptor expression
  ligand <- expr[LRpair[1], ]
  receptor <- expr[LRpair[2], ]
  LRexp <- rbind(ligand, receptor)
  
  # Compute neighborhood maximum expression
  neighexp <- apply(nnmatrix, 1, function(x) { apply(LRexp[, x[2:knn]], 1, max) })
  
  # Compute interaction score
  LRadd <- pmax(LRexp[1,] * neighexp[2,], LRexp[2,] * neighexp[1,])
  LRadd_max <- quantile(LRadd, probs = max.cut)
  LRadd[LRadd > LRadd_max] <- LRadd_max
  
  # Determine top cells based on ligand and receptor expression
  if (sum(ligand > 0) > topn) { 
    n1 <- order(ligand, sample(ncell, ncell), decreasing = TRUE)[1:topn]
  } else {
    n1 <- which(ligand > 0)
  }
  if (sum(receptor > 0) > topn) { 
    n2 <- order(receptor, sample(ncell, ncell), decreasing = TRUE)[1:topn]
  } else {
    n2 <- which(receptor > 0)
  }
  
  expcol <- rep(0, ncell)
  expcol[n1] <- 1
  expcol[n2] <- 2
  expcol[intersect(n1, n2)] <- 3
  
  tmp <- data.frame(x = location[, 1], y = location[, 2], Exp = as.factor(expcol))
  tmpLRadd <- data.frame(x = location[, 1], y = location[, 2], LR = LRadd)
  
  # Scale alpha based on interaction score
  alpha <- (LRadd - min(LRadd)) / (max(LRadd) - min(LRadd)) * (1 - alpha.min) + alpha.min
  
  # Add group information to ST object
  ST$LR_Exp <- tmp$Exp
  ST$LR_Exp <- factor(ST$LR_Exp, levels = c("0", "1", "2", "3"),
                      labels = c("Both Low", "Ligand High", "Receptor High", "Both High"))
  Idents(ST) <- "seurat_clusters"
  color_vector <- c("Both Low" = "grey", "Ligand High" = "red", "Receptor High" = "green", "Both High" = "blue")
  
  # Generate spatial dimension plot (p1)
  p1 <- SpatialDimPlot(ST, label.size = 13, cols = color_vector, group.by = "LR_Exp", pt.size.factor = 2) + 
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.key.size = unit(1, 'cm')) +
    guides(fill = guide_legend(override.aes = list(size = 7)))
  
  # Generate spatial feature plot (p2)
  ST$LR_score <- tmpLRadd$LR
  p2 <- SpatialFeaturePlot(ST, features = "LR_score", pt.size.factor = 2) +
    scale_fill_gradient(low = "white", high = "red",
                        labels = c("low", "high"),
                        breaks = c(min(ST$LR_score), max(ST$LR_score)),
                        guide = guide_colourbar(title = "Prob",
                                                title.position = 'top',
                                                title.hjust = 0,
                                                label.position = "right",
                                                title.theme = element_text(size = 16),
                                                label.theme = element_text(size = 14),
                                                barwidth = 1.2,
                                                barheight = 5)) +
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1, 'cm'),
          legend.position = "right")
  
  LRpair_title <- paste(LRpair, collapse = "-")
  
  # Combine the two plots and add a title
  combined_plot <- p1 + p2 + 
    plot_annotation(title = LRpair_title) +
    theme(plot.title = element_text(hjust = 0.5, size = 22))
  combined_plot
}


#' Plot All Ligand-Receptor Pairs
#'
#' This function takes a comma-separated string of ligand and receptor gene names, splits them into pairs,
#' and generates a plot for each pair using the plotLR function.
#'
#' @param ST A Seurat object containing spatial transcriptomics data.
#' @param LR_pairs A character string of ligand and receptor names separated by commas (e.g., "LIG1, REC1, LIG2, REC2").
#' @param normalize Logical; if TRUE, normalize the expression data (default: TRUE).
#' @param alpha.min Minimum alpha value for the plot (default: 0.1).
#' @param pt.size Numeric; point size for the plots (default: 1).
#'
#' @return A list of ggplot objects for each ligand-receptor pair.
#' @export
plot_all_LR_pairs <- function(ST, LR_pairs, normalize = TRUE, alpha.min = 0.1, pt.size = 1) {
  # Split the input string into a vector of gene names
  LR_pairs <- strsplit(LR_pairs, ",\\s*")[[1]]
  
  # Ensure that the number of elements is even (each pair has a ligand and receptor)
  if (length(LR_pairs) %% 2 != 0) {
    stop("Ligand_receptor must contain an even number of elements.")
  }
  
  # Loop through each ligand-receptor pair and generate a plot
  plots <- list()
  for (i in seq(1, length(LR_pairs), by = 2)) {
    ligand <- LR_pairs[i]
    receptor <- LR_pairs[i + 1]
    
    p <- plotLR(ST, LRpair = c(ligand, receptor), normalize = normalize, alpha.min = alpha.min, pt.size = pt.size)
    
    plots[[paste(ligand, receptor, sep = "-")]] <- p
  }
  
  return(plots)
}