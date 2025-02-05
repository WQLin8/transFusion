#' Run Colocalization Analysis
#'
#' This function performs a colocalization analysis using the "misty" framework.
#' It first creates an initial view from the deconvolution result, then (optionally)
#' adds a paraview using spatial positions, runs the misty algorithm, and collects results.
#'
#' @param decon_result A matrix or data.frame of deconvolution results with rows as spots and columns as cells or genes.
#'                     (Column names should be free of spaces; they will be converted to underscores.)
#' @param position A matrix or data.frame of spatial positions. The spot IDs must match those in decon_result.
#' @param outdir A character string specifying the output folder.
#' @param para_l Optional numeric value for the number of neighbors to use in the paraview. Default is NULL.
#'
#' @return A list of colocalization results (output of collect_results()).
#' @export
run_colocalization <- function(decon_result, position, outdir, para_l = NULL) {
  # Replace any spaces in column names with underscores
  colnames(decon_result) <- gsub(" ", "_", colnames(decon_result))
  
  # Create the initial view from the deconvolution result
  misty_intra <- create_initial_view(decon_result)
  
  # If the neighbor parameter is provided, add a paraview using the provided positions
  if (!is.null(para_l)) {
    misty_views <- misty_intra %>% add_paraview(position, l = para_l)
  } else {
    misty_views <- misty_intra
  }
  
  # Create the output directory if it does not exist
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    print(paste("outdir:", outdir, "is created"))
  }
  
  # Run the misty algorithm and collect the results
  results_folder <- run_misty(misty_views, results.folder = outdir)
  misty_results <- collect_results(results_folder)
  return(misty_results)
}


#' Get Default Image from an Object
#'
#' This function updates an object, retrieves the associated images based on the default assay,
#' and returns the first image.
#'
#' @param object An object that contains image data.
#'
#' @return The first image from the object.
#' @export
DefaultImage <- function(object) {
  object <- UpdateSlots(object = object)
  images <- Images(object = object, assay = DefaultAssay(object = object))
  if (length(images) < 1) {
    images <- Images(object = object)
  }
  return(images[[1]])
}


#' Downsample Seurat Clusters
#'
#' This function downsamples cells from each cluster in a Seurat object. For clusters with more cells
#' than the specified number, it randomly selects a subset of cells.
#'
#' @param seurat_obj A Seurat object.
#' @param cells_per_cluster An integer specifying the maximum number of cells to retain per cluster (default: 200).
#' @param seed An integer used as the random seed for reproducibility (default: 666).
#'
#' @return A downsampled Seurat object.
#' @export
downsample_seurat_clusters <- function(seurat_obj, cells_per_cluster = 200, seed = 666) {
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Get cluster identities
  clusters <- Idents(seurat_obj)
  
  # Initialize a vector to store cell names to keep
  cells_to_keep <- c()
  
  # Loop through each unique cluster and sample cells accordingly
  for (cluster in unique(clusters)) {
    cluster_cells <- names(clusters[clusters == cluster])
    if (length(cluster_cells) <= cells_per_cluster) {
      cells_to_keep <- c(cells_to_keep, cluster_cells)
    } else {
      cells_to_keep <- c(cells_to_keep, sample(cluster_cells, cells_per_cluster))
    }
  }
  
  # Subset the Seurat object with the selected cells
  downsampled_seurat <- subset(seurat_obj, cells = cells_to_keep)
  
  # Print summary information
  cat(paste("Original object cell count:", ncol(seurat_obj), "\n"))
  cat(paste("Downsampled object cell count:", ncol(downsampled_seurat), "\n"))
  cat("Cell counts per cluster after downsampling:\n")
  print(table(Idents(downsampled_seurat)))
  
  return(downsampled_seurat)
}


#' Plot Interaction Heatmap
#'
#' This function generates a heatmap that visualizes cell–cell interaction strengths based on
#' aggregated importance and improvement statistics from a misty analysis.
#'
#' @param misty.results A list containing the misty results (must include "importances.aggregated" and "improvements.stats").
#' @param view A character string specifying the view to use (must match a value in the 'view' column of importances.aggregated).
#' @param cutoff A numeric value that serves as the midpoint for the fill gradient (default: 1).
#' @param trim A numeric threshold for filtering targets (default: -Inf).
#' @param trim.measure A character vector indicating the metric to use for trimming. Options include "gain.R2", "multi.R2", "intra.R2",
#'                     "gain.RMSE", "multi.RMSE", "intra.RMSE" (default: "gain.R2").
#' @param clean Logical; if TRUE, only predictors and targets with a total importance above cutoff are kept.
#' @param title A character string for the plot title (default: "Cell-cell dependency").
#'
#' @return A ggplot object displaying the interaction heatmap.
#' @export
interaction_heatmap <- function(misty.results, view, cutoff = 1,
                                trim = -Inf, trim.measure = c(
                                  "gain.R2", "multi.R2", "intra.R2",
                                  "gain.RMSE", "multi.RMSE", "intra.RMSE"
                                ),
                                clean = FALSE,
                                title = "Cell-cell dependency") {
  trim.measure.type <- match.arg(trim.measure)
  
  # Check that necessary fields are present in the results list
  assertthat::assert_that(("importances.aggregated" %in% names(misty.results)),
                          msg = "The provided result list is malformed. Consider using collect_results().")
  assertthat::assert_that(("improvements.stats" %in% names(misty.results)),
                          msg = "The provided result list is malformed. Consider using collect_results().")
  assertthat::assert_that((view %in%
                             (misty.results$importances.aggregated %>% dplyr::pull(view))),
                          msg = "The selected view cannot be found in the results table.")
  
  # Determine the sign for filtering based on the chosen measure
  inv <- sign((stringr::str_detect(trim.measure.type, "gain") |
                  stringr::str_detect(trim.measure.type, "RMSE", negate = TRUE)) - 0.5)
  
  # Filter targets based on improvement statistics and the trim value
  targets <- misty.results$improvements.stats %>%
    dplyr::filter(
      measure == trim.measure.type,
      inv * mean >= inv * trim
    ) %>%
    dplyr::pull(target)
  
  # Filter the aggregated importances for the selected view and targets
  plot.data <- misty.results$importances.aggregated %>%
    dplyr::filter(view == !!view, Target %in% targets)
  
  # Optionally clean the data by removing predictors/targets with low total importance
  if (clean) {
    clean.predictors <- plot.data %>%
      dplyr::mutate(Importance = Importance * (Importance >= cutoff)) %>%
      dplyr::group_by(Predictor) %>%
      dplyr::summarize(total = sum(Importance, na.rm = TRUE)) %>%
      dplyr::filter(total > 0) %>%
      dplyr::pull(Predictor)
    clean.targets <- plot.data %>%
      dplyr::mutate(Importance = Importance * (Importance >= cutoff)) %>%
      dplyr::group_by(Target) %>%
      dplyr::summarize(total = sum(Importance, na.rm = TRUE)) %>%
      dplyr::filter(total > 0) %>%
      dplyr::pull(Target)
    plot.data.clean <- plot.data %>%
      dplyr::filter(Predictor %in% clean.predictors, Target %in% clean.targets)
  } else {
    plot.data.clean <- plot.data
  }
  
  # Define a custom blue color
  set2.blue <- "#008B00"
  plot.data.clean <- plot.data.clean %>% 
    dplyr::mutate(Importance = ifelse(Predictor == Target, NA, Importance))
  
  # Create the heatmap plot using ggplot2
  results.plot <- ggplot2::ggplot(
    plot.data.clean,
    ggplot2::aes(x = Predictor, y = Target)
  ) +
    ggplot2::geom_tile(ggplot2::aes(fill = Importance), colour = "grey", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(
      low = "white",
      mid = "white",
      high = set2.blue,
      midpoint = cutoff,
      limits = c(0, max(plot.data.clean$Importance, na.rm = TRUE)),
      oob = scales::squish,
      na.value = "white",
      breaks = c(0, max(plot.data.clean$Importance, na.rm = TRUE)),
      labels = c(">0", "max")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(face = "bold", size = 14),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 14, face = "bold", margin = margin(b = 15, unit = "pt")),
      legend.text = element_text(size = 12),
      legend.spacing.y = unit(10, "pt"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
    ) + 
    ggplot2::coord_equal() +
    ggplot2::ggtitle(title)
  
  print(results.plot)
}



#' Check if the Last Letter is Uppercase
#'
#' This function removes all non-letter characters from a string and checks if the last remaining character is an uppercase letter.
#'
#' @param input_string A character string.
#'
#' @return A logical value indicating whether the last letter (after removing non-letters) is uppercase.
#' @export
check_last_uppercase <- function(input_string) {
  # Remove all non-letter characters
  letters_only <- gsub("[^a-zA-Z]", "", input_string)
  
  # Return FALSE if no letters remain
  if (nchar(letters_only) == 0) {
    return(FALSE)
  }
  
  # Extract the last character and check if it is uppercase
  last_char <- substr(letters_only, nchar(letters_only), nchar(letters_only))
  return(last_char %in% LETTERS)
}


#' jjVolcano: Create a Volcano Plot for Differential Expression
#'
#' This function generates a volcano plot for differential expression data.
#' It filters the data based on fold-change and p-value cutoffs, assigns significance types,
#' selects top genes per cluster, and overlays labels on the plot.
#'
#' @param diffData A data.frame of differential expression results (must include columns like avg_log2FC, p_val, p_val_adj, cluster, and gene).
#' @param myMarkers Optional vector of gene names to display instead of automatically selected top genes.
#' @param order.by A character vector specifying the column by which to order the genes (default: "avg_log2FC").
#' @param log2FC.cutoff Numeric; fold-change cutoff (default: 0.25).
#' @param pvalue.cutoff Numeric; p-value cutoff (default: 0.05).
#' @param adjustP.cutoff Numeric; adjusted p-value cutoff (default: 0.01).
#' @param topGeneN Integer; number of top genes per cluster to label (default: 5).
#' @param col.type Character; determines how points are colored ("updown" or "adjustP", default: "updown").
#' @param back.col Background color for the tiles (default: "grey93").
#' @param pSize Numeric; point size for jittered points (default: 0.75).
#' @param aesCol A vector of colors used for significant up/down regulation (default: c("#0099CC", "#CC3333")).
#' @param legend.position A numeric vector specifying legend position (default: c(0.7, 0.9)).
#' @param base_size Base font size for the plot (default: 14).
#' @param tile.col A vector of colors used for the tile background (default obtained from jjAnno::useMyCol("paired", n = 9)).
#' @param cluster.order Optional vector specifying the order of clusters.
#' @param polar Logical; if TRUE, use polar coordinates for the plot (default: FALSE).
#' @param expand A numeric vector for expansion in the y-axis scale (default: c(-1, 1)).
#' @param flip Logical; if TRUE, flip the coordinate system (default: FALSE).
#' @param height Numeric; height of the tile (default: 0.8).
#' @param celltypeSize Numeric; text size for cluster labels (default: 3).
#' @param ... Additional parameters passed to ggrepel::geom_text_repel.
#'
#' @return A ggplot object representing the volcano plot.
#' @export
jjVolcano <- function(
    diffData = NULL,
    myMarkers = NULL,
    order.by = c("avg_log2FC"), # Options: "avg_log2FC" or "p_val"
    log2FC.cutoff = 0.25,
    pvalue.cutoff = 0.05,
    adjustP.cutoff = 0.01,
    topGeneN = 5,
    col.type = "updown",
    back.col = "grey93",
    pSize = 0.75,
    aesCol = c("#0099CC", "#CC3333"),
    legend.position = c(0.7, 0.9),
    base_size = 14,
    tile.col = jjAnno::useMyCol("paired", n = 9),
    cluster.order = NULL,
    polar = FALSE,
    expand = c(-1, 1),
    flip = FALSE,
    height = 0.8,
    celltypeSize = 3,
    ...) {
  
  # Filter data based on fold-change and p-value cutoffs
  diff.marker <- diffData %>%
    dplyr::filter(abs(avg_log2FC) >= log2FC.cutoff & p_val < pvalue.cutoff)
  
  # Assign significance type based on avg_log2FC and adjusted p-value
  diff.marker <- diff.marker %>%
    dplyr::mutate(type = ifelse(avg_log2FC >= log2FC.cutoff, "sigUp", "sigDown")) %>%
    dplyr::mutate(type2 = ifelse(p_val_adj < adjustP.cutoff,
                                 paste("adjust Pvalue < ", adjustP.cutoff, sep = ""),
                                 paste("adjust Pvalue >= ", adjustP.cutoff, sep = "")
    ))
  
  # Reorder clusters if a specific order is provided
  if (!is.null(cluster.order)) {
    diff.marker$cluster <- factor(diff.marker$cluster, levels = cluster.order)
  }
  
  # Get background data per cluster for tile plotting
  back.data <- purrr::map_df(unique(diff.marker$cluster), function(x) {
    tmp <- diff.marker %>% dplyr::filter(cluster == x)
    data.frame(
      cluster = x,
      min = min(tmp$avg_log2FC) - 0.2,
      max = max(tmp$avg_log2FC) + 0.2
    )
  })
  
  # Select top genes per cluster based on the ordering metric
  top.marker.tmp <- diff.marker %>% dplyr::group_by(cluster)
  top.marker.max <- top.marker.tmp %>% dplyr::slice_max(n = topGeneN, order_by = get(order.by))
  top.marker.min <- top.marker.tmp %>% dplyr::slice_min(n = topGeneN, order_by = get(order.by))
  
  top.marker <- rbind(top.marker.max, top.marker.min)
  
  # If custom markers are provided, use them instead
  if (!is.null(myMarkers)) {
    top.marker <- diff.marker %>% dplyr::filter(gene %in% myMarkers)
  }
  
  # Start constructing the volcano plot
  p1 <- ggplot2::ggplot(diff.marker, ggplot2::aes(x = cluster, y = avg_log2FC)) +
    ggplot2::geom_col(data = back.data, ggplot2::aes(x = cluster, y = min), fill = back.col) +
    ggplot2::geom_col(data = back.data, ggplot2::aes(x = cluster, y = max), fill = back.col)
  
  # Add jittered points colored by significance type
  if (col.type == "updown") {
    p2 <- p1 +
      ggplot2::geom_jitter(ggplot2::aes(color = type), size = pSize) +
      ggplot2::scale_color_manual(values = c("sigDown" = aesCol[1], "sigUp" = aesCol[2]))
  } else if (col.type == "adjustP") {
    p2 <- p1 +
      ggplot2::geom_jitter(ggplot2::aes(color = type2), size = pSize) +
      ggplot2::scale_color_manual(values = c(aesCol[2], aesCol[1]))
  }
  
  # Apply additional theme settings
  p3 <- p2 +
    ggplot2::scale_y_continuous(n.breaks = 6) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = legend.position,
                   legend.title = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank()) +
    ggplot2::xlab("Clusters") + ggplot2::ylab("Average log2FoldChange") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))
  
  # Add a tile layer and gene labels using ggrepel
  p4 <- p3 +
    ggplot2::geom_tile(ggplot2::aes(x = cluster, y = 0, fill = cluster),
                       color = "black",
                       height = height,
                       alpha = 0.3,
                       show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = tile.col) +
    ggrepel::geom_text_repel(data = top.marker,
                             ggplot2::aes(x = cluster, y = avg_log2FC, label = gene),
                             max.overlaps = 50,
                             ...)
  
  # Optionally use polar coordinates or flip the plot
  if (polar == TRUE) {
    p5 <- p4 +
      geomtextpath::geom_textpath(ggplot2::aes(x = cluster, y = 0, label = cluster)) +
      ggplot2::scale_y_continuous(n.breaks = 6, expand = ggplot2::expansion(mult = expand)) +
      ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(legend.position = legend.position, legend.title = ggplot2::element_blank()) +
      ggplot2::coord_polar(clip = "off", theta = "x")
  } else {
    if (flip == TRUE) {
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_label(ggplot2::aes(x = cluster, y = 0, label = cluster), size = celltypeSize) +
        ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank()) +
        ggplot2::coord_flip()
    } else {
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_text(ggplot2::aes(x = cluster, y = 0, label = cluster), size = celltypeSize) +
        ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    }
  }
  return(p5)
}
