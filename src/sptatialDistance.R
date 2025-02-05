#' Calculate Radial Distances for Spots in a Spatial Network
#'
#' @description
#' This function computes the radial distances for given spots within a spatial network.
#' The distance is defined as the nearest neighbor distance to the border of the region.
#' Spots outside the specified region receive positive distances, whereas spots inside
#' (but not on the border) receive negative distances. Optionally, distances can be
#' converted from pixels to microns and segmented by angle.
#'
#' @param object A data.frame, matrix, or tibble with exactly 4 columns: 
#'   - **barcode**: a character identifier for each spot,
#'   - **x**: numeric x-coordinate,
#'   - **y**: numeric y-coordinate, and
#'   - **sampleID**: a numeric identifier for the sample.
#' @param spots A character vector containing the barcodes of spots to analyze.
#' @param angles Optional numeric vector of length 2 specifying the lower and upper bounds (in degrees)
#'   for filtering spots by angle relative to the sample centroid. Default is `NULL`.
#' @param angles_nbreaks Optional single numeric value indicating the number of intervals to break the
#'   angle range into. Default is `NULL`.
#' @param remove_singletons Logical. If `TRUE`, spots with zero or one neighbor (singletons) will be removed
#'   from the analysis. Default is `TRUE`.
#' @param convert_to_microns Logical. If `TRUE`, convert the calculated radial distances from pixel units to microns.
#'   Default is `FALSE`.
#' @param verbose Logical. If `TRUE`, informative messages are printed during processing. Default is `TRUE`.
#' @param ... Additional arguments passed to helper functions (e.g., `GetSpatialNetwork`, `RegionNeighbors`).
#'
#' @return If no angle filtering is specified, a numeric vector of radial distances is returned (named by barcode).
#'   If angles are provided, a data.frame is returned with the columns: `barcode`, `sampleID`, `angle`, `r_dist`,
#'   and (if applicable) angle intervals.
#'
#' @details
#' The function first validates the input object and then computes a spatial network using
#' `GetSpatialNetwork()`. It filters the network to include only the specified `spots`
#' and optionally removes isolated spots. The border of the region is identified via `RegionNeighbors()`,
#' and the nearest neighbor distances from the border are calculated for both inside and outside spots.
#' Optionally, the distances are converted from pixels to microns. If an angle range is provided,
#' each spot's angle relative to the sample centroid is computed, and the data can be split into intervals.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'spatialData' is a data.frame with the required structure and
#'   # GetSpatialNetwork is defined, compute the radial distances as follows:
#'   distances <- RadialDistance.default(spatialData, spots = c("spot1", "spot2", "spot3"))
#' }
#'
#' @export
RadialDistance <- function(object,
                                   spots,
                                   angles = NULL,
                                   angles_nbreaks = NULL,
                                   remove_singletons = TRUE,
                                   convert_to_microns = FALSE,
                                   verbose = TRUE,
                                   ...) {
  # Validate input object class, dimensions, and required column names
  valid_classes <- c("data.frame", "matrix", "tbl")
  if (!any(class(object) %in% valid_classes)) {
    abort(glue("Invalid class '{class(object)}'."))
  }
  if (ncol(object) != 4) {
    abort(glue("Invalid number of columns '{ncol(object)}'. Expected 4."))
  }
  required_cols <- c("barcode", "x", "y", "sampleID")
  if (!all(colnames(object) == required_cols)) {
    abort("Required columns are: 'barcode', 'x', 'y', and 'sampleID'")
  }
  checks <- object |>
    summarize(
      check_barcode = is.character(barcode),
      check_x = is.numeric(x),
      check_y = is.numeric(y),
      check_sampleID = is.numeric(sampleID)
    ) |>
    unlist()
  if (!all(checks)) {
    abort("Invalid column class(es).")
  }
  
  # Validate 'spots' parameter
  stopifnot(
    inherits(spots, what = "character"),
    length(spots) > 0,
    all(spots %in% object$barcode)
  )
  
  # Get the spatial network and restrict to one tissue section
  spatnet <- GetSpatialNetwork(object, ...)
  if (length(spatnet) > 1) {
    abort(glue("Default method can only handle 1 tissue section at a time, got {length(spatnet)}"))
  }
  spatnet_region <- spatnet[[1]] |>
    filter(from %in% spots, to %in% spots) |>
    group_by(from) |>
    mutate(nn = n())
  
  # Optionally remove spots with zero or one neighbor (singletons)
  if (remove_singletons) {
    spatnet_region <- spatnet_region |> filter(nn > 1)
    filtered_spots <- unique(spatnet_region$from)
    if (verbose) {
      cli_alert_info("Removing {length(spots) - length(filtered_spots)} spots with 0 neighbors.")
    }
    spots <- filtered_spots
  }
  
  if (length(spots) == 0) {
    if (verbose) cli_alert_warning("Found no spots. Returning NA values")
    return(setNames(rep(NA_real_, nrow(object)), object$barcode))
  }
  if (verbose) cli_alert_info("Extracting border spots from a region with {length(spots)} spots")
  
  # Identify border spots within the specified region
  border_spots <- RegionNeighbors(spatnet, spots = spots, outer_border = FALSE, ...)
  inside_spots <- setdiff(spots, border_spots)
  outside_spots <- setdiff(object$barcode, spots)
  
  # Get indices for the various spot groups in the object
  border_idx  <- which(object$barcode %in% border_spots)
  inside_idx  <- which(object$barcode %in% inside_spots)
  outside_idx <- which(!(object$barcode %in% spots))
  
  if (verbose) {
    cli_alert("  Detected {length(border_idx)} spots on borders")
    cli_alert("  Detected {length(inside_idx)} spots inside borders")
    cli_alert("  Detected {length(outside_idx)} spots outside borders")
  }
  
  # Calculate nearest neighbor distances from border spots for outside and inside spots
  border_coords  <- as.matrix(object[border_idx, c("x", "y")])
  outside_coords <- as.matrix(object[outside_idx, c("x", "y")])
  inside_coords  <- as.matrix(object[inside_idx, c("x", "y")])
  
  knn_outside <- kNN(x = border_coords, query = outside_coords, k = 1)
  knn_inside  <- kNN(x = border_coords, query = inside_coords, k = 1)
  
  if (verbose) cli_alert_success("Returning radial distances")
  
  # Initialize the radial distances vector
  radial_dists <- setNames(numeric(nrow(object)), object$barcode)
  # Spots outside the region receive positive distances
  radial_dists[object$barcode[outside_idx]] <- knn_outside$dist[, 1]
  # Spots inside (but not on the border) receive negative distances
  radial_dists[object$barcode[inside_idx]] <- -knn_inside$dist[, 1]
  # Border spots remain at 0
  
  # Optionally convert pixel distances to microns (requires 'dbscan' package)
  if (convert_to_microns) {
    if (!requireNamespace("dbscan", quietly = TRUE)) {
      abort(glue("Package {cli::col_br_magenta('dbscan')} is required. Please install it with:\ninstall.packages('dbscan')"))
    }
    # Compute the scale for each sample based on the minimum nearest-neighbor distance
    pixel_scales <- object |>
      group_by(sampleID) |>
      group_split() |>
      sapply(function(df) {
        min(kNN(x = df |> select(x, y), k = 1)$dist)
      })
    # Convert distances (assumes scale is provided per sample; adjust if needed)
    radial_dists <- radial_dists / (pixel_scales / 100)
  }
  
  # If angle filtering is specified, compute each spot's angle relative to its sample's centroid
  if (!is.null(angles_nbreaks)) {
    stopifnot(is.numeric(angles_nbreaks), length(angles_nbreaks) == 1)
    angles <- angles %||% c(0, 360)
  }
  if (!is.null(angles)) {
    # Validate angles and retrieve centroids via a helper function
    c(angles, centroids) %<-% .check_angles(list(spatnet_region) |> setNames("1"),
                                             object, spots, angles, angles_nbreaks, verbose)
    
    # Calculate the angle for each spot relative to its sample centroid (in degrees)
    object_by_sample <- object |>
      group_by(sampleID) |>
      group_split() |>
      setNames(nm = sapply(., function(df) unique(df$sampleID)))
    angle_df <- do.call(rbind, lapply(names(object_by_sample), function(sid) {
      df <- object_by_sample[[sid]]
      cen <- centroids[[sid]]
      df <- df |>
        mutate(diff_x = x - cen[1],
               diff_y = y - cen[2],
               angle = atan2(diff_y, diff_x) * (180 / pi),
               angle = if_else(angle <= 0, angle + 360, angle)) |>
        filter(between(angle, angles[1], angles[2])) |>
        select(-diff_x, -diff_y)
      return(df)
    }))
    
    # If angle breaks are specified, segment the angle range into intervals
    if (!is.null(angles_nbreaks)) {
      breaks_seq <- seq(angles[1], angles[2], length.out = angles_nbreaks + 1)
      angle_df <- angle_df |>
        mutate(intervals = cut(angle, breaks = breaks_seq, include.lowest = TRUE)) |>
        arrange(intervals) |>
        mutate(intervals = factor(intervals, levels = unique(intervals)))
    }
    # Add the radial distances to the data frame
    angle_df <- angle_df |>
      mutate(r_dist = radial_dists[barcode])
    radial_dists <- angle_df |>
      select(barcode, sampleID, angle, r_dist, contains("intervals"))
  }
  
  return(radial_dists)
}

# ------------------------------------------------------------------------------
# Helper function: Create a distance gradient plot
# ------------------------------------------------------------------------------
create_distance_gradient_plot <- function(coords, image, ROI_df, radius) {
  # Compute the flipped y coordinates so that the image is oriented correctly
  flipped_y <- max(coords$imagerow) - coords$imagerow + min(coords$imagerow)
  
  p <- SingleSpatialPlot(data = coords, image = image, pt.size.factor = 0) +
    geom_point(data = ROI_df,
               aes(x = imagecol, y = flipped_y, fill = distance),
               size = radius, shape = 21, color = "transparent") +
    scale_fill_gradient2(low = "darkblue", high = "red") +
    labs(fill = "Distance") +
    ggtitle("Distance Gradient Heatmap") +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  return(p)
}

# ------------------------------------------------------------------------------
# Helper function: Create a pathway activity plot
# ------------------------------------------------------------------------------
create_pathway_activity_plot <- function(test_path, palette) {
  p <- ggplot(data = test_path,
              mapping = aes(x = Distance, y = Activity, color = Pathway, group = Pathway)) +
    geom_smooth(method = "loess", se = FALSE) +
    theme_minimal() +
    theme_dr(xlength = 1,
             ylength = 1,
             arrow = arrow(length = unit(0.2, "inches"))) +
    theme(panel.grid = element_blank(),
          axis.title.y = element_text(face = "bold", hjust = 0.5, size = 15),
          axis.title.x = element_text(face = "bold", hjust = 0.5, size = 15),
          axis.title = element_text(face = "bold", hjust = 0.03, size = 15),
          axis.line = element_line(linewidth = 1),
          axis.text.x = element_text(color = "black", size = 13),
          axis.text.y = element_blank(),
          legend.title = element_text(face = "bold", size = 15),
          legend.text = element_text(size = 13)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(min(test_path$Distance), 0, max(test_path$Distance)),
                       labels = c("INNER", "0", "OUTER")) +
    scale_y_continuous(labels = NULL) +
    scale_color_manual(values = palette) +
    ggtitle("Variation of pathway activities with distance gradient") +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  return(p)
}

# ------------------------------------------------------------------------------
# Helper function: Create a gene activity plot
# ------------------------------------------------------------------------------
create_gene_activity_plot <- function(test_gene, palette) {
  p <- ggplot(data = test_gene,
              mapping = aes(x = Distance, y = Expr, color = Gene, group = Gene)) +
    geom_smooth(method = "loess", se = FALSE) +
    theme_minimal() +
    theme_dr(xlength = 1,
             ylength = 1,
             arrow = arrow(length = unit(0.2, "inches"))) +
    theme(panel.grid = element_blank(),
          axis.title.y = element_text(face = "bold", hjust = 0.5, size = 15),
          axis.title.x = element_text(face = "bold", hjust = 0.5, size = 15),
          axis.title = element_text(face = "bold", hjust = 0.03, size = 15),
          axis.line = element_line(linewidth = 1),
          axis.text.x = element_text(color = "black", size = 13),
          axis.text.y = element_blank(),
          legend.title = element_text(face = "bold", size = 15),
          legend.text = element_text(size = 13)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(min(test_gene$Distance), 0, max(test_gene$Distance)),
                       labels = c("INNER", "0", "OUTER")) +
    scale_y_continuous(labels = NULL) +
    ggtitle("Variation of specified genes with distance gradient") +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  return(p)
}

# ------------------------------------------------------------------------------
# Main function: spatialDistance
# ------------------------------------------------------------------------------
#' Generate Spatial Distance, Pathway, and Gene Activity Plots
#'
#' @param ST A Seurat object (typically from global_data$ST).
#' @param input A list or environment containing input parameters such as:
#'   - roi_cluster: Cluster ID to use for ROI assignment.
#'   - loading_ST_data: A flag indicating whether ST data is being loaded.
#'   - sd_gene: A comma-separated string of gene names to plot.
#' @param selected_spots Optional. The result of selected_spots() (if available).
#' @param ST_assay A string indicating the assay name to be used in progeny analysis.
#'
#' @return A list of ggplot objects:
#'   - gradience: Distance gradient heatmap.
#'   - pathway_plot: Smoothed plot of pathway activities versus distance.
#'   - gene_plot: (Optional) Smoothed plot of specified gene expression versus distance.
#'

SpatialDistance <- function(ST, input, selected_spots = NULL, ST_assay) {
  # --- ROI assignment ---
  # If a set of selected spots is provided, mark them as ROI ("1")
  if (!is.null(selected_spots)) {
    ROI_ID <- rownames(selected_spots)
    ST$ROI <- "0"
    ST$ROI[names(ST$seurat_clusters) %in% ROI_ID] <- "1"
  }
  # If an ROI cluster is specified in input, mark those spots as ROI ("1")
  if (!is.null(input$roi_cluster)) {
    ST$ROI <- "0"
    ST$ROI[ST$seurat_clusters %in% input$roi_cluster] <- "1"
  }
  
  # --- Update the Seurat object and compute radial distances ---
  ST <- UpdateSeuratForSemla(ST)
  coords <- GetTissueCoordinates(ST)
  image <- ST[[DefaultImage(ST)]]
  ST <- RadialDistance(ST, column_name = "ROI", selected_groups = "1",
                       remove_singletons = FALSE, verbose = FALSE)
  distance <- ST$r_dist_1
  ROI_df <- cbind(coords, distance)
  
  # --- Determine plotting point size ---
  if (!is.null(input$loading_ST_data)) {
    radius <- 4
  } else {
    radius <- diff(range(ROI_df$imagerow)) * diff(range(ROI_df$imagecol))
    radius <- radius / nrow(ROI_df)
    radius <- radius / pi
    radius <- sqrt(radius) * 0.85
  }
  
  # --- Create the distance gradient plot ---
  gradience <- create_distance_gradient_plot(coords, image, ROI_df, radius)
  
  # --- Create a clusters data frame for downstream analysis ---
  Idents(ST) <- "r_dist_1"
  SpotsClusters <- data.frame(Spot = names(Idents(ST)),
                              Distance = ST$r_dist_1,
                              stringsAsFactors = FALSE)
  
  # --- Compute progeny scores for pathway analysis ---
  ST <- progeny(ST, scale = FALSE, organism = "Human", top = 500, perm = 1,
                return_assay = TRUE, assay = ST_assay)
  ST <- Seurat::ScaleData(ST, assay = "progeny")
  progeny_scores_df <-
    as.data.frame(t(GetAssayData(ST, slot = "scale.data", assay = "progeny"))) %>%
    rownames_to_column("Spot") %>%
    tidyr::gather(Pathway, Activity, -Spot)
  progeny_scores_df <- dplyr::inner_join(progeny_scores_df, SpotsClusters)
  
  # --- Create the pathway activity plot ---
  my_palette <- distinctColorPalette(14)
  p_path <- create_pathway_activity_plot(progeny_scores_df, my_palette)
  
  # --- Optionally, create the gene activity plot if specific genes are provided ---
  plots_list <- list(gradience = gradience, pathway_plot = p_path)
  
  if (!is.null(input$sd_gene) && nchar(trimws(input$sd_gene)) > 0) {
    gene_names <- strsplit(input$sd_gene, ",\\s*")[[1]]
    dat_df <-
      as.data.frame(t(GetAssayData(ST, slot = "data", assay = ST_assay))) %>%
      rownames_to_column("Spot") %>%
      tidyr::gather(Gene, Expr, -Spot)
    
    gene_expr_df <- dplyr::inner_join(dat_df, SpotsClusters)
    test_gene <- gene_expr_df[gene_expr_df$Gene %in% gene_names, ]
    
    p_gene <- create_gene_activity_plot(test_gene, my_palette)
    plots_list$gene_plot <- p_gene
  }
  
  return(plots_list)
}

