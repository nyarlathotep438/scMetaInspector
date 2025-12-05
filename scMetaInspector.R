# General Data Statistic####
# Check the number of NA and Empty in the meta.table
Check_NA_Empty = function(scdata){
  #require package
  require(readr)
  
  #Get meta data from Seurat object
  meta.table <- scdata@meta.data
  
  #create direction
  if (!dir.exists("./report/NAreport")) dir.create("./report/NAreport", recursive = TRUE)
  sample_name <- deparse(substitute(scdata))
  report_dir <- paste0("./report/NAreport/", sample_name, "_NAreport.csv")
  
  #Count NA and empty value
  na_counts <- colSums(is.na(meta.table))
  empty_counts <- apply(meta.table, 2, function(x){
    x_char = as.character(x)
    sum(x_char == "", na.rm = TRUE)
  })
  
  total_missing <- na_counts + empty_counts
  
  all_na <- na_counts == nrow(meta.table)
  all_empty <- empty_counts == (nrow(meta.table) - na_counts)
  
  #Make final result table for report the NA and empty and save it as CSV file
  result_df <- data.frame(
    Feature = names(na_counts),
    NA_count = na_counts,
    Empty_count = empty_counts,
    Total_missing = total_missing,
    All_NA = all_na,
    All_Empty = all_empty,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  write_csv(result_df, report_dir)
  return(result_df)
}

# Make a Frequency count(One Column Selection)
Frequency_Count = function(scdata, column, sample_name = NULL) {
  #require package
  require(dplyr)
  require(readr)
  
  # Safe column
  safe_column_name <- function(name) {
    gsub("[^[:alnum:]_]", "_", name)
  }
  
  # Get metadata
  meta.table <- scdata@meta.data
  
  # Get column names (raw)
  column_name <- if (is.character(column)) {
    column
  } else {
    deparse(substitute(column))
  }
  
  # Verify that the column exists
  if (!column_name %in% colnames(meta.table)) {
    stop("Column '", column_name, "' not found in metadata")
  }
  
  # Determine the sample name (the passed-in name is used first)
  if (is.null(sample_name)) {
    sample_name <- deparse(substitute(scdata))
  }
  
  # Create a directory (using a safe name)
  safe_sample_name <- safe_column_name(sample_name)
  sample_dir <- file.path("./report", "Frequencyreport", safe_sample_name)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Count using original column name
  counts_table <- meta.table %>%
    dplyr::count(dplyr::across(dplyr::all_of(column_name)), name = "Counts") %>%
    dplyr::mutate(Percentage = round(Counts / sum(Counts) * 100, 2)) %>%
    dplyr::arrange(dplyr::desc(Counts)) %>%
    dplyr::rename(!!column_name := dplyr::all_of(column_name))
  
  # Generate safe file names
  safe_name <- safe_column_name(column_name)
  table_path <- file.path(sample_dir, paste0(safe_name, "_freq.csv"))
  
  # Save CSV
  readr::write_csv(counts_table, table_path)
  
  cat("Frequency table saved to:", table_path, "\n")
  return(counts_table)
}

# Make a Frequency count table list(Manual Selection)
Frequency_Count_List = function(scdata, columns) {
  # Get the name of the object
  sample_name <- deparse(substitute(scdata))
  
  result_list <- lapply(columns, function(col) {
    # Pass the sample name to Frequency_Count
    Frequency_Count(scdata, column = col, sample_name = sample_name)
  })
  
  names(result_list) <- columns
  return(result_list)
}

# Make Frequency count table(Automatic Selection)
Frequency_Count_All = function(scdata, max_unique = 50) {
  # Safe Column
  safe_column_name <- function(name) {
    gsub("[^[:alnum:]_]", "_", name)
  }
  
  # Get the name of the passed object
  sample_name <- deparse(substitute(scdata))
  safe_sample_name <- safe_column_name(sample_name)
  
  # Create Directory Path
  sample_dir <- file.path("./report", "Frequencyreport", safe_sample_name)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  meta.table <- scdata@meta.data
  all_cols <- colnames(meta.table)
  
  result_list <- list()
  skipped_cols <- character(0)
  unique_list <- list()
  
  for (col in all_cols) {
    unique_values <- unique(meta.table[[col]])
    num_unique <- length(unique_values)
    
    if (num_unique > max_unique) {
      skipped_cols <- c(skipped_cols, col)
      next
    }
    
    if (num_unique == 1) {
      unique_list[[col]] <- data.frame(
        Name = safe_column_name(col),
        Value = as.character(unique_values),
        stringsAsFactors = FALSE
      )
    } else {
      tryCatch({
        # Pass the sample name to Frequency_Count
        counts_table <- Frequency_Count(scdata, col, sample_name = sample_name)
        result_list[[safe_column_name(col)]] <- counts_table
      }, error = function(e) {
        warning("Error processing column '", col, "': ", e$message)
        skipped_cols <<- c(skipped_cols, col)
      })
    }
  }
  
  # Convert a list of unique values
  unique_values_df <- if (length(unique_list) > 0) {
    do.call(rbind, unique_list)
  } else {
    data.frame(Name = character(), Value = character(), stringsAsFactors = FALSE)
  }
  
  # Save unique values data to CSV file
  if (nrow(unique_values_df) > 0) {
    unique_file_path <- file.path(sample_dir, "unique_value_summary.csv")
    readr::write_csv(unique_values_df, unique_file_path)
    cat("Unique value summary saved to:", unique_file_path, "\n")
  }
  
  return(list(
    frequency_tables = result_list,
    skipped_columns = skipped_cols,
    unique_value_summary = unique_values_df
  ))
}

# Cluster####
# Preparation before clustering
Preparation_folders <- function() {
  # Define the list of folders that need to be created
  folders <- c("plot", "report", "marker_gene")
  
  # Creat folders
  for (folder in folders) {
    if (!dir.exists(folder)) {
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      message("Folder created successfully: ", folder)
    } else {
      message("Folder already exists: ", folder)
    }
  }
  
  # Returns the creation result summary
  results <- sapply(folders, dir.exists)
  return(invisible(results))
}

# Run PCA and Plot Elbow to find suitable dimensions
RunPCA_PlotElbow <- function(scdata, 
                            sample_name = "Sample",
                            plot_dir = "./plot/",
                            ndims = 50,
                            width = 6, 
                            height = 6) {
  # Run PCA analysis
  scdata <- Seurat::RunPCA(scdata)
  
  # Create the output directory if it does not exist
  sample_dir <- file.path(plot_dir, sample_name)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    message("Created directory: ", sample_dir)
  }
  
  # Create a file path
  file_path <- file.path(sample_dir, paste0(sample_name, "_elbow.tiff"))
  
  # Generate and save the elbow plot
  elbow_plot <- Seurat::ElbowPlot(scdata, ndims = ndims)
  ggplot2::ggsave(
    filename = file_path,
    plot = elbow_plot,
    device = "tiff",
    width = width,
    height = height
  )
  
  message("Elbow plot saved to: ", file_path)
  return(scdata)
}

# Run cluster and umap repeat to find the best resolution
Try_different_resolutions <- function(scdata,
                                      resolutions,
                                      dims = 1:15,
                                      min.dist = 0.5,
                                      spread = 1.5,
                                      plot_dir = "./plot/") {
  # Load necessary packages
  require(Seurat)
  require(ggplot2)
  
  # Obtain the variable name of the passed object
  dataset_name <- deparse(substitute(scdata))
  
  # Clean special characters in dataset names
  clean_name <- gsub("[[:space:][:punct:]]", "_", dataset_name)
  
  # Create a directory structure: plot_dir/dataset_name/UMAP/
  dataset_dir <- file.path(plot_dir, clean_name)
  umap_dir <- file.path(dataset_dir, "UMAP")
  
  # Safely create directories
  if (!dir.exists(umap_dir)) {
    dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
    message("Creat PATH: ", umap_dir)
  }
  
  # Seurat analysis workflow
  scdata <- FindNeighbors(object = scdata, dims = dims)
  scdata <- RunUMAP(scdata, dims = dims, min.dist = min.dist, spread = spread)
  
  # Looping through different resolutions
  for (res in resolutions) {
    cluster_name <- paste0("TEST_res_", res)
    scdata <- FindClusters(object = scdata, resolution = res, cluster.name = cluster_name)
    
    # Creating UMAP Plots
    umap_plot <- DimPlot(scdata, reduction = "umap", label = TRUE, group.by = cluster_name) + 
      ggtitle(paste0("UMAP Plot (Resolution = ", res, ")"))
    
    # Constructing file paths
    plot_file <- paste0("umap_res_", res, ".png")
    plot_path <- file.path(umap_dir, plot_file)
    
    # Save the image
    ggsave(
      filename = plot_path,
      plot = umap_plot,
      device = "png",
      width = 8,
      height = 6,
      dpi = 300
    )
    
    message("Save UMAP plot: ", plot_path)
  }
  
  # Returns the Seurat object
  return(scdata)
}

# Run cluster and find marker gene repeat to make sure the best resolution
Get_marker_genes <- function(scdata,
                             resolutions,
                             genes_amount = 5,
                             table_dir = "./marker_gene/",
                             recompute_clusters = FALSE) {
  
  # Load necessary package
  require(dplyr)
  require(readr)
  
  # Obtain the variable name of the passed object
  dataset_name <- deparse(substitute(scdata))
  
  # Clean special characters in dataset names
  clean_name <- gsub("[[:space:][:punct:]]", "_", dataset_name)
  
  # Get the precision file path
  table_dir <- paste0(table_dir, clean_name)
  
  # Make sure the output directory exists
  if (!dir.exists(table_dir)) {
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    message("Created directory: ", table_dir)
  }
  
  # Traversing different resolutions
  for (res in resolutions) {
    cluster_col_name <- paste0("TEST_res_", res)
    
    # Conditional computation clustering
    if(recompute_clusters || !cluster_col_name %in% colnames(scdata@meta.data)) {
      scdata <- FindClusters(
        object = scdata,
        resolution = res,
        cluster.name = cluster_col_name
      )
    }
    
    # Set the current cluster as the active flag
    Idents(scdata) <- cluster_col_name
    
    # Finding marker genes
    cluster_markers <- FindAllMarkers(
      scdata,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25
    )
    
    # Extract the top genes for each cluster
    top_markers <- cluster_markers %>%
      group_by(cluster) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = genes_amount)  
    
    # Creating a gene list
    top_marks_table <- top_markers %>%
      group_by(cluster) %>%
      summarise(
        top_genes = paste(gene, collapse = "/"),
        .groups = "drop"
      )
    
    # Constructing file paths
    file_path <- file.path(table_dir, paste0("/cluster_genes_res", res, ".csv"))
    
    # Save the results
    write_csv(top_marks_table, file_path)
    message("Saved marker genes for resolution ", res, " to: ", file_path)
  }
  
  # Returns the modified Seurat object
  return(scdata)
}


Get_marker_genes_all <- function(scdata,
                             resolutions,
                             genes_amount = 5,
                             table_dir = "./marker_gene_all/",
                             recompute_clusters = FALSE) {
  
  # Load necessary packages
  require(dplyr)
  require(readr)
  require(Seurat)
  
  # Get the dataset name (clean special characters)
  dataset_name <- gsub("[[:space:][:punct:]]", "_", deparse(substitute(scdata)))
  
  # Create Output Directory
  dir.create(file.path(table_dir, dataset_name), 
             recursive = TRUE, showWarnings = FALSE)
  
  # Traversing different resolutions
  for (res in resolutions) {
    cluster_col <- paste0("TEST_res_", res)
    
    # Recalculate clusters (if necessary)
    if(recompute_clusters || !cluster_col %in% colnames(scdata@meta.data)) {
      scdata <- FindClusters(scdata, resolution = res, cluster.name = cluster_col)
    }
    
    # Set the active cluster ID
    Idents(scdata) <- cluster_col
    clusters <- levels(Idents(scdata))
    
    # Create a separate directory for each cluster
    res_dir <- file.path(table_dir, dataset_name, paste0("res_", res))
    dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Initialize the summary table
    top_marks_summary <- data.frame()
    
    # Traverse each cluster to find marker genes
    for (cluster_id in clusters) {
      # Find the marker gene of the current cluster (compared with all other clusters)
      markers <- FindMarkers(
        object = scdata,
        ident.1 = cluster_id,
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.25
      )
      
      # Add gene name and cluster information
      markers$gene <- rownames(markers)
      markers$cluster <- cluster_id
      
      # Save complete results to a separate CSV
      cluster_file <- file.path(res_dir, paste0("cluster_", cluster_id, "_markers.csv"))
      write_csv(markers, cluster_file)
      
      # directly merge all genes into the total table
      all_marks_summary <- bind_rows(all_marks_summary, markers)
    }
    
    # Save COMPLETE summary table for current resolution
    summary_file <- file.path(table_dir, dataset_name, paste0("res_", res, "_all_markers.csv"))
    write_csv(all_marks_summary, summary_file)
    message(sprintf("Saved ALL markers for res=%g to: %s", res, summary_file))
  }
  
  return(scdata)
}





# Marker_gene_dictionary
# Create a cell type-marker gene database (expandable S3 object)
create_marker_db <- function() {
  marker_db <- list(
    "T cell" = c("CD3D", "IL7R", "CD4", "CD8A"),
    "B cell" = c("CD19", "MS4A1", "CD79A"),
    "Macrophages" = c("LYZ", "CD14", "FCGR3A"),
    "NK cell" = c("GNLY", "NKG7"),
    "Melanocytes" = c("TYR", "MITF", "SOX10"),
    "Malignant cell" = c("BRAF", "CDKN2A", "TP53"),
    "Tumor-associated fibroblasts" = c("ACTA2","CTSK","CCND1","BCAN")
  )
  class(marker_db) <- "cell_marker_db"
  return(marker_db)
}

# Generic function for updating the database (for subsequent expansion)
add_markers <- function(db, cell_type, genes) {
  UseMethod("add_markers")
}

add_markers.cell_marker_db <- function(db, cell_type, genes) {
  db[[cell_type]] <- unique(c(db[[cell_type]], genes))
  return(db)
}

# Initialize the database
marker_db <- create_marker_db()
