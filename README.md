# scMetaInspector

An R toolkit for comprehensive single-cell RNA-seq data inspection, quality control, and analysis automation using Seurat objects.

## Overview

scMetaInspector provides a collection of functions to streamline the single-cell RNA-seq analysis workflow. It focuses on metadata inspection, frequency analysis, clustering exploration, and marker gene identification.

## Features

### 1. **Data Quality Control**
- **`Check_NA_Empty()`** - Identifies missing values (NA) and empty strings in metadata
  - Generates detailed CSV reports showing the count and distribution of missing values per feature
  - Useful for data cleaning and preprocessing decisions

### 2. **Frequency Analysis**
- **`Frequency_Count()`** - Analyzes value distribution for a single metadata column
  - Generates frequency tables with counts and percentages
  - Saves results sorted by frequency
  
- **`Frequency_Count_List()`** - Batch processes multiple selected columns
  - Apply frequency analysis to a list of specific metadata columns
  - Returns results as a named list
  
- **`Frequency_Count_All()`** - Automatically analyzes all metadata columns
  - Configurable maximum unique value threshold (default: 50)
  - Skips high-cardinality columns
  - Generates summary of columns with single unique values
  - Ideal for comprehensive dataset exploration

### 3. **Clustering & Dimensionality Reduction**
- **`Preparation_folders()`** - Sets up directory structure for analysis outputs
  - Creates `plot/`, `report/`, and `marker_gene/` directories
  
- **`RunPCA_PlotElbow()`** - Performs PCA and generates elbow plots
  - Helps determine optimal number of dimensions
  - Saves publication-quality TIFF plots
  
- **`Try_different_resolutions()`** - Tests multiple clustering resolutions
  - Runs clustering and UMAP for different Leiden resolutions
  - Generates comparative visualizations
  - Facilitates resolution selection

### 4. **Marker Gene Discovery**
- **`Get_marker_genes()`** - Identifies top marker genes per cluster
  - Returns top N genes (default: 5) for each cluster
  - Customizable statistics thresholds (min.pct = 0.25, logfc.threshold = 0.25)
  - Summarizes results with gene combinations per cluster
  
- **`Get_marker_genes_all()`** - Comprehensive marker gene analysis
  - Returns all identified marker genes, not just top candidates
  - Creates per-cluster and per-resolution output files
  - More detailed than `Get_marker_genes()`

### 5. **Cell Type Annotation Database**
- **`create_marker_db()`** - Initializes a cell type-marker gene database
  - Includes common immune and tumor-associated cell types
  - Pre-populated with canonical marker genes
  
- **`add_markers()`** - Generic S3 method for expanding the database
  - Allows users to add custom cell types and markers
  - Extensible design for project-specific applications

## Installation

```r
# Ensure you have the required packages installed
install.packages(c("dplyr", "readr", "ggplot2"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Seurat")

# Source the scMetaInspector functions
source("scMetaInspector.R")
```

## Dependencies

- **Seurat** - Single-cell analysis framework
- **dplyr** - Data manipulation
- **readr** - CSV I/O
- **ggplot2** - Visualization

## Usage Examples

### Quality Control
```r
# Check for missing values in your Seurat object
na_report <- Check_NA_Empty(seurat_object)
head(na_report)
```

### Metadata Exploration
```r
# Analyze frequency of cell types
freq_table <- Frequency_Count(seurat_object, column = "celltype")

# Analyze multiple columns at once
freq_list <- Frequency_Count_List(seurat_object, 
                                  columns = c("celltype", "batch", "treatment"))

# Automatically analyze all columns (skip high-cardinality columns)
all_freq <- Frequency_Count_All(seurat_object, max_unique = 50)
```

### Clustering Workflow
```r
# Prepare directories
Preparation_folders()

# Find optimal dimensions
seurat_object <- RunPCA_PlotElbow(seurat_object, sample_name = "MySample")

# Test different clustering resolutions
seurat_object <- Try_different_resolutions(seurat_object, 
                                           resolutions = c(0.2, 0.4, 0.6, 0.8),
                                           dims = 1:15)

# Get marker genes for each resolution
seurat_object <- Get_marker_genes(seurat_object, 
                                  resolutions = c(0.4, 0.6),
                                  genes_amount = 10)
```

### Cell Type Annotation
```r
# Create marker database
marker_db <- create_marker_db()

# Add custom cell types
marker_db <- add_markers(marker_db, "Fibroblast", c("COL1A1", "COL1A2", "PDGFRA"))

# View available cell types
names(marker_db)
```

## Output Structure

Functions automatically generate organized output:

```
├── plot/
│   └── {sample_name}/
│       ├── {sample_name}_elbow.tiff
│       └── UMAP/
│           ├── umap_res_0.2.png
│           ├── umap_res_0.4.png
│           └── ...
├── report/
│   ├── NAreport/
│   │   └── {sample_name}_NAreport.csv
│   └── Frequencyreport/
│       └── {sample_name}/
│           └── {column_name}_freq.csv
└── marker_gene/
    └── {sample_name}/
        └── cluster_genes_res{resolution}.csv
```

## Key Features

- **Automatic file naming**: Safe character conversion prevents path errors
- **Error handling**: Graceful handling of edge cases (single-value columns, high-cardinality features)
- **Flexible parameters**: Most functions accept customizable thresholds and parameters
- **Report generation**: All results saved as CSV for downstream analysis
- **Reproducible**: Consistent output formatting and structure

## Parameters Reference

| Function | Key Parameters | Notes |
|----------|----------------|-------|
| `Check_NA_Empty()` | `scdata` | Requires Seurat object |
| `Frequency_Count()` | `column`, `sample_name` | Column can be string or unquoted |
| `Frequency_Count_All()` | `max_unique` | Default: 50 unique values |
| `RunPCA_PlotElbow()` | `ndims`, `width`, `height` | Default: 50 dims, 6x6 plot |
| `Try_different_resolutions()` | `resolutions`, `dims`, `min.dist`, `spread` | Customizable UMAP parameters |
| `Get_marker_genes()` | `genes_amount`, `recompute_clusters` | Default: 5 top genes per cluster |
| `Get_marker_genes_all()` | `genes_amount` | Returns all markers (more detailed) |

## Notes

- All functions expect **Seurat objects** as input
- Output files are saved to relative paths (ensure working directory is set appropriately)
- Functions modify Seurat metadata with cluster assignments (e.g., `TEST_res_0.4`)
- Safe file naming handles special characters automatically

## Add new markers or cell type for marker database

To expand the marker gene database or add new cell type definitions, use `add_markers()`:

```r
# Example: Adding custom markers
marker_db <- add_markers(marker_db, 
                        cell_type = "Endothelial cell",
                        genes = c("PECAM1", "CDH5", "VWF"))
```

## License

MIT License

Copyright (c) 2026 Shixiang Yan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Author

Created by **Shixiang Yan**  
University of Liverpool  
ORCID: [0009-0006-8978-3786](https://orcid.org/0009-0006-8978-3786)

## Citation

If you use **scMetaInspector** in your research, please cite it as:

Yan S. (2026). *scMetaInspector: A lightweight R toolkit for Seurat metadata inspection*.  
GitHub repository: https://github.com/nyarlathotep438/scMetaInspector  
ORCID: 0009-0006-8978-3786

## Notes by human

*Most of the above content was generated by Github Copilot, with only minor modifications to some problematic parts. Copilot is very powerful and its judgments for my functions are generally accurate.*
```
