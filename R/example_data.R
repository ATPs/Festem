#' A toy example of scRNA-seq data in Seurat object
#'
#' A toy example of scRNA-seq data containing two cell types. Cell labels and counts are stored in Seurat object
#' 
#' @format A list containing expression counts and metadata of 100 genes across 485 cells. Can be convert to a Seurat object via \code{Seurat::CreateSeuratObject}
#' 
#' @examples
#' data("example_data")
#' example_data <- Seurat::CreateSeuratObject(example_data$counts,meta.data = example_data$metadata)
"example_data"