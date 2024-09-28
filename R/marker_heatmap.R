#' Heatmap visualization of markers identified by Festem
#' 
#' @description Heatmap of marker gene expressions with given group labels.
#' 
#' @param object A Seurat object, matrix or dgCMatrix containing the original count matrix. 
#' If a Seurat object is provided, further information on cells can be included in the metadata slot 
#' and specified with parameters \code{group_by} 
#' If a matrix or dgCMatrix is provided, its rows should represent different
#' genes and columns stand for cells. We recommend setting the rownames as the name of genes.
#' If row names are not specified, numerical indices will be created to represent genes
#' according to their order in the input matrix.
#' @param features a vector of genes to plot.
#' @param plot_cell_prop A float specifying the proportion of cells to plot. Only has an effect when the total
#' cell number is larger than 500.
#' @param color_bar_max A float specifying the maximum range of the color map.
#' @param ... other parameters
#' 
#' @examples 
#' # Load example data
#' data("example_data")
#' example_data <- Seurat::CreateSeuratObject(example_data$counts,meta.data = example_data$metadata)
#' MarkerHeatmap(example_data,features = rownames(example_data)[1:10],group_by = "labels")
#' @rdname MarkerHeatmap
#' @export

MarkerHeatmap <- function(object,features,plot_cell_prop = 0.1,color_bar_max = 5,...){
  UseMethod("MarkerHeatmap")
}

#' @rdname MarkerHeatmap
#' 
#' @param group_by A character string or \code{NULL} (default). If \code{NULL}, \code{active.ident} is used
#' If a string, then the corresponding column in \code{meta.data} is taken as labels for cells.
#' @param assay Name of Assay used
#' @export

MarkerHeatmap.Seurat <- function(object,features,
                                 plot_cell_prop = 0.1,color_bar_max = 5,
                                 group_by = NULL,assay = "RNA",...){
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running Festem on a Seurat object requires Seurat")
  }
  object@active.assay <- assay
  if (is.null(group_by)){
    cluster <- object@active.ident
  } else if (group_by %in% colnames(object@meta.data)){
    cluster <- object@meta.data[,group_by]
  } else{
    stop("group_by should be the name of a column in metadata.")
  }
  
  object <- Seurat::NormalizeData(object)
  object <- Seurat::ScaleData(object, features = features)
  set.seed(321)
  tmp <- data.frame(cell_id = colnames(object),
                    cluster = cluster)
  
  if (is.null(plot_cell_prop) || ncol(object) < 500){
    data <- SeuratObject::LayerData(object,assay = assay,layer = "scale.data", features = features)
  } else{
    tmp <- tmp %>%
      group_by(.data$cluster) %>%
      sample_frac(plot_cell_prop)
    data <- SeuratObject::LayerData(object,assay = assay,layer = "scale.data", features = features, cells = tmp$cell_id)
  }
  
  anno <- data.frame(cluster = tmp$cluster)
  rownames(anno) <- tmp$cell_id
  
  pheatmap::pheatmap(data, 
           color=colorRampPalette(c("navy", "white", "red"))(50),
           breaks = seq(-color_bar_max,color_bar_max,length.out = 51),
           cluster_rows = T, cluster_cols = F,
           border_color = F,
           show_colnames = F, show_rownames = F,
           annotation_col = anno,
           treeheight_row = 0,
           clustering_distance_cols = 'euclidean',
           clustering_distance_rows = 'euclidean',
           clustering_method = 'ward.D')
  
}

#' @rdname MarkerHeatmap
#' 
#' @param label A vector containing labels for cells.
#' 
#' @export
MarkerHeatmap.matrix <- function(object,features,
                                 plot_cell_prop = 0.1,color_bar_max = 5,
                                 label,...){
  
  set.seed(321)
  object <- t(scale(object))
  tmp <- data.frame(cell_id = colnames(object),
                    cluster = label)
  
  if (is.null(plot_cell_prop) || ncol(object) < 500){
    data <- object[features,]
  } else{
    tmp <- tmp %>%
      group_by(.data$cluster) %>%
      sample_frac(plot_cell_prop)
    data <- object[features,tmp$cell_id]
  }
  
  anno <- data.frame(cluster = tmp$cluster)
  rownames(anno) <- tmp$cell_id
  
  pheatmap::pheatmap(data, 
           color=colorRampPalette(c("navy", "white", "red"))(50),
           breaks = seq(-color_bar_max,color_bar_max,length.out = 51),
           cluster_rows = T, cluster_cols = F,
           border_color = F,
           show_colnames = F, show_rownames = F,
           annotation_col = anno,
           treeheight_row = 0,
           clustering_distance_cols = 'euclidean',
           clustering_distance_rows = 'euclidean',
           clustering_method = 'ward.D')
  
}

#' @rdname MarkerHeatmap
#' 
#' @param label A vector containing labels for cells.
#' 
#' @export
MarkerHeatmap.Matrix <- function(object,features,
                                 plot_cell_prop = 0.1,color_bar_max = 5,
                                 label,...){
  
  set.seed(321)
  object <- t(scale(object))
  tmp <- data.frame(cell_id = colnames(object),
                    cluster = label)
  
  if (is.null(plot_cell_prop) || ncol(object) < 500){
    data <- object[features,]
  } else{
    tmp <- tmp %>%
      group_by(.data$cluster) %>%
      sample_frac(plot_cell_prop)
    data <- object[features,tmp$cell_id]
  }
  
  anno <- data.frame(cluster = tmp$cluster)
  rownames(anno) <- tmp$cell_id
  
  pheatmap::pheatmap(data, 
           color=colorRampPalette(c("navy", "white", "red"))(50),
           breaks = seq(-color_bar_max,color_bar_max,length.out = 51),
           cluster_rows = T, cluster_cols = F,
           border_color = F,
           show_colnames = F, show_rownames = F,
           annotation_col = anno,
           treeheight_row = 0,
           clustering_distance_cols = 'euclidean',
           clustering_distance_rows = 'euclidean',
           clustering_method = 'ward.D')
  
}