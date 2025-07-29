\name{RunFestemData}
\alias{RunFestemData}
\alias{RunFestemData.Seurat}
\title{Run Festem with normalized data}
\description{
Run Festem algorithm with Seurat pipeline using normalized data instead of raw counts.
This version is specifically designed to work with Seurat objects that only have the "data" layer.
}
\usage{
RunFestemData(object, ...)

\method{RunFestemData}{Seurat}(
  object,
  G = NULL,
  prior = "HVG",
  batch = NULL,
  prior.weight = 0.05,
  prior.weight.filter = 0.9,
  earlystop = 1e-4,
  min.percent = 0.01,
  min.cell.num = 30,
  seed = 321,
  num.threads = 1,
  FDR_level = 0.05,
  assay = "RNA",
  prior_parameters = list(HVG_num = 2000, PC_dims = 50, resolution = 0.7),
  ...
)
}
\arguments{
\item{object}{A Seurat object with normalized data in the "data" layer}

\item{...}{Additional arguments passed to other methods}

\item{G}{Number of clusters. If NULL, determined automatically}

\item{prior}{Prior cluster labels. Can be "HVG" (default), "active.ident", a column name in meta.data, or a vector of labels}

\item{batch}{Batch information. Can be NULL (no batch), a column name in meta.data, or a vector of batch labels}

\item{prior.weight}{Prior weight for EM test}

\item{prior.weight.filter}{Prior weight for filtering}

\item{earlystop}{Early stopping threshold}

\item{min.percent}{Minimum percentage of cells expressing a gene}

\item{min.cell.num}{Minimum number of cells expressing a gene}

\item{seed}{Random seed}

\item{num.threads}{Number of threads for parallel processing}

\item{FDR_level}{FDR level for significant features}

\item{assay}{Name of assay to use}

\item{prior_parameters}{List of parameters for prior generation when prior="HVG"}
}
\details{
This function adapts the Festem algorithm to work with normalized data by:
\enumerate{
  \item Using a modified EM test suitable for continuous normalized data
  \item Skipping count-specific preprocessing steps
  \item Adjusting filtering criteria for normalized values
}
}
\value{
A Seurat object with Festem results added:
\itemize{
  \item Variable features are updated with significant markers
  \item Meta.features are added with Festem statistics:
    \itemize{
      \item Festem_p_adj: Adjusted p-values
      \item Festem_EM: EM test statistics
      \item Festem_rank: Rank in feature list
    }
}
}
\examples{
\dontrun{
# Create a Seurat object with normalized data
seurat_obj <- CreateSeuratObject(counts = your_data)
seurat_obj <- NormalizeData(seurat_obj)

# Run Festem on normalized data
result <- RunFestemData(seurat_obj)
}
} 