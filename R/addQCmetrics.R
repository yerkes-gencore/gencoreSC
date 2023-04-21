#' Add QC metrics
#'
#' Calculates qc metrics and reformats Seurat object metadata to be more readable
#'
#' @param obj Seurat object
#' @param mito.pattern pattern defining mitochondiral genes
#' @param ribo.pattern pattern defining ribosomal genes
#'
#' @return Seurat object with modified metadata
#'
#' @export
addQCmetrics <-
  function(obj, mito.pattern="^MT", ribo.pattern="^RP[SL]") {
    # Add number of genes per UMI for each cell to metadata
    obj$log10GenesPerUMI <- log10(obj$nFeature_RNA)/log10(obj$nCount_RNA)
    # Mitochondrial ratio
    obj$mitoRatio <- Seurat::PercentageFeatureSet(object=obj, pattern=mito.pattern)/100
    # Ribosomal ratio
    obj$riboRatio <- Seurat::PercentageFeatureSet(object=obj, pattern=ribo.pattern)/100
    # Create metadata dataframe
    metadata <- obj@meta.data
    # Add cell IDs to metadata
    metadata$cells <- rownames(metadata)
    # Rename columns to be more readable/intuitive
    metadata <- metadata %>%
      dplyr::rename(capID = .data$orig.ident,
                    nUMI = .data$nCount_RNA,
                    nGene = .data$nFeature_RNA)
    obj@meta.data <- metadata
    return(obj)
  }
