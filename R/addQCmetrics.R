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
    # Turn mitoRatio into categorical factor vector based on quartile values
    mito_qs <- summary(obj$mitoRatio)
    obj$mitoFr <- cut(obj$mitoRatio,
                      breaks=c(-Inf, mito_qs[2], mito_qs[3], mito_qs[4], Inf),
                      labels=c("Low","Medium","Med-high", "High"))

    # Ribosomal ratio
    obj$riboRatio <- Seurat::PercentageFeatureSet(object=obj, pattern=ribo.pattern)/100
    # Turn mitoRatio into categorical factor vector based on quartile values
    ribo_qs <- summary(obj$riboRatio)
    obj$riboFr <- cut(obj$riboRatio,
                      breaks=c(-Inf, ribo_qs[2], ribo_qs[3], ribo_qs[4], Inf),
                      labels=c("Low","Medium","Med-high", "High"))

    # Create metadata dataframe
    metadata <- obj@meta.data
    # Add cell IDs to metadata
    metadata$cellID <- rownames(metadata)
    # Rename columns to be more readable/intuitive
    metadata <- metadata %>%
      dplyr::rename(capID = "orig.ident",
                    nUMI = "nCount_RNA",
                    nGene = "nFeature_RNA")
    obj@meta.data <- metadata
    return(obj)
  }
