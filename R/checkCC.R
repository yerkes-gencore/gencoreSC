#' Perform cell cycle check analysis
#'
#' Calculates phase of cell cycle given a list of cell cycle genes
#'
#' @param obj Seurat object
#' @param what2return c("plot_list", "seurat_obj", "both")
#' @param cc.genes Vector of gene names for cell cycle genes. If you are starting with Ensembl IDs use gencoreSC::EnsDb2GeneName() first.
#' @param features2plot Vector of categorical features to plot PCA for e.g. c("Phase", "mitoRatio", "riboRatio"))
#'
#' @return Seurat object with Phase metadata column, a plot list with PCAs of features specified by features2plot, or both in a list object.
#'
#' @examples
#' \dontrun{
#'   # Inferred by orthology searches against Human list published in Tirosh, I, et al.:
#'   mus_url <- "https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv"
#'   cc_file <- getURL(mus_url)
#'   cc_cycle_genes <- read.csv(text = cc_file)
#'
#'   # Only return seurat object with new Phase metadata column
#'   obj <- checkCC(obj, what2return = "seurat_obj", cc.genes = cc_cycle_genes)
#'
#'   # Run cell cycle check and only return seurat object
#'   cc.plots <- checkCC(obj, what2return = "plot_list", cc.genes = cc_cycle_genes)
#'
#'   # Run cell cycle check and return both
#'   list[obj, cc.plots] %<-% checkCC(obj, what2return = "both", cc.genes = cc_cycle_genes)
#'
#' }
#'
#' @note
#'
#' This function draws on code from:
#'   - https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/cell_cycle_scoring.md
#'
#' To do:
#'  - Fix namespace warnings and notes from using methods defined by AnnotationHub
#'
#' @export
checkCC <- function(obj, what2return = c("plot_list", "seurat_obj", "both"),
                    cc.genes, features2plot = c("Phase", "mitoFr", "riboFr")) {
  # Note that this creates a normalized and scaled object with PCA, which may not want to keep depending on downstream analysis
  obj.cc <- scoreCC(obj, cc.genes)

  if (what2return == "plot_list") {
    p.list <- checkPCA(obj = obj.cc, features = features2plot)
    return(p.list)
  } else if (what2return == "seurat_obj") {
    # Only return obj with resulting metadata, to avoid returning normalized object
    obj.out <- obj
    obj.out@meta.data <- obj.cc@meta.data
    return(obj.out)
  } else if (what2return == "both") {
    p.list <- checkPCA(obj = obj.cc, features = features2plot)
    obj.out <- obj
    obj.out@meta.data <- obj.cc@meta.data
    return(list(p.list, obj.out))
  }
}

#' Add cell cycle metadata
#'
#' Calculates phase of cell cycle given a list of cell cycle genes
#'
#' @param obj Seurat object
#' @param cc.genes Vector of gene names for cell cycle genes. If you are starting with Ensembl IDs use gencoreSC::EnsDb2GeneName() first.
#'
#' @return Seurat object with cell cycle metadata (plus normalized and scaled data and PCA reduction)
#' @export
scoreCC <- function(obj, cc.genes) {
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s_genes <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes
  # Score cells for cell cycle
  obj.phase <- Seurat::NormalizeData(obj, verbose=FALSE)
  obj.phase <- Seurat::CellCycleScoring(obj.phase, verbose=FALSE,
                                   g2m.features = g2m_genes,
                                   s.features = s_genes)
  obj.phase <- Seurat::FindVariableFeatures(obj.phase, verbose=FALSE,
                                       selection.method = "vst",
                                       nfeatures = 2000)
  obj.phase <- Seurat::ScaleData(obj.phase, verbose=FALSE)
  obj.phase <- Seurat::RunPCA(obj.phase, npcs=20, verbose=FALSE)

  # Add CC.Difference
  ## i.e. signals separating non-cycling cells and cycling cells will be maintained,
  # but diffs in cell cycle phase among proliferating cells can be regressed out
  obj.phase$CC.Difference <- obj.phase$S.Score - obj.phase$G2M.Score

  return(obj.phase)
}

#' Plot PCAs for a list of categorical features
#'
#' Take a vector of features and plot PCA plots for each, returning a list.
#'
#' @param obj Seurat object
#' @param features Vector of features to to group.by and split.by in PCA plots
#'
#' @return List of plot objects
#' @export
checkPCA <- function(obj, features = c("Phase", "mitoRatio", "riboRatio")) {
  p.list <- list()
  for (feature in features) {
    p.list[[feature]] <- Seurat::DimPlot(obj,
                               reduction = "pca",
                               group.by= feature,
                               split.by = feature)
  }
  return(p.list)
}

#' Convert cell cycle Ensembl IDs to gene names
#'
#' Convert a vector of cell cycle Ensembl IDs to gene symbols
#'
#' @param cell_cycle_genes Vector of Ensemble IDs for cell cycle genes
#' @param species Binomial species epithet for organism (default = "Mus musculus")
#'
#' @return Vector of cell cycle gene names
#' @export
EnsDb2GeneName <- function(cell_cycle_genes, species = "Mus musculus") {
  ## Based on https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/cell_cycle_scoring.md
  ## These annotations are Ensemble IDs but we need gene names.
  # Connect to AnnotationHub
  ah <- AnnotationHub::AnnotationHub()
  # Access the Ensembl database for organism
  ahDb <- query(ah, pattern = c(species, "EnsDb"), ignore.case = TRUE)
  # Acquire the latest annotation files
  id <- ahDb %>%
    mcols() %>%
    rownames() %>%
    tail(n=1)
  # Download the appropriate Ensembldb databases
  edb <- ah[[id]]
  # Extract gene-level information from database
  annotations <- genes(edb, return.type = "data.frame")
  # Select annotaitons of interest
  annotations <- annotations %>% dplyr::select(.data$gene_id, .data$gene_name, .data$seq_name, .data$gene_biotype, .data$description)
  ## These annotations are Ensemble IDs but we need gene names.
  annotations <- EnsDb2GeneName(species = species)

  ## Now we can use these annotations to get the corresponding gene names for the Ensembl IDs of the cell cycle genes.
  # Get gene names for Ensemble IDs for each gene
  cell_cycle_markers <- cell_cycle_genes %>%
    dplyr::left_join(annotations, by =c("geneID" = "gene_id"))

  cc_genes <- list()
  # Acquire the S phase
  cc_genes$s_genes <- cell_cycle_markers %>%
    dplyr::filter(.data$phase == "S") %>%
    pull("gene_name")
  # Acquire the G2M phase genes
  cc_genes$g2m_genes <- cell_cycle_markers %>%
    dplyr::filter(.data$phase == "G2/M") %>%
    pull("gene_name")

  return(cc_genes)
}
