#' Add cell cycle metadata
#'
#' Calculates phase of cell cycle given a list of cell cycle genes
#'
#' @param obj Seurat object
#' @param cc.genes Named list of vectors of S and G2M gene names. The names must be `s.genes` and `g2m.genes`. If you are starting with Ensembl IDs see examples for code to convert them to gene names.
#'
#' @return Seurat object with cell cycle metadata (plus normalized and scaled data and PCA reduction)
#'
#' @examples
#' \dontrun{
#' ## For human data, it's easy because a list of cc genes, `cc.genes`, are loaded with Seurat automatically
#' load(Seurat)
#' seurat_objs <- lapply(seurat_objs, scoreCC, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes) # here seurat_objs is a list of seurat objects
#' qc_pca_plots <- lapply(seurat_objs, plotPCAs)
#' ggarrange(plotlist = lapply(qc_pca_plots, function(x) x$Phase))
#'
#' ## For another species, you may need to convert from ensembl ids to gene names, like so:
#' # Getting cell cycle genes for Mus musculus (see https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/cell_cycle_scoring.md)
#' library(RCurl)
#' library(AnnotationHub)
#' library(dplyr)
#'
#' # Inferred by orthology searches against Human list published in Tirosh, I, et al.:
#' cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv")
#' cc_cycle_genes <- read.csv(text = cc_file)
#'
#' # These annotations are Ensemble IDs but we need gene names.
#' # Connect to AnnotationHub
#' ah <- AnnotationHub()
#' # Access the Ensembl database for organism
#' ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"), ignore.case = TRUE)
#' # Acquire the latest annotation files
#' id <- ahDb %>% mcols() %>% rownames() %>% tail(n=1)
#' # Download the appropriate Ensembldb databases
#' edb <- ah[[id]]
#' # Extract gene-level information from database
#' annotations <- genes(edb, return.type = "data.frame")
#' # Select annotations of interest
#' annotations <- annotations %>% dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
#'
#' # Now we can use these annotaitons to get teh corresponding gene names for the Ensembl IDs of the cell cycle genes.
#' # Get gene names for Ensemble IDs for each gene
#' cell_cycle_markers <- cc_cycle_genes %>% dplyr::left_join(annotations, by =c("geneID" = "gene_id"))
#'
#' cc_genes <- list()
#' cc_genes$s.genes <- cell_cycle_markers %>% dplyr::filter(phase == "S") %>% pull("gene_name")
#' cc_genes$g2m.genes <- cell_cycle_markers %>% dplyr::filter(phase == "G2/M") %>% pull("gene_name")
#'
#' ## Now you can run scoreCC and check the PCAs as you would for humans
#' seurat_objs <- lapply(seurat_objs, scoreCC, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes) # here seurat_objs is a list of seurat objects
#' qc_pca_plots <- lapply(seurat_objs, plotPCAs)
#' ggarrange(plotlist = lapply(qc_pca_plots, function(x) x$Phase))
#' }
#'
#' @export
scoreCC <- function(obj, s.genes, g2m.genes) {
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
