#' Add cell cycle metadata
#'
#' Calculates phase of cell cycle given a list of cell cycle genes
#'
#' @param obj Seurat object
#' @param s.genes Vector of S phase genes as gene names.
#' @param g2m.genes Vector of G2M phase genes as gene names.
#'
#' @return Seurat object with cell cycle metadata (plus normalized and scaled data and PCA reduction)
#'
#' @note If you are starting with Ensembl IDs see examples for code to convert them to gene names (see examples below).
#'
#' @examples
#' \dontrun{
#' ## For human data, a list of cc genes, `cc.genes`, are loaded with Seurat automatically
#' load(Seurat)
#' # where seurat_objs list of seurat objects:
#' seurat_objs <- lapply(seurat_objs, scoreCC,
#'                       s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes)
#' qc_pca_plots <- lapply(seurat_objs, plotPCAs)
#' ggarrange(plotlist = lapply(qc_pca_plots, function(x) x$Phase))
#'
#' ## For another species, you may need to convert from ensembl ids to gene names, like so:
#' # Getting cell cycle genes for Mus musculus
#' # See https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/cell_cycle_scoring.md
#' library(RCurl)
#' library(AnnotationHub)
#' library(dplyr)
#'
#' # Inferred by orthology searches against Human list published in Tirosh, I, et al.:
#' url <- "https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv"
#' cc_file <- getURL(url)
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
#' annotations <- annotations %>%
#'   dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
#'
#' # Now we can use these annotations to get the corresponding gene names for the Ensembl IDs.
#' # Get gene names for Ensemble IDs for each gene
#' cell_cycle_markers <- cc_cycle_genes %>%
#'   dplyr::left_join(annotations, by =c("geneID" = "gene_id"))
#'
#' s.genes <- cell_cycle_markers %>%
#'   dplyr::filter(phase == "S") %>%
#'   pull("gene_name")
#' g2m.genes <- cell_cycle_markers %>%
#'   dplyr::filter(phase == "G2/M") %>%
#'   pull("gene_name")
#'
#' ## Now you can run scoreCC and check the PCAs as you would for humans
#' seurat_objs <- lapply(seurat_objs, scoreCC, s.genes = s.genes, g2m.genes = g2m.genes)
#' qc_pca_plots <- lapply(seurat_objs, plotPCAs)
#' ggarrange(plotlist = lapply(qc_pca_plots, function(x) x$Phase))
#' }
#'
#' @export
scoreCC <- function(obj, s.genes, g2m.genes) {
  # Score cells for cell cycle
  obj.phase <- Seurat::NormalizeData(obj, verbose=FALSE)
  obj.phase <- Seurat::CellCycleScoring(obj.phase, verbose=FALSE,
                                   g2m.features = g2m.genes,
                                   s.features = s.genes)
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
