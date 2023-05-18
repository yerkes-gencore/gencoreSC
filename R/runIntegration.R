#' Run integration via Seurat, Harmony, or STACAS
#'
#' A wrapper function for integration of RNA assays via Seurat, Harmony, or STACAS; normalizing, regressing out or blacklisting features etc. as specified.
#'
#' The intention is to streamline analyses where multiple integration methods or parameter combinations need to be compared. Thus, the input Seurat object need only contain raw RNA counts.
#'
#' @param s.split A split Seurat object containing all samples/captures to integrate
#' @param integration_method String specifying integration method; any of c("merge", "seurat", "harmony", "stacas")
#' @inheritParams NormFindVarFeatScaleData
#' @param npcs Total Number of PCs to compute and store (30 by default)
#' @inheritParams integrate_seurat
#' @inheritParams integrate_harmony
#' @param harmony.group.by.vars Variable(s) to remove (character vector).
#' @param harmony_norm_merged Whether to normalized merged object or normalize each capture separately (logical).
#' @inheritParams integrate_stacas
#' @param STACAS.supervision_labels Vector of cell labels to pass to supervised integration with STACAS. See \code{magrittr::\link[STACAS:Run.STACAS]{Run.STACAS}} `cell.labels` parameter for details.
#'
#' @returns Integrated Seurat object post Seurat::FindClusters() (with Seurat's default resolution = 0.8), ready for plotting UMAPs etc.
#' @export
runIntegration <- function(s.split,
                           integration_method = "merge",
                           nfeatures = 2000,
                           npcs = 30,
                           norm_method = "logNorm",
                           feature.blacklist = NULL,
                           harmony.group.by.vars,
                           regress.out = NULL,
                           STACAS.supervision_labels = NULL,
                           dim.reduct = "cca",
                           harmony_norm_merged = TRUE,
                           verbose = FALSE) {

  # Run integration
  print("Running integration.")
  if (integration_method == "merge") {
    s <- merge(x = s.split[[1]],
               y = s.split[2:length(s.split)],
               merge.data = FALSE)
    s <- NormFindVarFeatScaleData(s, norm_method = norm_method, nfeatures = nfeatures,
                                  regress.out = regress.out, feature.blacklist = feature.blacklist,
                                  verbose = verbose, scale_data = NULL)
    print("Running PCA, UMAP and Finding Neighbors")
    s <- s %>%
      Seurat::RunPCA(npcs = npcs, verbose = verbose) %>%
      Seurat::RunUMAP(reduction = "pca", dims = 1:npcs, verbose=verbose) %>%
      Seurat::FindNeighbors(verbose=verbose)
  } else if (integration_method == "stacas") {
    s <- integrate_stacas(s.split, nfeatures = nfeatures, npcs = npcs,
                          norm_method = norm_method, feature.blacklist = feature.blacklist,
                          regress.out = regress.out, cell.labels = STACAS.supervision_labels, verbose = verbose)
  } else if (integration_method == "seurat") {
    s <- integrate_seurat(s.split, dim.reduct = dim.reduct, nfeatures = nfeatures, npcs = npcs,
                          norm_method = norm_method, feature.blacklist = feature.blacklist,
                          regress.out = regress.out, verbose = verbose)
  } else if (integration_method == "harmony") {
    if (harmony_norm_merged) {
      print("Normalizing merged object before integrating with Harmony.")
      # call it s.split even though it's merged, just so we can recycle the same code in the integrate_harmony() call
      s.split <- merge(x = s.split[[1]],
                       y = s.split[2:length(s.split)],
                       merge.data = FALSE)
    } else if (harmony_norm_merged) {
      print("Normalizing split object before integrating with Harmony.")
    }
    s <- integrate_harmony(s.split, nfeatures = nfeatures, npcs = npcs,
                           norm_method = norm_method, feature.blacklist = feature.blacklist,
                           group.by.vars = harmony.group.by.vars, scale_data = NULL, reduction.save = "harmony",
                           regress.out = regress.out, verbose = verbose)
  } else {
    errorCondition("integration_method must be in c('merge', 'stacas', 'seurat', 'harmony')")
  }

  return(s)
}

#' Run integration via Seurat
#'
#' A wrapper function for integration of RNA assays via Seurat using CCA or RPCA; normalizing, regressing out or blacklisting features etc. as specified.
#'
#' The intention is to streamline analyses where multiple integration methods or parameter combinations need to be compared. Thus, the input Seurat object need only contain raw RNA counts.
#'
#' @param s.split A split Seurat object containing all samples/captures to integrate
#' @param dim.reduct Reduction to use for Seurat integration ("cca" or "rpca"). see \code{\link[Seurat:FindIntegrationAnchors]{Seurat::FindIntegrationAnchors}} reduction parameter for details.
#' @param npcs Total Number of PCs to compute and store (30 by default)
#' @inheritParams Seurat::FindIntegrationAnchors
#' @inheritParams NormFindVarFeatScaleData
#'
#' @returns Integrated Seurat object
#'
#' \dontrun{
#' # log normalize each capture separately and integrate using canonical correlation analysis dim reduction
#' s.Harmony <- integrate_seurat(s.split, norm_method = "logNorm", reduction = "cca")
#'
#' # SCTransform each capture separately and integrate using recripocal PCA dim reduction
#' s.Harmony <- integrate_seurat(s.split, norm_method = "SCT", dim.reduct = "rpca")
#'
#' # SCTransform each capture separately, blacklist TCR genes from variable features, and integrate using cca
#' library(scGate)
#' TCR_genes <- scGate::genes.blacklist.default$Mm$TCR
#' s.Harmony <- integrate_harmony(s, norm_method = "SCT", dim.reduct = "cca", feature.blacklist = "TCR_genes")
#' }
#'
#' @export
integrate_seurat <- function(s.split, dim.reduct = "cca", nfeatures = 2000, npcs = 30, verbose = FALSE,
                             norm_method = "logNorm", feature.blacklist = NULL, regress.out = NULL,
                             k.anchor = 5) {

  s.split <- lapply(X = s.split, FUN = function(x) {
    x <- x %>%
      NormFindVarFeatScaleData(
        norm_method = norm_method,
        nfeatures = nfeatures,
        regress.out = regress.out,
        feature.blacklist = feature.blacklist,
        scale_data = F,
        verbose = verbose)
  })

  # Select integration features
  # NormFindVarFeatScaleData() should have already blacklisted var features if desired, so they should not show up in integration anchors either
  anchor.features <- Seurat::SelectIntegrationFeatures(s.split, nfeatures = nfeatures, verbose = verbose)

  # Must treat SCTransformed data differently from log normalized data. see: https://satijalab.org/seurat/articles/integration_rpca.html and ?PrepSCTIntegration
  if (s.split[[1]]@active.assay == "SCT") {
    s.split <- Seurat::PrepSCTIntegration(object.list = s.split, anchor.features = anchor.features)
    if (dim.reduct == "rpca") {
      s.split <- lapply(X = s.split, FUN = function(x) {
        x <- Seurat::RunPCA(x, features = anchor.features, verbose = verbose)
      })
    }
  } else if (s.split[[1]]@active.assay == "RNA") {
    # Must calculate PCA first if using RPCA integration algorithm
    if (dim.reduct == "rpca") {
      s.split <- lapply(X = s.split, FUN = function(x) {
        x <- Seurat::ScaleData(x, features = anchor.features, verbose = verbose)
        x <- Seurat::RunPCA(x, features = anchor.features, verbose = verbose)
      })
    }
  }

  # Find integration anchors
  IntegAnchors <- Seurat::FindIntegrationAnchors(s.split, anchor.features = anchor.features,
                                         reduction = dim.reduct, dims = 1:npcs,
                                         verbose = verbose)
  # Integrate data
  normalization.method <- ifelse(s.split[[1]]@active.assay == "RNA", "LogNormalize", s.split[[1]]@active.assay)
  s.Integrated <- Seurat::IntegrateData(anchorset = IntegAnchors, dims = 1:npcs,
                                        normalization.method = normalization.method, verbose = verbose)
  # Setup UMAP, clustering
  s.Integrated <- s.Integrated %>%
    Seurat::ScaleData(verbose = verbose) %>%
    Seurat::RunPCA(npcs=npcs, verbose = verbose) %>%
    Seurat::RunUMAP(dims=1:npcs, verbose = verbose) %>%
    Seurat::FindNeighbors(reduction="pca", dims = 1:npcs, verbose = verbose) %>%
    Seurat::FindClusters(verbose = verbose)

  return(s.Integrated)
}

#' Run integration via Harmony
#'
#' A wrapper function for integration of RNA assays via \code{\link[Harmony:harmony]{Harmony}}; normalizing, regressing out or blacklisting features etc. as specified. The intention is to streamline analyses where multiple integration methods or parameter combinations need to be compared. Thus, the input Seurat object need only contain raw RNA counts.
#'
#' @param s A Seurat object. If it is a split object, will run normalize samples separately. If merged, will normalize samples together.
#' @param npcs Total Number of PCs to compute and store (30 by default)
#' @inheritParams harmony::RunHarmony
#' @inheritParams NormFindVarFeatScaleData
#'
#' @returns Integrated Seurat object
#'
#' \dontrun{
#' # log normalize each capture separately and integrate via Harmony
#' # s.split is a Seurat object split by "Sample_Name"
#' s.Harmony <- integrate_harmony(s.split, norm_method = "logNorm", group.by.vars = "Sample_Name)
#'
#' # SCTransform each capture separately and integrate via Harmony
#' s.Harmony <- integrate_harmony(s, norm_method = "SCT", group.by.vars = "Sample_Name)
#'
#' # SCTransform each capture separately, blacklist TCR genes from variable features, and integrate via Harmony
#' library(scGate)
#' TCR_genes <- scGate::genes.blacklist.default$Mm$TCR
#' s.Harmony <- integrate_harmony(s, norm_method = "SCT", group.by.vars = "Sample_Name, feature.blacklist = "TCR_genes")
#'
#' # SCTransform all captures together and integrate via Harmony
#' s <- merge(x = s.split[[1]], y = s.split[2:length(s.split)], merge.data = FALSE)
#' s.Harmony <- integrate_harmony(s, norm_method = "SCT", group.by.vars = "Sample_Name)
#' }
#'
#' @export
# Run Harmony; run v1 if given merged object, run v2 if given list object
integrate_harmony <- function(s, norm_method = "logNorm", nfeatures = 2000, npcs = 30,
                              group.by.vars = "Sample_Name", verbose = FALSE,
                              regress.out = NULL, feature.blacklist = NULL,
                              scale_data = NULL, reduction.save = "harmony") {

  # the native Seurat default for classic is to scale data, but for sctransform is to not scale data
  if (is.null(scale_data) & norm_method == "logNorm") {
    scale_data <- T
  } else if (is.null(scale_data) & norm_method == "SCT") {
    scale_data <- F
  }

  # normalize and scale, prep for Harmony
  if (!is.list(s)) {
    print("Normalizing merged object before integrating with Harmony.")
    DefaultAssay(s) <- "RNA"
    s <- s %>%
      NormFindVarFeatScaleData(
        norm_method = norm_method,
        nfeatures = nfeatures,
        regress.out = regress.out,
        feature.blacklist = feature.blacklist,
        scale_data = scale_data,
        verbose = verbose)
    s <- RunPCA(s, npcs = npcs, verbose = verbose)
  } else if (is.list(s)) {
    print("Normalizing split object before integrating with Harmony.")
    s.split <- s
    s.split <- lapply(X = s.split, FUN = function(x) {
      DefaultAssay(x) <- "RNA"
      x <- x %>%
        NormFindVarFeatScaleData(
          norm_method = norm_method,
          nfeatures = nfeatures,
          regress.out = regress.out,
          feature.blacklist = feature.blacklist,
          scale_data = scale_data,
          verbose = verbose)
    })
    # Find most variable features across samples to integrate
    integ_features <- SelectIntegrationFeatures(s.split, nfeatures = nfeatures, verbose = verbose)
    # Must treat SCTransformed data differently from log normalized data. see: https://satijalab.org/seurat/articles/integration_rpca.html and ?PrepSCTIntegration
    if (s.split[[1]]@active.assay == "SCT") {
      s.split <- PrepSCTIntegration(object.list = s.split, anchor.features =  integ_features)
    }
    # Merge normalized samples
    s <- merge(x = s.split[[1]], y = s.split[2:length(s.split)], merge.data = TRUE)
    Seurat::DefaultAssay(s) <- ifelse(norm_method == "SCT", "SCT", "RNA")
    # Manually set variable features of merged Seurat object
    Seurat::VariableFeatures(s) <- integ_features
    # Not sure why I need to scale data again after merging?
    s <- Seurat::ScaleData(s, verbose = verbose)
    # Calculate PCs using manually set variable features
    s <- Seurat::RunPCA(s, npcs = npcs, verbose = verbose)
  }
  # Run Harmony
  s.Harmony <- s %>%
    harmony::RunHarmony(group.by.vars = group.by.vars,
                        reduction = "pca",
                        reduction.save = reduction.save,
                        verbose = verbose) %>%
    Seurat::RunUMAP(reduction=reduction.save, dims=1:npcs, verbose = verbose) %>%
    Seurat::FindNeighbors(reduction=reduction.save, dims = 1:npcs, verbose = verbose) %>%
    Seurat::FindClusters(verbose = verbose)

  return(s.Harmony)
}

#' Run integration via STACAS
#'
#' A wrapper function for integration of RNA assays via \code{\link[STACAS:Run.STACAS]{STACAS}}; normalizing, regressing out or blacklisting features etc. as specified. The intention is to streamline analyses where multiple integration methods or parameter combinations need to be compared. Thus, the input Seurat object need only contain raw RNA counts.
#'
#' @param s.split A split Seurat object containing all samples/captures to integrate
#' @inheritParams STACAS::Run.STACAS
#' @inheritParams NormFindVarFeatScaleData
#' @param npcs Total Number of PCs to compute and store (30 by default)
#'
#' @returns Integrated Seurat object
#'
#' @examples
#' \dontrun{
#' # log normalize, find variable features blacklisting TCR genes (Mus musculus), integrate via unsupervised STACAS
#' library(scGate)
#' scGate::genes.blacklist.default$Mm %>% names()
#' TCR_genes <- scGate::genes.blacklist.default$Mm$TCR
#' s.STACAS <- integrate_stacas(s.split, norm_method = "logNorm", feature.blacklist = TCR_genes)
#'
#' # unsupervised STACAS, blacklisting cell cycle, heatshock, mito, ribo, IFn response and TCR genes
#' blacklist_all <- scGate::genes.blacklist.default$Mm %>% unlist() %>% unname()
#' s.STACAS <- integrate_stacas(s.split, norm_method = "logNorm", feature.blacklist = blacklist_all)
#'
#' # "fully-supervised" STACAS, using cell labels from previously running SingleR with the ImmGen dataset stored in the Seurat object metadata
#' # see https://carmonalab.github.io/STACAS.demo/STACAS.demo.html for details
#' s.STACAS <- integrate_stacas(s.split, norm_method = "logNorm", feature.blacklist = blacklist_all,
#'                              cell.labels = "ImmGen.main.pruned.labels")
#'
#' # "semi-supervised" STACAS, using a subset of cell labels from previously running SingleR with the ImmGen dataset stored in the Seurat object metadata
#' # see https://carmonalab.github.io/STACAS.demo/STACAS.demo.html for details
#' s.split$Tlike_cells <- ifelse(s.split$ImmGen.main.pruned.labels %in% c("T cells", "NKT", "NK cells", "ILC", "Tgd"),
#'                               s.split$ImmGen.main.pruned.labels,
#'                               NA)
#' s.STACAS <- integrate_stacas(s.split, norm_method = "logNorm", feature.blacklist = blacklist_all,
#'                              cell.labels = "Tlike_cells")
#' }
#'
#' @export
# Run Harmony; run v1 if given merged object, run v2 if given list object
integrate_stacas <- function(s.split, nfeatures = 2000, npcs = 30,
                             norm_method = "logNorm", feature.blacklist = NULL,
                             regress.out = NULL, cell.labels = NULL, verbose = FALSE) {

  s.split <- lapply(X = s.split, FUN = function(x) {
    x <- x %>%
      NormFindVarFeatScaleData(
        norm_method = norm_method,
        nfeatures = nfeatures,
        regress.out = regress.out,
        feature.blacklist = feature.blacklist,
        scale_data = F, # Run.STACAS
        verbose = verbose)
  })

  s.STACAS <- STACAS::Run.STACAS(s.split, dims = 1:npcs, anchor.features = nfeatures, verbose=verbose,
                                 genesBlockList = feature.blacklist,
                                 cell.labels = cell.labels) %>%
    Seurat::RunUMAP(dim = 1:npcs, verbose=verbose) %>%
    Seurat::FindNeighbors(reduction="pca", dims = 1:npcs, verbose = verbose) %>%
    Seurat::FindClusters(verbose = verbose)

  return(s.STACAS)
}
