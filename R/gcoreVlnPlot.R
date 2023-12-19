#' Violin plot of gene(s) expression grouped by metadata
#'
#' Similar to `Seurat::VlnPlot` but more customizable. Pass a seurat object and
#'  genes to plot. Optionally subset the data by a meta.data variable defined
#'  with `subset_vars`, selecting level(s) `subsets`. Useful to plot expression
#'  in one subset of your data without having to make separate objects.
#'  Group the data by meta.data variable `grouping_var`, selecting levels `groups`.
#'  Optionally filter zeros to focus on cells expressing gene.
#'
#'  If more than one value is provided to `subset_vars`, you can specify multiple
#'    subsetting criteria to `subsets` via a list. So if you wanted to subset
#'    your object on both treatment and condition, you could specify that with
#'    something like:
#'
#'    `subset_vars = c('timepoint', 'condition')` and
#'    `subsets = list(c('D0'), c('2ml', '4ml'))`,`
#'
#'  Requires `ggforce` to be installed
#'
#' @param obj Seurat object
#' @param genes Genes to plot
#' @param grouping_var Column of `obj@meta.data` to group data by.
#' @param groups Optional. Levels of `grouping_var` to include in the plot. Also used
#'  to specify order of levels.
#' @param subset_vars  Optional. Column of `obj@meta.data` to subset on. Default
#'  is `seurat_clusters` so you could `subset` on a specific cluster.
#' @param subsets Levels of `subset_vars` to subset data to
#' @param filter_zeros Remove 0s from plot: default `TRUE`.
#' @param assay Assay to pull expression data from, default `RNA`
#'
#' @returns A ggplot
#' @export
#'
#' @examples
#' \dontrun{
#' ## If you only have one value to subset for, the list convention isn't too picky
#' gcoreVlnPlot(obj = obj,
#'   genes = genes,
#'   assay = 'RNA',
#'   subset_var = 'cell_type',
#'   subset = 'ISG-rich',
#'   grouping_var = 'Group',
#'   groups = c("Control", "1D3"))
#'
#' ## If you want to subset for multiple values of a single subset_var, do this
#' gcoreVlnPlot(obj = obj,
#'   genes = genes,
#'   assay = 'RNA',
#'   subset_var = 'cell_type',
#'   subset = list(c('ISG-rich', 'Other')),
#'   grouping_var = 'Group',
#'   groups = c("Control", "1D3"))
#'
#' ## If you want to subset for multiple values of multiple subset_vars
#' gcoreVlnPlot(obj = obj,
#'   genes = genes,
#'   assay = 'RNA',
#'   subset_var = c('cell_type', 'timepoint'),
#'   subset = list(c('ISG-rich', 'Other'), c('pre')),
#'   grouping_var = 'Group',
#'   groups = c("Control", "1D3"))
#'
#' }
#'
gcoreVlnPlot <- function(obj,
                         genes,
                         grouping_var,
                         groups = NULL,
                         subsets = NULL,
                         subset_vars = NULL,
                         filter_zeros = TRUE,
                         assay = 'RNA'){
  if (!require(ggforce)) {
    stop('Error: Requires ggforce to be installed')
  }
  if (missing(grouping_var)) {
    stop('Use the "grouping_var" argument
                   to specify a metadata variable to group observations')
  }
  if (missing(genes)) {
    stop('Specify gene(s) to plot data for')
  }
  if (!missing(subsets) & missing(subset_vars)) {
    stop('Subset_vars not specified for subsetting')
  }
  if (missing(subsets) & !missing(subset_vars)) {
    warning('No subsets specified for subset_vars')
  }
  if (length(subsets) != length(subset_vars)) {
    stop('Error: subsets should be provided as a list of character vectors,
         where the number of vectors equals the number of entries in
         subset_vars. See examples for details')
  }

  if (!missing(subsets) & !missing(subset_vars)) {
    for (i in 1:length(subset_vars)) {
      subset_var <- subset_vars[i]
      if (subset_var %in% colnames(obj@meta.data)) {
        for (subset in subsets[i]) {
          if (subset %in% unique(obj@meta.data[[subset_var]])) {
            obj <- obj[,obj@meta.data[[subset_var]] %in% subset]
          } else {
            stop('Error: subset value :"',subsets,
                 '" not present as a value of ',
                 subset_var[i], 'column in object metadata.')
          }
        }
      } else {
        stop('Error: subset_var value :"',subset_var,
             '" not present as a column in the object metadata.')
      }
    }
  }

  good_genes <- c()
  bad_genes  <- c()
  for (gene in genes) {
    if (gene %in% rownames(obj@assays[[assay]]@data)) {
      good_genes <- c(good_genes, gene)
    } else {
      bad_genes <- c(bad_genes, gene)
    }
  }
  if (length(bad_genes) > 0) {
    warning('The following genes were not found in the object: ', genes)
  }
  if (length(good_genes) <= 0) {
    stop('No genes found to plot data for!')
  }
  genes <- good_genes
  mat_to_plot <- reshape2::melt(as.matrix(obj@assays[[assay]]@data)[genes,])
  if (length(genes) == 1) {
    mat_to_plot$Var2 <- rownames(mat_to_plot)
  }
  mat_to_plot <- merge(mat_to_plot,
                       obj@meta.data %>%
                         as.data.frame() %>%
                         dplyr::select(.data[[grouping_var]]),
                       by.x = 'Var2', by.y = 0)
  if (grouping_var %in% colnames(obj@meta.data)) {
    if (!is.null(groups)) {
      if (all(groups %in% unique(obj@meta.data[[grouping_var]]))) {
        mat_to_plot <- mat_to_plot %>%
          filter(.data[[grouping_var]] %in% groups)
        mat_to_plot[[grouping_var]] <- factor(mat_to_plot[[grouping_var]],
                                              levels = groups)
      } else {
        stop('Error: Specified groups are not all present in the ',
             grouping_var, ' column of the object metadata.')
      }
    }  else {
      ## Groups aren't specified
      mat_to_plot[[grouping_var]] <- factor(mat_to_plot[[grouping_var]])
    }
  } else {
    stop('Error: grouping_var value :"',grouping_var,
         '" not present as a column in the object metadata.')
  }



  if (filter_zeros) {mat_to_plot <- mat_to_plot %>% filter(.data[['value']] != 0)}

  out <- ggplot2::ggplot(mat_to_plot,
                         aes(x = .data[[grouping_var]],
                             y = .data[['value']],
                             color = .data[[grouping_var]])) +
    # geom_jitter() +
    ggplot2::geom_violin(draw_quantiles = 0.5) +
    ggforce::geom_sina(size = 0.01, alpha = 0.2) +
    ggplot2::theme_bw() +
    ggplot2::labs(caption = paste0(if(!is.null(subsets)) {paste0('Showing expression of ', paste0(subsets, collapse = ', '), ' cells\n')},
                                   if(!is.null(filter_zeros)) {'Only showing cells with non-zero expression'}),
                  y = 'Expression')
  if (length(genes) > 1) {
    out <- out + ggplot2::facet_wrap(~Var1, scale='free_y')
  }
  out
}

#' Extension of gcoreVlnPlot intended for facetting data for a single gene
#'
#' An extension of gcoreVlnPlot, focusing on a single gene but facetting
#' it by a metadata variable. The intended use case for this was to facet
#' by individual/capture to see if a single individual is driving an effect.
#' See `gencoreSC::gcoreVlnPlot()` for more details.
#'
#' Requires `ggforce` to be installed
#'
#' @param obj Seurat object
#' @param gene Gene to plot
#' @param facet_var Metadata variable to facet plot on
#' @param grouping_var Column of `obj@meta.data` to group data by.
#' @param groups Optional. Levels of `grouping_var` to include in the plot. Also used
#'  to specify order of levels.
#' @param subset_var  Optional. Column of `obj@meta.data` to subset on. Default
#'  is `seurat_clusters` so you could `subset` on a specific cluster.
#' @param subset Levels of `subset_var` to subset data to
#' @param filter_zeros Remove 0s from plot: default `TRUE`.
#' @param assay Assay to pull expression data from, default `RNA`
#'
#' @returns A ggplot
#' @export
#'
#' @examples
#' \dontrun{
#' ## Plots only the 'Intermediate' cells (as labeled in 'coarse_labels' column)
#' ## Groups results by 'Pre' and 'Post', as labeled in the 'stage' column'
#' ## Facets by individual
#'   gcoreVlnPlot_facetted(obj,
#'     gene = 'ISG15',
#'     facet_var = 'sample',
#'     subset = 'Intermediate', subset_var = 'coarse_labels',
#'     grouping_var = 'stage', groups = c('Pre', 'Post'))
#' }

gcoreVlnPlot_facetted <- function(obj,
                                  gene,
                                  facet_var,
                                  grouping_var,
                                  groups = NULL,
                                  subset = NULL,
                                  subset_var = NULL,
                                  filter_zeros = TRUE,
                                  assay = 'RNA') {
  if (!(facet_var %in% colnames(obj@meta.data))) {
    stop('Error: Value for facet_var not found as a object metadata column')
  }
  if (missing(grouping_var)) {
    stop('Use the "grouping_var" argument
         to specify a metadata variable to group observations')
    }
  if (missing(gene)) {
    stop('Specify gene(s) to plot data for')
  }
  if (!(gene %in% rownames(obj@assays[[assay]]@data))) {
    stop('Gene not found in object/assay')
  }

  if (!is.null(subset)) {
    if (!is.null(subset_var)) {
      if (subset_var %in% colnames(obj@meta.data)) {
        if (subset %in% unique(obj@meta.data[[subset_var]])) {
          obj <- obj[,obj@meta.data[[subset_var]] %in% subset]
        } else {
          stop('Error: subset value :"',subset,
               '" not present as a value of ',
               subset_var, 'column in object metadata.')
        }
      } else {
        stop('Error: subset_var value :"',subset_var,
             '" not present as a column in the object metadata.')
      }
    } else {
      stop('Error: Specify a variable to subset on via subset_var')
    }
  }

  mat_to_plot <- reshape2::melt(as.matrix(obj@assays[[assay]]@data)[gene,])
  mat_to_plot <- merge(mat_to_plot,
                       obj@meta.data %>%
                         as.data.frame() %>%
                         dplyr::select(.data[[grouping_var]], .data[[facet_var]]),
                       by = 0)
  if (grouping_var %in% colnames(obj@meta.data)) {
    if (!is.null(groups)) {
      if (all(groups %in% unique(obj@meta.data[[grouping_var]]))) {
        mat_to_plot <- mat_to_plot %>%
          filter(.data[[grouping_var]] %in% groups)
        mat_to_plot[[grouping_var]] <- factor(mat_to_plot[[grouping_var]],
                                              levels = groups)
      } else {
        stop('Error: Specified groups are not all present in the ',
             grouping_var, ' column of the object metadata.')
      }
    }  else {
      ## Groups aren't specified
      mat_to_plot[[grouping_var]] <- factor(mat_to_plot[[grouping_var]])
    }
  } else {
    stop('Error: grouping_var value :"',grouping_var,
         '" not present as a column in the object metadata.')
  }



  if (filter_zeros) {mat_to_plot <- mat_to_plot %>% filter(.data[['value']] != 0)}

  ggplot2::ggplot(mat_to_plot,
                         aes(x = .data[[grouping_var]],
                             y = .data[['value']],
                             color = .data[[grouping_var]])) +
    # geom_jitter() +
    ggplot2::geom_violin(draw_quantiles = 0.5) +
    ggforce::geom_sina(size = 0.01, alpha = 0.2) +
    ggplot2::theme_bw() +
    ggplot2::labs(caption = paste0(if(!is.null(subset)) {paste0('Showing expression of ', subset, ' cells\n')},
                                   if(!is.null(filter_zeros)) {'Only showing cells with non-zero expression'}),
                  y = 'Expression') +
    ggplot2::facet_wrap({{ facet_var }})
}
