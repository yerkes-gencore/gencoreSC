#' Plot genes with high degree of contamination correction
#'
#' Uses the results of `runSoupX()` to produce a plot of genes with a high
#'  degree of contamination correction through SoupX. This is defined as genes
#'  with a high portion of cell-level counts set to 0 after running SoupX.
#'  Ideally genes shown here should be strong cell-type markers,
#'  e.g. hemoglobins or immunoglobulins. For lightweight rendering
#'  and to focus on non-trivial signal, only genes found in a sufficient number
#'  of cells and with a sufficient degree of correction are plotted. You can
#'  also optionally just plot genes with a meaningful signal and exclude genes
#'  with only Ensembl IDs.
#'
#' @param sc SoupChannel object created from `runSoupX()`
#' @param ylim_min Minimum value of 'n'umber of original cells' to plot (y axis)
#' @param xlim_min Minimum value of 'portion of cells set to zero' to plot (x axis)
#' @param hide_ensembl Whether to include genes with `ens_pattern` in the plot
#' @param ens_pattern Regex pattern for identifying genes with no gene symbol.
#'  Only used to hide unknown genes from the plot if they aren't of interest
#'
#' @returns A plotly object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr drop_na
#' @importFrom plotly ggplotly
#' @importFrom Matrix rowSums
#'
#' @examples
#' \dontrun{
#' sc <- runSoupX(tod = tod,
#'   toc = toc,
#'   clus = clusters)
#' soupXGeneAdjustPlot(sc, hide_ensembl = TRUE)
#' }
plotSoupXGeneAdjustments <- function(sc,
                                     ylim_min = 10,
                                     xlim_min = 0.05,
                                     hide_ensembl = FALSE,
                                     ens_pattern = '^ENS'){
  cntSoggy = Matrix::rowSums(sc$toc > 0)
  cntStrained = Matrix::rowSums(sc$adjusted_counts > 0)
  ratio = (cntSoggy - cntStrained) / cntSoggy
  count_diff <- data.frame(proportion = ratio,
                           original = Matrix::rowSums(sc$toc > 0)) %>%
    tibble::rownames_to_column('Gene') %>%
    tidyr::drop_na()

  if (hide_ensembl){
    count_diff <- count_diff[!grepl(ens_pattern, count_diff$Gene),]
  }
  soup_plot <- ggplot2::ggplot(count_diff,
                              aes(x=.data$proportion,
                              y=.data$original,
                              text=.data$Gene,
                              alpha = 0.5)) +
    geom_point() +
    labs(x='Portion of cells set to zero',
         y='Original number of cells with expression',
         title = 'Genes most affected by SoupX correction') +
    theme_bw() +
    guides(alpha = 'none',
           label = 'none') +
    lims(x = c(xlim_min, 1)) +
    scale_y_log10(limits=c((ylim_min),NA))

  message(paste0('This plot focuses on genes found in at least ',
                 ylim_min, ' cells in the original dataset',
                 'and a minimum adjustment of ',
                 xlim_min, ' to reduce rendering burden.'))
  plotly::ggplotly(soup_plot, tooltip = 'text') %>%
    config(
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", 'lasso2d', 'pan2d', 'autoScale2d', 'zoom2d')
    )
}
