#' Read BD data
#'
#' Reads in BD counts and metadata and returns Seurat obj
#'
#' @param counts_csv path to csv of UMI count table
#' @param metadata_csv path to metadata for UMI count table
#'
#' @return Seurat object
#'
#' @note
#' To do:
#'   - add option to not remove multiplets and undetermined droplets
#'   - allow flexible assays; RNA only, or additional arbitrary assays
#'   - add working example
#'
#' @export
readCountsBD <- function(counts_csv, metadata_csv) {
  files <- counts_csv
  counts <- read.table(files, skip = 0, sep = ",", header = TRUE, row.names = 1)
  tcounts <- data.frame(t(counts), check.names = FALSE)

  #read in sample tags
  tag_meta <- read.table(metadata_csv, skip = 0, sep = ",", header = TRUE, row.names = 1)

  #remove multiplets
  remove_multiplets <- tag_meta[tag_meta$Sample_Tag != "Multiplet",]
  remove_un <- remove_multiplets[remove_multiplets$Sample_Tag != "Undetermined",]

  #remove barcodes not in final set (excluding multiplets and undetermined)
  final <- tcounts[,colnames(tcounts) %in% row.names(remove_un)]

  #now separate out abseqs
  #will need to change when all abseqs read in
  ADT <- final[grepl(".pAbO$", rownames(final)),]
  RNA <- final[!grepl(".pAbO$", rownames(final)),]

  #shorten abseq names
  row.names(ADT) <- sub(".AMM.*$", "", row.names(ADT))

  # Create the Seurat object using the RNA data as the first assay
  s <- Seurat::CreateSeuratObject(counts = RNA, meta.data = remove_un)
  SeuratObject::Assays(s)

  # And then add the AbSeq data as a second assay
  adt_assay <- SeuratObject::CreateAssayObject(counts = ADT)
  s[["ADT"]] <- adt_assay
  show(s)

  return(s)
}
