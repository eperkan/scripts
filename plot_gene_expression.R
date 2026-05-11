#' Plot gene expression across sample groups
#'
#' Create a boxplot with jittered sample-level points showing normalized
#' expression of a single gene across sample groups.
#'
#' @param count_matrix A numeric matrix or data frame with genes in rows and
#'   samples in columns. Data should be normalized before plotting. Row names
#'   must contain gene symbols or Ensembl gene IDs, and column names must
#'   contain sample IDs matching the row names of `sample_data`.
#' @param sample_data A data frame containing sample metadata. Row names must
#'   contain sample IDs matching the column names of `count_matrix`. Must
#'   include the column specified by `group`.
#' @param gene_of_interest A character string giving the gene symbol or Ensembl
#'   gene ID to plot. Must match a row name in `count_matrix`.
#' @param group A character string giving the column name in `sample_data` used
#'   for grouping samples on the x-axis.
#'
#' @return A `ggplot` object showing gene expression by group.
#'
#' @examples
#' PlotGeneExpression(
#'   count_matrix = counts_vst,
#'   sample_data = metadata,
#'   gene_of_interest = "TP53",
#'   group = "tumor_type"
#' )
PlotGeneExpression <- function(count_matrix, sample_data, gene_of_interest, group) {
        
        if (!group %in% colnames(sample_data)) {
                stop(paste0("'", group, "' is not a valid column name in sample_data."))
        }
        
        if (!gene_of_interest %in% rownames(count_matrix)) {
                stop(paste0("'", gene_of_interest, "' is not present in count_matrix."))
        }
        
        if (is.null(colnames(count_matrix))) {
                stop("count_matrix must have column names corresponding to sample IDs.")
        }
        
        if (is.null(rownames(sample_data))) {
                stop("sample_data must have row names corresponding to sample IDs.")
        }
        
        if (anyDuplicated(colnames(count_matrix))) {
                stop("count_matrix contains duplicated sample IDs.")
        }
        
        if (anyDuplicated(rownames(sample_data))) {
                stop("sample_data contains duplicated sample IDs.")
        }
        
        if (!setequal(colnames(count_matrix), rownames(sample_data))) {
                stop("Sample IDs in count_matrix columns and sample_data row names do not match.")
        }
        
        idx <- match(colnames(count_matrix), rownames(sample_data))
        sample_data <- sample_data[idx, , drop = FALSE]
        
        if (!identical(colnames(count_matrix), rownames(sample_data))) {
                stop("Sample order still does not match after reordering sample_data.")
        }
        
        gene_exp <- as.numeric(count_matrix[gene_of_interest, ])
        sample_data$expression <- gene_exp
        
        tmp_plot <- ggpubr::ggboxplot(
                data = sample_data,
                x = group,
                y = "expression",
                color = group,
                palette = "npg",
                title = gene_of_interest,
                add = "jitter",
                xlab = group,
                ylab = "Normalized mRNA expression",
                ggtheme = ggplot2::theme_bw()
        )
        
        return(tmp_plot)
}