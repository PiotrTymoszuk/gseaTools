# Utilities.

#' Filter out empty strings.
#'
#' @description Removes empty strings from a character vector.
#' @return a character vector.
#' @param x a character vector.

  empty_null_ <- function(x) {

    x <- unlist(x)

    x <- x[x != '']

    unique(unname(x))

  }

#' Create an expression matrix.
#'
#' @description Creates an matrix with the gene expression values in a format
#' accepted by the \code{\link[GSVA]{gsva}} function.
#' @return an expression matrix: gene names are presented in row names,
#' columns represent the samples.
#' @param data a data frame, column names should contain gene names or
#' identifiers.
#' @param container a dbsig database object or a list with the
#' gene name vectors.
#' @export

  as_exprs_matrix <- function(data, container) {

    ## entry control

    if(!is.data.frame(data)) {

      stop('The data argument needs to be a data frame.', call. = FALSE)

    }

    if(!gseaTools::is_dbsig(container)) {

      if(!is.list(container)) {

        stop('Please provide a valid dbsig object or a list with gene names.',
             call. = FALSE)

      }

    }

    ## unique gene names/identifiers

    if(gseaTools::is_dbsig(container)) {

      all_genes <- purrr::reduce(container[['genes']], union)

    } else {

      all_genes <- purrr::reduce(container, union)

    }

    missing_genes <- all_genes[!all_genes %in% names(data)]

    if(length(missing_genes) > 0) {

      if(length(missing_genes) > 10) {

        miss_txt <- paste0(paste(missing_genes[1:10], collapse = ', '), ', ...')

      } else {

        miss_txt <- paste(missing_genes, collapse = ', ')

      }

      warning(paste(length(missing_genes),
                    'genes missing from data:',
                    miss_txt),
              call. = FALSE)

    }

    ## output

    t(data[names(data) %in% all_genes])

  }

# END ----
