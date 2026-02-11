# Utilities.

#' Filter out empty strings.
#'
#' @description Removes empty strings from a character vector.
#' @return a character vector.
#' @param x a character vector.

  empty_null_ <- function(x) {

    x <- unlist(x)

    x <- x[!is.na(x)]

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
#' gene name vectors.
#' @export

  as_exprs_matrix <- function(data) {

    ## entry control

    if(!is.data.frame(data)) {

      stop('The data argument needs to be a data frame.', call. = FALSE)

    }

    ## output

    t(data)

  }

# END ----
