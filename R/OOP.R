# S3 OOP interface.

# Constructor ------

#' Create a dbsig object from a data tibble.
#'
#' @description Converts a data tibble into a dbsig object.
#' @details Technically, the object inherits most of the properties and methods
#' from the tibble and data.frame classes.
#' @param x a tibble.
#' @return an object of class dbsig.
#' @export

  dbsig <- function(x) {

    if(!is_tibble(x)) {

      stop('Please provide a valid data tibble.',
           call. = FALSE)

    }

    structure(x, class = c('dbsig', 'tbl_df', 'tbl', 'data.frame'))

  }

# Class checking -----

#' Check for the dbsig class.
#'
#' @description Checks is an object is an instance of the dbsig class.
#' @param x an object.
#' @return a logical value.
#' @export

  is_dbsig <- function(x) any(class(x) == 'dbsig')

# Generics -----

#' Search for a database entry.
#'
#' @description Filters entries of a dbsig database with a regular expression.
#' The fields available for searching are 'genes', 'sign_name'
#' (official signature name) and 'sign_link' (signature URL).
#'
#' @param x a dbsig object.
#' @param key the dbsig database field, see description. 'sign_name' by default.
#' @param value a regular expression used to filter the database.
#' @param ... arguments passed to the methods.
#'
#' @return a dbsig object.
#'
#' @export

  search <- function(x, ...) UseMethod('search')

#' @rdname search
#' @export search.dbsig
#' @export

  search.dbsig <- function(x,
                           key = c('sign_name', 'genes', 'sign_link'),
                           value = '.*', ...) {

    ## entry control

    stopifnot(is_dbsig(x))

    key <- match.arg(key[1], c('sign_name', 'genes', 'sign_link'))

    ## filtering

    if(key %in% c('sign_name', 'sign_link')) {

      filter(x, stri_detect(.data[[key]], regex = value))

    } else if(key == 'genes') {

      pos_records <- map_lgl(x[['genes']],
                             ~any(stri_detect(.x, regex = value)))

      filter(x, pos_records)

    }

  }

# Score calculation ------

#' Calculate GSVA scores.
#'
#' @description Calculates multi-gene scores via \code{\link[GSVA]{gsva}}
#' function, given a character vector with gene names and a data frame with
#' expression values.
#' @details See: \code{\link[GSVA]{gsva}}.
#'
#' @param x a list with gene name/identifier vectors or a dbsig object.
#' @param data a data frame with the expression values.
#' @param id_col name of the column to store identifiers/names of the samples
#' in the output data frames.
#' @param BPPARAM An object of class BiocParallelParam
#' (https://www.bioconductor.org/packages//release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.html)
#' with parameters of paralellization used by \code{\link[GSVA]{gsva}}.
#' @param ... extra arguments passed to \code{\link[GSVA]{gsvaParam}}.
#' @return A data frame with GSVA enrichment scores:
#' rows are samples, columns are scores.
#'
#'
#'
#' @export calculate.default
#' @export

  calculate.default <- function(x,
                                data,
                                id_col = "sample_id",
                                BPPARAM = SerialParam(progressbar = TRUE),
                                ...) {

    ## expression

    expr_mtx <- as_exprs_matrix(data)

    ## calculation

    gsva_input <- gsvaParam(exprData = expr_mtx,
                            geneSets = x, ...)

    scores <- gsva(gsva_input, BPPARAM = BPPARAM)

    ## formatting

    scores <- as.data.frame(t(scores))

    if(!is.null(rownames(scores))) {

      scores <- rownames_to_column(scores, id_col)

    }

    as_tibble(scores, .name_repair = 'universal')

  }

#' @rdname calculate.default
#' @export calculate.dbsig
#' @export

  calculate.dbsig <- function(x, data, id_col = "sample_id", ...) {

    stopifnot(is_dbsig(x))

    gene_container <- set_names(x[['genes']], x[['sign_name']])

    calculate.default(x = gene_container,
                      data = data,
                      id_col = id_col, ...)

  }

# END -----
