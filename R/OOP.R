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

    if(!tibble::is_tibble(x)) {

      stop('Please provide a valid data tibble.',
           call. = FALSE)

    }

    structure(x,
              class = c('dbsig', 'tbl_df', 'tbl', 'data.frame'))

  }

# Class checking -----

#' Check for the dbsig class.
#'
#' @description Checks is an object is an instance of the dbsig class.
#' @param x an object.
#' @return a logical value.
#' @export

  is_dbsig <- function(x) {

    any(class(x) == 'dbsig')

  }

# Generics -----

#' Calculate object statistic.
#'
#' @description Calculate object statistics.
#' @param x an object.
#' @param ... arguments passed to the methods.
#' @export

  calculate <- function(x, ...) {

    UseMethod('calculate', x)

  }

#' Search for a database entry.
#'
#' @description Search for a database entry.
#' @param x an object.
#' @param ... arguments passed to the methods.
#' @export

  search <- function(x, ...) {

    UseMethod('search', x)

  }

# Extraction methods ------

#' Search in the dbsig database with a text expression.
#'
#' @description Filters entries of a dbsig database with a regular expression.
#' The fields available for searching are 'genes', 'sign_name'
#' (official signature name) and 'sign_link' (signature URL).
#' @return a dbsig object.
#' @param x a dbsig object.
#' @param key the dbsig database field, see description. 'sign_name' by default.
#' @param value a regular expression used to filter the database.
#' @export search.dbsig
#' @export

  search.dbsig <- function(x,
                           key = c('sign_name', 'genes', 'sign_link'),
                           value = '.*') {

    ## entry control

    stopifnot(gseaTools::is_dbsig(x))

    key <- match.arg(key[1], c('sign_name', 'genes', 'sign_link'))

    ## filtering

    if(key %in% c('sign_name', 'sign_link')) {

      dplyr::filter(x, stringi::stri_detect(.data[[key]], regex = value))

    } else if(key == 'genes') {

      pos_records <- purrr::map_lgl(x[['genes']],
                                    ~any(stringi::stri_detect(.x, regex = value)))

      dplyr::filter(x, pos_records)

    }

  }

# Score calculation ------

#' Calculate GSVA scores.
#'
#' @description Calculates multi-gene scores via \code{\link[GSVA]{gsva}}
#' function, given a character vector with gene names and a data frame with
#' expression values.
#' @details See: \code{\link[GSVA]{gsva}}.
#' @param x a list with gene name/identifier vectors.
#' @param data a data frame with the expression values.
#' @param ... extra arguments passed to \code{\link[GSVA]{gsva}}.
#' @return A data frame with GSVA enrichment scores:
#' rows are samples, columns are scores.
#' @export calculate.default
#' @export

  calculate.default <- function(x, data, ...) {

    ## expression

    expr_mtx <- gseaTools::as_exprs_matrix(data)

    ## calculation

    scores <- GSVA::gsva(expr = expr_mtx,
                         gset.idx.list = x, ...)

    ## formatting

    scores <- as.data.frame(t(scores))

    tibble::as_tibble(scores, .name_repair = 'universal')

  }

#' Calculate GSVA scores.
#'
#' @description Calculates multi-gene scores via \code{\link[GSVA]{gsva}}
#' function, given a character vector with gene names and a data frame with
#' expression values.
#' @details See: \code{\link[GSVA]{gsva}}.
#' @param x a dbsig object.
#' @inheritParams calculate.default
#' @return A data frame with GSVA enrichment scores:
#' rows are samples, columns are scores.
#' @export calculate.dbsig
#' @export

  calculate.dbsig <- function(x, data, ...) {

    stopifnot(gseaTools::is_dbsig(x))

    gene_container <- rlang::set_names(x[['genes']],
                                       x[['sign_name']])

    gseaTools::calculate.default(x = gene_container, data = data, ...)

  }

# END -----
