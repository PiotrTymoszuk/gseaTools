# Functions for loading and management of the dbsig databases.

# Database loading from the disc -----

#' Load a GSEA database.
#'
#' @description Loads a GSEA database from a local file and converts it
#' to a data frame of the dbsig class.
#' @return a dbsig object.
#' @param path a path to the .gmt file containing GSEA signatures.
#' @export

  load_dbsig <- function(path) {

    ## entry control

    if(!stringi::stri_detect(path, fixed = '.gmt')) {

      stop('Please provide a path to a valid .gmt file', call. = FALSE)

    }

    ## benchmarking

    start_time <- Sys.time()
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    ## reading and formatting

    message('Reading the signature database')

    db <- tibble::as_tibble(utils::read.csv(path,
                                            header = FALSE,
                                            sep = '\t'))

    message(paste('Clearing', nrow(db), 'database records'))

    gene_lst <- purrr::map(1:nrow(db),
                           ~gseaTools:::empty_null_(db[.x, 3:ncol(db)]))

    db <- rlang::set_names(db[1:2],
                           c('sign_name', 'sign_link'))

    gseaTools::dbsig(dplyr::mutate(db, genes = gene_lst))

  }
