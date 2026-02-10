# Package tests.

  library(tidyverse)
  library(stringi)

  library(gseaTools)
  library(microViz)

  library(GSVA)
  library(GSVAdata)

  select <- dplyr::select

# tests --------

  test_data <- tibble(CXCR1 = rnorm(100, 20),
                      CXCR2 = rnorm(100, 20, 15),
                      CXCR3 = rnorm(100, 17, 5),
                      CXCR4 = rnorm(100, 45, 5),
                      CXCR9 = rchisq(100, 2),
                      GPX4 = rnorm(100))

  test_container <- list(sign1 = c('CXCR4', 'CXCR1', 'CXCR3'),
                         sign2 = c('CXCR2', 'CXCR5', 'CXCR4'),
                         sign3 = c('CXCR2', 'CXCR9'))

  as_exprs_matrix(test_data)

  test_db <- load_dbsig(path = './inst/examples/msigdb.v7.4.symbols.gmt')

  test_res <- filter(test_db, stri_detect(sign_name, fixed = 'chr'))

  test_res <- search(test_db, key = 'gene', value = 'GPX4')

  is_dbsig(select(test_db, sign_name))

  is_dbsig(test_db)

  is_dbsig(test_res)

  calculate(test_container, data = test_data)

  calculate(test_res[1:100, ], data = test_data)

# calculation of ssGSEA scores for the whole transcriptome -------

  reactome_db <- test_db %>%
    search(key = "sign_name", value = "^REACTOM")

  brca_expression <- microViz::brca %>%
    column_to_rownames("sample_id") %>%
    select(-patient_id, -timepoint, -metastasis, -histology, -er_status)

  reactome_scores <- calculate(reactome_db[1:50, ], brca_expression)

# END -----
