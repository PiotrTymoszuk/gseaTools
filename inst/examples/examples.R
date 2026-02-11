# Package tests.

  library(tidyverse)

  library(gseaTools)

  library(microViz) ### breast cancer gene expression

  library(BiocParallel) ### paralellization

  select <- dplyr::select ### resolution of possible function conflicts

# expression data for testing -----------

  ## synthetic gene expression data

  chemokine_expression <- tibble(CXCR1 = rnorm(100, 20),
                                 CXCR2 = rnorm(100, 20, 15),
                                 CXCR3 = rnorm(100, 17, 5),
                                 CXCR4 = rnorm(100, 45, 5),
                                 CXCR9 = rchisq(100, 2),
                                 GPX4 = rnorm(100))

  ## whole transcriptome expression in breast cancer:
  ## only gene expression information without clinical/pathological data

  brca_expression <- microViz::brca %>%
    column_to_rownames("sample_id") %>%
    select(-patient_id, -timepoint, -metastasis, -histology, -er_status)

# gene set containers ----------

  ## synthetic signatures of chemokines

  chemokine_signatures <- list(sign1 = c("CXCR4", "CXCR1", "CXCR3"),
                               sign2 = c("CXCR2", "CXCR5", "CXCR4"),
                               sign3 = c("CXCR2", "CXCR9"))

  ## databases from GMT files, extraction of Reactome Pathway signatures:
  ## a legacy and a newer version

  msig_db <- c(legacy = "./inst/examples/msigdb.v7.4.symbols.gmt",
               recent = "./inst/examples/msigdb.v2026.1.Hs.symbols.gmt") %>%
    map(load_dbsig)

  reactome_db <- msig_db %>%
    map(search, key = "sign_name", value = "^REACTOME")

# calculation of ssGSEA scores for the whole transcriptome -------

  ## scores of chemokine signatures, synthetic expression data

  chemokine_scores <- chemokine_signatures %>%
    calculate(data = chemokine_expression)

  ## scores of Reactome Pathway signatures, breast cancer expression
  ## parallel run via `snow`

  reactome_scores <- reactome_db %>%
    map(calculate,
        data = brca_expression,
        BPPARAM = SnowParam(workers = 7))

# END -----
