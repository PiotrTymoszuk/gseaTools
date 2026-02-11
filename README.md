# gseaTools

_Tools for management of gene signatures_

A development package providing tools for reading locally stored 
MSig Data Base files (https://www.gsea-msigdb.org/gsea/index.jsp) and converting 
them to a _tidyverse_ compatible format for easy search of the signatures by gene identifiers or signature names. 
Signature scores can be calculated with the user-provided gene expression data 
with the Gene Set Variance Analysis (GSVA) technique via the interface of 
[GSVA](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html) package.

## Terms of use

The package is available under a [GPL-3 license](https://github.com/PiotrTymoszuk/gseaTools/blob/main/LICENSE).

## Contact

The package maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).

## Acknowledgements

`gseaTools` uses tools provided by the 
[tidyverse](https://www.tidyverse.org/) package bundle, and R packages 
[stringi](https://stringi.gagolewski.com/), 
[utils](https://www.rdocumentation.org/packages/utils/versions/3.6.2), 
[GSVA](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html), 
[BiocParallel](https://www.bioconductor.org/packages/release/bioc/html/BiocParallel.html). 
Many thanks to their developers, maintainers, and contributors.

## Installation

You may install the package using `devtools`:

```r

devtools::install_github('PiotrTymoszuk/gseaTools')

```

## Basic usage

In the following example, we will compute single sample gene set enrichment analysis scores 
for synthetic chemokine signatures and simulated gene expression data, and Reactome Pathways 
gene signatures and breast cancer mRNA sequencing data. 
The Reactome Pathways gene signatures will be derived from .GMT files available from 
[MSig Database](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp).  
The pre-formatted breast cancer gene expression data come with 
the [microViz package](https://github.com/PiotrTymoszuk/microViz) package.

```r 

  ## required packages

  library(tidyverse)

  library(gseaTools)

  library(microViz) ### breast cancer gene expression

  library(BiocParallel) ### paralellization

  select <- dplyr::select ### resolution of possible function conflicts

```

```r

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

```
The chemokine gene signatures are defined as a simple list of vectors with 
official gene symbols. 
Database objects are created from .GMT files with `load_dbsig()` function; 
Reactome Pathways are selected with `search()`, which takes a regular expression #
as a search term. 

```r

  ## synthetic signatures of chemokines

  chemokine_signatures <- list(sign1 = c("CXCR4", "CXCR1", "CXCR3"),
                               sign2 = c("CXCR2", "CXCR5", "CXCR4"),
                               sign3 = c("CXCR2", "CXCR9"))

  ## databases from GMT files, extraction of Reactome Pathway signatures:
  ## a legacy and a newer version. Please use paths to your local .GMT files

  msig_db <- c(legacy = "msigdb.v7.4.symbols.gmt",
               recent = "msigdb.v2026.1.Hs.symbols.gmt") %>%
    map(load_dbsig)

  reactome_db <- msig_db %>%
    map(search, key = "sign_name", value = "^REACTOME")

```

```r

> head(reactome_db$legacy)

# A tibble: 6 × 3
  sign_name                                                                                   sign_link                       genes
  <chr>                                                                                       <chr>                           <lis>
1 REACTOME_INTERLEUKIN_6_SIGNALING                                                            http://www.gsea-msigdb.org/gse… <chr>
2 REACTOME_APOPTOSIS                                                                          http://www.gsea-msigdb.org/gse… <chr>
3 REACTOME_HEMOSTASIS                                                                         http://www.gsea-msigdb.org/gse… <chr>
4 REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS                                                    http://www.gsea-msigdb.org/gse… <chr>
5 REACTOME_MAPK3_ERK1_ACTIVATION                                                              http://www.gsea-msigdb.org/gse… <chr>
6 REACTOME_TRANSLESION_SYNTHESIS_BY_Y_FAMILY_DNA_POLYMERASES_BYPASSES_LESIONS_ON_DNA_TEMPLATE http://www.gsea-msigdb.org/gse… <chr>


```
The signature scores are computed with `calculate()`, which takes a list of gene identifier vectors or 
a database object as its first argument and gene/protein expression data as the second argument. 
Please note, that the expression data is provided in data frames, where genes are represented by columns and 
samples are represented by rows. 
If the input expression data frame has row names, they will be stored in the 
first column of the output data with the signature scores, named by default `sample_id`. 
By specifying `BPPARAM`, you may define and configure a parallel backend; by default the function runs in a serial manner.

```r

  ## scores of chemokine signatures, synthetic expression data

  chemokine_scores <- chemokine_signatures %>%
    calculate(data = chemokine_expression)

  ## scores of Reactome Pathway signatures, breast cancer expression
  ## parallel run via `snow`

  reactome_scores <- reactome_db %>%
    map(calculate,
        data = brca_expression,
        BPPARAM = SnowParam(workers = 7))

```

```r

> head(chemokine_scores)

# A tibble: 6 × 4
  sample_id   sign1  sign2 sign3
  <chr>       <dbl>  <dbl> <dbl>
1 1          0.167  -1     -0.5 
2 2          0.75    0.167  0.5 
3 3          0.0833 -0.667 -1   
4 4         -0.133   0.2    0.75
5 5          1       0.5   -1   
6 6          0      -0.5    0  

```

```r

> head(reactome_scores$legacy)

# A tibble: 6 × 1,605
  sample_id             REACTOME_INTERLEUKIN…¹ REACTOME_APOPTOSIS REACTOME_HEMOSTASIS REACTOME_INTRINSIC_P…² REACTOME_MAPK3_ERK1_…³
  <chr>                                  <dbl>              <dbl>               <dbl>                  <dbl>                  <dbl>
1 MBC-MBCProject_GvHkH…                  0.142            -0.226              -0.0959                -0.0765                -0.115 
2 MBC-MBCProject_N4srs…                 -0.465             0.0434             -0.167                  0.321                 -0.301 
3 MBC-MBCProject_N4srs…                 -0.366             0.134              -0.110                  0.150                  0.0523
4 MBC-MBCProject_4MF1F…                 -0.543            -0.0216             -0.241                  0.0973                -0.192 
5 MBC-MBCProject_4MF1F…                  0.142             0.107               0.149                  0.243                  0.222 
6 MBC-MBCProject_wAiri…                  0.216            -0.135               0.0211                -0.0195                -0.0840
# ℹ abbreviated names: ¹​REACTOME_INTERLEUKIN_6_SIGNALING, ²​REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS,
#   ³​REACTOME_MAPK3_ERK1_ACTIVATION
# ℹ 1,599 more variables: REACTOME_TRANSLESION_SYNTHESIS_BY_Y_FAMILY_DNA_POLYMERASES_BYPASSES_LESIONS_ON_DNA_TEMPLATE <dbl>,
#   REACTOME_RECOGNITION_OF_DNA_DAMAGE_BY_PCNA_CONTAINING_REPLICATION_COMPLEX <dbl>, REACTOME_TRANSLESION_SYNTHESIS_BY_POLH <dbl>,
#   REACTOME_RECOGNITION_AND_ASSOCIATION_OF_DNA_GLYCOSYLASE_WITH_SITE_CONTAINING_AN_AFFECTED_PURINE <dbl>,
#   REACTOME_DISPLACEMENT_OF_DNA_GLYCOSYLASE_BY_APEX1 <dbl>, REACTOME_POLB_DEPENDENT_LONG_PATCH_BASE_EXCISION_REPAIR <dbl>,
#   REACTOME_RESOLUTION_OF_AP_SITES_VIA_THE_MULTIPLE_NUCLEOTIDE_PATCH_REPLACEMENT_PATHWAY <dbl>, …
# ℹ Use `colnames()` to see all variable names

```
```r

> head(reactome_scores$recent)

# A tibble: 6 × 1,840
  sample_id      REACTOME_2_LTR_CIRCL…¹ REACTOME_ABACAVIR_ADME REACTOME_ABACAVIR_TR…² REACTOME_ABC_FAMILY_…³ REACTOME_ABC_TRANSPO…⁴
  <chr>                           <dbl>                  <dbl>                  <dbl>                  <dbl>                  <dbl>
1 MBC-MBCProjec…                 -0.749                -0.192                -0.192                -0.347                    0.0852
2 MBC-MBCProjec…                  0.487                -0.0348                0.00454               0.000824                 0.456 
3 MBC-MBCProjec…                  0.395                -0.302                -0.160                 0.268                    0.522 
4 MBC-MBCProjec…                  0.562                -0.403                -0.618                -0.0190                  -0.211 
5 MBC-MBCProjec…                 -0.424                 0.336                 0.520                 0.153                    0.273 
6 MBC-MBCProjec…                 -0.749                 0.488                 0.627                -0.218                    0.0312
# ℹ abbreviated names: ¹​REACTOME_2_LTR_CIRCLE_FORMATION, ²​REACTOME_ABACAVIR_TRANSMEMBRANE_TRANSPORT,
#   ³​REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT, ⁴​REACTOME_ABC_TRANSPORTERS_IN_LIPID_HOMEOSTASIS
# ℹ 1,834 more variables: REACTOME_ABC_TRANSPORTER_DISORDERS <dbl>,
#   REACTOME_ABERRANT_REGULATION_OF_MITOTIC_EXIT_IN_CANCER_DUE_TO_RB1_DEFECTS <dbl>,
#   REACTOME_ABERRANT_REGULATION_OF_MITOTIC_G1_S_TRANSITION_IN_CANCER_DUE_TO_RB1_DEFECTS <dbl>,
#   REACTOME_ABORTIVE_ELONGATION_OF_HIV_1_TRANSCRIPT_IN_THE_ABSENCE_OF_TAT <dbl>,
#   REACTOME_ACETYLCHOLINE_BINDING_AND_DOWNSTREAM_EVENTS <dbl>, …
# ℹ Use `colnames()` to see all variable names

```
