Table 1: The Conservation Plane
================

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#notebook-resources" id="toc-notebook-resources">Notebook
  resources</a>
  - <a href="#libraries" id="toc-libraries">Libraries</a>
  - <a href="#dataset" id="toc-dataset">Dataset</a>
- <a href="#table-1-the-conservation-plane-table-data"
  id="toc-table-1-the-conservation-plane-table-data">Table 1: The
  Conservation Plane (table data)</a>

# Introduction

This R markdown document tabulates summary data for the different
conservation plane categories. The dataset used for this analysis is
obtained from our pre-calculated MSA column statistics for Pfam,
included in this repository. The output includes data used in Table 1 in
the manuscript.

# Notebook resources

## Libraries

``` r
suppressPackageStartupMessages(library(tidyverse))
```

## Dataset

``` r
# Read and process dataset
pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv <-
  read_csv(
    "data/pfam-gnomAD-clinvar-pdb-colstats_c7c3e19.csv.gz",
    comment = "#",
    show_col_types = FALSE
  )

# Processing
avs <-
  pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv  # avs: alignments variants structure

# Filter unnecessary variables
whitelist <- c(
  "family",
  "column",
  "splice_region_variant&synonymous_variant",
  "synonymous_variant",
  "occupancy",
  "shenkin",
  "_missense_all",
  "mes_or",
  "mes_p",
  "occupancy_minus_threshold",
  "occupancy_gteq_threshold",
  "shenkin_nrank",
  "mes_or_nrank",
  "mes_or_log",
  "CV_Pathogenic",
  "PDB_sequences_with_contacts",
  "PDB_jury_column_rsa_unb"
)
blacklist <- !names(avs) %in% whitelist
avs[, blacklist] <- NULL

# Better variable name format
names(avs)[names(avs) == '_missense_all'] <- 'missense_all'

# Calculate synonymous all
avs[, 'synonymous_all'] <-
  avs[, 'splice_region_variant&synonymous_variant'] +
  avs[, 'synonymous_variant']

# Clean-up
rm(pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv)
```

Below is where I define the conservation plane categories.

``` r
# Apply thresholds to categorise MES and Shenkin and the combined class
avs <- (
  avs
  |> mutate(mes_class = factor(
    ifelse(mes_p < 0.1, ifelse(mes_or < 1, 'md', 'me'), 'mn'), levels = c('md', 'mn', 'me')
  ))
  |> mutate(cons_class = factor(
    ifelse(shenkin_nrank > 0.5, 'u', 'c'), levels = c('u', 'c')
  ))
  |> mutate(xclass = factor(
    interaction(cons_class, mes_class, sep = ''),
    levels = c('cmd', 'cmn', 'cme', 'umd', 'umn', 'ume')
  ))
)
```

``` r
# Filters
apply_filters <- function(df,
                          occupancy = FALSE,
                          has_RSA = FALSE,
                          has_depleted = FALSE,
                          has_enriched = FALSE,
                          has_PDB = FALSE,
                          has_Pathogenic = FALSE) {
  if (occupancy)
    df <- df |> filter(occupancy > 0)  # must have human sequences
  if (has_RSA)
    df <-
      df |> filter(PDB_jury_column_rsa_unb %in% c('Core', 'Part. Exposed', 'Surface'))
  
  # Heuristic for well-powered Pfams
  if (has_depleted)
    df <-
      df |>
      group_by(family) |>
      filter(any(mes_or < 1 & mes_p < 0.1)) |> ungroup()  # >1 depleted
  if (has_enriched)
    df <-
      df |>
      group_by(family) |>
      filter(any(mes_or > 1 & mes_p < 0.1)) |> ungroup()  # >1 enriched
  
  # Has PDB data
  if (has_PDB)
    df <-
      df |>
      group_by(family) |>
      filter(any(PDB_sequences_with_contacts >= 1)) |> ungroup()
  
  # Has ClinVar Pathogenic
  if (has_Pathogenic)
    df <- df |> group_by(family) |> filter(any(CV_Pathogenic > 0)) |> ungroup()
  return(df)
}
```

# Table 1: The Conservation Plane (table data)

This table reports summary statistics for each category defined by the
classification of Pfam domains by population missense constraint and
evolutionary divergence.

``` r
avs |>
  apply_filters(occupancy = TRUE,
                has_RSA = FALSE,
                has_depleted = TRUE,
                has_enriched = FALSE,
                has_PDB = FALSE,
                has_Pathogenic = FALSE) |>
  group_by(xclass) |>
  summarise(Families = n_distinct(family),
            Columns = n(),
            Residues = sum(occupancy),
            Missense = sum(missense_all)) ->
  avs_summary
avs_summary
```

    ## # A tibble: 6 × 5
    ##   xclass Families Columns Residues Missense
    ##   <fct>     <int>   <int>    <dbl>    <dbl>
    ## 1 cmd         683    3774   304635    85512
    ## 2 cmn         746  100754  1193815   529116
    ## 3 cme         581    2476    84073    65998
    ## 4 umd         518    1312    60665    17249
    ## 5 umn         746   60762  1405971   674276
    ## 6 ume         630    2763   233112   146168

``` r
avs_summary |>
  summarise(Families = max(Families),
            across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) ->
  avs_summary_totals
avs_summary_totals
```

    ## # A tibble: 1 × 4
    ##   Families Columns Residues Missense
    ##      <int>   <int>    <dbl>    <dbl>
    ## 1      746  171841  3282271  1518319
