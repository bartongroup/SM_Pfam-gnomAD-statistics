Supplementary Table 4: Feature enrichment in UMDs vs. UMEs
================

``` r
library(broom)
suppressPackageStartupMessages(library(dplyr))
library(readr)
library(tidyr)


# Load and process the data
pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv <-
  read_csv("data/pfam-gnomAD-clinvar-pdb-colstats_c7c3e19.csv.gz",
           comment = "#")
```

    ## Rows: 1466478 Columns: 92
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (3): family, PDB_jury_column_rsa_unb, PDB_jury_column_rsa_bound
    ## dbl (88): column, frameshift_variant, frameshift_variant&splice_region_varia...
    ## lgl  (1): occupancy_gteq_threshold
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
avs <-
  pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv  # avs: alignments variants structure

# Drop unneccesary columns
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
  "PDB_protein_ligand_interactions",
  "PDB_total_ligand_interactions",
  "PDB_protein_protein_interactions",
  "PDB_total_protein_interactions",
  "PDB_biolip",
  "PDB_total_biolip",
  "PDB_jury_column_rsa_unb"
)
blacklist <- !names(avs) %in% whitelist
avs[, blacklist] <- NULL

# Rename variable for convenience
names(avs)[names(avs) == '_missense_all'] <- 'missense_all'

# Categories columns with MES and conservation thresholds
(
  avs <- avs
  
  # Categorise MES and Shenkin
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
  
  
  
  # Minimal filters for this notebook
  # Column filters
  |> filter(occupancy > 0)  # must have human sequences
  # # |> filter(PDB_sequences_with_contacts > 0)  # will also drop Pfams with no PDB data
  #
  # Filter Pfams that don't have at least one missense enriched and missense depleted position
  # Leaving this out attenuates the results but doesn't change any overall conclusion
  |> group_by(family) |> filter(any(mes_or < 1 & mes_p < 0.1)) |> ungroup()
  |> group_by(family) |> filter(any(mes_or > 1 & mes_p < 0.1)) |> ungroup()
  
  # Minor data cleaning steps
  |> mutate(PDB_jury_column_rsa_unb = as.factor(PDB_jury_column_rsa_unb))
  |> replace_na(
    list(
      PDB_protein_ligand_interactions = 0,
      PDB_biolip = 0,
      PDB_sequences_with_contacts = 0
    )
  )
)
```

    ## # A tibble: 167,171 × 26
    ##    family  column splice_…¹ synon…² occup…³ shenkin misse…⁴ mes_or mes_p occup…⁵
    ##    <chr>    <dbl>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl> <dbl>   <dbl>
    ##  1 PF00001     41         0       0       1    6.30       1  1.89  1         -70
    ##  2 PF00001     42         0       1       1    6.33       1  1.89  1         -70
    ##  3 PF00001     43         0       0       2    6.37       1  0.945 1         -69
    ##  4 PF00001     44         0       1       3    6.46       0  0     0.556     -68
    ##  5 PF00001     45         0       1       5    6.59       2  0.756 1         -66
    ##  6 PF00001     46         0       1       5    6.82       4  1.51  0.507     -66
    ##  7 PF00001     47         0       3       6    6.84       7  2.20  0.154     -65
    ##  8 PF00001     48         0       2       7    6.96       6  1.62  0.393     -64
    ##  9 PF00001     49         0       4       7    7.04       4  1.08  1         -64
    ## 10 PF00001     50         0       4       7    7.09       5  1.35  0.763     -64
    ## # … with 167,161 more rows, 16 more variables: occupancy_gteq_threshold <lgl>,
    ## #   shenkin_nrank <dbl>, mes_or_nrank <dbl>, mes_or_log <dbl>,
    ## #   CV_Pathogenic <dbl>, PDB_sequences_with_contacts <dbl>,
    ## #   PDB_protein_ligand_interactions <dbl>, PDB_total_ligand_interactions <dbl>,
    ## #   PDB_protein_protein_interactions <dbl>,
    ## #   PDB_total_protein_interactions <dbl>, PDB_biolip <dbl>,
    ## #   PDB_total_biolip <dbl>, PDB_jury_column_rsa_unb <fct>, mes_class <fct>, …

``` r
# Unload unprocessed dataset
rm(pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv)

# Print some basic statistics
print(paste('N Pfam domains =', n_distinct(avs$family)))
```

    ## [1] "N Pfam domains = 697"

``` r
print(paste('N positions =', avs |> nrow()))
```

    ## [1] "N positions = 167171"

``` r
print(paste('N positions with human =', avs |> filter(occupancy > 0) |> nrow()))
```

    ## [1] "N positions with human = 167171"

``` r
(
  avs
  |> filter(xclass %in% c('umd', 'ume'))
  |> group_by(xclass)
  |> summarise(
    n_biolip = sum(PDB_biolip),
    seq_pdb_cov = sum(PDB_sequences_with_contacts),
    n_positions = n(),
    n_human_residues = sum(occupancy)  # for additional context
  )
  |> as.data.frame()
  -> umd_ume_summary
)
```

    ##   xclass n_biolip seq_pdb_cov n_positions n_human_residues
    ## 1    umd      655       37311        1278            60352
    ## 2    ume     1203      116520        2763           233112

``` r
rownames(umd_ume_summary) <- umd_ume_summary[['xclass']]
test_df <- umd_ume_summary[c('umd', 'ume'), c('n_biolip', 'seq_pdb_cov')]
test_df
```

    ##     n_biolip seq_pdb_cov
    ## umd      655       37311
    ## ume     1203      116520

``` r
fisher.test(test_df)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  test_df
    ## p-value < 2.2e-16
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  1.542572 1.873078
    ## sample estimates:
    ## odds ratio 
    ##   1.700345

``` r
(
  avs
  |> filter(xclass %in% c('umd', 'ume'))
  |> group_by(xclass)
  |> summarise(
    n_biolip = sum(PDB_biolip),
    seq_pdb_cov = sum(PDB_sequences_with_contacts)
  )
  |> ungroup()
  |> arrange(xclass)  # umd > ume
  |> select(n_biolip, seq_pdb_cov)
  |> summarise(broom::tidy(fisher.test(matrix(
    c(n_biolip, seq_pdb_cov), ncol = 2
  ))))
)
```

    ## # A tibble: 1 × 6
    ##   estimate  p.value conf.low conf.high method                            alter…¹
    ##      <dbl>    <dbl>    <dbl>     <dbl> <chr>                             <chr>  
    ## 1     1.70 6.46e-26     1.54      1.87 Fisher's Exact Test for Count Da… two.si…
    ## # … with abbreviated variable name ¹​alternative

``` r
test_umd_vs_ume <- function(.data, feature, background) {
  (
    .data
    |> filter(xclass %in% c('umd', 'ume'))
    |> group_by(xclass)
    |> summarise(
      n_feature = sum({
        {
          feature
        }
      }, na.rm = TRUE),
      n_background = sum({
        {
          background
        }
      }, na.rm = TRUE)
    )
    |> ungroup()
    |> arrange(xclass)  # umd > ume
    |> select(n_feature, n_background)
    |> summarise(broom::tidy(fisher.test(
      matrix(c(n_feature, n_background - n_feature), ncol = 2)
    )))
  )
}

# BioLiP
test_umd_vs_ume(avs, PDB_biolip, PDB_sequences_with_contacts)
```

    ## # A tibble: 1 × 6
    ##   estimate  p.value conf.low conf.high method                            alter…¹
    ##      <dbl>    <dbl>    <dbl>     <dbl> <chr>                             <chr>  
    ## 1     1.71 1.58e-26     1.55      1.89 Fisher's Exact Test for Count Da… two.si…
    ## # … with abbreviated variable name ¹​alternative

``` r
# All protein-ligand
test_umd_vs_ume(avs,
                PDB_protein_ligand_interactions,
                PDB_sequences_with_contacts)
```

    ## # A tibble: 1 × 6
    ##   estimate  p.value conf.low conf.high method                            alter…¹
    ##      <dbl>    <dbl>    <dbl>     <dbl> <chr>                             <chr>  
    ## 1     1.27 2.67e-26     1.21      1.32 Fisher's Exact Test for Count Da… two.si…
    ## # … with abbreviated variable name ¹​alternative

``` r
# Protein-protein interfaces
test_umd_vs_ume(avs,
                PDB_protein_protein_interactions,
                PDB_sequences_with_contacts)
```

    ## # A tibble: 1 × 6
    ##   estimate     p.value conf.low conf.high method                         alter…¹
    ##      <dbl>       <dbl>    <dbl>     <dbl> <chr>                          <chr>  
    ## 1     1.10 0.000000106     1.06      1.14 Fisher's Exact Test for Count… two.si…
    ## # … with abbreviated variable name ¹​alternative

``` r
# ClinVar pathogenic enrichment relative to gnomAD missense
test_umd_vs_ume(avs, CV_Pathogenic, missense_all)
```

    ## # A tibble: 1 × 6
    ##   estimate  p.value conf.low conf.high method                            alter…¹
    ##      <dbl>    <dbl>    <dbl>     <dbl> <chr>                             <chr>  
    ## 1     4.75 3.72e-37     3.79      5.91 Fisher's Exact Test for Count Da… two.si…
    ## # … with abbreviated variable name ¹​alternative

``` r
apply_test_over_rsa_classes <-
  function(.data, feature, background) {
    rsa_classes <- c('Core', 'Part. Exposed', 'Surface')
    do.call(rbind,
            lapply(list(1:3, 1, 2, 3), function(rsa_select) {
              test_umd_vs_ume(.data |> filter(PDB_jury_column_rsa_unb %in% rsa_classes[rsa_select]),
                              {{feature}},
                              {{background}})
            }))
  }

by_rsa_biolip <- apply_test_over_rsa_classes(avs, PDB_biolip, PDB_sequences_with_contacts)
by_rsa_protein <- apply_test_over_rsa_classes(avs, PDB_protein_protein_interactions, PDB_sequences_with_contacts)
by_rsa_pathog <- apply_test_over_rsa_classes(avs, CV_Pathogenic, missense_all)

by_rsa_biolip
```

    ## # A tibble: 4 × 6
    ##   estimate  p.value conf.low conf.high method                            alter…¹
    ##      <dbl>    <dbl>    <dbl>     <dbl> <chr>                             <chr>  
    ## 1    1.59  6.41e-16    1.42       1.77 Fisher's Exact Test for Count Da… two.si…
    ## 2    0.939 7.78e- 1    0.641      1.38 Fisher's Exact Test for Count Da… two.si…
    ## 3    2.20  2.12e-12    1.75       2.79 Fisher's Exact Test for Count Da… two.si…
    ## 4    1.37  1.01e- 4    1.17       1.60 Fisher's Exact Test for Count Da… two.si…
    ## # … with abbreviated variable name ¹​alternative

``` r
by_rsa_protein
```

    ## # A tibble: 4 × 6
    ##   estimate  p.value conf.low conf.high method                            alter…¹
    ##      <dbl>    <dbl>    <dbl>     <dbl> <chr>                             <chr>  
    ## 1     1.14 4.54e-11    1.09       1.18 Fisher's Exact Test for Count Da… two.si…
    ## 2     1.14 1.92e- 1    0.937      1.39 Fisher's Exact Test for Count Da… two.si…
    ## 3     1.27 7.25e- 8    1.16       1.39 Fisher's Exact Test for Count Da… two.si…
    ## 4     1.50 7.89e-65    1.44       1.57 Fisher's Exact Test for Count Da… two.si…
    ## # … with abbreviated variable name ¹​alternative

``` r
by_rsa_pathog
```

    ## # A tibble: 4 × 6
    ##   estimate  p.value conf.low conf.high method                            alter…¹
    ##      <dbl>    <dbl>    <dbl>     <dbl> <chr>                             <chr>  
    ## 1     6.41 4.98e-47     5.08      8.08 Fisher's Exact Test for Count Da… two.si…
    ## 2     5.05 8.59e- 5     2.06     13.6  Fisher's Exact Test for Count Da… two.si…
    ## 3     4.90 3.46e-10     2.88      8.55 Fisher's Exact Test for Count Da… two.si…
    ## 4     6.48 1.48e-27     4.79      8.69 Fisher's Exact Test for Count Da… two.si…
    ## # … with abbreviated variable name ¹​alternative
