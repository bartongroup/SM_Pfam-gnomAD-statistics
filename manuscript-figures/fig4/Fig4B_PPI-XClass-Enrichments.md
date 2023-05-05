Figure 4B: Enrichment of protein-protein interactions
================

- <a href="#notebook-resources" id="toc-notebook-resources">Notebook
  resources</a>
  - <a href="#libraries" id="toc-libraries">Libraries</a>
  - <a href="#dataset" id="toc-dataset">Dataset</a>
- <a href="#calculate-and-plot-odds-ratios"
  id="toc-calculate-and-plot-odds-ratios">Calculate and plot odds
  ratios</a>
  - <a
    href="#odds-ratios-for-protein-protein-interaction-enrichment-in-x-classes"
    id="toc-odds-ratios-for-protein-protein-interaction-enrichment-in-x-classes">Odds
    ratios for protein-protein interaction enrichment in X-Classes</a>
  - <a
    href="#odds-ratios-for-protein-protein-interaction-enrichment-in-x-class--rsa"
    id="toc-odds-ratios-for-protein-protein-interaction-enrichment-in-x-class--rsa">Odds
    ratios for protein-protein interaction enrichment in X-Class * RSA</a>
- <a href="#figure-4b-odds-ratio-plot"
  id="toc-figure-4b-odds-ratio-plot">Figure 4B: Odds ratio plot</a>

# Notebook resources

## Libraries

``` r
suppressPackageStartupMessages(library(tidyverse))
library(patchwork)
```

## Dataset

``` r
# Read and process dataset
pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv <- read_csv("data/pfam-gnomAD-clinvar-pdb-colstats_c7c3e19.csv.gz",
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
avs <- pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv  # avs: alignments variants structure

# Filter unnecessary variables
whitelist <- c(
  "family",
  "column",
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
  "PDB_sequences_with_contacts",
  "PDB_protein_protein_interactions",
  "PDB_total_protein_interactions",
  "PDB_non_self_ppis",
  "PDB_total_non_self_ppis",
  "PDB_surface_ppis",
  "PDB_total_surface_ppis",
  "PDB_surface_homo_ppis",
  "PDB_total_surface_homo_ppis",
  "PDB_homo_ppis",
  "PDB_total_homo_ppis",
  "PDB_intra_domain_x_interactions",
  "PDB_total_intra_domain_x_interactions",
  "PDB_intra_domain_domain_interactions",
  "PDB_total_intra_domain_domain_interactions",
  "PDB_jury_column_rsa_unb",
  "PDB_jury_column_rsa_bound"
)
blacklist <- !names(avs) %in% whitelist
avs[, blacklist] <- NULL

# Better variable name format
names(avs)[names(avs) == '_missense_all'] <- 'missense_all'

avs <- (
  avs
  
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
  
  # Other
  |> mutate(PDB_jury_column_rsa_unb = as.factor(PDB_jury_column_rsa_unb))  #include all RSA/MES intersections
  |> replace_na(list(PDB_protein_protein_interactions = 0,
                     PDB_surface_ppis = 0,
                     PDB_sequences_with_contacts = 0))
)


# Clean-up
rm(pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv)


COLOUR_SCHEME <- c(
  'cmd' = '#CC79A7',
  'ume' = '#56B4E9',
  'umn' = '#999999',
  'cmn' = '#999999',
  'umd' = '#E69F00',
  'cme' = '#009E73'
)
```

# Calculate and plot odds ratios

1.  Filter Pfams without any protein-protein interactions
2.  Filter low-powered Pfams
3.  Calculate summary statistics for each X-Class
4.  Calculate odds ratios for each X-class
5.  Repeat steps 3-4 for individual RSA classes
6.  Compose plot

First we apply the filters:

``` r
# Summarise counts
avs_test_set <- (
  avs
  
  # Column filters
  |> filter(occupancy > 0)  # must have human sequences
  |> filter(PDB_sequences_with_contacts > 0)  # drop columns with no PDB data
  
  # Pfam filters
  |> group_by(family) |> filter(any(PDB_protein_protein_interactions > 0)) |> ungroup()
  
  # Heuristic to select well-powered Pfams
  |> group_by(family) |> filter(any(mes_or < 1 & mes_p < 0.1)) |> ungroup()  # at least one depleted site
  |> group_by(family) |> filter(any(mes_or > 1 & mes_p < 0.1)) |> ungroup()  # at least one enriched site
)
```

Then summarise the data over the X-Classes:

``` r
summarize_data <- function(data, group_vars) {
  data |>
    group_by(across(all_of(group_vars))) |>
    summarise(
      n_pfams = n_distinct(family),
      n_columns = n(),
      n_residues = sum(occupancy),
      n_missense = sum(missense_all),
      across(c(PDB_protein_protein_interactions, PDB_sequences_with_contacts), sum),
      cons_class = first(cons_class),
      mes_class = first(mes_class)
    )
}

add_ratios <-
  function(data) {
    data %>%
      mutate(
        missense_per_res = n_missense / n_residues,
        ppi_per_pdb = PDB_protein_protein_interactions / PDB_sequences_with_contacts,
        PDB_per_column = PDB_sequences_with_contacts / n_columns
      )
  }

add_test_backgrounds <- function(data, group_vars) {
  data |>
    ungroup() |>
    group_by(across(all_of(group_vars))) |>
    mutate(across(
      c(n_missense, n_residues, PDB_protein_protein_interactions, PDB_sequences_with_contacts),
      sum,
      .names = "test_background_{.col}"
    ))
}


# Calculate group statistics for different X-Classes
avs_xclass_summary <- avs_test_set |>
  summarize_data(group_vars = "xclass") |>
  add_ratios() |>
  add_test_backgrounds(group_vars = NULL)
avs_xclass_summary
```

    ## # A tibble: 6 × 16
    ##   xclass n_pfams n_col…¹ n_res…² n_mis…³ PDB_p…⁴ PDB_s…⁵ cons_…⁶ mes_c…⁷ misse…⁸
    ##   <fct>    <int>   <int>   <dbl>   <dbl>   <dbl>   <dbl> <fct>   <fct>     <dbl>
    ## 1 cmd        535    3082  266547   72783   13238  164750 c       md        0.273
    ## 2 cmn        581   51712  954148  417352   93094  963031 c       mn        0.437
    ## 3 cme        464    1660   68164   50900    4872   47816 c       me        0.747
    ## 4 umd        401     991   48059   13113    4766   37061 u       md        0.273
    ## 5 umn        581   47032 1190891  563948  148812 1206996 u       mn        0.474
    ## 6 ume        524    2287  211780  129342   13690  116136 u       me        0.611
    ## # … with 6 more variables: ppi_per_pdb <dbl>, PDB_per_column <dbl>,
    ## #   test_background_n_missense <dbl>, test_background_n_residues <dbl>,
    ## #   test_background_PDB_protein_protein_interactions <dbl>,
    ## #   test_background_PDB_sequences_with_contacts <dbl>, and abbreviated variable
    ## #   names ¹​n_columns, ²​n_residues, ³​n_missense,
    ## #   ⁴​PDB_protein_protein_interactions, ⁵​PDB_sequences_with_contacts,
    ## #   ⁶​cons_class, ⁷​mes_class, ⁸​missense_per_res

``` r
# Totals
avs_xclass_summary |> ungroup() |> select(where(is.numeric)) |> summarise_all(sum)
```

    ## # A tibble: 1 × 13
    ##   n_pfams n_columns n_residues n_misse…¹ PDB_p…² PDB_s…³ misse…⁴ ppi_p…⁵ PDB_p…⁶
    ##     <int>     <int>      <dbl>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1    3086    106764    2739589   1247438  278472 2535790    2.81   0.649    215.
    ## # … with 4 more variables: test_background_n_missense <dbl>,
    ## #   test_background_n_residues <dbl>,
    ## #   test_background_PDB_protein_protein_interactions <dbl>,
    ## #   test_background_PDB_sequences_with_contacts <dbl>, and abbreviated variable
    ## #   names ¹​n_missense, ²​PDB_protein_protein_interactions,
    ## #   ³​PDB_sequences_with_contacts, ⁴​missense_per_res, ⁵​ppi_per_pdb,
    ## #   ⁶​PDB_per_column

``` r
# Chi-square test
chi <- chisq.test(avs_xclass_summary[, c("PDB_protein_protein_interactions", "PDB_sequences_with_contacts")])
chi
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  avs_xclass_summary[, c("PDB_protein_protein_interactions", "PDB_sequences_with_contacts")]
    ## X-squared = 4556.5, df = 5, p-value < 2.2e-16

## Odds ratios for protein-protein interaction enrichment in X-Classes

Now apply Fishers exact test for each X-Class (CMD, UMD, etc…).

``` r
perform_fisher_test <- function(df, group_vars) {
  df %>%
    group_by(!!!syms(group_vars)) %>%
    summarise(
      test = broom::tidy(
        fisher.test(matrix(
          c(
            PDB_protein_protein_interactions,
            PDB_sequences_with_contacts,
            test_background_PDB_protein_protein_interactions - PDB_protein_protein_interactions,
            test_background_PDB_sequences_with_contacts - PDB_sequences_with_contacts
          ),
          nrow = 2
        ), alternative = 'two.sided')
      ),
      cons_class = first(cons_class),
      mes_class = first(mes_class)
    ) %>%
    unnest(test)
}

avs_xclass_fisher <- perform_fisher_test(avs_xclass_summary, "xclass")
avs_xclass_fisher
```

    ## # A tibble: 6 × 9
    ##   xclass estimate   p.value conf.low conf.high method    alter…¹ cons_…² mes_c…³
    ##   <fct>     <dbl>     <dbl>    <dbl>     <dbl> <chr>     <chr>   <fct>   <fct>  
    ## 1 cmd       0.718 7.42e-304    0.705     0.731 Fisher's… two.si… c       md     
    ## 2 cmn       0.820 0            0.813     0.827 Fisher's… two.si… c       mn     
    ## 3 cme       0.927 4.04e-  7    0.899     0.955 Fisher's… two.si… c       me     
    ## 4 umd       1.17  3.88e- 24    1.14      1.21  Fisher's… two.si… u       md     
    ## 5 umn       1.26  0            1.25      1.27  Fisher's… two.si… u       mn     
    ## 6 ume       1.08  1.77e- 15    1.06      1.10  Fisher's… two.si… u       me     
    ## # … with abbreviated variable names ¹​alternative, ²​cons_class, ³​mes_class

Now do the same while stratifying for RSA class as well.

``` r
avs_xclass_summary_byRSA <- avs_test_set |>
  filter(PDB_jury_column_rsa_unb %in% c('Core', 'Part. Exposed', 'Surface')) |>
  summarize_data(group_vars = c("PDB_jury_column_rsa_unb", "xclass")) |>
  add_ratios() |>
  add_test_backgrounds(group_vars = "PDB_jury_column_rsa_unb")
```

    ## `summarise()` has grouped output by 'PDB_jury_column_rsa_unb'. You can override
    ## using the `.groups` argument.

``` r
avs_xclass_summary_byRSA
```

    ## # A tibble: 18 × 17
    ## # Groups:   PDB_jury_column_rsa_unb [3]
    ##    PDB_jury_col…¹ xclass n_pfams n_col…² n_res…³ n_mis…⁴ PDB_p…⁵ PDB_s…⁶ cons_…⁷
    ##    <fct>          <fct>    <int>   <int>   <dbl>   <dbl>   <dbl>   <dbl> <fct>  
    ##  1 Core           cmd        392    1358   98170   26007    2216   70217 c      
    ##  2 Core           cmn        542   15626  311009  129243   11857  358914 c      
    ##  3 Core           cme        181     288    8005    6714     338    8127 c      
    ##  4 Core           umd        117     184    5536    1319     282    6243 u      
    ##  5 Core           umn        464    6555  114201   51071    6336  153554 u      
    ##  6 Core           ume         92     127    3877    3025     178    4439 u      
    ##  7 Part. Exposed  cmd        340     850   70785   19187    3946   38274 c      
    ##  8 Part. Exposed  cmn        555   11678  215027   91658   21641  221443 c      
    ##  9 Part. Exposed  cme        218     394   14615   11793     972   11374 c      
    ## 10 Part. Exposed  umd        173     257   14773    4232    1198   11221 u      
    ## 11 Part. Exposed  umn        533    9255  207489   95553   21435  228315 u      
    ## 12 Part. Exposed  ume        208     312   16456   11091    1022   11930 u      
    ## 13 Surface        cmd        294     727   68726   20495    5323   23245 c      
    ## 14 Surface        cmn        569   21854  376037  175013   55695  330258 c      
    ## 15 Surface        cme        355     875   38041   27829    2843   21021 c      
    ## 16 Surface        umd        292     521   24384    6665    2819   16246 u      
    ## 17 Surface        umn        569   29776  801297  387794  112746  761380 u      
    ## 18 Surface        ume        484    1747  172668  105626    9686   79391 u      
    ## # … with 8 more variables: mes_class <fct>, missense_per_res <dbl>,
    ## #   ppi_per_pdb <dbl>, PDB_per_column <dbl>, test_background_n_missense <dbl>,
    ## #   test_background_n_residues <dbl>,
    ## #   test_background_PDB_protein_protein_interactions <dbl>,
    ## #   test_background_PDB_sequences_with_contacts <dbl>, and abbreviated variable
    ## #   names ¹​PDB_jury_column_rsa_unb, ²​n_columns, ³​n_residues, ⁴​n_missense,
    ## #   ⁵​PDB_protein_protein_interactions, ⁶​PDB_sequences_with_contacts, …

``` r
# Totals
avs_xclass_summary_byRSA |> ungroup() |> select(where(is.numeric)) |> summarise_all(sum)
```

    ## # A tibble: 1 × 13
    ##   n_pfams n_columns n_residues n_misse…¹ PDB_p…² PDB_s…³ misse…⁴ ppi_p…⁵ PDB_p…⁶
    ##     <int>     <int>      <dbl>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1    6378    102384    2561096   1174315  260533 2355592    8.77    1.78    568.
    ## # … with 4 more variables: test_background_n_missense <dbl>,
    ## #   test_background_n_residues <dbl>,
    ## #   test_background_PDB_protein_protein_interactions <dbl>,
    ## #   test_background_PDB_sequences_with_contacts <dbl>, and abbreviated variable
    ## #   names ¹​n_missense, ²​PDB_protein_protein_interactions,
    ## #   ³​PDB_sequences_with_contacts, ⁴​missense_per_res, ⁵​ppi_per_pdb,
    ## #   ⁶​PDB_per_column

## Odds ratios for protein-protein interaction enrichment in X-Class \* RSA

``` r
avs_xclass_fisher_byRSA <- perform_fisher_test(avs_xclass_summary_byRSA, c("PDB_jury_column_rsa_unb", "xclass"))
```

    ## `summarise()` has grouped output by 'PDB_jury_column_rsa_unb'. You can override
    ## using the `.groups` argument.

``` r
avs_xclass_fisher_byRSA
```

    ## # A tibble: 18 × 10
    ## # Groups:   PDB_jury_column_rsa_unb [3]
    ##    PDB_jury_co…¹ xclass estim…²   p.value conf.…³ conf.…⁴ method alter…⁵ cons_…⁶
    ##    <fct>         <fct>    <dbl>     <dbl>   <dbl>   <dbl> <chr>  <chr>   <fct>  
    ##  1 Core          cmd      0.883 3.09e-  8   0.844   0.923 Fishe… two.si… c      
    ##  2 Core          cmn      0.857 9.64e- 28   0.834   0.881 Fishe… two.si… c      
    ##  3 Core          cme      1.18  3.41e-  3   1.06    1.32  Fishe… two.si… c      
    ##  4 Core          umd      1.28  7.79e-  5   1.14    1.45  Fishe… two.si… u      
    ##  5 Core          umn      1.24  1.13e- 44   1.21    1.28  Fishe… two.si… u      
    ##  6 Core          ume      1.14  9.48e-  2   0.974   1.32  Fishe… two.si… u      
    ##  7 Part. Exposed cmd      1.08  1.45e-  5   1.04    1.12  Fishe… two.si… c      
    ##  8 Part. Exposed cmn      1.03  1.81e-  3   1.01    1.05  Fishe… two.si… c      
    ##  9 Part. Exposed cme      0.887 3.33e-  4   0.830   0.948 Fishe… two.si… c      
    ## 10 Part. Exposed umd      1.11  5.30e-  4   1.05    1.18  Fishe… two.si… u      
    ## 11 Part. Exposed umn      0.960 1.46e-  5   0.942   0.978 Fishe… two.si… u      
    ## 12 Part. Exposed ume      0.889 3.19e-  4   0.833   0.949 Fishe… two.si… u      
    ## 13 Surface       cmd      1.51  2.90e-143   1.46    1.55  Fishe… two.si… c      
    ## 14 Surface       cmn      1.14  2.90e-125   1.13    1.15  Fishe… two.si… c      
    ## 15 Surface       cme      0.879 8.17e- 11   0.845   0.914 Fishe… two.si… c      
    ## 16 Surface       umd      1.13  2.73e-  9   1.09    1.18  Fishe… two.si… u      
    ## 17 Surface       umn      0.912 7.90e- 75   0.903   0.921 Fishe… two.si… u      
    ## 18 Surface       ume      0.783 4.69e-114   0.767   0.801 Fishe… two.si… u      
    ## # … with 1 more variable: mes_class <fct>, and abbreviated variable names
    ## #   ¹​PDB_jury_column_rsa_unb, ²​estimate, ³​conf.low, ⁴​conf.high, ⁵​alternative,
    ## #   ⁶​cons_class

# Figure 4B: Odds ratio plot

``` r
plot_odds_ratio <- function(df, facet_var = NULL) {
  dodge_width <- 0.6
  
  plot <- df |>
    ggplot() +
    aes(x = xclass, y = estimate, colour = xclass) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high),
      width = 0.5,
      position = position_dodge(dodge_width)
    ) +
    geom_point(position = position_dodge(dodge_width)) +
    scale_y_log10() +
    ylab('Odds Ratio') +
    xlab('') +
    theme_minimal() +
    theme(legend.position = "none", aspect.ratio = 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
    scale_color_manual(values = COLOUR_SCHEME)
  
  if (!is.null(facet_var)) {
    plot <- plot + facet_grid(paste0(" ~ ", facet_var))
  }
  
  return(plot)
}


feature_oddsratio_plot <- plot_odds_ratio(avs_xclass_fisher)
feature_oddsratio_byRSA_plot <- plot_odds_ratio(avs_xclass_fisher_byRSA, "PDB_jury_column_rsa_unb")

patchworked_plot <- (feature_oddsratio_plot + feature_oddsratio_byRSA_plot
  + plot_layout(guides = 'collect', widths = c(1, 3))
  )


p_ranges_y <- 10^c(ggplot_build(patchworked_plot[[1]])$layout$panel_scales_y[[1]]$range$range,
                   ggplot_build(patchworked_plot[[2]])$layout$panel_scales_y[[1]]$range$range)

combined_plot <- (
    patchworked_plot
  & scale_y_log10(limits = c(1/1.6, 1.6),
                  breaks = c(0.7, 1, 1.5),
                  minor_breaks = c(0.6, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.6)
  )
  & scale_x_discrete(limits = c('', 'cmd', 'cmn', 'cme', '', 'umd', 'umn', 'ume', ''))
  & theme(
    axis.line = element_line(color='black'),
    panel.grid.minor = element_line(size = 0.1),
    panel.grid.major = element_line(size = 0.1),
    panel.grid.major.x = element_blank()
    )
)
```

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.
    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

``` r
combined_plot
```

![](Fig4B_PPI-XClass-Enrichments_files/figure-gfm/Combined%20plot-1.png)<!-- -->

``` r
ggsave('plot4b.svg', combined_plot)
```

    ## Saving 7 x 5 in image
