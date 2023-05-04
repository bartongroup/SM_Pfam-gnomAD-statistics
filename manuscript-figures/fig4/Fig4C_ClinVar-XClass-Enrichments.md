Figure 4C: Enrichment of Pathogenic variants in Conservation Plane
================

- <a href="#notebook-resources" id="toc-notebook-resources">Notebook
  resources</a>
  - <a href="#libraries" id="toc-libraries">Libraries</a>
  - <a href="#dataset" id="toc-dataset">Dataset</a>
- <a href="#calculate-and-plot-odds-ratios"
  id="toc-calculate-and-plot-odds-ratios">Calculate and plot odds
  ratios</a>
  - <a href="#odds-ratios-for-pathogenic-enrichment-in-x-classes"
    id="toc-odds-ratios-for-pathogenic-enrichment-in-x-classes">Odds ratios
    for Pathogenic enrichment in X-Classes</a>
  - <a href="#odds-ratios-for-pathogenic-enrichment-in-x-class--rsa"
    id="toc-odds-ratios-for-pathogenic-enrichment-in-x-class--rsa">Odds
    ratios for Pathogenic enrichment in X-Class * RSA</a>
- <a href="#figure-4c-odds-ratio-plot"
  id="toc-figure-4c-odds-ratio-plot">Figure 4C: Odds ratio plot</a>
- <a href="#supporting-statistics"
  id="toc-supporting-statistics">Supporting statistics</a>
- <a href="#future-work" id="toc-future-work">Future work</a>

This notebook presents code that generates Figure 4C and supporting
statistics, which describe the enrichment of ClinVar Pathogenic variants
within different regions of the conservation plane.

The notebook works with our pre-calculated MSA column statistics for
population and clinical variants, and structural features, for relevant
Pfams. If you want to go straight to the plot code then follow the link
in the Table of Contents above.

# Notebook resources

## Libraries

``` r
suppressPackageStartupMessages(library(tidyverse))
library(patchwork)
```

## Dataset

``` r
# Read and process dataset
pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv <-
  read_csv("./data/pfam-gnomAD-clinvar-pdb-colstats_c7c3e19.csv.gz",
           comment = "#", show_col_types = FALSE)

avs <-
  pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv  # avs: alignments variants structure

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
  "CV_Benign",
  "CV_Benign/Likely_benign",
  "CV_Conflicting_interpretations_of_pathogenicity",
  "CV_Likely_benign",
  "CV_Likely_pathogenic",
  "CV_Pathogenic",
  "CV_Uncertain_significance",
  "CV_Pathogenic/Likely_pathogenic",
  "PDB_jury_column_rsa_unb"
)
blacklist <- !names(avs) %in% whitelist
avs[, blacklist] <- NULL

# FYI: There are several other ClinVar entries I could look at in future
# "CV_not_provided",
# "CV_drug_response",
# "CV_association",
# "CV_risk_factor",
# "CV_other",
# "CV_Affects",
# "CV_protective",

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
  |> replace_na(list(CV_Pathogenic = 0))
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

1.  Filter Pfams without any Pathogenic variants
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
  |> filter(occupancy > 0)
  
  # Pfam filters
  |> group_by(family) |> filter(any(CV_Pathogenic > 0)) |> ungroup()
  
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
      across(c(CV_Pathogenic, CV_Benign), sum),
      cons_class = first(cons_class),
      mes_class = first(mes_class)
    )
}

add_ratios <-
  function(data, pathogenic_col, missense_col, residues_col) {
    data %>%
      mutate(
        missense_per_res = !!sym(missense_col) / !!sym(residues_col),
        patho_per_missense = !!sym(pathogenic_col) / !!sym(missense_col),
        patho_per_column = !!sym(pathogenic_col) / !!sym(residues_col)
      )
  }

add_test_backgrounds <- function(data, group_vars) {
  data |>
    ungroup() |>
    group_by(across(all_of(group_vars))) |>
    mutate(across(
      c(n_missense, n_residues, CV_Pathogenic, CV_Benign),
      sum,
      .names = "test_background_{.col}"
    ))
}


# Calculate group statistics for different X-Classes
avs_xclass_summary <- avs_test_set |>
  summarize_data(group_vars = "xclass") |>
  add_ratios(
    pathogenic_col = "CV_Pathogenic",
    missense_col = "n_missense",
    residues_col = "n_residues"
  ) |>
  add_test_backgrounds(group_vars = NULL)
avs_xclass_summary
```

    ## # A tibble: 6 × 16
    ##   xclass n_pfams n_col…¹ n_res…² n_mis…³ CV_Pa…⁴ CV_Be…⁵ cons_…⁶ mes_c…⁷ misse…⁸
    ##   <fct>    <int>   <int>   <dbl>   <dbl>   <dbl>   <dbl> <fct>   <fct>     <dbl>
    ## 1 cmd        382    2482  198790   52109     819      24 c       md        0.262
    ## 2 cmn        406   67269  776297  341192    2618     264 c       mn        0.440
    ## 3 cme        354    1718   51479   40057     437      22 c       me        0.778
    ## 4 umd        301     848   40514   11069     129      11 u       md        0.273
    ## 5 umn        406   38376  889607  421071    2027     379 u       mn        0.473
    ## 6 ume        372    1798  165464  100376     232      72 u       me        0.607
    ## # … with 6 more variables: patho_per_missense <dbl>, patho_per_column <dbl>,
    ## #   test_background_n_missense <dbl>, test_background_n_residues <dbl>,
    ## #   test_background_CV_Pathogenic <dbl>, test_background_CV_Benign <dbl>, and
    ## #   abbreviated variable names ¹​n_columns, ²​n_residues, ³​n_missense,
    ## #   ⁴​CV_Pathogenic, ⁵​CV_Benign, ⁶​cons_class, ⁷​mes_class, ⁸​missense_per_res

``` r
# Totals
avs_xclass_summary |> ungroup() |> select(where(is.numeric)) |> summarise_all(sum)
```

    ## # A tibble: 1 × 13
    ##   n_pfams n_columns n_residues n_misse…¹ CV_Pa…² CV_Be…³ misse…⁴ patho…⁵ patho…⁶
    ##     <int>     <int>      <dbl>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1    2221    112491    2122151    965874    6262     772    2.83  0.0531  0.0228
    ## # … with 4 more variables: test_background_n_missense <dbl>,
    ## #   test_background_n_residues <dbl>, test_background_CV_Pathogenic <dbl>,
    ## #   test_background_CV_Benign <dbl>, and abbreviated variable names
    ## #   ¹​n_missense, ²​CV_Pathogenic, ³​CV_Benign, ⁴​missense_per_res,
    ## #   ⁵​patho_per_missense, ⁶​patho_per_column

## Odds ratios for Pathogenic enrichment in X-Classes

Now apply Fishers exact test for each X-Class (CMD, UMD, etc…).

``` r
perform_fisher_test <- function(df, group_vars) {
  df %>%
    group_by(!!!syms(group_vars)) %>%
    summarise(
      test = broom::tidy(
        fisher.test(matrix(
          c(
            CV_Pathogenic,
            n_missense,
            test_background_CV_Pathogenic - CV_Pathogenic,
            test_background_n_missense - n_missense
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
    ## 1 cmd       2.64  5.16e-116    2.45      2.84  Fisher's… two.si… c       md     
    ## 2 cmn       1.32  4.12e- 26    1.25      1.38  Fisher's… two.si… c       mn     
    ## 3 cme       1.73  1.17e- 24    1.57      1.91  Fisher's… two.si… c       me     
    ## 4 umd       1.81  9.80e- 10    1.51      2.16  Fisher's… two.si… u       md     
    ## 5 umn       0.619 2.82e- 73    0.587     0.653 Fisher's… two.si… u       mn     
    ## 6 ume       0.332 6.54e- 86    0.290     0.378 Fisher's… two.si… u       me     
    ## # … with abbreviated variable names ¹​alternative, ²​cons_class, ³​mes_class

Now do the same while stratifying for RSA class as well.

``` r
avs_xclass_summary_byRSA <- avs_test_set |>
  filter(PDB_jury_column_rsa_unb %in% c('Core', 'Part. Exposed', 'Surface')) |>
  summarize_data(group_vars = c("PDB_jury_column_rsa_unb", "xclass")) |>
  add_ratios(
    pathogenic_col = "CV_Pathogenic",
    missense_col = "n_missense",
    residues_col = "n_columns"
  ) |>
  add_test_backgrounds(group_vars = "PDB_jury_column_rsa_unb")
```

    ## `summarise()` has grouped output by 'PDB_jury_column_rsa_unb'. You can override
    ## using the `.groups` argument.

``` r
avs_xclass_summary_byRSA
```

    ## # A tibble: 18 × 17
    ## # Groups:   PDB_jury_column_rsa_unb [3]
    ##    PDB_jury_col…¹ xclass n_pfams n_col…² n_res…³ n_mis…⁴ CV_Pa…⁵ CV_Be…⁶ cons_…⁷
    ##    <chr>          <fct>    <int>   <int>   <dbl>   <dbl>   <dbl>   <dbl> <fct>  
    ##  1 Core           cmd        285    1020   68827   17299     138       8 c      
    ##  2 Core           cmn        364   11823  225006   92775     888      55 c      
    ##  3 Core           cme        139     233    6240    5197     112       3 c      
    ##  4 Core           umd         97     153    4851    1184      17       2 u      
    ##  5 Core           umn        315    5268   91587   40879     320      33 u      
    ##  6 Core           ume         78     110    3342    2619       8       1 u      
    ##  7 Part. Exposed  cmd        245     633   48834   12508     219       4 c      
    ##  8 Part. Exposed  cmn        372    8905  157099   66683     617      52 c      
    ##  9 Part. Exposed  cme        159     296    8379    7141     115       3 c      
    ## 10 Part. Exposed  umd        122     195   11385    3194      42       3 u      
    ## 11 Part. Exposed  umn        362    6954  151762   69598     411      61 u      
    ## 12 Part. Exposed  ume        157     240   10912    7459      23       5 u      
    ## 13 Surface        cmd        201     516   49535   14716     322       7 c      
    ## 14 Surface        cmn        381   16535  267618  124569     849     101 c      
    ## 15 Surface        cme        245     663   27905   20543     172      15 c      
    ## 16 Surface        umd        211     402   20450    5811      66       5 u      
    ## 17 Surface        umn        380   21715  550035  264221    1109     243 u      
    ## 18 Surface        ume        324    1242  130771   78908     162      54 u      
    ## # … with 8 more variables: mes_class <fct>, missense_per_res <dbl>,
    ## #   patho_per_missense <dbl>, patho_per_column <dbl>,
    ## #   test_background_n_missense <dbl>, test_background_n_residues <dbl>,
    ## #   test_background_CV_Pathogenic <dbl>, test_background_CV_Benign <dbl>, and
    ## #   abbreviated variable names ¹​PDB_jury_column_rsa_unb, ²​n_columns,
    ## #   ³​n_residues, ⁴​n_missense, ⁵​CV_Pathogenic, ⁶​CV_Benign, ⁷​cons_class

``` r
# Totals
avs_xclass_summary_byRSA |> ungroup() |> select(where(is.numeric)) |> summarise_all(sum)
```

    ## # A tibble: 1 × 13
    ##   n_pfams n_columns n_residues n_misse…¹ CV_Pa…² CV_Be…³ misse…⁴ patho…⁵ patho…⁶
    ##     <int>     <int>      <dbl>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1    4437     76903    1834538    835304    5590     655    352.   0.184    3.39
    ## # … with 4 more variables: test_background_n_missense <dbl>,
    ## #   test_background_n_residues <dbl>, test_background_CV_Pathogenic <dbl>,
    ## #   test_background_CV_Benign <dbl>, and abbreviated variable names
    ## #   ¹​n_missense, ²​CV_Pathogenic, ³​CV_Benign, ⁴​missense_per_res,
    ## #   ⁵​patho_per_missense, ⁶​patho_per_column

## Odds ratios for Pathogenic enrichment in X-Class \* RSA

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
    ##    PDB_jury_col…¹ xclass estim…²  p.value conf.…³ conf.…⁴ method alter…⁵ cons_…⁶
    ##    <chr>          <fct>    <dbl>    <dbl>   <dbl>   <dbl> <chr>  <chr>   <fct>  
    ##  1 Core           cmd      0.846 6.43e- 2   0.704   1.01  Fishe… two.si… c      
    ##  2 Core           cmn      1.08  1.46e- 1   0.973   1.20  Fishe… two.si… c      
    ##  3 Core           cme      2.43  1.86e-15   1.98    2.96  Fishe… two.si… c      
    ##  4 Core           umd      1.55  9.13e- 2   0.900   2.51  Fishe… two.si… u      
    ##  5 Core           umn      0.801 4.11e- 4   0.706   0.908 Fishe… two.si… u      
    ##  6 Core           ume      0.326 2.64e- 4   0.140   0.645 Fishe… two.si… u      
    ##  7 Part. Exposed  cmd      2.23  3.46e-23   1.92    2.58  Fishe… two.si… c      
    ##  8 Part. Exposed  cmn      1.14  1.46e- 2   1.03    1.27  Fishe… two.si… c      
    ##  9 Part. Exposed  cme      1.96  3.26e-10   1.60    2.37  Fishe… two.si… c      
    ## 10 Part. Exposed  umd      1.55  8.58e- 3   1.11    2.11  Fishe… two.si… u      
    ## 11 Part. Exposed  umn      0.564 6.14e-24   0.501   0.633 Fishe… two.si… u      
    ## 12 Part. Exposed  ume      0.349 4.70e- 9   0.221   0.527 Fishe… two.si… u      
    ## 13 Surface        cmd      4.58  1.86e-99   4.06    5.16  Fishe… two.si… c      
    ## 14 Surface        cmn      1.43  4.66e-17   1.32    1.55  Fishe… two.si… c      
    ## 15 Surface        cme      1.63  7.52e- 9   1.39    1.90  Fishe… two.si… c      
    ## 16 Surface        umd      2.19  2.72e- 8   1.68    2.79  Fishe… two.si… u      
    ## 17 Surface        umn      0.653 1.02e-27   0.604   0.706 Fishe… two.si… u      
    ## 18 Surface        ume      0.350 2.51e-51   0.297   0.411 Fishe… two.si… u      
    ## # … with 1 more variable: mes_class <fct>, and abbreviated variable names
    ## #   ¹​PDB_jury_column_rsa_unb, ²​estimate, ³​conf.low, ⁴​conf.high, ⁵​alternative,
    ## #   ⁶​cons_class

# Figure 4C: Odds ratio plot

``` r
plot_odds_ratio <- function(df, facet_var = NULL) {
  dodge_width <- 0.6
  
  plot <- df |>
    ggplot() +
    aes(x = xclass, y = estimate, colour = xclass) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_errorbar(
      aes(ymin = sapply(conf.low, function(y) max(y, 0.2)),  # clip
          ymax = conf.high),
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
  & scale_y_log10(limits = c(0.2, 6),  # limit in conjunction with clip
                  breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5, 6),
                  labels = c('', '', 0.3, '', 0.5, 1, 2, 3, 4, 5, ''),
                  minor_breaks = c(seq(0.1, 1, 0.1), seq(1, 5, 1)),
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

![](Fig4C_ClinVar-XClass-Enrichments_files/figure-gfm/Combined%20plot-1.png)<!-- -->

``` r
ggsave('plot4c.svg', combined_plot)
```

    ## Saving 7 x 5 in image

# Supporting statistics

We quoted the maximum discrimination observed with respect to X-Class
and ClinVar Pathogenic variants. This is found by comparing Surface CMD
to UME positions:

``` r
(
  avs
  |> filter(PDB_jury_column_rsa_unb == 'Surface')
  |> filter(xclass %in% c('cmd', 'ume'))
  |> group_by(xclass)
  |> summarise(
    n_patho = sum(CV_Pathogenic),
    n_missense = sum(missense_all)
  )
  -> cmd_vs_ume
)
```

    ## # A tibble: 2 × 3
    ##   xclass n_patho n_missense
    ##   <fct>    <dbl>      <dbl>
    ## 1 cmd        322      20644
    ## 2 ume        182     115498

``` r
(
  cmd_vs_ume
  |> pivot_wider(names_from = xclass, values_from = starts_with("n_"))
  |> summarise(
      test =  broom::tidy(
        fisher.test(matrix(c(n_patho_cmd, n_missense_cmd, n_patho_ume, n_missense_ume),
                           nrow = 2), alternative = 'two.sided'))
    )
  |> unnest(test)
)
```

    ## # A tibble: 1 × 6
    ##   estimate   p.value conf.low conf.high method                           alter…¹
    ##      <dbl>     <dbl>    <dbl>     <dbl> <chr>                            <chr>  
    ## 1     9.90 5.52e-135     8.22      11.9 Fisher's Exact Test for Count D… two.si…
    ## # … with abbreviated variable name ¹​alternative

It’s sensible to check that we observe the widely reported enrichment of
Pathogenic variants in buried positions:

``` r
avs |>
  group_by(PDB_jury_column_rsa_unb) |>
  summarise(n_missense = sum(missense_all), CV_Pathogenic = sum(CV_Pathogenic)) |>
  mutate(across(c(n_missense, CV_Pathogenic), sum, .names = "test_background_{.col}")) |>
  group_by(PDB_jury_column_rsa_unb, .drop = TRUE) |>
  summarise(
    # Biolip enrichment
    test =  broom::tidy(
      fisher.test(matrix(c(CV_Pathogenic, n_missense,
                           test_background_CV_Pathogenic - CV_Pathogenic,
                           test_background_n_missense - n_missense),
                         nrow = 2), alternative = 'two.sided'))
  ) |>
  unnest(test)
```

    ## # A tibble: 5 × 7
    ##   PDB_jury_column_rsa_unb estimate   p.value conf.low conf.high method   alter…¹
    ##   <chr>                      <dbl>     <dbl>    <dbl>     <dbl> <chr>    <chr>  
    ## 1 -                          0     1   e+  0    0        10.8   Fisher'… two.si…
    ## 2 Core                       2.05  3.35e-183    1.95      2.14  Fisher'… two.si…
    ## 3 Part. Exposed              1.63  1.40e- 81    1.55      1.71  Fisher'… two.si…
    ## 4 Surface                    0.796 3.50e- 29    0.764     0.828 Fisher'… two.si…
    ## 5 <NA>                       0.530 2.04e-164    0.505     0.556 Fisher'… two.si…
    ## # … with abbreviated variable name ¹​alternative

# Future work

- ClinVar enrichments in different positions with various structural
  features
