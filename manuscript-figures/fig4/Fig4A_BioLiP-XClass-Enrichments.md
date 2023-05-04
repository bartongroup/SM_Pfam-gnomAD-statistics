Figure 4A: Enrichment of protein-ligand interactions
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
- <a href="#figure-4a-odds-ratio-plot"
  id="toc-figure-4a-odds-ratio-plot">Figure 4A: Odds ratio plot</a>

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
  "PDB_biolip",
  "PDB_protein_ligand_interactions",
  "PDB_sequences_with_contacts",
  "PDB_jury_column_rsa_unb"
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
  |> replace_na(list(PDB_protein_ligand_interactions = 0,
                     PDB_biolip = 0,
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

1.  Filter Pfams without any BioLiP anotations
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
  |> group_by(family) |> filter(any(PDB_biolip > 0)) |> ungroup()
  
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
      across(c(PDB_protein_ligand_interactions, PDB_biolip, PDB_sequences_with_contacts), sum),
      cons_class = first(cons_class),
      mes_class = first(mes_class)
    )
}

add_ratios <-
  function(data) {
    data %>%
      mutate(
        missense_per_res = n_missense / n_residues,
        biolip_per_pdb = PDB_biolip / PDB_sequences_with_contacts,
        PDB_per_column = PDB_sequences_with_contacts / n_columns
      )
  }

add_test_backgrounds <- function(data, group_vars) {
  data |>
    ungroup() |>
    group_by(across(all_of(group_vars))) |>
    mutate(across(
      c(n_missense, n_residues, PDB_biolip, PDB_sequences_with_contacts),
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

    ## # A tibble: 6 × 17
    ##   xclass n_pfams n_col…¹ n_res…² n_mis…³ PDB_p…⁴ PDB_b…⁵ PDB_s…⁶ cons_…⁷ mes_c…⁸
    ##   <fct>    <int>   <int>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl> <fct>   <fct>  
    ## 1 cmd        307    1910  206997   58711   15839    6248  139002 c       md     
    ## 2 cmn        330   32821  653607  289913   68574   15758  760356 c       mn     
    ## 3 cme        268    1084   51445   36862    3677     641   40335 c       me     
    ## 4 umd        236     605   34637   10138    2779     655   30559 u       md     
    ## 5 umn        330   28029  808232  383723   75079   12621  958228 u       mn     
    ## 6 ume        307    1419  170283  100663    7172    1203   98450 u       me     
    ## # … with 7 more variables: missense_per_res <dbl>, biolip_per_pdb <dbl>,
    ## #   PDB_per_column <dbl>, test_background_n_missense <dbl>,
    ## #   test_background_n_residues <dbl>, test_background_PDB_biolip <dbl>,
    ## #   test_background_PDB_sequences_with_contacts <dbl>, and abbreviated variable
    ## #   names ¹​n_columns, ²​n_residues, ³​n_missense,
    ## #   ⁴​PDB_protein_ligand_interactions, ⁵​PDB_biolip,
    ## #   ⁶​PDB_sequences_with_contacts, ⁷​cons_class, ⁸​mes_class

``` r
# Totals
avs_xclass_summary |> ungroup() |> select(where(is.numeric)) |> summarise_all(sum)
```

    ## # A tibble: 1 × 14
    ##   n_pfams n_columns n_residues n_misse…¹ PDB_p…² PDB_b…³ PDB_s…⁴ misse…⁵ bioli…⁶
    ##     <int>     <int>      <dbl>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1    1778     65868    1925201    880010  173120   37126 2026930    2.80   0.128
    ## # … with 5 more variables: PDB_per_column <dbl>,
    ## #   test_background_n_missense <dbl>, test_background_n_residues <dbl>,
    ## #   test_background_PDB_biolip <dbl>,
    ## #   test_background_PDB_sequences_with_contacts <dbl>, and abbreviated variable
    ## #   names ¹​n_missense, ²​PDB_protein_ligand_interactions, ³​PDB_biolip,
    ## #   ⁴​PDB_sequences_with_contacts, ⁵​missense_per_res, ⁶​biolip_per_pdb

``` r
# Chi-square test
chi <- chisq.test(avs_xclass_summary[, c("PDB_biolip", "PDB_sequences_with_contacts")])
chi
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  avs_xclass_summary[, c("PDB_biolip", "PDB_sequences_with_contacts")]
    ## X-squared = 6980.1, df = 5, p-value < 2.2e-16

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
            PDB_biolip,
            PDB_sequences_with_contacts,
            test_background_PDB_biolip - PDB_biolip,
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
    ##   xclass estimate  p.value conf.low conf.high method     alter…¹ cons_…² mes_c…³
    ##   <fct>     <dbl>    <dbl>    <dbl>     <dbl> <chr>      <chr>   <fct>   <fct>  
    ## 1 cmd       2.75  0           2.67      2.83  Fisher's … two.si… c       md     
    ## 2 cmn       1.23  5.99e-83    1.20      1.25  Fisher's … two.si… c       mn     
    ## 3 cme       0.865 2.51e- 4    0.799     0.936 Fisher's … two.si… c       me     
    ## 4 umd       1.17  8.56e- 5    1.08      1.27  Fisher's … two.si… u       md     
    ## 5 umn       0.574 0           0.562     0.587 Fisher's … two.si… u       mn     
    ## 6 ume       0.656 1.85e-52    0.619     0.695 Fisher's … two.si… u       me     
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

    ## # A tibble: 18 × 18
    ## # Groups:   PDB_jury_column_rsa_unb [3]
    ##    PDB_jury_col…¹ xclass n_pfams n_col…² n_res…³ n_mis…⁴ PDB_p…⁵ PDB_b…⁶ PDB_s…⁷
    ##    <chr>          <fct>    <int>   <int>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 Core           cmd        229     889   76381   21285    3610    1170   57774
    ##  2 Core           cmn        305    9903  210625   89337   18341    4421  283316
    ##  3 Core           cme        114     183    5459    4394     476     106    6707
    ##  4 Core           umd         75     124    4230    1120     385      67    5441
    ##  5 Core           umn        259    4193   74701   34038    8270    1349  125127
    ##  6 Core           ume         66      90    3240    2418     263      51    3963
    ##  7 Part. Exposed  cmd        198     520   56699   16211    3784    1128   31065
    ##  8 Part. Exposed  cmn        314    7039  143401   62608   21610    5178  172574
    ##  9 Part. Exposed  cme        130     248   10490    8476    1005     199    9153
    ## 10 Part. Exposed  umd         94     149   10664    3254     929     232    9012
    ## 11 Part. Exposed  umn        297    5363  134848   62719   16040    2954  179336
    ## 12 Part. Exposed  ume        130     196   11087    7524     699     113    9314
    ## 13 Surface        cmd        164     362   45134   14127    2888    1145   17164
    ## 14 Surface        cmn        321   13679  249827  117625   23619    4624  253341
    ## 15 Surface        cme        206     563   28848   20085    1530     249   17244
    ## 16 Surface        umd        175     308   16431    4873     993     210   12765
    ## 17 Surface        umn        321   17331  534209  259243   44011    6595  590879
    ## 18 Surface        ume        278    1042  137279   81256    4243     749   64839
    ## # … with 9 more variables: cons_class <fct>, mes_class <fct>,
    ## #   missense_per_res <dbl>, biolip_per_pdb <dbl>, PDB_per_column <dbl>,
    ## #   test_background_n_missense <dbl>, test_background_n_residues <dbl>,
    ## #   test_background_PDB_biolip <dbl>,
    ## #   test_background_PDB_sequences_with_contacts <dbl>, and abbreviated variable
    ## #   names ¹​PDB_jury_column_rsa_unb, ²​n_columns, ³​n_residues, ⁴​n_missense,
    ## #   ⁵​PDB_protein_ligand_interactions, ⁶​PDB_biolip, …

``` r
# Totals
avs_xclass_summary_byRSA |> ungroup() |> select(where(is.numeric)) |> summarise_all(sum)
```

    ## # A tibble: 1 × 14
    ##   n_pfams n_columns n_residues n_misse…¹ PDB_p…² PDB_b…³ PDB_s…⁴ misse…⁵ bioli…⁶
    ##     <int>     <int>      <dbl>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1    3676     62182    1753553    810593  152696   30540 1849014    8.81   0.369
    ## # … with 5 more variables: PDB_per_column <dbl>,
    ## #   test_background_n_missense <dbl>, test_background_n_residues <dbl>,
    ## #   test_background_PDB_biolip <dbl>,
    ## #   test_background_PDB_sequences_with_contacts <dbl>, and abbreviated variable
    ## #   names ¹​n_missense, ²​PDB_protein_ligand_interactions, ³​PDB_biolip,
    ## #   ⁴​PDB_sequences_with_contacts, ⁵​missense_per_res, ⁶​biolip_per_pdb

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
    ##    PDB_jury_co…¹ xclass estim…²   p.value conf.…³ conf.…⁴ method alter…⁵ cons_…⁶
    ##    <chr>         <fct>    <dbl>     <dbl>   <dbl>   <dbl> <chr>  <chr>   <fct>  
    ##  1 Core          cmd      1.43  6.34e- 27   1.35    1.53  Fishe… two.si… c      
    ##  2 Core          cmn      1.13  3.79e-  7   1.08    1.19  Fishe… two.si… c      
    ##  3 Core          cme      1.06  5.09e-  1   0.869   1.29  Fishe… two.si… c      
    ##  4 Core          umd      0.827 1.28e-  1   0.640   1.05  Fishe… two.si… u      
    ##  5 Core          umn      0.662 4.03e- 45   0.623   0.703 Fishe… two.si… u      
    ##  6 Core          ume      0.865 3.55e-  1   0.643   1.14  Fishe… two.si… u      
    ##  7 Part. Exposed cmd      1.59  4.83e- 42   1.49    1.69  Fishe… two.si… c      
    ##  8 Part. Exposed cmn      1.54  1.52e- 99   1.48    1.61  Fishe… two.si… c      
    ##  9 Part. Exposed cme      0.908 2.00e-  1   0.784   1.05  Fishe… two.si… c      
    ## 10 Part. Exposed umd      1.08  2.50e-  1   0.942   1.23  Fishe… two.si… u      
    ## 11 Part. Exposed umn      0.556 1.17e-163   0.532   0.581 Fishe… two.si… u      
    ## 12 Part. Exposed ume      0.502 8.58e- 16   0.413   0.605 Fishe… two.si… u      
    ## 13 Surface       cmd      5.04  0           4.73    5.37  Fishe… two.si… c      
    ## 14 Surface       cmn      1.43  2.91e- 83   1.38    1.49  Fishe… two.si… c      
    ## 15 Surface       cme      1.02  7.70e-  1   0.893   1.15  Fishe… two.si… c      
    ## 16 Surface       umd      1.16  3.50e-  2   1.01    1.33  Fishe… two.si… u      
    ## 17 Surface       umn      0.584 4.03e-210   0.565   0.605 Fishe… two.si… u      
    ## 18 Surface       ume      0.803 2.52e-  9   0.745   0.865 Fishe… two.si… u      
    ## # … with 1 more variable: mes_class <fct>, and abbreviated variable names
    ## #   ¹​PDB_jury_column_rsa_unb, ²​estimate, ³​conf.low, ⁴​conf.high, ⁵​alternative,
    ## #   ⁶​cons_class

# Figure 4A: Odds ratio plot

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
  & scale_y_log10(limits = c(min(p_ranges_y), max(p_ranges_y)),
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

![](Fig4A_BioLiP-XClass-Enrichments_files/figure-gfm/Combined%20plot-1.png)<!-- -->

``` r
ggsave('plot4a.svg', combined_plot)
```

    ## Saving 7 x 5 in image
