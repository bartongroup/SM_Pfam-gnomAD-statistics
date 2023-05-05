Figure 3A: The Conservation Plane
================

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#notebook-resources" id="toc-notebook-resources">Notebook
  resources</a>
  - <a href="#libraries" id="toc-libraries">Libraries</a>
  - <a href="#dataset" id="toc-dataset">Dataset</a>
- <a href="#plot-sample" id="toc-plot-sample">Plot sample</a>
- <a href="#figure-3a-the-conservation-plane"
  id="toc-figure-3a-the-conservation-plane">Figure 3A: The Conservation
  Plane</a>

# Introduction

This R markdown document visualizes the relationship between
evolutionary conservation and population constraint in Pfam domains. The
dataset used for this analysis is obtained from our pre-calculated MSA
column statistics for Pfam, included in this repository. The resulting
plot corresponds to Figure 3A in the manuscript.

# Notebook resources

## Libraries

``` r
suppressPackageStartupMessages(library(tidyverse))


COLOUR_SCHEME <- c(
  'cmd' = '#CC79A7',
  'ume' = '#56B4E9',
  'umn' = '#999999',
  'cmn' = '#999999',
  'umd' = '#E69F00',
  'cme' = '#009E73'
)
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
  # Minimal filters for notebook
avs <- (
  avs
  # Column filters
  |> filter(occupancy > 0)  # must have human sequences
)
```

# Plot sample

``` r
# We will plot a sample of the data so setting a seed to ensure reproducibility
set.seed(777)
(mes_vs_cons_data <- 
  avs
  |> filter(occupancy_gteq_threshold)
  |> filter(occupancy > 20)
  |> sample_n(10000)  # good balance for processing, completeness and legibility
  |> mutate(
    mes_class = factor(
      ifelse(mes_p < 0.1, ifelse(mes_or < 1, 'md', 'me'), 'mn'),
      levels = c('md', 'mn', 'me')
    ),
    cons_class = factor(
      ifelse(shenkin_nrank > 0.5, 'u', 'c'),
      levels = c('u', 'c')
    ),
    xclass = factor(
      interaction(cons_class, mes_class, sep = ''),
      levels = c('cmd', 'cmn', 'cme', 'umd', 'umn', 'ume')
    )
  )
)
```

    ## # A tibble: 10,000 × 19
    ##    family  column splice…¹ synon…² occup…³ shenkin misse…⁴ mes_or  mes_p occup…⁵
    ##    <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>
    ##  1 PF00271    104        1       9      49    24.5      16  0.894 0.780       23
    ##  2 PF01839     24        0       8      35    55.4      15  0.783 0.459       27
    ##  3 PF00337    102        0       3      22    22.5      10  0.911 1           17
    ##  4 PF00003    400        0       5      22    30.2       4  0.369 0.0609      17
    ##  5 PF01839     16        0      12      35    37.0      12  0.624 0.168       27
    ##  6 PF00225    605        0       9      40    62.2      18  1.02  1           29
    ##  7 PF15779     90        0       0      23    10.6       6  1.36  0.452       18
    ##  8 PF00615    161        0       7      38    73.2      10  0.583 0.158       29
    ##  9 PF00531     29        1       9      31    82.7      17  1.13  0.757       24
    ## 10 PF00201    385        0       4      21    22.2       9  0.649 0.351       16
    ## # … with 9,990 more rows, 9 more variables: occupancy_gteq_threshold <lgl>,
    ## #   shenkin_nrank <dbl>, mes_or_nrank <dbl>, mes_or_log <dbl>,
    ## #   PDB_jury_column_rsa_unb <chr>, synonymous_all <dbl>, mes_class <fct>,
    ## #   cons_class <fct>, xclass <fct>, and abbreviated variable names
    ## #   ¹​`splice_region_variant&synonymous_variant`, ²​synonymous_variant,
    ## #   ³​occupancy, ⁴​missense_all, ⁵​occupancy_minus_threshold

``` r
mes_vs_cons_data |>
  summarise(n(), n_distinct(family))
```

    ## # A tibble: 1 × 2
    ##   `n()` `n_distinct(family)`
    ##   <int>                <int>
    ## 1 10000                  283

``` r
mes_vs_cons_data |>
  group_by(PDB_jury_column_rsa_unb) |>
  summarise(n(), n_distinct(family))
```

    ## # A tibble: 4 × 3
    ##   PDB_jury_column_rsa_unb `n()` `n_distinct(family)`
    ##   <chr>                   <int>                <int>
    ## 1 Core                     1996                  225
    ## 2 Part. Exposed            1900                  232
    ## 3 Surface                  4771                  251
    ## 4 <NA>                     1333                   46

# Figure 3A: The Conservation Plane

``` r
mes_vs_cons_data |>
  arrange(mes_p) |>  # points with lower p will be plotted on top
  ggplot() +
  aes(x = shenkin_nrank, y = mes_or) +
  geom_point(aes(colour = xclass, alpha = ifelse(mes_p < 0.1, 1-mes_p, 0.01))) +
  scale_alpha_continuous(range = c(0.1, 0.7)) +
  scale_color_manual(values = COLOUR_SCHEME) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_vline(xintercept = 0.5, lty = 2) +
  
  # Annotation
  # https://stackoverflow.com/questions/28709331/logarithmic-grid-for-plot-with-ggplot2
  scale_y_log10(breaks = 10^(-10:10),
                minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9)),
                labels = 10^(-10:10)) +
  annotation_logticks(sides = 'l') +
  xlab('Normalised Shenkin') +
  ylab('MES') +
  annotate("text", x = 0.4, y = 1.20, label = "CME") +
  annotate("text", x = 0.4, y = 0.83, label = "CMD") +
  annotate("text", x = 0.6, y = 1.20, label = "UME") +
  annotate("text", x = 0.6, y = 0.83, label = "UMD") +

  # Plot theme
  theme_minimal() +
  theme(legend.position = "none", aspect.ratio = 1) +
  theme(
    axis.line = element_line(color='black'),  # heavy line for axes
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    # Hide facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    ) ->
  mes_vs_cons_plot

mes_vs_cons_plot
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](Fig3A_Conservation-Plane_files/figure-gfm/Plot%20of%20the%20conservation%20plane-1.png)<!-- -->
