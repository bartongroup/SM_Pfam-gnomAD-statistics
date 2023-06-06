Figure 2A: RSA class stratified MES and Shenkin residue ECDFs
================

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#notebook-resources" id="toc-notebook-resources">Notebook
  resources</a>
  - <a href="#libraries" id="toc-libraries">Libraries</a>
  - <a href="#dataset" id="toc-dataset">Dataset</a>
- <a
  href="#figure-2a-residue-shenkingmes-ecdfs-for-each-exposure-class-manuscript-figure"
  id="toc-figure-2a-residue-shenkingmes-ecdfs-for-each-exposure-class-manuscript-figure">Figure
  2A: Residue Shenking/MES ECDFs for each exposure class (manuscript
  figure)</a>
- <a href="#calculating-some-quoted-figures-in-the-manuscript"
  id="toc-calculating-some-quoted-figures-in-the-manuscript">Calculating
  some quoted figures in the manuscript</a>

# Introduction

This R markdown document visualizes relationships between solvent
accessibility with both evolutionary conservation and population
constraint in Pfam domains.

The dataset used for this analysis is obtained from our pre-calculated
MSA column statistics for Pfam, included in this repository. The
resulting plots corresponds to Figure 2A in the manuscript.

# Notebook resources

## Libraries

``` r
suppressPackageStartupMessages(library(tidyverse))
```

## Dataset

``` r
# Load pfam-gnomAD dataset}
pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv <-
  read_csv("data/pfam-gnomAD-clinvar-pdb-colstats_c7c3e19.csv.gz", comment = "#")
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
# Processing
avs <- pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv  # avs: alignments variants structure

# Filter columns
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
  "PDB_jury_column_rsa_unb",
  "PDB_jury_column_rsa_bound"
)

blacklist <- !names(avs) %in% whitelist
avs[, blacklist] <- NULL

# Simplify column names
names(avs)[names(avs) == '_missense_all'] <- 'missense_all'
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
(avs <- avs
  # Column filters
  |> filter(occupancy > 0)  # must have human sequences
  # |> filter(PDB_sequences_with_contacts > 0)  # drop Pfams with no PDB data?

  # Powered Pfams
  |> group_by(family) |> filter(any(mes_or < 1 & mes_p < 0.1)) |> ungroup()
  |> group_by(family) |> filter(any(mes_or > 1 & mes_p < 0.1)) |> ungroup()
  
  # Other
  |> mutate(PDB_jury_column_rsa_unb = as.factor(PDB_jury_column_rsa_unb))  # count all RSA*MES intersections
  |> replace_na(list(PDB_protein_ligand_interactions = 0, PDB_biolip = 0,
                     PDB_sequences_with_contacts = 0))
)
```

    ## # A tibble: 167,171 × 18
    ##    family  column occupancy shenkin misse…¹ mes_or mes_p occup…² occup…³ shenk…⁴
    ##    <chr>    <dbl>     <dbl>   <dbl>   <dbl>  <dbl> <dbl>   <dbl> <lgl>     <dbl>
    ##  1 PF00001     41         1    6.30       1  1.89  1         -70 FALSE         0
    ##  2 PF00001     42         1    6.33       1  1.89  1         -70 FALSE         0
    ##  3 PF00001     43         2    6.37       1  0.945 1         -69 FALSE         0
    ##  4 PF00001     44         3    6.46       0  0     0.556     -68 FALSE         0
    ##  5 PF00001     45         5    6.59       2  0.756 1         -66 FALSE         0
    ##  6 PF00001     46         5    6.82       4  1.51  0.507     -66 FALSE         0
    ##  7 PF00001     47         6    6.84       7  2.20  0.154     -65 FALSE         0
    ##  8 PF00001     48         7    6.96       6  1.62  0.393     -64 FALSE         0
    ##  9 PF00001     49         7    7.04       4  1.08  1         -64 FALSE         0
    ## 10 PF00001     50         7    7.09       5  1.35  0.763     -64 FALSE         0
    ## # … with 167,161 more rows, 8 more variables: mes_or_nrank <dbl>,
    ## #   mes_or_log <dbl>, PDB_sequences_with_contacts <dbl>,
    ## #   PDB_jury_column_rsa_unb <fct>, PDB_jury_column_rsa_bound <chr>,
    ## #   mes_class <fct>, cons_class <fct>, xclass <fct>, and abbreviated variable
    ## #   names ¹​missense_all, ²​occupancy_minus_threshold, ³​occupancy_gteq_threshold,
    ## #   ⁴​shenkin_nrank

``` r
rm(pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv)
```

``` r
 RSA_LEVELS <- c('Core', 'Part. Exposed', 'Surface')

# N residues per exposure with high/low conservation
(avs
 |> filter(PDB_jury_column_rsa_unb %in% c('Core', 'Part. Exposed', 'Surface'))
 |> group_by(shenkin_nrank < 0.5, PDB_jury_column_rsa_unb)
 |> summarise(n())
 )
```

    ## `summarise()` has grouped output by 'shenkin_nrank < 0.5'. You can override
    ## using the `.groups` argument.

    ## # A tibble: 6 × 3
    ## # Groups:   shenkin_nrank < 0.5 [2]
    ##   `shenkin_nrank < 0.5` PDB_jury_column_rsa_unb `n()`
    ##   <lgl>                 <fct>                   <int>
    ## 1 FALSE                 Core                     7320
    ## 2 FALSE                 Part. Exposed           10403
    ## 3 FALSE                 Surface                 33745
    ## 4 TRUE                  Core                    17896
    ## 5 TRUE                  Part. Exposed           13645
    ## 6 TRUE                  Surface                 24723

``` r
# N residues per exposure with high/low MES
(avs
 |> filter(PDB_jury_column_rsa_unb %in% c('Core', 'Part. Exposed', 'Surface'))
 |> group_by(mes_or < 1, PDB_jury_column_rsa_unb)
 |> summarise(n())
 )
```

    ## `summarise()` has grouped output by 'mes_or < 1'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 6 × 3
    ## # Groups:   mes_or < 1 [2]
    ##   `mes_or < 1` PDB_jury_column_rsa_unb `n()`
    ##   <lgl>        <fct>                   <int>
    ## 1 FALSE        Core                     9574
    ## 2 FALSE        Part. Exposed           10669
    ## 3 FALSE        Surface                 31593
    ## 4 TRUE         Core                    15642
    ## 5 TRUE         Part. Exposed           13379
    ## 6 TRUE         Surface                 26875

# Figure 2A: Residue Shenking/MES ECDFs for each exposure class (manuscript figure)

``` r
COLOUR_SCHEME <- c(
  'partex_cum_frac' = '#56B4E9',
  'core_cum_frac' = '#E69F00',
  'surface_cum_frac' = '#009E73'
)

(avs
  |> filter(PDB_jury_column_rsa_unb %in% RSA_LEVELS)
  |> filter(occupancy_gteq_threshold)
  |> filter(occupancy > 0)
  |> filter(occupancy >= 20)
  # |> filter(occupancy < 50). # variation with bounded occupancy window
  |> arrange(shenkin_nrank)
  |> select(family, column, shenkin_nrank, PDB_jury_column_rsa_unb)
  |> pivot_wider(names_from = PDB_jury_column_rsa_unb,
                 values_from = PDB_jury_column_rsa_unb,
                 values_fn = function(x) return(1),
                 values_fill = 0)
  |> mutate(core_cum_frac = cumsum(Core) / sum(Core),
            partex_cum_frac = cumsum(`Part. Exposed`) / sum(`Part. Exposed`),
            surface_cum_frac = cumsum(Surface) / sum(Surface))
  |> pivot_longer(core_cum_frac:surface_cum_frac)
  |> ggplot()
  + aes(x = shenkin_nrank, y = value*100, colour = name)
  + geom_abline(slope = 100, intercept = 0, lty = 2)
  + geom_line()
  
  # annotation
  + xlab('Normalised Shenkin')
  + ylab('Proportion of positions (%)')

  # plot theme
  + scale_colour_manual(values = COLOUR_SCHEME)
  + theme_minimal()
  + theme(aspect.ratio = 1)
  + theme(
    axis.line = element_line(color='black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.1),
    legend.position = "none"
  )
)
```

![](Fig2A_RSA-Class-Distribution_files/figure-gfm/Figure%202A%20(right%20panel):%20ECDF%20of%20residue%20Shenkin%20for%20different%20exposure%20classes-1.png)<!-- -->

``` r
ggsave('plot2a_right.svg')
```

    ## Saving 7 x 5 in image

``` r
# ggsave('plot2a_right_supp-occ200.svg')
```

``` r
(avs
  |> filter(PDB_jury_column_rsa_unb %in% RSA_LEVELS)
  |> filter(occupancy_gteq_threshold)
  |> filter(occupancy > 0)
  |> filter(occupancy >= 20)
  |> arrange(mes_or_nrank)
  |> select(family, column, mes_or_nrank, PDB_jury_column_rsa_unb)
  |> pivot_wider(names_from = PDB_jury_column_rsa_unb,
                 values_from = PDB_jury_column_rsa_unb,
                 values_fn = function(x) return(1),
                 values_fill = 0)
  |> mutate(core_cum_frac = cumsum(Core) / sum(Core),
            partex_cum_frac = cumsum(`Part. Exposed`) / sum(`Part. Exposed`),
            surface_cum_frac = cumsum(Surface) / sum(Surface))
  |> pivot_longer(core_cum_frac:surface_cum_frac)
  |> ggplot()
  + aes(x = mes_or_nrank, y = value*100, colour = name)
  + geom_abline(slope = 100, intercept = 0, lty = 2)
  + geom_line()
 
   # annotation
  + xlab('Normalised MES')
  + ylab('Proportion of positions (%)')

  # plot theme
  + scale_colour_manual(values = COLOUR_SCHEME)
  + theme_minimal()
  + theme(aspect.ratio = 1)
  + theme(
    axis.line = element_line(color='black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.1),
    legend.position = "none"
  )
)
```

![](Fig2A_RSA-Class-Distribution_files/figure-gfm/Figure%202A%20(left%20panel):%20ECDF%20of%20residue%20MES%20for%20different%20exposure%20classes-1.png)<!-- -->

``` r
ggsave('plot2a_left.svg')
```

    ## Saving 7 x 5 in image

``` r
# ggsave('plot2a_left_supp-occ200.svg')
```

# Calculating some quoted figures in the manuscript

``` r
(
avs
  |> filter(PDB_jury_column_rsa_unb %in% RSA_LEVELS)
  |> filter(occupancy_gteq_threshold)
  |> filter(occupancy > 0)
  |> filter(occupancy > 20)
  |> arrange(mes_or_nrank)
  |> select(family, column, mes_or_nrank, PDB_jury_column_rsa_unb)
  |> pivot_wider(names_from = PDB_jury_column_rsa_unb,
                 values_from = PDB_jury_column_rsa_unb,
                 values_fn = function(x) return(1),
                 values_fill = 0)
  |> mutate(core_cum_frac = cumsum(Core) / sum(Core),
            partex_cum_frac = cumsum(`Part. Exposed`) / sum(`Part. Exposed`),
            surface_cum_frac = cumsum(Surface) / sum(Surface))
  |> pivot_longer(core_cum_frac:surface_cum_frac)
  |> filter(mes_or_nrank >= 0.5)
  |> group_by(name)
  |> summarise_all(min)
)
```

    ## # A tibble: 3 × 8
    ##   name             family  column mes_or_nrank Part. Expos…¹ Surface  Core value
    ##   <chr>            <chr>    <dbl>        <dbl>         <dbl>   <dbl> <dbl> <dbl>
    ## 1 core_cum_frac    PF00001      2          0.5             0       0     0 0.671
    ## 2 partex_cum_frac  PF00001      2          0.5             0       0     0 0.581
    ## 3 surface_cum_frac PF00001      2          0.5             0       0     0 0.407
    ## # … with abbreviated variable name ¹​`Part. Exposed`

``` r
(avs
  |> filter(PDB_jury_column_rsa_unb %in% RSA_LEVELS)
  |> filter(occupancy_gteq_threshold)
  |> filter(occupancy > 0)
  |> filter(occupancy >= 20)
  # |> filter(occupancy < 50)
  |> arrange(shenkin_nrank)
  |> select(family, column, shenkin_nrank, PDB_jury_column_rsa_unb)
  |> pivot_wider(names_from = PDB_jury_column_rsa_unb,
                 values_from = PDB_jury_column_rsa_unb,
                 values_fn = function(x) return(1),
                 values_fill = 0)
  |> mutate(core_cum_frac = cumsum(Core) / sum(Core),
            partex_cum_frac = cumsum(`Part. Exposed`) / sum(`Part. Exposed`),
            surface_cum_frac = cumsum(Surface) / sum(Surface))
  |> pivot_longer(core_cum_frac:surface_cum_frac)
  |> filter(shenkin_nrank >= 0.5)
  |> group_by(name)
  |> summarise_all(min)
)
```

    ## # A tibble: 3 × 8
    ##   name             family  column shenkin_nrank Part. Expo…¹  Core Surface value
    ##   <chr>            <chr>    <dbl>         <dbl>        <dbl> <dbl>   <dbl> <dbl>
    ## 1 core_cum_frac    PF00001      2           0.5            0     0       0 0.746
    ## 2 partex_cum_frac  PF00001      2           0.5            0     0       0 0.527
    ## 3 surface_cum_frac PF00001      2           0.5            0     0       0 0.292
    ## # … with abbreviated variable name ¹​`Part. Exposed`

``` r
sessionInfo()
```

    ## R version 4.1.3 (2022-03-10)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur/Monterey 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] forcats_0.5.2   stringr_1.4.1   dplyr_1.0.9     purrr_0.3.4    
    ## [5] readr_2.1.2     tidyr_1.2.0     tibble_3.1.8    ggplot2_3.3.6  
    ## [9] tidyverse_1.3.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] svglite_2.1.0       lubridate_1.8.0     assertthat_0.2.1   
    ##  [4] digest_0.6.29       utf8_1.2.2          R6_2.5.1           
    ##  [7] cellranger_1.1.0    backports_1.4.1     reprex_2.0.2       
    ## [10] evaluate_0.16       highr_0.9           httr_1.4.4         
    ## [13] pillar_1.8.1        rlang_1.0.4         googlesheets4_1.0.1
    ## [16] readxl_1.4.1        rstudioapi_0.14     rmarkdown_2.16     
    ## [19] textshaping_0.3.6   labeling_0.4.2      googledrive_2.0.0  
    ## [22] bit_4.0.4           munsell_0.5.0       broom_1.0.0        
    ## [25] compiler_4.1.3      modelr_0.1.9        xfun_0.32          
    ## [28] systemfonts_1.0.4   pkgconfig_2.0.3     htmltools_0.5.3    
    ## [31] tidyselect_1.1.2    fansi_1.0.3         crayon_1.5.1       
    ## [34] tzdb_0.3.0          dbplyr_2.2.1        withr_2.5.0        
    ## [37] grid_4.1.3          jsonlite_1.8.0      gtable_0.3.0       
    ## [40] lifecycle_1.0.1     DBI_1.1.3           magrittr_2.0.3     
    ## [43] scales_1.2.1        cli_3.3.0           stringi_1.7.8      
    ## [46] vroom_1.5.7         farver_2.1.1        fs_1.5.2           
    ## [49] xml2_1.3.3          ellipsis_0.3.2      ragg_1.2.2         
    ## [52] generics_0.1.3      vctrs_0.4.1         tools_4.1.3        
    ## [55] bit64_4.0.5         glue_1.6.2          hms_1.1.2          
    ## [58] parallel_4.1.3      fastmap_1.1.0       yaml_2.3.5         
    ## [61] colorspace_2.0-3    gargle_1.2.0        rvest_1.0.3        
    ## [64] knitr_1.40          haven_2.5.1
