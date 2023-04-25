Figures 1C-D: Frequency distribution of variants in MSAs
================

This notebook creates plots from pre-computed missense variant
frequencies over MSA columns in the Pfam-gnomAD dataset.

The first plot shows the frequency distributions of the number of
missense variants in Pfam MSA columns in Pfams with different numbers of
human sequences.

The second plot shows the relationship between the number of missense
and synonymous variants in an MSA column and the Shenkin diversity in
SH2 domains.

``` r
suppressPackageStartupMessages(library(tidyverse))


# Load and process dataset
pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv <- read_csv("~/Library/CloudStorage/OneDrive-UniversityofDundee/Projects/VarAlign-2022/data/pfam-gnomAD-clinvar-pdb-colstats_c7c3e19.csv.gz",
comment = "#", show_col_types = FALSE)
avs <- pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv  # avs: alignments variants structure

# Filter unneccesary columns
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
  "shenkin_nrank",
  "mes_or_log"
)
blacklist <- !names(avs) %in% whitelist
avs[, blacklist] <- NULL

# Fix variable name
names(avs)[names(avs) == '_missense_all'] <- 'missense_all'

# Calculate total synonymous count
avs[, 'synonymous_all'] <- avs[, 'splice_region_variant&synonymous_variant'] +
  avs[, 'synonymous_variant']

# Filters
avs <- avs |> filter(occupancy > 0)  # must have human sequences

# Clear raw dataset from memory
rm(pfam_gnomAD_clinvar_pdb_colstats_c7c3e19_csv)
```

``` r
(avs
 |> filter(occupancy %in% c(1, 5, 9))  # N human sequences to be plotted, filters columns
 # |> sample_n(10000)  # for prototyping
 |> group_by(occupancy, missense_all)
 |> summarise(n_cols = n())
 -> avs_missense_freqs_by_nhuman
 )
```

    ## `summarise()` has grouped output by 'occupancy'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 48 × 3
    ## # Groups:   occupancy [3]
    ##    occupancy missense_all n_cols
    ##        <dbl>        <dbl>  <int>
    ##  1         1            0 391034
    ##  2         1            1 183476
    ##  3         1            2  51035
    ##  4         1            3   9143
    ##  5         1            4   1163
    ##  6         1            5     96
    ##  7         1            6     14
    ##  8         1            7      3
    ##  9         1            8      1
    ## 10         1            9      1
    ## # … with 38 more rows

``` r
(avs_missense_freqs_by_nhuman
  # data
  |> filter(occupancy %in% c(1, 5, 9))
  
  # plot
  |> ggplot()
  + aes(x = missense_all, y = n_cols)
  + geom_bar(stat = 'identity')
  + facet_wrap(~ occupancy, scales = 'free_y')
  
  # annotation
  + xlab('# Missense')
  # + ylab('# Pfam Columns')
  + scale_y_continuous(name = '# Pfam Columns', labels = scales::comma)
  
  # plot theme
  + theme_minimal()
  + theme(aspect.ratio = 1)
  + theme(
    axis.line = element_line(color='black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.1),
    panel.grid.major.x = element_blank(),
    
    # hide facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    )
  )
```

![](Figs1C-D_variant-totals-in-MSA-columns_files/figure-gfm/Figure%201C:%20Column%20missense%20frequency%20histograms%20by%20N%20human%20sequences-1.png)<!-- -->

``` r
ggsave('plot1c.svg')
```

    ## Saving 7 x 5 in image

``` r
# Graph to help pick occupancy threshold for the figure
# This is required because Shenkin is strongly dependent on occupancy in the
# range from 1 - ~100, which applies to extremely gappy regions found in Pfams.
# n.b. The regression model handles this specifically but here I think its best
# to filter these artefacts to focus on the informative region of the data.
(avs
 |> filter(family == 'PF00017')
 |> ggplot()
 + aes(x = occupancy, y = shenkin)
 + geom_point()
 + geom_vline(xintercept = 90, lty = 2)
 + theme_minimal()
 + theme(aspect.ratio = 1)
 + theme(axis.line = element_line(color='black'))
)
```

![](Figs1C-D_variant-totals-in-MSA-columns_files/figure-gfm/Figure%201D:%20Column%20variant%20totals%20vs.%20Shenkin%20in%20SH2%20domains-1.png)<!-- -->

``` r
(avs
 |> filter(family == 'PF00017')
 |> filter(occupancy >= 90)
 |> select(column, occupancy, shenkin, missense_all, synonymous_all)
 |> pivot_longer(missense_all:synonymous_all, names_to = 'variant', values_to = 'n')
  
  # Plot
 |> ggplot()
 + aes(x = shenkin, y = n)
 + facet_grid(~ variant)
 + geom_smooth(formula = y ~ x, method = 'lm', colour = 'black')
 + geom_point()
  
  # Annotation
  + scale_x_continuous(name = 'Shenkin', breaks = c(0, 25, 50, 75, 100),
                       labels = c(0, 25, 50, 75, 100), limits = c(0, 100))
  + ylab('# Missense')
 
  # Theme
  + theme_minimal()
  + theme(aspect.ratio = 1)
  + theme(
    axis.line = element_line(color='black'),  # heavy line for axes
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    
    # hide facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    )
 )
```

![](Figs1C-D_variant-totals-in-MSA-columns_files/figure-gfm/Figure%201D:%20Column%20variant%20totals%20vs.%20Shenkin%20in%20SH2%20domains-2.png)<!-- -->

``` r
ggsave('plot1d.svg')
```

    ## Saving 7 x 5 in image

``` r
(avs
 |> filter(family == 'PF00017')
 |> group_by(occupancy >= 90, mes_or < 1, mes_p < 0.1)
 |> summarise(max(mes_or_log))
 )
```

    ## `summarise()` has grouped output by 'occupancy >= 90', 'mes_or < 1'. You can
    ## override using the `.groups` argument.

    ## # A tibble: 7 × 4
    ## # Groups:   occupancy >= 90, mes_or < 1 [4]
    ##   `occupancy >= 90` `mes_or < 1` `mes_p < 0.1` `max(mes_or_log)`
    ##   <lgl>             <lgl>        <lgl>                     <dbl>
    ## 1 FALSE             FALSE        FALSE                    1.60  
    ## 2 FALSE             FALSE        TRUE                     2.00  
    ## 3 FALSE             TRUE         FALSE                   -0.0780
    ## 4 TRUE              FALSE        FALSE                    0.270 
    ## 5 TRUE              FALSE        TRUE                     0.467 
    ## 6 TRUE              TRUE         FALSE                   -0.0136
    ## 7 TRUE              TRUE         TRUE                    -0.435
