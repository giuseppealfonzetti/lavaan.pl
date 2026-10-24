<!-- README.md is generated from README.Rmd. Please edit that file -->



# lavaan.pl

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/giuseppealfonzetti/lavaan.pl/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/giuseppealfonzetti/lavaan.pl/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/giuseppealfonzetti/lavaan.pl/graph/badge.svg)](https://app.codecov.io/gh/giuseppealfonzetti/lavaan.pl)
<!-- badges: end -->

The `{lavaan.pl}` package allows the estimation of confirmatory factor models for ordinal data using stochastic and numeric pairwise likelihood optimisation. 

## Installation

You can install the development version of lavaan.pl from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("giuseppealfonzetti/lavaan.pl")
```

When the package is installed, you can load it with:

``` r
library(lavaan.pl)
```

The startup message indicates that the `{lavaan.pl}` package provides its own version of the `cfa()` function, which can be used to estimate lavaan.pl models using the `lavaan` syntax.
An example is given below.

## Fitting an ordinal CFA model

Consider the data containing 25 personality items measured on 6-point scale from the International Personality Item Pool (IPIP) dataset.
The data set is available in the `{lavaan.pl}` package and can be loaded with the following code:

``` r
data("bfi")
head(bfi)
```
|A1|A2|A3|A4|A5|C1|C2|C3|C4|C5|E1|E2|E3|E4|E5|N1|N2|N3|N4|N5|O1|O2|O3|O4|O5|gender|education|age|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|--:|--:|--:|
|6|6|5|6|5|6|6|6|1|3|2|1|6|5|6|3|5|2|2|3|4|3|5|6|1|2|3|21|
|4|3|1|5|1|3|2|4|2|4|3|6|4|2|1|6|3|2|6|4|3|2|4|5|3|1|2|19|
|4|4|5|6|5|4|3|5|3|2|1|3|2|5|4|3|3|4|2|3|5|3|5|6|3|1|1|21|
|4|5|2|2|1|5|5|5|2|2|3|4|3|6|5|2|4|2|2|3|5|2|5|5|5|1|1|17|
|1|5|6|5|6|4|3|2|4|5|2|1|2|5|2|2|2|2|2|2|6|1|5|5|2|1|5|68|
|2|6|5|6|5|3|5|6|3|6|2|2|4|6|6|4|4|4|6|6|6|1|5|6|1|2|2|27|


The `{lavaan}` model syntax to fit a confirmatory factor analysis model with 5 latent factors, where the first 5 items load on the first factor and so on, is:

``` r
mod <- "
  opn =~ O1 + O2 + O3 + O4 + O5  # Openness
  con =~ C1 + C2 + C3 + C4 + C5  # Conscientiousness
  ext =~ E1 + E2 + E3 + E4 + E5  # Extraversion
  agr =~ A1 + A2 + A3 + A4 + A5  # Agreeableness
  neu =~ N1 + N2 + N3 + N4 + N5  # Neuroticism
"
```

Users of `{lavaan}` will be familiar with this syntax. 
The call to the `cfa()` function is made as usual:

``` r
# Fit the model (note, this cfa() function is from the lavaan.pl package)
fit <- cfa(model = mod, data = bfi, std.lv = TRUE, test = "mean.var.adjusted")
summary(fit)
#> lavaan 0.7-2 ended normally after 89 iterations
#> 
#>   Estimator                                        PML
#>   Optimization method                           NLMINB
#>   Number of model parameters                       160
#> 
#>   Number of observations                          2236
#> 
#> Model Test User Model:
#>                                               Standard      Scaled
#>   Test Statistic                              4584.031    7206.841
#>   Degrees of freedom                               265     333.112
#>   P-value (Unknown)                                 NA       0.000
#>   Scaling correction factor                                  0.636
#>     mean+var adjusted correction (PML)                            
#> 
#> Parameter Estimates:
#> 
#>   Parameterization                               Delta
#>   Standard errors                             Sandwich
#>   Information bread                           Observed
#>   Observed information based on                Hessian
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   opn =~                                              
#>     O1                0.647    0.025   25.909    0.000
#>     O2               -0.444    0.032  -14.065    0.000
#>     O3                0.808    0.027   30.169    0.000
#>     O4                0.136    0.042    3.253    0.001
#>     O5               -0.486    0.029  -16.593    0.000
#>   con =~                                              
#>     C1                0.564    0.026   21.442    0.000
#>     C2                0.571    0.026   21.929    0.000
#>     C3                0.543    0.024   22.851    0.000
#>     C4               -0.760    0.020  -38.438    0.000
#>     C5               -0.706    0.020  -34.670    0.000
#>   ext =~                                              
#>     E1                0.540    0.023   23.359    0.000
#>     E2                0.727    0.018   40.150    0.000
#>     E3               -0.690    0.018  -38.311    0.000
#>     E4               -0.742    0.018  -42.384    0.000
#>     E5               -0.620    0.020  -31.383    0.000
#>   agr =~                                              
#>     A1                0.335    0.027   12.363    0.000
#>     A2               -0.663    0.021  -31.722    0.000
#>     A3               -0.771    0.017  -44.760    0.000
#>     A4               -0.545    0.022  -24.412    0.000
#>     A5               -0.774    0.017  -45.158    0.000
#>   neu =~                                              
#>     N1                0.837    0.017   50.232    0.000
#>     N2                0.801    0.018   43.435    0.000
#>     N3                0.755    0.015   51.794    0.000
#>     N4                0.666    0.019   35.196    0.000
#>     N5                0.570    0.021   27.636    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   opn ~~                                              
#>     con               0.322    0.032   10.156    0.000
#>     ext              -0.466    0.031  -15.076    0.000
#>     agr              -0.269    0.032   -8.358    0.000
#>     neu              -0.146    0.034   -4.274    0.000
#>   con ~~                                              
#>     ext              -0.390    0.029  -13.477    0.000
#>     agr              -0.370    0.029  -12.717    0.000
#>     neu              -0.301    0.030  -10.199    0.000
#>   ext ~~                                              
#>     agr               0.702    0.020   34.779    0.000
#>     neu               0.267    0.028    9.453    0.000
#>   agr ~~                                              
#>     neu               0.232    0.028    8.208    0.000
#> 
#> Thresholds:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     O1|t1            -2.447    0.090  -27.339    0.000
#>     O1|t2            -1.695    0.046  -36.776    0.000
#>     O1|t3            -1.190    0.035  -34.449    0.000
#>     O1|t4            -0.430    0.027  -15.689    0.000
#>     O1|t5             0.455    0.028   16.517    0.000
#>     O2|t1            -0.565    0.028  -20.081    0.000
#>     O2|t2             0.131    0.027    4.925    0.000
#>     O2|t3             0.495    0.028   17.872    0.000
#>     O2|t4             1.011    0.032   31.622    0.000
#>     O2|t5             1.557    0.042   36.956    0.000
#>     O3|t1            -1.990    0.058  -34.468    0.000
#>     O3|t2            -1.451    0.039  -36.738    0.000
#>     O3|t3            -0.956    0.031  -30.481    0.000
#>     O3|t4            -0.122    0.027   -4.596    0.000
#>     O3|t5             0.839    0.030   27.763    0.000
#>     O4|t1            -2.187    0.069  -31.746    0.000
#>     O4|t2            -1.567    0.042  -36.916    0.000
#>     O4|t3            -1.246    0.035  -35.101    0.000
#>     O4|t4            -0.607    0.028  -21.419    0.000
#>     O4|t5             0.255    0.027    9.502    0.000
#>     O5|t1            -0.584    0.028  -20.678    0.000
#>     O5|t2             0.261    0.027    9.740    0.000
#>     O5|t3             0.771    0.029   26.137    0.000
#>     O5|t4             1.341    0.037   36.128    0.000
#>     O5|t5             1.931    0.055   35.093    0.000
#>     C1|t1            -2.011    0.059  -34.160    0.000
#>     C1|t2            -1.432    0.039  -36.693    0.000
#>     C1|t3            -0.960    0.031  -30.619    0.000
#>     C1|t4            -0.269    0.027  -10.017    0.000
#>     C1|t5             0.736    0.029   25.105    0.000
#>     C2|t1            -1.880    0.053  -35.553    0.000
#>     C2|t2            -1.203    0.035  -34.609    0.000
#>     C2|t3            -0.785    0.030  -26.467    0.000
#>     C2|t4            -0.146    0.027   -5.501    0.000
#>     C2|t5             0.828    0.030   27.509    0.000
#>     C3|t1            -1.887    0.053  -35.472    0.000
#>     C3|t2            -1.191    0.035  -34.479    0.000
#>     C3|t3            -0.778    0.030  -26.293    0.000
#>     C3|t4            -0.022    0.026   -0.850    0.396
#>     C3|t5             0.921    0.031   29.700    0.000
#>     C4|t1            -0.550    0.028  -19.594    0.000
#>     C4|t2             0.208    0.027    7.822    0.000
#>     C4|t3             0.664    0.029   23.227    0.000
#>     C4|t4             1.289    0.036   35.748    0.000
#>     C4|t5             2.033    0.060   33.932    0.000
#>     C5|t1            -0.886    0.031  -28.854    0.000
#>     C5|t2            -0.254    0.027   -9.472    0.000
#>     C5|t3             0.047    0.026    1.791    0.073
#>     C5|t4             0.630    0.028   22.187    0.000
#>     C5|t5             1.303    0.036   35.771    0.000
#>     E1|t1            -0.728    0.029  -24.896    0.000
#>     E1|t2            -0.063    0.026   -2.371    0.018
#>     E1|t3             0.311    0.027   11.577    0.000
#>     E1|t4             0.789    0.030   26.609    0.000
#>     E1|t5             1.373    0.038   36.231    0.000
#>     E2|t1            -0.867    0.031  -28.386    0.000
#>     E2|t2            -0.153    0.027   -5.757    0.000
#>     E2|t3             0.169    0.027    6.357    0.000
#>     E2|t4             0.754    0.029   25.754    0.000
#>     E2|t5             1.320    0.037   35.904    0.000
#>     E3|t1            -1.625    0.044  -36.954    0.000
#>     E3|t2            -1.010    0.032  -31.607    0.000
#>     E3|t3            -0.518    0.028  -18.684    0.000
#>     E3|t4             0.267    0.027    9.966    0.000
#>     E3|t5             1.148    0.034   33.814    0.000
#>     E4|t1            -1.628    0.044  -36.986    0.000
#>     E4|t2            -1.061    0.033  -32.560    0.000
#>     E4|t3            -0.716    0.029  -24.681    0.000
#>     E4|t4            -0.263    0.027   -9.824    0.000
#>     E4|t5             0.640    0.029   22.399    0.000
#>     E5|t1            -1.813    0.050  -36.110    0.000
#>     E5|t2            -1.199    0.035  -34.583    0.000
#>     E5|t3            -0.798    0.030  -26.829    0.000
#>     E5|t4            -0.168    0.027   -6.320    0.000
#>     E5|t5             0.790    0.030   26.546    0.000
#>     A1|t1            -0.403    0.027  -14.747    0.000
#>     A1|t2             0.363    0.027   13.379    0.000
#>     A1|t3             0.780    0.030   26.386    0.000
#>     A1|t4             1.270    0.036   35.515    0.000
#>     A1|t5             1.888    0.053   35.548    0.000
#>     A2|t1            -2.156    0.067  -32.227    0.000
#>     A2|t2            -1.558    0.042  -37.111    0.000
#>     A2|t3            -1.223    0.035  -34.998    0.000
#>     A2|t4            -0.513    0.028  -18.491    0.000
#>     A2|t5             0.456    0.028   16.529    0.000
#>     A3|t1            -1.849    0.051  -35.946    0.000
#>     A3|t2            -1.332    0.037  -36.154    0.000
#>     A3|t3            -0.995    0.032  -31.405    0.000
#>     A3|t4            -0.350    0.027  -12.948    0.000
#>     A3|t5             0.593    0.028   20.950    0.000
#>     A4|t1            -1.716    0.047  -36.670    0.000
#>     A4|t2            -1.206    0.035  -34.690    0.000
#>     A4|t3            -0.927    0.031  -29.909    0.000
#>     A4|t4            -0.423    0.027  -15.474    0.000
#>     A4|t5             0.202    0.027    7.556    0.000
#>     A5|t1            -2.015    0.059  -34.227    0.000
#>     A5|t2            -1.364    0.037  -36.442    0.000
#>     A5|t3            -0.946    0.031  -30.442    0.000
#>     A5|t4            -0.277    0.027  -10.339    0.000
#>     A5|t5             0.652    0.029   22.691    0.000
#>     N1|t1            -0.711    0.029  -24.451    0.000
#>     N1|t2            -0.060    0.026   -2.281    0.023
#>     N1|t3             0.335    0.027   12.435    0.000
#>     N1|t4             0.892    0.031   29.085    0.000
#>     N1|t5             1.489    0.040   36.776    0.000
#>     N2|t1            -1.163    0.034  -34.002    0.000
#>     N2|t2            -0.473    0.028  -17.154    0.000
#>     N2|t3            -0.094    0.026   -3.555    0.000
#>     N2|t4             0.570    0.028   20.350    0.000
#>     N2|t5             1.259    0.036   35.233    0.000
#>     N3|t1            -0.924    0.031  -29.730    0.000
#>     N3|t2            -0.210    0.027   -7.896    0.000
#>     N3|t3             0.106    0.026    4.000    0.000
#>     N3|t4             0.685    0.029   23.769    0.000
#>     N3|t5             1.351    0.037   36.094    0.000
#>     N4|t1            -0.961    0.032  -30.494    0.000
#>     N4|t2            -0.235    0.027   -8.788    0.000
#>     N4|t3             0.155    0.027    5.819    0.000
#>     N4|t4             0.768    0.029   26.054    0.000
#>     N4|t5             1.341    0.037   36.051    0.000
#>     N5|t1            -0.714    0.029  -24.521    0.000
#>     N5|t2            -0.035    0.026   -1.337    0.181
#>     N5|t3             0.297    0.027   11.058    0.000
#>     N5|t4             0.825    0.030   27.481    0.000
#>     N5|t5             1.358    0.038   36.138    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .O1                0.581                           
#>    .O2                0.803                           
#>    .O3                0.347                           
#>    .O4                0.982                           
#>    .O5                0.764                           
#>    .C1                0.681                           
#>    .C2                0.674                           
#>    .C3                0.705                           
#>    .C4                0.422                           
#>    .C5                0.502                           
#>    .E1                0.708                           
#>    .E2                0.471                           
#>    .E3                0.524                           
#>    .E4                0.450                           
#>    .E5                0.615                           
#>    .A1                0.888                           
#>    .A2                0.561                           
#>    .A3                0.406                           
#>    .A4                0.703                           
#>    .A5                0.401                           
#>    .N1                0.300                           
#>    .N2                0.359                           
#>    .N3                0.429                           
#>    .N4                0.556                           
#>    .N5                0.675                           
#>     opn               1.000                           
#>     con               1.000                           
#>     ext               1.000                           
#>     agr               1.000                           
#>     neu               1.000                           
#> 
```


All^[Most :)--we're working on it!] `lavaan` methods should carry over, such as:

``` r
head(coef(fit), 10)
#>    opn=~O1    opn=~O2    opn=~O3    opn=~O4    opn=~O5    con=~C1    con=~C2    con=~C3 
#>  0.6474239 -0.4443221  0.8077954  0.1355670 -0.4862491  0.5644095  0.5707376  0.5430038 
#>    con=~C4    con=~C5 
#> -0.7602739 -0.7055655 
head(parameterestimates(fit, se = TRUE))
```
|lhs|op|rhs|est|se|z|pvalue|ci.lower|ci.upper|
|---|---|---|--:|--:|--:|--:|--:|--:|
|opn|=~|O1| 0.647|0.025| 25.909|0.000| 0.598| 0.696|
|opn|=~|O2|-0.444|0.032|-14.065|0.000|-0.506|-0.382|
|opn|=~|O3| 0.808|0.027| 30.169|0.000| 0.755| 0.860|
|opn|=~|O4| 0.136|0.042|  3.253|0.001| 0.054| 0.217|
|opn|=~|O5|-0.486|0.029|-16.593|0.000|-0.544|-0.429|
|con|=~|C1| 0.564|0.026| 21.442|0.000| 0.513| 0.616|


``` r
head(fitted(fit)$cov)  # model-implied covariance matrix for the UVs
```
| |O1|O2|O3|O4|O5|C1|C2|C3|C4|C5|E1|E2|E3|E4|E5|A1|A2|A3|A4|A5|N1|N2|N3|N4|N5|
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|O1| 1.000|-0.288| 0.523| 0.088|-0.315| 0.118| 0.119| 0.113|-0.159|-0.147|-0.163|-0.219| 0.208| 0.224| 0.187|-0.058| 0.116| 0.134| 0.095| 0.135|-0.079|-0.076|-0.072|-0.063|-0.054|
|O2|-0.288| 1.000|-0.359|-0.060| 0.216|-0.081|-0.082|-0.078| 0.109| 0.101| 0.112| 0.151|-0.143|-0.154|-0.128| 0.040|-0.079|-0.092|-0.065|-0.093| 0.054| 0.052| 0.049| 0.043| 0.037|
|O3| 0.523|-0.359| 1.000| 0.110|-0.393| 0.147| 0.149| 0.141|-0.198|-0.184|-0.203|-0.274| 0.260| 0.279| 0.233|-0.073| 0.144| 0.168| 0.119| 0.168|-0.099|-0.095|-0.089|-0.079|-0.067|
|O4| 0.088|-0.060| 0.110| 1.000|-0.066| 0.025| 0.025| 0.024|-0.033|-0.031|-0.034|-0.046| 0.044| 0.047| 0.039|-0.012| 0.024| 0.028| 0.020| 0.028|-0.017|-0.016|-0.015|-0.013|-0.011|
|O5|-0.315| 0.216|-0.393|-0.066| 1.000|-0.088|-0.089|-0.085| 0.119| 0.111| 0.122| 0.165|-0.156|-0.168|-0.141| 0.044|-0.087|-0.101|-0.071|-0.101| 0.060| 0.057| 0.054| 0.047| 0.041|
|C1| 0.118|-0.081| 0.147| 0.025|-0.088| 1.000| 0.322| 0.306|-0.429|-0.398|-0.119|-0.160| 0.152| 0.163| 0.137|-0.070| 0.138| 0.161| 0.114| 0.161|-0.142|-0.136|-0.128|-0.113|-0.097|


And also plot from the `{semPlot}` package:

``` r
library(semPlot)

grouplist <- list(
  lat = c("opn", "con", "ext", "agr", "neu"),
  man = rep(c("opn", "con", "ext", "agr", "neu"), each = 5)
)

# Define a named vector of colors for each group (latent factor)
groupcols <- RColorBrewer::brewer.pal(5, "Set2") |>
  rep(each = 5) |>
  adjustcolor(alpha.f = 0.5)

# Now call semPaths
semPaths(fit, whatLabels = "est", layout = "circle", residuals = FALSE, 
         intercept = FALSE, groups = grouplist, color = groupcols, 
         legend = FALSE)
```
![](<man/figures/README-fig-bfi-1.png>)

## Session information

``` r
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Golden Gate 27.0.1
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Europe/Rome
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] semPlot_1.2.0        lavaan_0.7-2         lavaan.pl_0.1.0.9002
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        xfun_0.60           ggplot2_4.0.3       htmlwidgets_1.6.4  
#>  [5] psych_2.6.5         processx_3.9.0      lattice_0.23-1      callr_3.8.0        
#>  [9] numDeriv_2016.8-1.1 quadprog_1.5-8      vctrs_0.7.3         tools_4.6.1        
#> [13] pak_0.11.1          ps_1.9.3            generics_0.1.4      stats4_4.6.1       
#> [17] parallel_4.6.1      tibble_3.3.1        cluster_2.1.8.2     ucminf_1.2.3       
#> [21] pkgconfig_2.0.3     qgraph_1.10.1       Matrix_1.7-6        data.table_1.18.6.1
#> [25] checkmate_2.3.4     RColorBrewer_1.1-3  S7_0.2.2            RcppParallel_6.2.1 
#> [29] lifecycle_1.0.5     compiler_4.6.1      farver_2.1.2        stringr_1.6.0      
#> [33] mnormt_2.1.2        codetools_0.2-20    litedown_0.11       glasso_1.11        
#> [37] htmltools_0.5.9     fdrtool_1.2.18      htmlTable_2.5.0     Formula_1.2-6      
#> [41] pillar_1.11.1       MASS_7.3-65         Hmisc_5.2-6         abind_1.4-8        
#> [45] rpart_4.1.27        nlme_3.1-169        gtools_3.9.5        commonmark_2.0.0   
#> [49] tidyselect_1.2.1    digest_0.6.39       stringi_1.8.9       reshape2_1.4.5     
#> [53] dplyr_1.2.1         fastmap_1.2.0       grid_4.6.1          colorspace_2.1-3   
#> [57] cli_3.6.6           magrittr_2.0.5      base64enc_0.1-6     pbivnorm_0.6.0     
#> [61] foreign_0.8-91      corpcor_1.6.10      scales_1.4.0        backports_1.5.1    
#> [65] rmarkdown_2.31      jpeg_0.1-11         igraph_2.3.3        otel_0.2.0         
#> [69] nnet_7.3-20         gridExtra_2.3.1     RcppEigen_0.3.4.0.2 png_0.1-9          
#> [73] pbapply_1.7-4       evaluate_1.0.5      knitr_1.51          rlang_1.3.0        
#> [77] Rcpp_1.1.2          glue_1.8.1          rstudioapi_0.19.0   R6_2.6.1           
#> [81] plyr_1.8.9         
```

<!-- To fit a latent factor analysis model with 5 factors, where the first 5 items load on the first factor and so on, we can use the `fit_plFA()` function. -->
<!-- In preparation, we must build the following objects to be supplied as a list to the `CONSTR_LIST` argument: -->

<!-- - `CONSTRMAT`: A matrix of constraints on the factor loadings. -->
<!-- - `CORRFLAG`: A flag indicating whether the factor correlations are estimated. -->
<!-- - `CONSTRVAR`: A vector of constraints on the factor variances, if any. -->
<!-- - `STDLV`: A flag indicating whether the factor loadings are standardised. -->
<!-- - `LLC`: A list indicating linear constraints on factor loading estimations, if any. -->

<!-- The above ingredients are prepared using the code below. -->

<!-- ```{r} -->
<!-- p <- 25         # number of items -->
<!-- q <- 5          # number of factors -->
<!-- n <- nrow(bfi)  # number of observations -->

<!-- # Convert data frame to matrix, with each item taking values 0,1,2,...,5 -->
<!-- D <- as.matrix(data.frame(lapply(bfi[, 1:p], as.numeric))) - 1 -->
<!-- head(D) -->

<!-- # Constraint matrix: NA indicate free loadings -->
<!-- constrmat <- build_constrMat(P = p, Q = q, STRUCT = "simple") -->
<!-- head(constrmat, 12) -->

<!-- # Build constraint list -->
<!-- constr_list <- list( -->
<!--   CONSTRMAT = constrmat, -->
<!--   CORRFLAG = TRUE, -->
<!--   CONSTRVAR = rep(1, q),  # these options standardise the   -->
<!--   STDLV = TRUE,           # factor variances -->
<!--   LLC = NULL -->
<!-- ) -->
<!-- ``` -->

<!-- We can now fit the factor model: -->

<!-- ```{r} -->
<!-- fit <- fit_plFA( -->
<!--   DATA = D,  -->
<!--   CONSTR_LIST = constr_list, -->
<!--   METHOD = "ucminf", -->
<!--   VERBOSE = TRUE -->
<!-- ) -->
<!-- print(fit) -->
<!-- ``` -->

<!-- The default method is using the optimiser from the `{ucminf}` package, which works well for small data sizes. -->
<!-- Alternatively, the user may choose `METHOD = "SA"` to turn on the stochastic approximation algorithm. -->
<!-- Use the function `getPar()` to extract the parameter estimates: -->

<!-- ```{r} -->
<!-- getPar(fit) -->
<!-- ``` -->

<!-- Computing the standard errors is done separately using the `computeVar()` function: -->

<!-- ```{r} -->
<!-- var <- computeVar(OBJ = fit, DATA = D, VERBOSE = TRUE) -->
<!-- str(var) -->

<!-- # Standard errors -->
<!-- head(sqrt(diag(var$vcov) / n)) -->
<!-- ``` -->

<!-- ## `{lavaan}` wrapper -->

<!-- This package also provides a user-friendly `{lavaan}`-style interface to estimate factor analysis models using the `cfa()` function. -->
<!-- Let's rerun the above example, but this time using the stochastic approximation algorithm. -->
