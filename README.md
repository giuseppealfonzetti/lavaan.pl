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
# install.packages("pak")
pak::pak("giuseppealfonzetti/lavaan.pl")
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
#> lavaan 0.7-2 ended normally after 1292 iterations
#> 
#>   Estimator                                        PML
#>   Optimization method                               SA
#>   Number of model parameters                       160
#> 
#>   Number of observations                          2236
#> 
#> Model Test User Model:
#>                                               Standard      Scaled
#>   Test Statistic                              4596.077    7480.325
#>   Degrees of freedom                               265     343.773
#>   P-value (Unknown)                                 NA       0.000
#>   Scaling correction factor                                  0.614
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
#>     O1                0.647    0.027   24.302    0.000
#>     O2               -0.445    0.033  -13.345    0.000
#>     O3                0.804    0.029   27.835    0.000
#>     O4                0.142    0.045    3.140    0.002
#>     O5               -0.480    0.031  -15.455    0.000
#>   con =~                                              
#>     C1                0.562    0.028   20.332    0.000
#>     C2                0.587    0.027   21.551    0.000
#>     C3                0.554    0.024   23.107    0.000
#>     C4               -0.759    0.020  -37.325    0.000
#>     C5               -0.698    0.021  -32.668    0.000
#>   ext =~                                              
#>     E1                0.540    0.024   22.140    0.000
#>     E2                0.718    0.020   36.590    0.000
#>     E3               -0.692    0.019  -36.434    0.000
#>     E4               -0.737    0.019  -38.615    0.000
#>     E5               -0.617    0.021  -29.929    0.000
#>   agr =~                                              
#>     A1                0.321    0.029   11.230    0.000
#>     A2               -0.656    0.022  -29.455    0.000
#>     A3               -0.756    0.018  -41.400    0.000
#>     A4               -0.547    0.023  -23.998    0.000
#>     A5               -0.786    0.018  -43.629    0.000
#>   neu =~                                              
#>     N1                0.843    0.022   38.776    0.000
#>     N2                0.797    0.023   34.417    0.000
#>     N3                0.762    0.017   45.758    0.000
#>     N4                0.665    0.022   30.436    0.000
#>     N5                0.577    0.021   26.932    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   opn ~~                                              
#>     con               0.328    0.033   10.009    0.000
#>     ext              -0.466    0.034  -13.819    0.000
#>     agr              -0.269    0.033   -8.082    0.000
#>     neu              -0.143    0.035   -4.096    0.000
#>   con ~~                                              
#>     ext              -0.400    0.030  -13.250    0.000
#>     agr              -0.370    0.029  -12.558    0.000
#>     neu              -0.304    0.031   -9.924    0.000
#>   ext ~~                                              
#>     agr               0.720    0.022   32.475    0.000
#>     neu               0.272    0.030    9.210    0.000
#>   agr ~~                                              
#>     neu               0.235    0.029    8.201    0.000
#> 
#> Thresholds:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     O1|t1            -2.447    0.090  -27.345    0.000
#>     O1|t2            -1.695    0.046  -36.779    0.000
#>     O1|t3            -1.190    0.035  -34.450    0.000
#>     O1|t4            -0.430    0.027  -15.689    0.000
#>     O1|t5             0.455    0.028   16.517    0.000
#>     O2|t1            -0.565    0.028  -20.081    0.000
#>     O2|t2             0.131    0.027    4.924    0.000
#>     O2|t3             0.495    0.028   17.872    0.000
#>     O2|t4             1.011    0.032   31.622    0.000
#>     O2|t5             1.557    0.042   36.957    0.000
#>     O3|t1            -1.990    0.058  -34.485    0.000
#>     O3|t2            -1.451    0.039  -36.746    0.000
#>     O3|t3            -0.956    0.031  -30.483    0.000
#>     O3|t4            -0.122    0.027   -4.597    0.000
#>     O3|t5             0.839    0.030   27.763    0.000
#>     O4|t1            -2.187    0.069  -31.747    0.000
#>     O4|t2            -1.567    0.042  -36.918    0.000
#>     O4|t3            -1.246    0.035  -35.103    0.000
#>     O4|t4            -0.607    0.028  -21.420    0.000
#>     O4|t5             0.255    0.027    9.502    0.000
#>     O5|t1            -0.584    0.028  -20.677    0.000
#>     O5|t2             0.261    0.027    9.742    0.000
#>     O5|t3             0.771    0.029   26.138    0.000
#>     O5|t4             1.341    0.037   36.135    0.000
#>     O5|t5             1.931    0.055   35.116    0.000
#>     C1|t1            -2.011    0.059  -34.165    0.000
#>     C1|t2            -1.432    0.039  -36.693    0.000
#>     C1|t3            -0.960    0.031  -30.619    0.000
#>     C1|t4            -0.269    0.027  -10.017    0.000
#>     C1|t5             0.736    0.029   25.104    0.000
#>     C2|t1            -1.880    0.053  -35.533    0.000
#>     C2|t2            -1.203    0.035  -34.604    0.000
#>     C2|t3            -0.785    0.030  -26.469    0.000
#>     C2|t4            -0.146    0.027   -5.503    0.000
#>     C2|t5             0.828    0.030   27.505    0.000
#>     C3|t1            -1.886    0.053  -35.466    0.000
#>     C3|t2            -1.191    0.035  -34.479    0.000
#>     C3|t3            -0.778    0.030  -26.297    0.000
#>     C3|t4            -0.023    0.026   -0.854    0.393
#>     C3|t5             0.921    0.031   29.695    0.000
#>     C4|t1            -0.551    0.028  -19.598    0.000
#>     C4|t2             0.208    0.027    7.813    0.000
#>     C4|t3             0.664    0.029   23.220    0.000
#>     C4|t4             1.289    0.036   35.732    0.000
#>     C4|t5             2.034    0.060   33.870    0.000
#>     C5|t1            -0.886    0.031  -28.854    0.000
#>     C5|t2            -0.254    0.027   -9.472    0.000
#>     C5|t3             0.047    0.026    1.790    0.073
#>     C5|t4             0.630    0.028   22.184    0.000
#>     C5|t5             1.303    0.036   35.768    0.000
#>     E1|t1            -0.728    0.029  -24.894    0.000
#>     E1|t2            -0.063    0.026   -2.369    0.018
#>     E1|t3             0.312    0.027   11.581    0.000
#>     E1|t4             0.790    0.030   26.611    0.000
#>     E1|t5             1.373    0.038   36.230    0.000
#>     E2|t1            -0.867    0.031  -28.385    0.000
#>     E2|t2            -0.153    0.027   -5.753    0.000
#>     E2|t3             0.169    0.027    6.358    0.000
#>     E2|t4             0.754    0.029   25.755    0.000
#>     E2|t5             1.320    0.037   35.910    0.000
#>     E3|t1            -1.625    0.044  -36.948    0.000
#>     E3|t2            -1.010    0.032  -31.607    0.000
#>     E3|t3            -0.518    0.028  -18.686    0.000
#>     E3|t4             0.267    0.027    9.966    0.000
#>     E3|t5             1.148    0.034   33.809    0.000
#>     E4|t1            -1.628    0.044  -36.996    0.000
#>     E4|t2            -1.061    0.033  -32.563    0.000
#>     E4|t3            -0.716    0.029  -24.683    0.000
#>     E4|t4            -0.263    0.027   -9.826    0.000
#>     E4|t5             0.640    0.029   22.398    0.000
#>     E5|t1            -1.813    0.050  -36.099    0.000
#>     E5|t2            -1.199    0.035  -34.576    0.000
#>     E5|t3            -0.798    0.030  -26.826    0.000
#>     E5|t4            -0.168    0.027   -6.314    0.000
#>     E5|t5             0.790    0.030   26.548    0.000
#>     A1|t1            -0.404    0.027  -14.751    0.000
#>     A1|t2             0.363    0.027   13.369    0.000
#>     A1|t3             0.780    0.030   26.374    0.000
#>     A1|t4             1.270    0.036   35.500    0.000
#>     A1|t5             1.889    0.053   35.529    0.000
#>     A2|t1            -2.157    0.067  -32.187    0.000
#>     A2|t2            -1.559    0.042  -37.092    0.000
#>     A2|t3            -1.223    0.035  -34.990    0.000
#>     A2|t4            -0.513    0.028  -18.485    0.000
#>     A2|t5             0.456    0.028   16.531    0.000
#>     A3|t1            -1.849    0.051  -35.965    0.000
#>     A3|t2            -1.332    0.037  -36.157    0.000
#>     A3|t3            -0.995    0.032  -31.404    0.000
#>     A3|t4            -0.350    0.027  -12.947    0.000
#>     A3|t5             0.593    0.028   20.950    0.000
#>     A4|t1            -1.716    0.047  -36.656    0.000
#>     A4|t2            -1.206    0.035  -34.686    0.000
#>     A4|t3            -0.927    0.031  -29.908    0.000
#>     A4|t4            -0.423    0.027  -15.473    0.000
#>     A4|t5             0.202    0.027    7.557    0.000
#>     A5|t1            -2.014    0.059  -34.222    0.000
#>     A5|t2            -1.363    0.037  -36.447    0.000
#>     A5|t3            -0.946    0.031  -30.449    0.000
#>     A5|t4            -0.277    0.027  -10.347    0.000
#>     A5|t5             0.652    0.029   22.688    0.000
#>     N1|t1            -0.711    0.029  -24.448    0.000
#>     N1|t2            -0.060    0.026   -2.277    0.023
#>     N1|t3             0.335    0.027   12.441    0.000
#>     N1|t4             0.892    0.031   29.088    0.000
#>     N1|t5             1.489    0.040   36.777    0.000
#>     N2|t1            -1.163    0.034  -34.001    0.000
#>     N2|t2            -0.473    0.028  -17.154    0.000
#>     N2|t3            -0.094    0.026   -3.555    0.000
#>     N2|t4             0.570    0.028   20.349    0.000
#>     N2|t5             1.259    0.036   35.232    0.000
#>     N3|t1            -0.924    0.031  -29.727    0.000
#>     N3|t2            -0.210    0.027   -7.898    0.000
#>     N3|t3             0.106    0.026    4.000    0.000
#>     N3|t4             0.685    0.029   23.772    0.000
#>     N3|t5             1.351    0.037   36.090    0.000
#>     N4|t1            -0.961    0.032  -30.493    0.000
#>     N4|t2            -0.235    0.027   -8.787    0.000
#>     N4|t3             0.155    0.027    5.821    0.000
#>     N4|t4             0.768    0.029   26.055    0.000
#>     N4|t5             1.341    0.037   36.052    0.000
#>     N5|t1            -0.714    0.029  -24.520    0.000
#>     N5|t2            -0.035    0.026   -1.336    0.182
#>     N5|t3             0.297    0.027   11.061    0.000
#>     N5|t4             0.825    0.030   27.482    0.000
#>     N5|t5             1.358    0.038   36.131    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .O1                0.581                           
#>    .O2                0.802                           
#>    .O3                0.353                           
#>    .O4                0.980                           
#>    .O5                0.770                           
#>    .C1                0.684                           
#>    .C2                0.655                           
#>    .C3                0.693                           
#>    .C4                0.425                           
#>    .C5                0.512                           
#>    .E1                0.708                           
#>    .E2                0.484                           
#>    .E3                0.521                           
#>    .E4                0.457                           
#>    .E5                0.619                           
#>    .A1                0.897                           
#>    .A2                0.570                           
#>    .A3                0.428                           
#>    .A4                0.701                           
#>    .A5                0.383                           
#>    .N1                0.289                           
#>    .N2                0.365                           
#>    .N3                0.420                           
#>    .N4                0.558                           
#>    .N5                0.668                           
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
#>    opn=~O1    opn=~O2    opn=~O3    opn=~O4    opn=~O5    con=~C1    con=~C2 
#>  0.6469172 -0.4450754  0.8043163  0.1415558 -0.4797192  0.5623444  0.5869601 
#>    con=~C3    con=~C4    con=~C5 
#>  0.5537472 -0.7585446 -0.6984811 
head(parameterestimates(fit, se = TRUE))
```
|lhs|op|rhs|est|se|z|pvalue|ci.lower|ci.upper|
|---|---|---|--:|--:|--:|--:|--:|--:|
|opn|=~|O1| 0.647|0.027| 24.302|0.000| 0.595| 0.699|
|opn|=~|O2|-0.445|0.033|-13.345|0.000|-0.510|-0.380|
|opn|=~|O3| 0.804|0.029| 27.835|0.000| 0.748| 0.861|
|opn|=~|O4| 0.142|0.045|  3.140|0.002| 0.053| 0.230|
|opn|=~|O5|-0.480|0.031|-15.455|0.000|-0.541|-0.419|
|con|=~|C1| 0.562|0.028| 20.332|0.000| 0.508| 0.617|


``` r
head(fitted(fit)$cov)  # model-implied covariance matrix for the UVs
```
| |O1|O2|O3|O4|O5|C1|C2|C3|C4|C5|E1|E2|E3|E4|E5|A1|A2|A3|A4|A5|N1|N2|N3|N4|N5|
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|O1| 1.000|-0.288| 0.520| 0.092|-0.310| 0.119| 0.125| 0.118|-0.161|-0.148|-0.163|-0.216| 0.208| 0.222| 0.186|-0.056| 0.114| 0.132| 0.095| 0.137|-0.078|-0.074|-0.070|-0.061|-0.053|
|O2|-0.288| 1.000|-0.358|-0.063| 0.214|-0.082|-0.086|-0.081| 0.111| 0.102| 0.112| 0.149|-0.143|-0.153|-0.128| 0.038|-0.079|-0.091|-0.066|-0.094| 0.054| 0.051| 0.048| 0.042| 0.037|
|O3| 0.520|-0.358| 1.000| 0.114|-0.386| 0.149| 0.155| 0.146|-0.200|-0.185|-0.202|-0.269| 0.259| 0.276| 0.231|-0.069| 0.142| 0.164| 0.118| 0.170|-0.097|-0.091|-0.087|-0.076|-0.066|
|O4| 0.092|-0.063| 0.114| 1.000|-0.068| 0.026| 0.027| 0.026|-0.035|-0.032|-0.036|-0.047| 0.046| 0.049| 0.041|-0.012| 0.025| 0.029| 0.021| 0.030|-0.017|-0.016|-0.015|-0.013|-0.012|
|O5|-0.310| 0.214|-0.386|-0.068| 1.000|-0.089|-0.092|-0.087| 0.120| 0.110| 0.121| 0.160|-0.155|-0.165|-0.138| 0.041|-0.085|-0.098|-0.071|-0.101| 0.058| 0.055| 0.052| 0.046| 0.039|
|C1| 0.119|-0.082| 0.149| 0.026|-0.089| 1.000| 0.330| 0.311|-0.427|-0.393|-0.121|-0.162| 0.156| 0.166| 0.139|-0.067| 0.137| 0.158| 0.114| 0.164|-0.144|-0.136|-0.130|-0.114|-0.098|


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

## Checking the stochastic approximation

The estimates come from pairwise stochastic approximation by default. `sa_diagnostics()` reports how the run behaved, how far the estimates are predicted to be from the exact pairwise likelihood estimate, and how much noise the sampling of pairs adds to the standard errors. See also `plot_sa()`, `sa_compare()` and `sa_stability()`. Use `control = set_ucminf_control()` for the exact estimate.

``` r
sa_diagnostics(fit)
#> Stochastic approximation: converged 
#> 
#>  Pairs per iteration: 16 of 300 | epoch: 19 iterations | gain: hessian 
#>           phase epochs iterations seconds
#>            Adam      8        152  0.0145
#>  Newton burn-in     10        190  0.0935
#>        Averaged     50        950  0.0940
#> 
#>  Burn-in plateau: TRUE 
#>  Skipped steps in the averaged phase: 0 
#>  Rejected steps: 101 
#>  Observed information positive definite: TRUE 
#>  Largest predicted distance from exact PML: 0.9 standard errors
#> 
#>  Optimisation noise relative to the standard error:
#>   block parameters median_ratio max_ratio max_distance
#>     tau        125       0.0023     0.011        0.017
#>  lambda         25       0.3402     0.880        0.817
#>    corr         10       0.2311     0.399        0.897
#> 
#>  EPOCHS = 3868 would bring the largest ratio to about 0.1 
```

## Session information

``` r
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Golden Gate 27.0.1
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8
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
#>  [1] gtable_0.3.6        xfun_0.60           ggplot2_4.0.3      
#>  [4] htmlwidgets_1.6.4   psych_2.6.5         lattice_0.23-1     
#>  [7] numDeriv_2016.8-1.1 quadprog_1.5-8      vctrs_0.7.3        
#> [10] tools_4.6.1         generics_0.1.4      stats4_4.6.1       
#> [13] parallel_4.6.1      tibble_3.3.1        cluster_2.1.8.2    
#> [16] pkgconfig_2.0.3     qgraph_1.10.1       Matrix_1.7-6       
#> [19] data.table_1.18.6.1 checkmate_2.3.4     RColorBrewer_1.1-3 
#> [22] S7_0.2.2            RcppParallel_6.2.1  lifecycle_1.0.5    
#> [25] compiler_4.6.1      farver_2.1.2        stringr_1.6.0      
#> [28] mnormt_2.1.2        codetools_0.2-20    litedown_0.11      
#> [31] glasso_1.11         htmltools_0.5.9     fdrtool_1.2.18     
#> [34] htmlTable_2.5.0     Formula_1.2-6       pillar_1.11.1      
#> [37] Hmisc_5.2-6         abind_1.4-8         rpart_4.1.27       
#> [40] nlme_3.1-169        gtools_3.9.5        commonmark_2.0.0   
#> [43] tidyselect_1.2.1    digest_0.6.39       stringi_1.8.9      
#> [46] dplyr_1.2.1         reshape2_1.4.5      fastmap_1.2.0      
#> [49] grid_4.6.1          colorspace_2.1-4    cli_3.6.6          
#> [52] magrittr_2.0.5      base64enc_0.1-6     pbivnorm_0.6.0     
#> [55] foreign_0.8-91      corpcor_1.6.10      scales_1.4.0       
#> [58] backports_1.5.1     rmarkdown_2.31      jpeg_0.1-11        
#> [61] igraph_2.3.3        otel_0.2.0          nnet_7.3-20        
#> [64] gridExtra_2.3.1     RcppEigen_0.3.4.0.2 png_0.1-9          
#> [67] pbapply_1.7-4       evaluate_1.0.5      knitr_1.51         
#> [70] rlang_1.3.0         Rcpp_1.1.2          glue_1.8.1         
#> [73] rstudioapi_0.19.0   R6_2.6.1            plyr_1.8.9         
```

<!-- To fit a latent factor analysis model with 5 factors, where the first 5 items load on the first factor and so on, we can use the `fit_latpl()` function. -->
<!-- In preparation, we must build the following objects to be supplied as a list to the `CONSTR_LIST` argument: -->

<!-- - `CONSTRMAT`: A matrix of constraints on the factor loadings. -->
<!-- - `CORRFLAG`: A flag indicating whether the factor correlations are estimated. -->
<!-- - `CONSTRVAR`: A vector of constraints on the factor variances, if any. -->
<!-- - `STDLV`: A flag indicating whether the factor loadings are standardised. -->
<!-- - `CONSTRAINTS`: Linear equality constraints `list(A, b)` on the parameters, if any. -->

<!-- The above ingredients are prepared using the code below. -->

<!-- ```{r} -->
<!-- p <- 25         # number of items -->
<!-- q <- 5          # number of factors -->
<!-- n <- nrow(bfi)  # number of observations -->

<!-- # Convert data frame to matrix, with each item taking values 0,1,2,...,5 -->
<!-- D <- as.matrix(data.frame(lapply(bfi[, 1:p], as.numeric))) - 1 -->
<!-- head(D) -->

<!-- # Constraint matrix: NA indicate free loadings -->
<!-- constrmat <- build_constr_mat(P = p, Q = q, STRUCT = "simple") -->
<!-- head(constrmat, 12) -->

<!-- # Build constraint list -->
<!-- constr_list <- list( -->
<!--   CONSTRMAT = constrmat, -->
<!--   CORRFLAG = TRUE, -->
<!--   CONSTRVAR = rep(1, q),  # these options standardise the   -->
<!--   STDLV = TRUE,           # factor variances -->
<!--   CONSTRAINTS = NULL -->
<!-- ) -->
<!-- ``` -->

<!-- We can now fit the factor model: -->

<!-- ```{r} -->
<!-- fit <- fit_latpl( -->
<!--   DATA = D,  -->
<!--   CONSTR_LIST = constr_list, -->
<!--   METHOD = "ucminf", -->
<!--   VERBOSE = TRUE -->
<!-- ) -->
<!-- print(fit) -->
<!-- ``` -->

<!-- The default method is using the optimiser from the `{ucminf}` package, which works well for small data sizes. -->
<!-- Alternatively, the user may choose `METHOD = "SA"` to turn on the stochastic approximation algorithm. -->
<!-- Use the function `get_par()` to extract the parameter estimates: -->

<!-- ```{r} -->
<!-- get_par(fit) -->
<!-- ``` -->

<!-- Computing the standard errors is done separately using the `compute_var()` function: -->

<!-- ```{r} -->
<!-- var <- compute_var(OBJ = fit, DATA = D, VERBOSE = TRUE) -->
<!-- str(var) -->

<!-- # Standard errors -->
<!-- head(sqrt(diag(var$vcov) / n)) -->
<!-- ``` -->

<!-- ## `{lavaan}` wrapper -->

<!-- This package also provides a user-friendly `{lavaan}`-style interface to estimate factor analysis models using the `cfa()` function. -->
<!-- Let's rerun the above example, but this time using the stochastic approximation algorithm. -->
