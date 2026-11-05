library(lavaan.pl)
data("bfi")
mod <- "
  opn =~ O1 + O2 + O3 + O4 + O5  # Openness
  con =~ C1 + C2 + C3 + C4 + C5  # Conscientiousness
  ext =~ E1 + E2 + E3 + E4 + E5  # Extraversion
  agr =~ A1 + A2 + A3 + A4 + A5  # Agreeableness
  neu =~ N1 + N2 + N3 + N4 + N5  # Neuroticism

  # opn ~~ 0*con
  # opn ~~ 0*ext
  # opn ~~ 0*agr
  # opn ~~ 0*neu
  # con ~~ 0*ext
  # con ~~ 0*agr
  # con ~~ 0*neu
  # ext ~~ 0*agr
  # ext ~~ 0*neu
  # agr ~~ 0*neu
"
fit <- pmlcfa(model = mod, data = bfi, std.lv = TRUE, test = "none", verbose = TRUE,
           estimator.args = list(computevar_numderiv = TRUE))
summary(fit)

# tictoc::tic()
# Hinv <- lavaan:::lav_model_information_observed(
#   lavmodel =fit@Model,
#   lavsamplestats = fit@SampleStats,
#   lavdata = fit@Data,
#   lavoptions = fit@Options,
#   lavcache = fit@Cache,
#   inverted = TRUE
# )
# tictoc::toc()
# tictoc::tic()
# J <- lavaan:::lav_model_information_firstorder(
#   lavmodel =fit@Model,
#   lavsamplestats = fit@SampleStats,
#   lavdata = fit@Data,
#   lavoptions = fit@Options,
#   lavcache = fit@Cache
# )
# tictoc::toc()
