library(lavaan.pl)
library(tidyverse)
data("bfi")
mod <- "agr =~ A1 + A2 + A3 + A4 + A5"

# Convert to binary
bfiord <-
  bfi |>
  mutate(
    across(A1:O5, function(x) ordered(as.numeric(x >= 4))),
    gender = ifelse(gender == 1, "Male", "Female")
  )

# Unconstrained model (aka configural invariance)
fit1 <- lavaan::cfa(mod, bfiord, std.lv = TRUE, group = "gender")
parameterestimates(fit1)
inspect(fit1, "free")  # Now a list of 2, one for each group

# (Weak) Invariance model. Note: Only group 1 LV is set to 1
fit2 <- lavaan::cfa(mod, bfiord, std.lv = TRUE, group = "gender",
                    group.equal = c("loadings"))
inspect(fit2, "free")
lavTables(fit2)  # freq counts also accounts for multi groups

# (Strong) Invariance model. Note: because of the restriction on the thresholds,
# the total variance of the y* cannot be 1 in both groups now, so the var(y*)
# needs to be estimated in the second group.
fit3 <- lavaan::cfa(mod, bfiord, std.lv = TRUE, group = "gender",
                    group.equal = c("loadings", "thresholds"))
inspect(fit3, "free")

# Possible to set a specific constraint, e.g. same loading for one item in each
# group. Here we set loading for A5 to be the same in both groups.
fit4 <- lavaan::cfa("agr =~ A1 + A2 + A3 + A4 + c(a,a)*A5", bfiord, std.lv = TRUE,
                    group = "gender")
partable(fit4)[partable(fit4)$label == "a", ]

## ---- Constraints in lavaan --------------------------------------------------

# Setting equality constraints is like doing a linear transformation of the full
# or unpacked version of the coefficients. The optimiser will work with the
# *packed* version, exploring a smaller parameter space. The matrix
# transformation is provided by lavaan.

x <- coef(fit2)  # "unpacked" coefficients. note the ordering
if (fit2@Model@eq.constraints) {
  # To make it smaller (input for fit function and optimiser)
  x_pack <- as.numeric(
    (x - fit2@Model@eq.constraints.k0) %*% fit2@Model@eq.constraints.K
  )
  # To make it bigger (converting the result of the optimiser)
  x_unpack <- as.numeric(
    fit2@Model@eq.constraints.K %*% x_pack + fit2@Model@eq.constraints.k0
  )
}

length(x)
length(x_pack)
all(x == x_unpack)
