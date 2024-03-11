# data for Annalisa on liver function
library(mvtnorm)
library(tidyverse)
library(tourr)

liver_stats <- tribble(
  ~ALT, ~AST, ~ALP, ~Albumin, ~Bilirubin, ~GGT, ~LD, ~PT, ~sex, ~stat,
  29,     14,   44,      3.5,        0.1,    5, 140,   0,   "m", "mn",
  33,     20,  147,      5.5,        1.2,   40, 280, 1.1,   "m", "mx",
  19,     10,   44,      3.5,        0.1,    5, 140,   0,   "f", "mn",
  25,     36,  147,      5.5,        1.2,   40, 280, 1.1,   "f", "mx"
)

liver_stats |>
  pivot_longer(ALT:PT, names_to = "test", values_to = "value") |>
  group_by(sex, test) |>
  pivot_wider(names_from = "stat", values_from = "value") |>
  mutate(dif = mx - mn)

# ALT, ALP IU/L
# AST U/L
# Albumin g/dL
# Bilirubin mg/dL
# GGT, LD U/L

# Using standardised measurements
set.seed(1257)
vc_norm <- diag(1, 5, 5)
vc_norm[upper.tri(vc_norm) == TRUE] <- sample(c(0.3, 0.4, 0.5), 10, replace=T)
for (i in 1:4)
  for (j in (i+1):5)
    vc_norm[j, i] <- vc_norm[i, j]
vc_norm[1,2] <- 0.6
vc_norm[2,1] <- 0.6

liver_norm <- rmvnorm(193, mean=rep(0, 5), sigma=vc_norm)
animate_xy(liver_norm, axes = "off", 
           ellipse=vc_norm, 
           ellsize=15.5, half_range = 7)
animate_xy(liver_norm, 
           guided_anomaly_tour(anomaly_index(),
             ellipse=vc_norm, ellsize=15.5), 
           ellipse=vc_norm, ellsize=15.5, half_range = 7,
           axes="off")
colnames(liver_norm) <- colnames(liver_stats)[1:5]
write_csv(as.data.frame(liver_norm), file="data/liver_norm.csv")
colnames(vc_norm) <- colnames(liver_stats)[1:5]
write_csv(as.data.frame(vc_norm), file="data/liver_norm_vc.csv")

vc_f <- vc_norm
vc_f[1,2] <- 0.9
vc_f[2,1] <- 0.9

liver_f_means <- c(-0.3, 1.5, 0, 0, 0)
liver_f <- rmvnorm(267, mean=liver_f_means, sigma=vc_f)

colnames(liver_f) <- colnames(liver_stats)[1:5]
write_csv(as.data.frame(liver_f), file="data/liver_f.csv")

colnames(vc_f) <- colnames(liver_stats)[1:5]
write_csv(as.data.frame(vc_f), file="data/liver_f_vc.csv")


# Try it out
vc_norm <- read_csv("data/liver_norm_vc.csv")
liver_f <- read_csv("data/liver_f.csv")
animate_xy(liver_f, axes = "off", ellipse=vc_norm, ellsize=3)
animate_xy(liver_f, 
           guided_anomaly_tour(anomaly_index(),
             ellipse=vc_norm, ellsize=20), 
           ellipse=vc_norm, ellsize=20, 
           half_range=7, axes="off")

d <- bind_rows(
  bind_cols(liver_norm, set="norm"),
  bind_cols(liver_f, set="females"))
d <- d |>
  mutate(set = factor(set))

animate_xy(d[,1:5], axes = "off", col=d$set)
animate_xy(d[,1:5], guided_tour(lda_pp(d$set)), col=d$set)


