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
norm_vc <- diag(1, 5, 5)
norm_vc[upper.tri(norm_vc) == TRUE] <- sample(c(0.3, 0.4, 0.5), 10, replace=T)
for (i in 1:4)
  for (j in (i+1):5)
    norm_vc[j, i] <- norm_vc[i, j]
norm_vc[1,2] <- 0.8
norm_vc[2,1] <- 0.8

liver_norm <- rmvnorm(193, mean=rep(0, 5), sigma=as.matrix(norm_vc))
animate_xy(liver_norm, axes = "off", 
           ellipse=as.matrix(norm_vc))
animate_xy(liver_norm, 
           guided_anomaly_tour(anomaly_index(),
             ellipse=as.matrix(norm_vc)), 
           ellipse=as.matrix(norm_vc),
           axes="off")
colnames(liver_norm) <- colnames(liver_stats)[1:5]
write_csv(as.data.frame(liver_norm), file="data/liver_norm.csv")
colnames(norm_vc) <- colnames(liver_stats)[1:5]
write_csv(as.data.frame(norm_vc), file="data/liver_norm_vc.csv")

vc_f <- norm_vc
vc_f[1,2] <- 0.5
vc_f[2,1] <- 0.5

liver_f_means <- c(-0.3, 0.5, 0, 0, 0)
liver_f <- rmvnorm(267, mean=liver_f_means, sigma=as.matrix(vc_f))

colnames(liver_f) <- colnames(liver_stats)[1:5]
write_csv(as.data.frame(liver_f), file="data/liver_f.csv")

colnames(vc_f) <- colnames(liver_stats)[1:5]
write_csv(as.data.frame(vc_f), file="data/liver_f_vc.csv")

norm_mu <- liver_f_means*(-1)
norm_mu <- matrix(norm_mu, ncol=5)
colnames(norm_mu) <- colnames(liver_stats)[1:5]
write_csv(as.data.frame(norm_mu), file="data/liver_norm_means.csv")

# Try it out
norm_vc <- read_csv("data/liver_norm_vc.csv")
liver_f <- read_csv("data/liver_f.csv")
norm_mu <- read_csv("data/liver_norm_means.csv")
animate_xy(liver_f, axes = "off", ellipse=as.matrix(norm_vc), 
           ellc = 23, ellmu = t(norm_mu))
animate_xy(liver_f, 
           guided_anomaly_tour(anomaly_index(),
             ellipse=as.matrix(norm_vc), 
             ellc = 23,
             ellmu = t(norm_mu)), 
           start = basis_random(5, 2),
           ellipse=as.matrix(norm_vc), 
           ellc = 23,
           ellmu = t(norm_mu), 
           axes="bottomleft")

d <- bind_rows(
  bind_cols(liver_norm, set="norm"),
  bind_cols(liver_f, set="females"))
d <- d |>
  mutate(set = factor(set))

animate_xy(d[,1:5], axes = "off", col=d$set)
animate_xy(d[,1:5], guided_tour(lda_pp(d$set)), col=d$set)


