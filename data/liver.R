# data for Annalisa on liver function
library(mvtnorm)
library(tidyverse)
library(tourr)
#library(GGally)

liver_stats <- tribble(
  ~ALT, ~AST, ~ALP, ~albumin, ~bilirubin, ~GGT, ~protein, ~sex, ~stat,
  5,     5,   30,      35,        0,    2, 60,   "m", "mn",
  40,     30,  120,      50,        20,   44, 80,   "m", "mx",
  5,     5,   30,      35,        0,    2, 60,   "f", "mn",
  40,     30,  120,      50,        20,   70, 80,    "f", "mx"
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
# Protein
# https://www.racgp.org.au/getattachment/36d1c5e0-9c1d-43fc-a8a0-b323e3ed8fbe/Liver-function-tests.aspx
# https://www.mayoclinic.org/tests-procedures/liver-function-tests/about/pac-20394595
# https://www.thepathologycentre.org/test/liver-panel-lft-albumin-bilirubin-alt-alp-ggt-and-protein-2/

# Aging women have elevated ALP
# ALT decreases with age, for all demog
# https://link.springer.com/referenceworkentry/10.1007/978-3-319-90761-1_49-1

# Using standardised measurements
set.seed(1257)
norm_vc <- diag(1, 4, 4)
norm_vc[upper.tri(norm_vc) == TRUE] <- sample(c(0.3, 0.4, 0.5), 6, replace=T)
for (i in 1:3)
  for (j in (i+1):4)
    norm_vc[j, i] <- norm_vc[i, j]
norm_vc[1,2] <- 0.8
norm_vc[2,1] <- 0.8

liver_norm <- rmvnorm(193, mean=rep(0, 4), sigma=as.matrix(norm_vc))
animate_xy(liver_norm, axes = "off", 
           ellipse=as.matrix(norm_vc))
animate_xy(liver_norm, axes = "off", 
           ellipse=as.matrix(norm_vc), ellc=qchisq(0.99, 4))
animate_xy(liver_norm, 
           guided_anomaly_tour(anomaly_index(),
             ellipse=as.matrix(norm_vc), ellc=qchisq(0.99, 4)), 
           ellipse=as.matrix(norm_vc), ellc=qchisq(0.99, 4),
           axes="off")
colnames(liver_norm) <- colnames(liver_stats)[c(1,2,4,5)]
write_csv(as.data.frame(liver_norm), file="data/liver_norm.csv")
colnames(norm_vc) <- colnames(liver_stats)[c(1,2,4,5)]
write_csv(as.data.frame(norm_vc), file="data/liver_norm_vc.csv")

vc_f <- norm_vc
vc_f[1,2] <- 0.1
vc_f[2,1] <- 0.1
vc_f <- vc_f/4

liver_f_means <- c(-0.3, 0.5, -0.1, 0.1)
set.seed(523)
liver_f <- rmvnorm(54, mean=liver_f_means, sigma=as.matrix(vc_f))
colnames(liver_f) <- colnames(liver_stats)[c(1,2,4,5)]

x <- rbind(liver_f, liver_norm)
x <- data.frame(x)
x$type <- factor(c(rep("f", 54), rep("norm", 193)))
ggscatmat(x, columns = 1:4, color="type", alpha=0.5)

norm_mu <- liver_f_means*(-1)
norm_mu <- matrix(norm_mu, ncol=4)

animate_xy(liver_f, axes = "off", ellipse=as.matrix(norm_vc), 
           ellc = qchisq(0.99, 4), ellmu = t(norm_mu), half_range=6)
animate_xy(liver_f, 
           guided_anomaly_tour(anomaly_index(),
                               ellipse=as.matrix(norm_vc), 
                               ellc = qchisq(0.99, 4),
                               ellmu = t(norm_mu)), 
           start = basis_random(4, 2),
           ellipse=as.matrix(norm_vc), 
           ellc = qchisq(0.99, 4),
           ellmu = t(norm_mu), 
           axes="bottomleft", half_range=6)
write_csv(as.data.frame(liver_f), file="data/liver_f.csv")

colnames(vc_f) <- colnames(liver_stats)[c(1,2,4,5)]
write_csv(as.data.frame(vc_f), file="data/liver_f_vc.csv")

colnames(norm_mu) <- colnames(liver_stats)[c(1,2,4,5)]
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

# Try 2 with actual values
liver_stats |>
  filter(sex == "m")

liver_stats |>
  pivot_longer(ALT:PT, names_to = "test", values_to = "value") |>
  group_by(sex, test) |>
  pivot_wider(names_from = "stat", values_from = "value") |>
  mutate(rng = mx - mn,
         mid = (mn + mx)/2)

set.seed(1158)
norm_vc <- diag(1, 4, 4)
norm_vc[upper.tri(norm_vc) == TRUE] <- sample(c(0.3, 0.4, 0.5), 6, replace=T)
for (i in 1:3)
  for (j in (i+1):4)
    norm_vc[j, i] <- norm_vc[i, j]
colnames(norm_vc) <- c("GGT", "AST", "ALP", "ALT")
norm_samp <- rmvnorm(500, sigma=as.matrix(norm_vc))
colnames(norm_samp) <- c("GGT", "AST", "ALP", "ALT")

set.seed(217)
women_vc <- norm_vc
women_vc[1,4] <- 0.5
women_vc[4,1] <- 0.5
women <- rmvnorm(16, mean=c(0.3, -0.3, -1.5, -1.7),
                 sigma = women_vc)
colnames(women) <- c("GGT", "AST", "ALP", "ALT")
women <- as.data.frame(women)

  #tibble(GGT = rnorm(16, mean=0.3, sd=1),
  #              ALP = rnorm(16, mean=-0.3, sd=0.8),
  #              AST = rnorm(16, mean=-1, sd=1.5), 
  #              ALT = rnorm(16, mean=-1.5, sd=1.1))

animate_xy(women,            
           ellipse=norm_vc, 
           axes = "bottomleft", 
           half_range=5, 
           center=FALSE)

animate_xy(women, guided_anomaly_tour(anomaly_index(),
                                    ellipse=norm_vc), 
           ellipse=norm_vc, 
           axes = "bottomleft", 
           half_range=5, 
           center=FALSE)

set.seed(1208)
pt1 <- tibble(GGT = rnorm(5, mean=1.3, sd=0.15),
              AST = rnorm(5, mean=1, sd=0.2),
              ALP = c(1, 1.3, 1.6, 2.8, 2.9), 
              ALT = c(-1, -1.3, -1.6, -2.9, -3.1))

animate_xy(pt1,            
           ellipse=norm_vc, 
           axes = "bottomleft", 
           half_range=5, 
           center=FALSE,
           edges = matrix(c(1:4, 2:5), byrow=FALSE, ncol=2))

animate_xy(pt1, guided_anomaly_tour(anomaly_index(),
                                     ellipse=norm_vc), 
           ellipse=norm_vc, 
           axes = "bottomleft", 
           half_range=5, 
           center=FALSE,
           edges = matrix(c(1:4, 2:5), byrow=FALSE, ncol=2),
           obs_labels = c("45", "50", "55", "65", "70"))


