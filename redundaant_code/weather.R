# script following the vignette of ShapleyOutlier
# for weather data example

library(ShapleyOutlier)
library(robustHD)
library(tidyverse)
library(tourr)
library(fpc)

data("WeatherVienna")

weather_summer <- WeatherVienna %>% dplyr::select(-c(`t`, t_max, t_min, p_max, p_min)) %>%
  drop_na() %>%
  filter(month %in% c("JUN", "JUL", "AUG")) %>%
  filter(year >= 1955) %>%
  group_by(year) %>%
  dplyr::select(-month) %>%
  summarise(across(.cols = everything(), function(x) mean(x)))

X <- weather_summer %>% dplyr::select(-c(num_frost, num_ice, year))
rownames(X) <- weather_summer$year
#> Warning: Setting row names on a tibble is deprecated.
X <- robStandardize(X)
as.data.frame(X) |>
  rownames_to_column("year") |>
  select(-year) |>
  write_csv("data_weather.csv")

set.seed(1)
MCD <- covMcd(X, alpha = 0.5, nsamp = "best")
#> Warning in .fastmcd(x, h, nsamp, nmini, kmini, trace = as.integer(trace)): 'nsamp = "best"' allows maximally 100000 subsets;
#> computing these subsets of size 17 out of 68
mu <-MCD$center
Sigma <- MCD$cov
Sigma_inv <- solve(MCD$cov)

size_weather <- qchisq(2*pnorm(5)-1, ncol(X))

set.seed(129)
animate_xy(X, 
           guided_anomaly_tour(anomaly_index(),
                               ellipse=as.matrix(Sigma),
                               ellc = size_weather), 
           ellipse=as.matrix(Sigma), ellc = size_weather,
           axes="off")

## now look for angular clusters
dist_pd <- mahalanobis(X, mu, Sigma)
X_outside <- X[dist_pd > size_weather,]
X_outside_scaled <- t(apply(X_outside,
                            1,
                            function(x) (x)/sqrt(sum(x^2))))
#check that it worked
rowSums(X_outside_scaled^2)
for(k in 2:6){
  set.seed(150)
  X_km <- kmeans(X_outside_scaled, k)
  cs <- cluster.stats(dist(X_outside_scaled), X_km$cluster)
  print(cs$dunn2)
}

#preferred solution 4-5 clusters using dunn2/dunn index
#will work with k=4

set.seed(150)
X_km <- kmeans(X_outside_scaled, 4)
X_km$cluster

# cluster 1: high temperature
# cluster 2: 
# cluster 3: temperature + sun_h, num_clear
# cluster 4: avg_t_min small for the outlying points
X_outside_2 <- X_outside[X_km$cluster == 4,]
set.seed(129)
animate_xy(X_outside_2, 
           guided_anomaly_tour(anomaly_index(),
                               ellipse=as.matrix(Sigma),
                               ellmu = mu,
                               ellc = size_weather), 
           ellipse=as.matrix(Sigma), ellc = size_weather,
           ellmu = mu, center = FALSE,
           axes="off")


# steps:
# 1 - find outliers, split them into angular clusters
# 1a - use some cluster metric to decide on k
# 2 - with all observations run anomaly index, show the
# clusters in different colors
# 2a - define an axis display that only labels "long" axes
# 3 - for some subset(s) of interest run anomaly index
# only on the subset

