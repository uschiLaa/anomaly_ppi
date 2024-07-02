# script following the vignette of ShapleyOutlier
# to extract the data examples into csv files

library(ShapleyOutlier)
library(robustHD)
library(tidyverse)
library(tourr)

data(TopGear)

rownames(TopGear) = paste(TopGear[,1],TopGear[,2])
myTopGear <- TopGear[,-31] #removing the verdict variable
myTopGear <- myTopGear[,sapply(myTopGear,function(x)any(is.numeric(x)))]
myTopGear <- myTopGear[!apply(myTopGear,1, function(x)any(is.na(x))),]
myTopGear <- myTopGear[,-2]
# Transform some variables to get roughly gaussianity in the center:
transTG = myTopGear
transTG$Price = log(myTopGear$Price)
transTG$Displacement = log(myTopGear$Displacement)
transTG$BHP = log(myTopGear$BHP)
transTG$Torque = log(myTopGear$Torque)
transTG$TopSpeed = log(myTopGear$TopSpeed)

transTG <- transTG %>% rename("log(Price)" = Price,
                              "log(Displacement)" = Displacement,
                              "log(BHP)" = BHP,
                              "log(Torque)" = Torque,
                              "log(TopSpeed)" = TopSpeed)

X <- as.matrix(transTG)
X <- robStandardize(X)

data_topgear <- as.data.frame(X) |>
  rownames_to_column("ID") |>
  select(-ID)

# for comparison we also export robust estimates
# for mean and covariance
set.seed(1)
MCD <- covMcd(X, nsamp = "best")
#> Warning in .fastmcd(x, h, nsamp, nmini, kmini, trace = as.integer(trace)): 'nsamp = "best"' allows maximally 100000 subsets;
#> computing these subsets of size 12 out of 245
mu <-MCD$center
Sigma_topgear <- MCD$cov
Sigma_inv_topgear <- solve(MCD$cov)

mu_topgear <- colMeans(data_topgear)

animate_xy(data_topgear, axes = "bottomleft", 
           ellipse=Sigma_topgear, #ellc = 100,
           # the default ellipse seems much too small
           # what would be the better way to center the ellipse??
           # data has median zero, but mean value is quite different
           #ellmu = mu_topgear
           )

# size at 5 sigma
size_topgear <- qchisq(2*pnorm(5)-1, ncol(data_topgear))

animate_xy(data_topgear, 
           guided_anomaly_tour(anomaly_index(),
                               ellipse=as.matrix(Sigma_topgear),
                               ellc = size_topgear), 
           ellipse=as.matrix(Sigma_topgear), ellc = size_topgear,
           axes="bottomleft")

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

animate_xy(X_outside_2, 
           grand_tour(),
           ellipse=as.matrix(Sigma), ellc = size_weather,
           ellmu = mu, center = FALSE,
           axes="off", half_range = 7)

#### try again but with a bit smaller p

X <- weather_summer %>% dplyr::select(-c(num_frost, num_ice, year,
                                        wind_v, num_precp_01))
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
## to get angular clusters normalize the lengths
X_outside_scaled <- t(apply(X_outside,
                          1,
                          function(x) (x)/sum(x^2)))
#check that it worked
rowSums(X_outside_scaled)
set.seed(150)
X_km <- kmeans(X_outside_scaled, 3)
