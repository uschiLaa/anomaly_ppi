# Jenny Wadsworth's river data
library(tidyverse)
river <- read_csv("data/riverFlow1.csv") |>
  rename(V1 = `1`, V2 = `2`, V3 = `3`, V4 = `4`) |>
  mutate(type = "river")

surface <- read_csv("data/SurfacePoints.csv") |>
  select(-`...1`) |>
  mutate(type = "surface")

river_all <- bind_rows(river, surface) |>
  mutate(type = factor(type))

animate_xy(river_all[,1:4], col=river_all$type, scale=TRUE)
# Data needs to be scaled to 0-1 to match the surface range

river_m <- apply(river[,1:4], 2, mean)
river_s <- apply(river[,1:4], 2, sd)

river_std <- river |>
  mutate_at(vars(V1:V4), function(x) (x-min(x))/(max(x)-min(x)))

river_all <- bind_rows(river_std, surface) |>
  mutate(type = factor(type, levels = c("surface", "river")))

animate_xy(river_all[,1:4], 
           col=river_all$type, 
           cex = c(rep(1, nrow(river_std)), 
                   rep(0.2, nrow(surface))),
           scale=TRUE, center=FALSE, 
           axes="bottomleft")

ptsize <- c(rep(1, nrow(river_std)), 
            rep(0.2, nrow(surface)))
render_gif(river_all[,1:4], 
           grand_tour(),
           display_xy(col=river_all$type, 
                      cex = ptsize, 
                      scale=TRUE, center=FALSE, 
                      axes="bottomleft",
                      half_range=1.5),
           gif_file = "river.gif",
           frames = 500,
           width = 500,
           height = 500
)

