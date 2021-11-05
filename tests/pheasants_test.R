library(raster)
load("./data/Chrysolophus.rda")
load("./data/env_stack.rda")
all_source <- list.files("./R","R",full.names = T)
lapply(all_source, source)
source("./R/multimaxent.R")
library(spatstat)
background_mask <- rasterToPolygons(env_stack[[1]])
Chrysolophus$species <- as.factor(Chrysolophus$species)

res <- multimaxent(Chrysolophus, env_stack,
                   feature_classes = "h",
                   verbos = T, interaction = MultiStrauss(matrix(50000,2,2)),
                   background_mask = background_mask,
                   n_background = 1000, addsamplestobackground = F)

