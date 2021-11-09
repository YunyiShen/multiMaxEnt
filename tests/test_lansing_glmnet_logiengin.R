library(spatstat)
library(glmnet)
source("./R/glmnet_logiengine.R")



lansing <- unique.ppp( spatstat.data::lansing)
unique_marks <- levels(marks(lansing))

radii1 <- matrix(1e-1,length(unique_marks),length(unique_marks))
radii2 <- matrix(5e-2,length(unique_marks),length(unique_marks))
interaction1 <- MultiStrauss(radii1)
interaction2 <- MultiStrauss(radii2)
interaction <- Hybrid(interaction1, interaction2)

n_p <- 2250
n_d <- 1000

dummy <- data.frame(x = runif(n_d, 0,1), 
            y = runif(n_d, 0,1), 
            mark = factor( sample(unique_marks, n_d, T), 
            levels = unique_marks))

dummy_pp <- ppp(x = dummy[,1],y = dummy[,2], marks = dummy[,3] ,
             window = lansing$window)

Q <- quadscheme.logi(lansing, dummy_pp)

env_p <- data.frame(e1 = rnorm(n_p), e2 = rnorm(n_p))
env_d <- data.frame(e1 = rnorm(n_d), e2 = rnorm(n_d))

trend <- ~marks + e1:marks+e2:marks -1
penalty.factor <- c(rep(0,6) ,rep(1, 12))
lambda <- NULL
covariates <- rbind(env_p,env_d)
res_obj <- glmnet.logi.engine(Q, trend, interaction = interaction, 
            penalty.factor = penalty.factor, lambda = lambda, 
            covariates = covariates)

res_obj_logi <- ppm(Q, trend = trend, 
                    interaction = interaction, 
                    covariates  = covariates)
