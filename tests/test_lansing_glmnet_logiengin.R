library(spatstat)
library(glmnet)
source("./R/glmnet_logiengine.R")



lansing <- unique.ppp( spatstat.data::lansing)
unique_marks <- levels(marks(lansing))


radii3 <- matrix(5e-3,length(unique_marks),length(unique_marks))
interaction3 <- MultiStrauss(radii3)
interaction <- Hybrid(interaction3)

n_p <- 2250
n_d <- 3000

set.seed(42)
dummy <- data.frame(x = runif(n_d, 0,1), 
            y = runif(n_d, 0,1), 
            mark = factor( sample(unique_marks, n_d, T), 
            levels = unique_marks))

dummy_pp <- ppp(x = dummy[,1],y = dummy[,2], marks = dummy[,3] ,
             window = lansing$window)

Q <- quadscheme.logi(lansing, dummy_pp)

trend <- ~marks - 1
penalty.factor <- c(rep(0,6) )
lambda <- NULL

res_obj <- glmnet.logi.engine(Q, trend, interaction = interaction, 
            penalty.factor = penalty.factor, lambda = lambda, 
            covariates = NULL, cv = F)

res_obj_logi <- ppm(Q, trend = trend, 
                    interaction = interaction, 
                    covariates  = NULL)
