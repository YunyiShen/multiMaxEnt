library(spatstat)
library(glmnet)
library(Matrix)
source("./R/glmnet_logiengine.R")

set.seed(42)
n_spp <- 5
spp_names <- paste0("spp_",1:n_spp)
radii <- matrix(5e-3,n_spp,n_spp)
log_gammas <- rsparsematrix(n_spp,n_spp, 
            density = 0.2,
            rand.x = function(n){runif(n,-3,-1)})
log_gammas <- (log_gammas + t(log_gammas))/2
rownames(log_gammas) <- colnames(log_gammas) <- spp_names

beta <- exp(rnorm(n_spp, 6,0.5))

pp <- rmh(model = rmhmodel(cif = "straussm", 
                            par = list(beta = beta, 
                            gamma = exp(as.matrix(log_gammas)),
                            radii = radii),
                            trend = NULL,
                            w = owin(xrange = c(0,1),
                                yrange = c(0,1)),
                            types = spp_names))
pp
n_d <- 2000
dummy <- data.frame(x = runif(n_d, 0,1), 
            y = runif(n_d, 0,1), 
            mark = factor( sample(spp_names, n_d, T), 
            levels = spp_names))

dummy_pp <- ppp(x = dummy[,1],y = dummy[,2], marks = dummy[,3] ,
             window = owin(xrange = c(0,1),
                           yrange = c(0,1)))

Q <- quadscheme.logi(pp, dummy_pp)


trend <- ~marks  -1
interaction <- MultiStrauss(1*radii)

penalty.factor <- c(rep(0,n_spp))
lambda <- NULL
res_obj <- glmnet.logi.engine(Q, trend, interaction = interaction, 
            penalty.factor = penalty.factor, lambda = lambda, 
            covariates = NULL, cv = T)

res_obj_logi <- ppm(Q, trend = trend, 
                    interaction = interaction, 
                    covariates  = NULL)
