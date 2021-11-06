library(spatstat)
library(glmnet)
source("./R/glmnet_logiengine.R")


radii <- matrix(10,2,2)
diag(radii) <- 10
interaction <- MultiStrauss(radii)
interaction <- Poisson()

## some default options there

subsetexpr=NULL
clipwin=NULL
correction="border"
rbord=reach(interaction)
covfunargs=list()
allcovar=FALSE
vnamebase=c("Interaction", "Interact.")
vnameprefix=NULL
justQ = FALSE
savecomputed = FALSE
precomputed = NULL




n_p <- 300 # presence
n_d <- 50000 # dummy


presence <- data.frame(x = runif(n_p, -100,100),
            y = runif(n_p, -100,100),
            mark = as.factor( rnorm(n_p)>0))

dummy <- data.frame(x = runif(n_d, -100,100),
            y = runif(n_d, -100,100),
            mark = as.factor( rnorm(n_d)>0))


env_p <- data.frame(e1 = rnorm(n_p), e2 = rnorm(n_p))
env_d <- data.frame(e1 = rnorm(n_d), e2 = rnorm(n_d))

pp <- ppp(x = presence[,1],y = presence[,2], marks = presence[,3] ,
             window = owin(xrange = c(-100,100),
                           yrange = c(-100,100)))

dummy_pp <- ppp(x = dummy[,1],y = dummy[,2], marks = dummy[,3] ,
             window = owin(xrange = c(-100,100),
                           yrange = c(-100,100)))

Q <- quadscheme.logi(pp, dummy_pp)

trend <- ~e1:marks+e2:marks+marks-1
penalty.factor <- c(1,1,1,1,1,1)
lambda <- c(10, 1, 0.1)
covariates <- rbind(env_p,env_d)
res_obj <- maxnet.logi.engine(Q, trend, interaction = interaction,
            penalty.factor = penalty.factor, lambda = lambda,
            covariates = covariates)

