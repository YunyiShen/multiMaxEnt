#' Fitting a multispecies MaxEnt model via Gibbs point process
#' The main algorithm to fit a multivarate extension of Phillip's MaxEnt model based on Gibbs point process using spatstat's ppm routine
#' @param p a data.frame of presence record, should have three columns, namely x,y and species (the exact same order).
#' @param env_layers a Raster.stack object from package raster, the environmental layers
#' @param background_mask a spatial polygon from package sp, the place to take background points, should be smaller than raster
#' @param f a function that make formula's for the fitting, default is the maxent's method, or a formula
#' @param feature_classes feature classes used if f is maxent's feature generator
#' @param interaction the interaction model for the species, default the multitype Strauss process with radii 1000 units (usually means 1km), see ?ppm for this, usual selection is HirStrauss or MultiStrauss depends on whether the interaction is symmetric
#' @param regmult a constant to adjust regularization.
#' @param regfun a function to compute regularization constant for each feature.
#' @param addsamplestobackground bool whether to add the presence point to the background
#' @param gitter noise added to the sample when adding to background, will be uniform [-gitter, gitter].
#' @param n_background number of background points
#' @param type_background method to sample the background passed to sp::spsample, argument 'type', recommand to use random
#' @param verbos whether to be verbos
#' @param ... other arguments to be used in ppm

# we need to use spatstat.geom::quadscheme to specify quadrature scheme with our own dummy points, which will be used later to calculate the entropy
multimaxent <-
function(p, # presence points, the third column should be species (or what ever taxa)
         env_layers, # environmental layers
         background_mask = NULL, # a spatial polygon
         f = maxent.formula,
         feature_classes = "default",
         interaction = NULL,
         regmult = 1.0,
         regfun = maxnet.default.regularization,
         addsamplestobackground=FALSE,
         gitter = 0.5,
         n_background = 10000,
         type_background = "random",
         verbos = TRUE,...)
{
   if(!all(colnames(p[,1:3])==c("x","y","species"))){
      stop("p should have the first three columns being x,y and species\n")
   }
   #browser()
   if(!"factor" %in% class(p[,3])){
     p[,3] <- as.factor(p[,3])
   }
   spp_names <- as.character(levels(p[,3])) # save spp name
   n_spp <- length(spp_names)
   #p[,3] <- as.numeric(p[,3])
   colnames(p) <- c("x","y","marks")

   if(is.null(interaction)){
      radii <- matrix(1000,n_spp, n_spp)
      diag(radii) <- 0
      interaction <- MultiStrauss(radii)
   }
   if(is.null(background_mask)){
      if(verbos) cat("no mask provided, will use environmental layers to generate background points...\n")
      background_mask <- rasterToPolygons(env_layers[[1]])
   }

   if(verbos) cat("generating background points...")
   dummy_points <- spsample(background_mask, n = n_background, type = type_background)
   rm(background_mask)
   dummy_points <- as.data.frame(dummy_points)
   dummy_points <- dummy_points[!duplicated(dummy_points),]
   colnames(dummy_points) <- c("x","y")
   if (addsamplestobackground) {
       dummy_points <- rbind(dummy_points, p[,1:2]+runif(2*nrow(p),-gitter,gitter))
   }

   env_presence <- extract(env_layers, p[,1:2])
   env_presence <- as.data.frame(env_presence)
   env_dummy <- extract(env_layers, dummy_points)
   dummy_points <- dummy_points[!is.na(rowSums(env_dummy)),]
   env_dummy <- na.omit(env_dummy)
   env_dummy <- as.data.frame(env_dummy)
   if (anyNA(env_presence)|anyNA(env_dummy)) stop("NA values in data table. This is either because presence points outside the rasters or rasters are not properly alligned.")

   if(verbos) cat("done\nmodel fitting...")
   if(class(f)=="function"){
     maxent_f <- f(p,env_presence,feature_classes)
   }
   else if(class(f)!="formula"){
     stop("f must be either a formula or a function that generated formula from presence point and environmental variables")
   }


   pp <- ppp(x = p[,1],y = p[,2], marks = p[,3] ,
             window = owin(xrange = c(env_stack@extent[1],
                                      env_stack@extent[2]),
                           yrange = c(env_stack@extent[3],
                                      env_stack@extent[4])))# point pattern for presence
   dummy_marks <- factor( sample(spp_names, nrow(dummy_points), TRUE ), levels = spp_names)
   dummy_pp <- ppp(x = dummy_points[,1],
                   y = dummy_points[,2],
                   marks = dummy_marks,
                   window = owin(xrange = c(env_stack@extent[1],
                                                  env_stack@extent[2]),
                                       yrange = c(env_stack@extent[3],
                                                  env_stack@extent[4])))# point pattern for presence
   Q <- quadscheme.logi(pp, dummy_pp)

   #browser()
   mm <- as.data.frame( model.matrix(maxent_f, rbind(env_presence)))
   #reg <- regfun(mm) * regmult
   reg <- lapply(1:n_spp, function(w,mm,regmult){regfun(mm) * regmult},mm,regmult )
   reg <- Reduce(c, reg)
   reg <- c(reg, rep(0,n_spp))
   mm <- rbind(mm,as.data.frame( model.matrix(maxent_f, rbind(env_dummy))))
   feature_name <- colnames(mm)
   colnames(mm) <- paste0("x",1:ncol(mm))
   lambda=10^(seq(4,0,length.out=200))*sum(reg)/length(reg)*nrow(env_presence)/nrow(mm)
   
   # TODO: when too many features exits, we will have infinite coefficients in usual logistics, should move to elstic nets, go logistic.R and change it, together with predict.ppm
   forla <- paste0("marks:x",1:ncol(mm))
   forla <- paste(forla, collapse = "+")
   forla <- paste0(c("~",forla,"+marks-1"),collapse = "")
   forla <- as.formula(forla)
   ppm_fit <- maxnet.logi.engine(Q = Q,trend = forla,
                  interaction = interaction,
                  covariates = mm, penalty.factor = reg, lambda = lambda,...)

   model <- list()
   model$ppm_fit <- ppm_fit


   if(verbos) cat("done\ncalculating entropy...")
   preds <- predict.ppm(object = ppm_fit, locations = dummy_pp, covariates = mm[1:nrow(env_dummy)+nrow(env_presence),], type = "cif")
   raw <- lapply(spp_names,function(spp, preds, dummy_marks){
     preds[dummy_marks==spp]
   }, preds, dummy_marks)
   model$alpha <- lapply(raw,function(w){-log(sum(w))})
   raw <- lapply(raw, function(w){w/(sum(w))})
   model$entropy <- lapply(raw, function(rawi){-sum(rawi * log(rawi))})

   if(verbos) cat("finished\n")
   full_env <- rbind(env_presence, env_dummy)
   vv <- (sapply(full_env, class)!="factor")
   model$varmin <- apply(full_env[,vv, drop = FALSE], 2, min)
   model$varmax <- apply(full_env[,vv, drop = FALSE], 2, max)
   model$featuremins <- apply(mm, 2, min)
   model$featuremaxs <- apply(mm, 2, max)
   means <- apply(env_presence[,vv, drop = FALSE], 2, mean)
   majorities <- sapply(names(full_env)[!vv],
      function(n) which.max(table(full_env[p==1,n, drop = FALSE])))
   names(majorities) <- names(full_env)[!vv]
   model$samplemeans <- unlist(c(means, majorities))
   model$levels <- lapply(full_env, levels)
   model$spp_names <- spp_names
   model$n_spp <- n_spp
   model$maxent_f <- maxent_f
   model$feature_name <- feature_name
   cat("done\n")
   model
}
