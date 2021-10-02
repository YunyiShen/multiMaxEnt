#' @import stats, spatstat.core, raster, sp
#' @export

# we need to use spatstat.geom::quadscheme to specify quadrature scheme with our own dummy points, which will be used later to calculate the entropy
multimaxent <-
function(p, # presence points, the third column should be species (or what ever taxa)
         env_layers, # environmental layers 
         background_mask = NULL, # a spatial polygon
         f = maxent.formula, 
         interaction = MultiHard(), 
         addsamplestobackground=TRUE, 
         n_background = 10000, 
         type_background = "random", 
         verbos = TRUE,...)
{
   spp_names <- as.character(levels(p[,3])) # save spp name
   n_spp <- length(spp_names)
   p[,3] <- as.numeric(p[,3])
   colnames(p) <- c("x","y","marks")
   if(is.null(background_mask)){
      if(verbos) cat("no mask provided, will use environmental layers...")
      background_mask = rasterToPolygons(env_layers)
   }

   if(verbos) cat("generating background points...")
   dummy_points <- spsample(background_mask, n = n_background, type = type_background)
   dummy_points <- as.data.frame(dummy_points)
   colnames(dummy_points) <- c("x","y")
   if (addsamplestobackground) {
       dummy_points <- rbind(dummy_points, p[,1:2])
   }

   env_presence <- extract(env_layers, p[,1:2])
   env_dummy <- extract(env_layers, dummy_points)
   if (anyNA(env_presence)|anyNA(env_dummy)) stop("NA values in data table. This is either because presence points outside the rasters or rasters are not properly alligned.")
      
   
   if(verbos) cat("model fitting...")
   maxent_f <- f(p,env_presence) 
   mm <- model.matrix(maxent_f, rbind(env_presence, env_dummy))
   
   pp <- ppp(x = p[,1],y = p[,2], marks = p[,3] )# point pattern for presence
   Q <- quadscheme(pp, dummy_points)

   ppm_fit <- ppm(Q,trend=~(.-1):marks-1, 
                  interaction = interaction, 
                  covariates = mm, method="logi")
   
   model <- list()
   model$ppm_fit <- ppm_fit

   pred_lists <- lapply(1:n_spp, function(i, dummy, feature_dummy, ppm_fit){
      loci_df <- dummy
      loci_df$marks <- i
      predict(ppm_fit, location = loci_df, covariates = feature_dummy, type = "intensity")
   }, dummy_points, model.matrix(maxent_f,env_dummy), ppm_fit)


   
   raw <- lapply(pred_lists,function(w){w/sum(w)})
   model$entropy <- lapply(raw, function(rawi){-sum(rawi * log(rawi))})
   model$alpha <- lapply(pred_lists,function(w){-log(sum(w))})

   model$featuremins <- apply(mm, 2, min)
   model$featuremaxs <- apply(mm, 2, max)
   full_env <- rbind(env_presence, env_dummy)
   vv <- (sapply(full_env, class)!="factor")
   model$varmin <- apply(full_env[,vv, drop = FALSE], 2, min)
   model$varmax <- apply(full_env[,vv, drop = FALSE], 2, max)
   means <- apply(env_presence[,vv, drop = FALSE], 2, mean)
   majorities <- sapply(names(full_env)[!vv], 
      function(n) which.max(table(full_env[p==1,n, drop = FALSE])))
   names(majorities) <- names(full_env)[!vv]
   model$samplemeans <- unlist(c(means, majorities))
   model$levels <- lapply(full_env, levels)
   model$spp_names <- spp_names
   model$n_spp <- n_spp
   model$maxent_f <- maxent_f
   model
}
