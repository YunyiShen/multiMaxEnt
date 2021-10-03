#' @export
#' Predict the multimaxent model
#' @param obj the fitted multimaxent
#' @param location x,y location to be predicted
#' @param spp_list the species to be predicted
#' @param newdata environmental variables on the location
#' @param clamp whether to feature and variable clampping
#' @param type which type to predict? can be exponential for the intensity of the underlying Gibbs process, otherwise cloglog or logistic for the probability output
predict.multimaxent <-
function(object, locations, spp_list = object$spp_names, newdata, clamp=T, type=c("exponential","cloglog","logistic"))
{
   spp_in <- intersect(object$spp_names, spp_list)

   if (clamp) {
     for (v in intersect(names(object$varmax), names(newdata))) {
       newdata[,v] <- pmin(pmax(newdata[,v], object$varmin[v]), object$varmax[v])
     }
   }
   terms <- sub("hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)", object$feature_name)
   terms <- sub("categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\"\\2\")", terms)
   terms <- sub("thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)", terms)
   f <- formula(paste("~", paste(terms, collapse=" + "), "-1"))
   mm <- model.matrix(f, data.frame(newdata))
   if (clamp) mm <- t(pmin(pmax(t(mm), object$featuremins[object$feature_name]),
                           object$featuremaxs[object$feature_name]))
   mm <- as.data.frame(mm)
   colnames(mm) <- paste0("x",1:ncol(mm))

   res <- lapply(spp_in, function(spp, spp_full,locations, newdata, ppm_fit){
      loci_df <- locations
      loci_df$marks <- factor(spp, levels = spp_list)
      predict(ppm_fit, location = loci_df, covariates = mm, type = "cif")
   }, object$spp_names,locations, mm, object$ppm_fit)
   type <- match.arg(type)
   if (type=="exponential") return(res)
   if (type=="cloglog"){
      res <- lapply(spp_in, function(i, expon, entropy){
         1-exp(0-exp(object$entropy[[i]]+log(res[[i]])))
      }, res, object$entropy)
      return(res)
   }
   if (type=="logistic"){
      res <- lapply(spp_in, function(i, expon, entropy){
         1/(1+exp(-object$entropy[[i]]-log(res[[i]])))
      }, res, object$entropy)

      return(res)
   }
}
