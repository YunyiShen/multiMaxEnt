#' @export
predict.multimaxent <-
function(object, locations, spp_list = object$spp_names, newdata, clamp=T, type=c("exponential","cloglog","logistic"), ...)
{
   spp_in <- intersect(object$spp_names, spp_list)
   spp_num <- sapply(spp_in, function(w, thelist){
      which(thelist==w)
   }, object$spp_names)
   if (clamp) {
      for (v in intersect(names(object$varmax), names(newdata))) {
         newdata[,v] <- pmin(pmax(newdata[,v], object$varmin[v]), object$varmax[v])
      }
   }

   mm <- model.matrix(object$maxent_f, data.frame(newdata))
   if (clamp) mm <- t(pmin(pmax(t(mm), object$featuremins), 
                 object$featuremaxs))
   res <- lapply(spp_num, function(i, locations, mm, ppm_fit){
      loci_df <- locations
      loci_df$marks <- i
      predict(ppm_fit, location = loci_df, covariates = mm, type = "intensity")
   }, locations, mm, object$ppm_fit)
   type <- match.arg(type)
   if (type=="exponential") return(res)
   if (type=="cloglog"){
      res <- lapply(1:object$n_spp, function(i, expon, entropy){
         1-exp(0-exp(object$entropy[[i]]+log(res[[i]])))
      }, res, object$entropy)
      return(res)
   } 
   if (type=="logistic"){
      res <- lapply(1:object$n_spp, function(i, expon, entropy){
         1/(1+exp(-object$entropy[[i]]-log(res[[i]])))
      }, res, object$entropy)

      return(res)
   }
}
