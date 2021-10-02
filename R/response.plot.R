#' @export
response.plot <-
function(mod, spp_list, v, loci,type, mm=mod$samplemeans, min=mod$varmin[v], max=mod$varmax[v], levels=unlist(mod$levels[v]), plot=T, xlab=v, ylab=tools::toTitleCase(type), ...) {
   spp <- c(spp)
   nr <- if (is.null(levels)) 100 else length(levels)
   m <- data.frame(matrix(mm,nr,length(mm),byrow=T))
   colnames(m) <- names(mm)
   m[,v] <- if (!is.null(levels)) levels else 
      seq(min - 0.1*(max-min), max+0.1*(max-min), length=100)
   locations <- data.frame(x = rep(loci[1], nr), y = loci[2])
   preds <- predict.multimaxent(mod, locations, spp_list, m, type=type)
   if (plot) {
      par(mfrow = c(1,length(spp_list)))
      for(i in 1:length(spp_list)){
         if (is.null(levels)) {
         plot(m[,v], preds[[i]], xlab=xlab, ylab=ylab, type="l", main = spp_list[i],...) 
      } else {
         graphics::barplot(as.vector(preds[[i]]), names.arg=levels, xlab=xlab, ylab=ylab, main = spp_list[i], ...)
      }
      }
      
   }
   else{
      preds <- lapply(preds, 
      function(w, mv,v){setNames(data.frame(mv, w), c(v, 'pred'))},
      m[, v],v)
      names(preds) <- spp_list
      return(preds)
   }
}
