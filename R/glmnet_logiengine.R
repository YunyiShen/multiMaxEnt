## use glmnet to fit the Gibbs process, used some internal function in spatstat, developed under spatstat version 2.2.2

maxnet.logi.engine <- function(Q,
                        trend = ~1,
                        interaction,
                        penalty.factor,
                        lambda,
                        ...,
                        covariates=NULL,
                        subsetexpr=NULL,
                        clipwin=NULL,
                        correction="border",
                        rbord=reach(interaction),
                        covfunargs=list(),
                        allcovar=FALSE,
                        vnamebase=c("Interaction", "Interact."),
                        vnameprefix=NULL,
                        justQ = FALSE,
                        savecomputed = FALSE,
                        precomputed = NULL#,
                        #VB=FALSE
                        ){
  if(is.null(trend)) trend <- ~1 
  if(is.null(interaction)) interaction <- spatstat.core:::Poisson()
  want.trend <- !spatstat.utils:::identical.formulae(trend, ~1)
  want.inter <- !spatstat.core:::is.poisson(interaction)
  want.subset <- !is.null(subsetexpr)
  # validate choice of edge correction
  correction <- pickoption("correction", correction,
                           c(border="border",
                             periodic="periodic",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             none="none"))
  # rbord applies only to border correction
  if(correction == "border") {
    spatstat.utils::check.1.real(rbord, "In ppm")
    spatstat.utils:::explain.ifnot(rbord >= 0, "In ppm")
  } else rbord <- 0
  # backdoor stuff
  if(!missing(vnamebase)) {
    if(length(vnamebase) == 1)
      vnamebase <- rep.int(vnamebase, 2)
    if(!is.character(vnamebase) || length(vnamebase) != 2)
      stop("Internal error: illegal format of vnamebase")
  }
  if(!is.null(vnameprefix)) {
    if(!is.character(vnameprefix) || length(vnameprefix) != 1)
      stop("Internal error: illegal format of vnameprefix")
  }
  # create dummy points
  #if(inherits(Q, "ppp")){
  #Xplus <- Q
  #Q <- spatstat.geom::quadscheme.logi(Xplus, ...) # ignored, we will always give quadscheme
  D <- Q$dummy
  Dinfo <- Q$param
  #} 
  #else if(checkfields(Q, c("data", "dummy"))) {
  #  Xplus <- Q$data
  #  D <- Q$dummy
  #  Dinfo <- Q$param
  #  if(is.null(Dinfo)){
  #    Dinfo <- list(how="given", rho=npoints(D)/(area(D)*markspace.integral(D)))
  #  }
  #  Q <- quadscheme.logi(Xplus, D)
  #} else stop("Format of object Q is not understood")
  ## clip to subset?
  if(!is.null(clipwin)) {
    if(is.data.frame(covariates)) {
      ok <- spatstat.geom::inside.owin(spatstat.geom::union.quad(Q), w=clipwin)
      covariates <- covariates[ok, , drop=FALSE]
    }
    Q <- Q[clipwin]
    Xplus <- Q$data
    D     <- Q$dummy
  }
  if (justQ) 
    return(Q)
  ### Dirty way of recording arguments so that the model can be refitted later (should probably be done using call, eval, envir, etc.):
  extraargs <- list(covfunargs = covfunargs, allcovar = allcovar, vnamebase = vnamebase, vnameprefix = vnameprefix)
  extraargs <- append(extraargs, list(...))
  ## Dummy intensity
  if(correction == "border" && Dinfo$how=="grid"){
    Dbord <- D[spatstat.geom::bdist.points(D)>=rbord]
    Dinfo$rho <- spatstat.geom::npoints(Dbord)/(spatstat.geom::eroded.areas(spatstat.geom::as.owin(Dbord), rbord)*spatstat.geom:::markspace.integral(Dbord))
  }
  rho <- Dinfo$rho
  ##Setting the B from Barker dynamics (relative to dummy intensity)
  B <- list(...)$Barker
  if(is.null(B))
    B <- 1
  B <- B*rho
  Dinfo <- append(Dinfo, list(B=B))
  Dinfo <- append(Dinfo, list(extraargs=extraargs))
  # 
  Wplus <- spatstat.geom::as.owin(Xplus)
  nXplus <- spatstat.geom::npoints(Xplus)
  U <- spatstat.geom::superimpose(Xplus, D, W=Wplus, check=FALSE)
#  E <- equalpairs(U, Xplus, marked = is.marked(Xplus))
  E <- cbind(1:nXplus, 1:nXplus)
#  
  computed <- if (savecomputed) list(X = Xplus, Q = Q, U = U) else list()
  # assemble covariate data frame
  if(want.trend || want.subset) {
    tvars <- spatstat.utils::variablesinformula(trend)
    if(want.subset)
      tvars <- union(tvars, all.vars(subsetexpr))
    if(!is.data.frame(covariates)) {
      ## resolve 'external' covariates
      externalvars <- setdiff(tvars, c("x", "y", "marks"))
      tenv <- environment(trend)
      covariates <- spatstat.utils:::getdataobjects(externalvars, tenv, covariates, fatal=TRUE)
    }
    wantxy <- c("x", "y") %in% tvars
    wantxy <- wantxy | rep.int(allcovar, 2)
    cvdf <- data.frame(x=U$x, y=U$y)[, wantxy, drop=FALSE]
    if(!is.null(covariates)) {
      df <- spatstat.core:::mpl.get.covariates(covariates, U, "quadrature points", covfunargs)
      cvdf <- cbind(cvdf, df)
    }
    wantmarks <- "marks" %in% tvars
    if(wantmarks) cvdf <- cbind(cvdf, marks = spatstat.geom::marks(U))
  } else cvdf <- NULL
  # evaluate interaction sufficient statistics
  if (!is.null(ss <- interaction$selfstart)) 
    interaction <- ss(Xplus, interaction)
  V <- spatstat.core:::evalInteraction(Xplus, U, E, interaction, correction, precomputed = precomputed, savecomputed = savecomputed)
  if(!is.matrix(V))
    stop("evalInteraction did not return a matrix")
  if (savecomputed) 
    computed <- append(computed, attr(V, "computed"))
  IsOffset <- attr(V, "IsOffset")
  if(is.null(IsOffset)) IsOffset <- rep.int(FALSE, ncol(V))
  # determine names
  if(ncol(V) > 0) {
    Vnames <- colnames(V)
    if(is.null(Vnames)) {
      nc <- ncol(V)
      Vnames <- if(nc == 1) vnamebase[1L] else paste(vnamebase[2L], 1:nc, sep="")
      colnames(V) <- Vnames
    } else if(!is.null(vnameprefix)) {
      Vnames <- paste(vnameprefix, Vnames, sep="")
      colnames(V) <- Vnames
    }
  } else Vnames <- character(0)
  # combine all data
  glmdata <- as.data.frame(V)
  if(!is.null(cvdf)) glmdata <- cbind(glmdata, cvdf)
  # construct response and weights
  ok <- if(correction == "border") (spatstat.geom::bdist.points(U) >= rbord) else rep.int(TRUE, spatstat.geom::npoints(U))
  # Keep only those quadrature points for which the
  # conditional intensity is nonzero.
  KEEP  <- if(ncol(V)>0) spatstat.utils:::matrowall(V != -Inf) else rep.int(TRUE, spatstat.geom::npoints(U))
  ok <- ok & KEEP
  wei <- c(rep.int(1,spatstat.geom::npoints(Xplus)),rep.int(B/rho,spatstat.geom::npoints(D)))
  resp <- c(rep.int(1,spatstat.geom::npoints(Xplus)),rep.int(0,spatstat.geom::npoints(D)))
  ## User-defined subset:
  if(!is.null(subsetexpr)) {
    USERSUBSET <- eval(subsetexpr, glmdata, environment(trend))
    ok <- ok & USERSUBSET
  }
  # add offset, subset and weights to data frame
  # using reserved names beginning with ".logi."
  glmdata <- cbind(glmdata,
                   .logi.Y = resp,
                   .logi.B = B,
                   .logi.w = wei,
                   .logi.ok =ok)
  # build glm formula 
  # (reserved names begin with ".logi.")
  trendpart <- paste(as.character(trend), collapse=" ")
  fmla <- paste(".logi.Y ", trendpart)
  # Interaction terms
  if(want.inter) {
    VN <- Vnames
    # enclose offset potentials in 'offset(.)'
    if(any(IsOffset))
      VN[IsOffset] <- paste("offset(", VN[IsOffset], ")", sep="")
    offset_vec <- rowSums(as.matrix(glmdata[,VN[IsOffset]])) # get this part of offset
    #glmdata <- glmdata[,!colnames(glmdata)%in% VN[IsOffset]] # remove offset
    #fmla <- paste(c(fmla, VN), collapse="+")
  }
  # add offset intrinsic to logistic technique
  #fmla <- paste(fmla, "offset(-log(.logi.B))", sep="+")
  offset_vec <- offset_vec - log(glmdata[,".logi.B"]) # offset due to the Gibbs part
  fmla <- as.formula(fmla)
  # to satisfy package checker: 
  .logi.B <- B
  .logi.w <- wei
  .logi.ok  <- ok
  .logi.Y   <- resp
  # suppress warnings from code checkers
  spatstat.utils:::dont.complain.about(.logi.B, .logi.w, .logi.ok, .logi.Y)
  Xmat <- model.matrix(fmla, glmdata)
  

  # go
  glmnet::glmnet.control(pmin=1.0e-8, fdev=0) 
  fit <- glmnet::glmnet(x=Xmat[.logi.ok,], y=as.factor(.logi.Y[.logi.ok]), 
    family="binomial", standardize=F, 
    penalty.factor=penalty.factor, 
    lambda=lambda, weights=weights)
  #fit <- if(VB) 
  #         vblogit.fmla(fmla, data = glmdata, 
  #                      subset = .logi.ok, weights = .logi.w, ...)
  #       else 
  #         glm(fmla, data = glmdata, 
  #             family = binomial(), subset = .logi.ok, weights = .logi.w)
  environment(fit$terms) <- sys.frame(sys.nframe())
  ## Fitted coeffs
  co <- coef(fit)
  fitin <- spatstat.core:::fii(interaction, co, Vnames, IsOffset)

  ## Saturated log-likelihood:
  satlogpl <- sum(ok*resp*log(B))
  ## Max. value of log-likelihood:
  maxlogpl <- logLik(fit) + satlogpl

  # Stamp with spatstat version number
  spv <- package_version(spatstat.geom:::versionstring.spatstat())
  the.version <- list(major=spv$major,
                      minor=spv$minor,
                      release=spv$patchlevel,
                      date="$Date: 2020/11/27 03:04:30 $")

  ## Compile results
  fit <- list(method      = "logi",
              fitter      = "glmnet",
              projected   = FALSE,
              coef        = co,
              trend       = trend,
              interaction = interaction,
              fitin       = fitin,
              Q           = Q,
              maxlogpl    = maxlogpl,
              satlogpl    = satlogpl,
              internal    = list(Vnames  = Vnames,
                                 IsOffset=IsOffset,
                                 glmdata = glmdata,
                                 glmfit = fit,
                                 logistic = Dinfo,
                                 computed = computed,
                                 vnamebase=vnamebase,
                                 vnameprefix=vnameprefix,
                                 ),
              covariates  = spatstat.core:::mpl.usable(covariates),
              covfunargs= covfunargs,
              subsetexpr = subsetexpr,
              correction  = correction,
              rbord       = rbord,
              fisher      = NULL,
              varcov      = NULL, # if(VB) fit$S else NULL,
              terms       = terms(trend),
              version     = the.version,
              problems    = list()
              )
  class(fit) <- "ppm"
  return(fit)
}
