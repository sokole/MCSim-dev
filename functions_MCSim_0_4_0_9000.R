# ---------------------------------------------------------------------------------------
# -- v0.4.0.9000 dev functions
# ---------------------------------------------------------------------------------------
#' fn.make.landscape
#' 
#' @title make a simulation landscape

fn.make.landscape<-function(
  # -------------------------------
  # -------------------------------
  # -------------------------------
  # -- data frame inputs
  # -- need one or the other, or they have to match. Priority given to dist if they don't
  site.coords = data.frame(),
  dist.mat = data.frame(),
  
  # -- will fill in info if none given
  site.info = data.frame(),
  
  # -- default metacommunity parameters if none given
  JM = 1000,
  I.rate.m2 = 1,
  area.m2 = 1,
  Ef.specificity = 0, # 0 is point specificity along env. gradient
  Ef = .5,
  guess.site.coords = TRUE,
  list.of.stuff = NA
){
  # -------------------------------
  # -------------------------------
  # -------------------------------
  
  # -- ESCAPE VAR
  get.the.f.out<-FALSE
  
  
  # -------------------------------
  # -- Make geo landscape data.frame
  # -------------------------------
  dist.mat<-as.data.frame(as.matrix(dist.mat)) #make into a data.frame
  
  # --  check info
  if(nrow(dist.mat)>0){
    dat.geo.dist.out<-dist.mat
    if(nrow(site.coords)!=nrow(dist.mat) & guess.site.coords){
      dat.geo.out<-data.frame(cmdscale(dat.geo.dist.out))
      print('I assigned dat.geo for you')
    }else if(nrow(site.coords)!=nrow(dist.mat)){
      dat.geo.out<-data.frame(site.label=c(1:nrow(dat.geo.dist.out)))
      print('Rock on')
    }else if(nrow(site.coords)==nrow(dist.mat)){
      dat.geo.out<-site.coords
    }
  }else if(nrow(site.coords)>0){
    dat.geo.out<-site.coords
    dat.geo.dist.out<-data.frame(as.matrix(dist(site.coords)))
    print('gangsta')
  }else if(nrow(site.info)>0){
    get.the.f.out<-TRUE
    print('no geo info!')
  }else{
    get.the.f.out<-TRUE
    print('no landscape info!')
  }
  
  n.obs<-nrow(dat.geo.dist.out)
  
  if(!get.the.f.out){
    # -- turn area into a vector if it is not a vector, otherwise, it remains the same, get's fed up if it's the wrong length
    area.m2<-data.frame(
      dummy=c(1:n.obs),
      area.m2=area.m2)$area.m2
    
    # -------------------------------
    # -- calculate assemblage sizes at sites, JL 
    # -- JL influenced by management
    # -------------------------------
    JL.wts <- area.m2 / sum(area.m2)
    JL.wts <- JL.wts/sum(JL.wts)
    JL <- round(JL.wts * JM,0)
    
    # -------------------------------
    # -- calculate immigration at sites, IL 
    # -------------------------------
    I.site <- I.rate.m2 * area.m2
    I.site <- round(I.site,0)
    m.site <- I.site/(I.site + JL - 1)
    
    # -------------------------------
    # -- dat with info
    # -------------------------------
    dat.info.default<-data.frame(
      site.ID = c(1:nrow(dat.geo.out)),
      area.m2 = area.m2,
      JL = JL,
      Ef = Ef,
      Ef.specificity = Ef.specificity,
      IL = I.site,
      m = m.site
    )
    
    if(nrow(site.info)>0){
      dat.info.out<-site.info  
      # -- check to see if specific vars need to be filled in
      for(i.var in c('site.ID','area.m2','JL','Ef','Ef.specificity','IL','m')){
        if(!i.var%in%names(dat.info.out)) dat.info.out[,i.var] <- dat.info.default[,i.var]
      }
    }else{
      dat.info.out<-dat.info.default
    }
    
    return(
      list(site.info=dat.info.out,
           site.coords=dat.geo.out,
           dist.mat=dat.geo.dist.out,
           list.of.stuff=list.of.stuff)
    )
  }else{
    print('no landscape for you!')
  }
}

# ---------------------------------------------------------------------------------------
#' fn.recruit.Jt
#' 
#' @title Recruitment in a metacommunity with context
#' 
#' @description Internal function called by \link{fn.metaSIM} to calculate local 
#' recruitment pools for all sites in a metacommunity, and call \link{fn.lottery.recruit}
#'to create a new assemblages (J) at time (t) after accounting for metacommunity 
#'composition at time (t-1), dispersal dynamics, landscape topology, and environmental 
#'filtering.
#' 
#' @usage 
#' fn.recruit.Jt(landscape.site.coords = landscape$dat[, c("x", "y")], nu = nu,
#'               SWM.slope = SWM.slope, J.t.minus.1 = J.t.minus.1, 
#'               taxa.list = taxa.list, traits.Ef = dat.gamma.t0$trait.Ef, 
#'               trait.Ef.sd = trait.Ef.sd, traits.dispersal = dat.gamma.t0$trait.dispersal, 
#'               landscape.m = landscape$dat$m, landscape.Ef = landscape$dat$Ef, 
#'               landscape.Ef.specificity = landscape$dat$Ef.specificity, 
#'               landscape.JL = landscape$dat$JL)
#'               
#' @param landscape.site.coords A data.frame with xy-coordinates for each site.  Each row is a 
#' site in the metacommunity landscape.
#' @param nu Numeric, Hubbell's "speciation rate", but can be interpreted as the 
#' probability of the appearance of a novel species.  If set to 0, no novel taxa 
#' will appear during the simulation.
#' @param SWM.slope Slope of dispersal kernel, see \link{fn.metaSIM}.  This value 
#' is used for \code{w} in eqn. 6 from Gravel et al. (2006) to model dispersal 
#' limitation.
#' @param J.t.minus.1 Species composition of sites in a metacommunity at previous 
#' time step.
#' @param taxa.list A vector of character strings of species' names.
#' @param traits.Ef A vector of species' trait scores, numeric.
#' @param trait.Ef.sd A scalar value representing species' niche widths (numeric, 
#' currently one value is used for all species).
#' @param traits.dispersal A vector of species' dispersal trait scores.
#' @param landscape.m A vector of values for \code{m} for each site in the landscape.
#' @param landscape.Ef A vector of environmental filter \code{Ef} values for each site 
#' in the landscape.
#' @param landscape.Ef.specificity Environmental specificity at each site (numeric). 
#' @param landscape.JL A vector of assemblage sizes (positive integers) for all sites 
#' in the landscape.
#' 
#' @seealso \link{fn.lottery.recruit}, \link{fn.make.landscape}, \link{fn.metaSIM}
#' 
#' @export
#' 
fn.recruit.Jt <- function(
  # calls calculates R.t (expected pool) to calculate extant pool
  par.list<-list(mat.geodist = as.matrix(dist(1:3)),
                 nu = .001,
                 SWM.slope = 0,
                 J.t.minus.1 = matrix(JL/length(taxa.list), 
                                      nrow=nrow(as.matrix(mat.geodist)), 
                                      ncol=length(taxa.list)),
                 taxa.list = letters[1:5],
                 traits.Ef = array(.5, length(taxa.list)),
                 trait.Ef.sd = 0,
                 traits.dispersal = array(1, length(taxa.list)),
                 m = 1,
                 Ef = array(.5, nrow(mat.geodist)),
                 Ef.specificity = 0,
                 JL = 10){
  
  # -- Calculate R.t, estimates of locally available pools at each site based on 
  # local relative abundances at time t-1 (or t0)
  # regional relative abundances at time t-1 (or t0), which is calculated from neighbors weighted by the SWM
  # the balance of local and regional pools is determined by m (where m = 0 is only local)
  
  # deal with sites with 0s
  J.t.minus.1.RAs<-as.matrix(J.t.minus.1/rowSums(J.t.minus.1))
  J.t.minus.1.RAs[!is.finite(J.t.minus.1.RAs)]<-0
  J.t.minus.1.RAs<-as.data.frame(J.t.minus.1.RAs)
  
  # ----------------------------
  # -- calculate recruitment probabilities from I for each site
  # create SWM to weight all neighboring sites contributions to "Immigrant" pool
  # use eq 6 from Gravel et al. 2006, but set diags to 0
  #   mat.geodist<-as.matrix(dist(landscape.site.coords))
  
  max.val<-max(mat.geodist[mat.geodist!=Inf]) #max finite value
  mat.geodist.scaled<-mat.geodist/max.val
  SWM<-exp(-1*SWM.slope*mat.geodist.scaled^2)
  SWM<-SWM/rowSums(SWM)
  diag(SWM)<-0
  
  
  # -- calculate I for each site based on composition of neighboring sites
  I.RAs<-SWM%*%as.matrix(J.t.minus.1.RAs) #distance decay without species bias
  
  # -- Alter RAs for I based on dispersal traits -- create bias toward species with better dispersal
  I.RAs.dispersal.biased<-I.RAs*(t(traits.dispersal%o%rep(1,nrow(I.RAs))))
  I.RAs.dispersal.biased<-I.RAs.dispersal.biased/rowSums(I.RAs.dispersal.biased)
  I.RAs.dispersal.biased[is.nan(I.RAs.dispersal.biased)]<-0
  # -- calculate recruitment pool by combining local RAs and regional RAs, weighted by m
  R.t<-m*I.RAs.dispersal.biased+(1-m)*J.t.minus.1.RAs
  
  # -- all species that have 0 regional RA at time t-1 are assigned a non-zero probability of recruitment,
  # which is nu / (count of unobserved species in the list).  A "speciation event" in this simulation is 
  # recruitment of a previously unobserved species. 
  unobserved.spp.count <- sum(colSums(R.t)==0)
  speciation.recruitment.prob <- ifelse(
    unobserved.spp.count == 0,
    0,
    nu / unobserved.spp.count)
  
  # -- alter all the probabilities accordingly, all rows must add to 1
  mat.nu<-R.t
  mat.nu[,colSums(mat.nu)>0]<-mat.nu[,colSums(mat.nu)>0]*(1-nu)
  mat.nu[,colSums(mat.nu)==0]<-speciation.recruitment.prob
  R.t.probs<-mat.nu
  
  # -- Local environmental filtering
  d.temp<-expand.grid(Ef=Ef,stringsAsFactors=FALSE,
                      trait.Ef=traits.Ef)
  d.temp$Ef.specificity<-Ef.specificity
  
  lambda.Ef.siteBYspp<-matrix(
    data=mapply(FUN=fn.lambda,
                trait.optimum=d.temp$trait.Ef,
                Ef=d.temp$Ef,
                Ef.specificity=d.temp$Ef.specificity,
                MoreArgs=list(niche.breadth=trait.Ef.sd)
    ),
    nrow=length(Ef), #number sites
    ncol=length(traits.Ef), #number spp
    byrow=FALSE
  )
  lambda.Ef.siteBYspp<-lambda.Ef.siteBYspp/rowSums(lambda.Ef.siteBYspp) #rescale so row sums are 1
  
  # reweight R.t.probs based on local env. filters
  R.t.envfiltered<-lambda.Ef.siteBYspp*R.t.probs
  R.t.envfiltered<-R.t.envfiltered/rowSums(R.t.envfiltered)
  
  # make lists of recruitment weights for each site for lottery call below
  R.t.list<-as.list(data.frame(t(R.t.envfiltered)))
  
  # -- Recruit extant community for time t from R.t
  J.t1<-data.frame(
    row.names=NULL,
    t(mapply(
      FUN=fn.lottery.recruit,
      vect.recruitment.weights=R.t.list,
      scalar.JL=as.list(JL),
      MoreArgs=list(vect.taxa.list=taxa.list)
    )))
  return(J.t1)
}
