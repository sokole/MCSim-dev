#' fn.metaSIM
#' 
#' @title A metacommunity simulation for ecologists
#' 
#' @description This function is a lottery-based, zero-sum, spatially explicit 
#' simulation that can include neutral and/or niche-based dynamics.  Results are 
#' written to a .csv file in a SIM_OUTPUT directory.  Output includes parameter 
#' settings, diversity partitioning outcomes, and variation partitioning outcomes.
#' 
#' @usage 
#' fn.metaSIM3(landscape)
#' 
#' @param landscape landscape list from \link{fn.make.landscape3}
#' @param scenario.ID see \link{fn.metaSIM}
#' @param alpha.fisher see \link{fn.metaSIM}
#' @param nu see \link{fn.metaSIM}
#' @param speciation.limit see \link{fn.metaSIM}
#' @param n.timestep see \link{fn.metaSIM}
#' @param SWM.slope see \link{fn.metaSIM}
#' @param trait.dispersal.median see \link{fn.metaSIM}
#' @param trait.dispersal.range see \link{fn.metaSIM}
#' @param trait.Ef.sd see \link{fn.metaSIM}
#' @param save.sim see \link{fn.metaSIM}
#' @param output.dir.path see \link{fn.metaSIM}
#' @param \dots other parameters to be passed to internal functions
#' 
#' @seealso \link{fn.metaSIM}
#' 
#' @references 
#' Chao, A., C. H. Chiu, and T. C. Hsieh. 2012. Proposing a resolution to debates on diversity 
#' partitioning. Ecology 93:2037--2051.
#' 
#' Etienne, R. S. 2005. A new sampling formula for neutral biodiversity. Ecology Letters 8:253--260.
#' 
#' Gravel, D., C. D. Canham, M. Beaudet, and C. Messier. 2006. Reconciling niche and neutrality: 
#' the continuum hypothesis. Ecology Letters 9:399--409.
#' 
#' Hubbell, S. P. 2001. A unified theory of biodiversity and biogeography. Princeton University 
#' Press.
#' 
#' Jost, L. 2007. Partitioning diversity into independent alpha and beta components. 
#' Ecology 88:2427--2439.
#' 
#' @export
#' 
fn.metaSIM_0.3.1<-function(
 landscape=landscape,
  scenario.ID = NA, 
  alpha.fisher = 2,
  nu = 1e-04,
  speciation.limit = NA,
  JM.limit = 1e5, # may need to put a limit on how many individuals can be in a simulation
  n.timestep = 10,
  SWM.slope = 0,
  trait.dispersal.median = 1 ,
  trait.dispersal.range = 0,
  trait.Ef.sd = 0.3,
  save.sim = TRUE, 
  output.dir.path = "SIM_OUTPUT", 
  ...
){
  try({  
    # --------------------------------------------------------
    # -- calculate JM from landscape
    # --------------------------------------------------------
    attach(landscape$site.info)
    JM <- sum(JL)
    n.sites <- length(JL)
    
    #' -- create regional pool
    dat.gamma.t0 <- fn.set.regional.species.pool(n.timestep = n.timestep, 
                                                 nu = nu, 
                                                 speciation.limit = speciation.limit, 
                                                 JM = ifelse(JM>JM.limit,JM.limit,JM), 
                                                 alpha.fisher = alpha.fisher, 
                                                 trait.dispersal.median = trait.dispersal.median, 
                                                 trait.dispersal.range = trait.dispersal.range)
    taxa.list <- as.character(dat.gamma.t0$taxa.list)
    d.temp <- expand.grid(Ef = Ef, stringsAsFactors = FALSE, 
                          trait.Ef = dat.gamma.t0$trait.Ef)
    d.temp$Ef.specificity <- Ef.specificity
    lambda.Ef.siteBYspp <- matrix(data = mapply(FUN = fn.lambda, 
                                                trait.optimum = d.temp$trait.Ef, Ef = d.temp$Ef, Ef.specificity = d.temp$Ef.specificity, 
                                                MoreArgs = list(niche.breadth = trait.Ef.sd)), nrow = n.sites, 
                                  ncol = length(taxa.list), byrow = FALSE)
    lambda.Ef.siteBYspp <- lambda.Ef.siteBYspp/rowSums(lambda.Ef.siteBYspp)
    
    R.probs.t0 <- lambda.Ef.siteBYspp * (rep(1, n.sites) %o% 
                                           dat.gamma.t0$regional.RA)
    R.probs.t0 <- R.probs.t0/rowSums(R.probs.t0)
    R.probs.list <- as.list(data.frame(t(R.probs.t0)))
    
    J.t0 <- data.frame(row.names = NULL, t(mapply(FUN = fn.lottery.recruit, 
                                                  vect.recruitment.weights = R.probs.list, scalar.JL = as.list(JL), 
                                                  MoreArgs = list(vect.taxa.list = taxa.list))))
    J <- list()
    J[[1]] <- J.t0
    J.t.minus.1 <- J.t0
    
    # ----------------------------------------------------------------------
    # loop for generational turnover in metacommunity
    # ----------------------------------------------------------------------
    for (t.index in 2:n.timestep) {
      J.t <- fn.recruit.Jt_0.3.1(mat.geodist=landscape$dist.mat,
                                 nu=nu,
                                 SWM.slope=SWM.slope,
                                 J.t.minus.1=J.t.minus.1,
                                 taxa.list=taxa.list,
                                 traits.Ef=dat.gamma.t0$trait.Ef,
                                 trait.Ef.sd=trait.Ef.sd,
                                 traits.dispersal=dat.gamma.t0$trait.dispersal,
                                 m=m,
                                 Ef=Ef,
                                 Ef.specificity=Ef.specificity,
                                 JL=JL)
      J[[t.index]] <- J.t
      J.t.minus.1 <- J.t
      print(paste("Timestep:", t.index))
    }
    
    # ----------------------------------------------------------------------
    sim.result.name <- paste("SIM_", scenario.ID, "_", format(Sys.time(), 
                                                              "%Y%m%d_%H%M%S"), "_", trunc(runif(1, 1e+05, 999999)), 
                             sep = "")
    
    sim.result <- list(scenario.ID = scenario.ID, 
                       sim.result.name = sim.result.name, 
                       alpha.fisher = alpha.fisher, 
                       nu.sim = nu, 
                       trait.Ef.sd = trait.Ef.sd, 
                       trait.dispersal = dat.gamma.t0$trait.dispersal, 
                       trait.Ef = dat.gamma.t0$trait.Ef, 
                       landscape = landscape, 
                       dat.gamma.t0 = dat.gamma.t0, 
                       SWM.slope = SWM.slope, 
                       J = J, 
                       n.timestep = n.timestep, 
                       taxa.list = taxa.list)
    
    #'   fn.sim.metadata.archive4(sim.result = sim.result, 
    #'                            save.sim = save.sim, var.dir = output.dir.path, 
    #'                            keep.timesteps = keep.timesteps,
    #'                            q.order=NA,
    #'                            ...)
    sim.result.filename<-paste(output.dir.path,"/",sim.result.name,".rda",sep="")
    sim.result.metadata <- data.frame(row.names = sim.result.name, 
                                      scenario.ID = scenario.ID, 
                                      sim.ID = sim.result.name, 
                                      sim.result.filename = sim.result.filename, 
                                      n.sites = n.sites, 
                                      n.timestep = n.timestep, 
                                      alpha.fisher = alpha.fisher, 
                                      nu.sim = nu, 
                                      JM = JM, 
                                      JL.mean = mean(JL), 
                                      JL.sd = sd(JL), 
                                      m.mean = mean(m), 
                                      m.sd = sd(m), 
                                      IL.mean = mean(IL), 
                                      IL.sd = sd(IL), 
                                      SWM.slope = SWM.slope, 
                                      Ef.mean = mean(Ef), 
                                      Ef.sd = sd(Ef), 
                                      Ef.specificity.mean = mean(Ef.specificity), 
                                      Ef.specificity.sd = sd(Ef.specificity), 
                                      Tr.disp.mean = mean(dat.gamma.t0$trait.dispersal), 
                                      Tr.disp.sd = sd(dat.gamma.t0$trait.dispersal), 
                                      Niche.breadth = sim.result$trait.Ef.sd, 
                                      stringsAsFactors = FALSE)
    
    #' -- check for director
    if(!output.dir.path%in%list.files())  dir.create(output.dir.path)
    
    #' -- save sim
    if(save.sim) save(sim.result,file=sim.result.filename)
    
    #' -- check to see if data for other reps from this scenario exist
    filname.sim.metadata <- paste(output.dir.path, "/sim.metadata_", 
                                  sim.result$scenario.ID, ".csv", sep = "")
    
    #' -- create new file if none exists, write results to file, otherwise append to existing file
    if (file.exists(filname.sim.metadata)) {
      dat.sim.metadata <- read.csv(filname.sim.metadata, row.names = 1, 
                                   header = TRUE)
      if (!sim.result.name %in% row.names(dat.sim.metadata)) {
        dat.sim.metadata <- rbind(dat.sim.metadata, sim.result.metadata)
        write.csv(dat.sim.metadata, filname.sim.metadata)
      }
    }else {
      write.csv(sim.result.metadata, filname.sim.metadata)
    }
    
    try(detach(landscape$site.info), silent = TRUE)
    try(detach(landscape), silent = TRUE)
    
    #' -- return results
    return(sim.result)
  })
}
