fn.make.landscape_0.3.1<-function(
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
  area_m2 = 1,
  Ef.specificity = 0, # 0 is point specificity along env. gradient
  Ef = .5,
  guess.site.coords = FALSE,
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
    area_m2<-data.frame(
      dummy=c(1:n.obs),
      area_m2=area_m2)$area_m2
    
    # -------------------------------
    # -- calculate assemblage sizes at sites, JL 
    # -- JL influenced by management
    # -------------------------------
    JL.wts <- area_m2 / sum(area_m2)
    JL.wts<-JL.wts/sum(JL.wts)
    JL <- round(JL.wts * JM,0)
    
    # -------------------------------
    # -- calculate immigration at sites, IL 
    # -------------------------------
    I.site <- I.rate.m2 * area_m2
    I.site <- round(I.site,0)
    m.site <- I.site/(I.site + JL - 1)
    
    # -------------------------------
    # -- dat with info
    # -------------------------------
    dat.info.default<-data.frame(
      site.ID = c(1:nrow(dat.geo.out)),
      area_m2 = area_m2,
      JL = JL,
      Ef = Ef,
      Ef.specificity = Ef.specificity,
      IL = I.site,
      m = m.site
    )
    
    if(nrow(site.info)>0){
      dat.info.out<-site.info  
      # -- check to see if specific vars need to be filled in
      for(i.var in c('site.ID','area_m2','JL','Ef','Ef.specificity','IL','m')){
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
