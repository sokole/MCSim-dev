source("RCODE_1_load_scripts.R")

# --------------------------------------------------------
# --------------------------------------------------------
# --------------------------------------------------------
# --------------------------------------------------------
# -- site some sim characteristics
# --------------------------------------------------------
scenario.ID<-'DISPERSAL_EXAMPLES'

n.sites<-16
m<-.01
SWM.slope<-70
n.generations<-100
nu<-1e-3

# --------------------------------------------------------
# -- make landscape using igraph
# --------------------------------------------------------
require(igraph)

# ----------
# -- graph a tree
edge.size<-round(sqrt(n.sites),0)
i.graph<-make_lattice(c(edge.size,edge.size), directed=TRUE)
dist.mat<-distances(i.graph, mode='in')
# -- watershed size
n.parents<-apply(
  dist.mat,
  1,
  function(x){sum(x[x!=Inf]>0)})
# hist(degree(g.tree))
graphics.off()
par(mfrow=c(1,2), mar=c(4,2,2,0))
plot(i.graph)
hist(log(n.parents+1))


# ----------
# -- landscape metadata
d.site.info<-data.frame(site.ID=c(1:n.sites),
                        size=runif(n.sites,1,25),
                        watershed.size=n.parents,
                        site.degree=degree(i.graph, mode='all'),
                        Ef=runif(n.sites,0,1),
                        m=m)


# ----------
# -- make landscape object
landscape<-fn.make.landscape(site.info=d.site.info,
                             dist.mat=as.matrix(dist.mat),
                             list.of.stuff=list(i.graph=i.graph))


# --------------------------------------------------------
# -- Run a simulation

sim.result<-fn.metaSIM(
  landscape=landscape,
  scenario.ID = scenario.ID, 
  alpha.fisher = .7,
  nu = nu,
  speciation.limit = NA,
  JM.limit = 1e5, # may need to put a limit on how many individuals can be in a simulation
  n.timestep = n.generations,
  SWM.slope = SWM.slope,
  trait.dispersal.median = 1 ,
  trait.dispersal.range = 0,
  trait.Ef.sd = 0.3,
  save.sim = TRUE, 
  output.dir.path = paste('SIM',scenario.ID,sep='_'))




