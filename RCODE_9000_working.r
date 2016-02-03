# devtools::load_all("MCSim")
# devtools::document("MCSim")
source("RCODE_1_load_scripts.R")

# -------------------------------------
# testing fn.lambda()
# -- a function that calculates a species recruitment prob given
# -- local environmental conditions (Ef)
# -------------------
# -- WHAT NEW? -- set default parameter values

trait.optimum.vect<-c(.2,.3,.01)
env.filter.vect<-c(.5,.5,.9)
fn.lambda(trait.optimum = trait.optimum.vect,
          Ef = env.filter.vect)

# -------------------------------------
# testing fn.lottery.recruit ()
# -- a function that does lottery recruitment, given a species list
# -- and a vector of lambdas
# -------------------

# WHAT NEW? -- set default par values
fn.lottery.recruit(vect.taxa.list = letters[1:10])

fn.lottery.recruit(
  vect.recruitment.weights = c(.3, .4, .5, .9),
  vect.taxa.list = letters[1:4],
  scalar.JL = 10L) 

# -------------------------------------
# testing fn.set.regional.species.pool()
# -- a function that sets the initial regional species pool
# -------------------
set.seed(10)
fn.set.regional.species.pool(alpha.fisher = 1)
