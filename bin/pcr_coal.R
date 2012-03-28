#!/usr/bin/env Rscript

#
# Script simulating PCR amplifications and dilutions. Makes use of the pcrcoal package.
#

library(ape)
library(pcrcoal)

args <- commandArgs(trailingOnly = TRUE)

# Process command line arguments:
name            <- args[1]
init.popsize    <- as.numeric(args[2])
pcr.eff         <- as.numeric(args[3])
nr.cycles.mut   <- as.numeric(args[4])
dilf.after.mut  <- as.numeric(args[5])
nr.cycles.cln   <- as.numeric(args[6])
dilf.after.cln  <- as.numeric(args[7])
nr.cycles.cov   <- as.numeric(args[8])
sample.size.mut <- as.numeric(args[9])
mut.only        <- as.numeric(args[10])

#print(dilf.after.mut)
#print(dilf.after.cln)

nwk.out         <- paste(name,".nwk",sep="")
cov.out         <- paste(name,".cov",sep="")

if(length(args) == 0) {
    stop("Not enough arguments!")
}

# Mutagenic PCR:
eff.mut         <- rep(pcr.eff, nr.cycles.mut)

sample.size.cln <- 1   # Sample size does not matter in this case.

# Cleanup PCR:
eff.cln         <- rep(pcr.eff, nr.cycles.cln)


##
## Begin simulation:
##

# Construct pcrcoal object for the mutagenic PCR:
pcoal.mut<-PCRcoal(
    initial.size    = init.popsize,     # Number of template molecules.
    sample.size     = sample.size.mut,  # Number of molecules sampled.
    nr.cycles       = nr.cycles.mut,    # Number of PCR cycles.
    efficiencies    = eff.mut           # A vector with the per-cycle PCR efficiencies.
)


# Simulate mutagenic PCR:
mut.res    <-sample.tnt(pcoal.mut)
tree       <- mut.res$phylo
sm         <- c()

if(mut.only == 0){

# Dilute the PCR product:
size.after.mut   <- sum(mut.res$trajectories[,ncol(mut.res$trajectories)]) 

drate.after.mut  <- size.after.mut/dilf.after.mut

# Sample the number of molecues after dilution from a Poisson distribution:
init.size.cln   <- rpois(1, drate.after.mut)

# Finish simulation if got no molecules after dilution:
if(init.size.cln == 0){
    cat('',file=cov.out)
    cat('',file=nwk.out)
    q(status=0)
}

# Simulate cleanup PCR:
pcoal.cln<-PCRcoal(
    initial.size    = init.size.cln,
    sample.size     = sample.size.cln,
    nr.cycles       = nr.cycles.cln,
    efficiencies    = eff.cln
)

cln.res <- sample.trs(pcoal.cln)

# Dilute after cleanup PCR:
size.after.cln      <- sum(cln.res$trajectories[,ncol(cln.res$trajectories)])
drate.after.cln     <- size.after.cln/dilf.after.cln

final.size          <- rpois(1, drate.after.cln)
# Finish simulation if got no molecules after dilution:
if(final.size == 0){
    cat('',file=cov.out)
    cat('',file=nwk.out)
    q(status=0)
}

tips                <- tree$tip.label
sm                  <- sample((1:length(tips)),final.size,replace=FALSE)

if(final.size > sample.size.mut){
    q(status=1)
}

# Simulate coverage PCR:
eff.cov         <- rep(pcr.eff, nr.cycles.cov)
pcoal.cov<-PCRcoal(
    initial.size    = length(sm),
    sample.size     = 1,        # We are not interested in the tree
    nr.cycles       = nr.cycles.cov,
    efficiencies    = eff.cov
)

final.traj <- sample.trs(pcoal.cov)$trajectories

cov.perc<-as.numeric(final.traj[,ncol(final.traj)])
cov.perc <- cov.perc/sum(cov.perc)

}

ii <- 1
cat('',file=cov.out)
for(i in sm){
    s <- tips[[ i ]]
    cat(paste(s,"\t",cov.perc[ii],"\n",sep=""), file=cov.out, append=TRUE)
    ii <- ii + 1
}

write.tree(tree, file=nwk.out)

