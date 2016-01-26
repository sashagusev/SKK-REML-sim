source('func_reml.R')
load("sim.rbin")

seed = as.numeric(commandArgs(trailingOnly=T))[1]
set.seed(seed)

P.sample = 5000

K = list()

## regenerate kinship with subsampled SNPs
keep = sample(1:P,P.sample)
K[[1]] = Z[,keep] %*% t(Z[,keep]) / P.sample
reml = aiML(K,y,c(0.5,0.5),verbose=F)
cat( seed , "down-sampled" , reml$h2 , reml$se[1] , '\n' )
