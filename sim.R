source('func_reml.R')

# reproducibility!
set.seed(1234)

N = 2000
P = 50000
h2 = 0.75

## generate genotype matrix with MAF = 0.5
Z = matrix(rbinom(N*P,2,0.5),nrow=N,ncol=P)

# standardize
for ( i in 1:P ) {
Z[,i] = (Z[,i] - mean(Z[,i]))/sd(Z[,i])
if ( i %% 1000 == 0 ) cat(i,'\n',file=stderr())
}

## generate kinship
A = Z %*% t(Z) / P

## generate a heritable phenotype
# 1. SNP effects
u = rnorm(P,0,1)
# 2. genetic value
g = Z %*% u
g = (g - mean(g))/sd(g)
# 3. add environmental noise
y = sqrt(h2) * g + rnorm(N,0,sqrt(1-h2))
y = (y - mean(y))/sd(y)

## save data
save(N,P,h2,y,g,u,A,Z,file="sim.rbin")

## estimate heritability
K = list()
K[[1]] = A
reml = aiML(K,y,c(0.5,0.5))
cat( "estimate" , reml$h2 , reml$se[1] , '\n' )
