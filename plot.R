truth = read.table("sim.out",as.is=T)
truth.mu = unlist(truth[2])
truth.se = unlist(truth[3])
h2 = 0.75

clr = c("black","#a6cee3","#1f78b4")

rand = read.table("run_random.out",as.is=T)
P = 50000
P.samp = 5000

svg("plot_krishna_kumar_pnas.svg",height=4)

## plot down-sampled results
keep = rand[,2] == "down-sampled"
vals = sort(rand[keep,3] / P.samp)
xlim = 2*range(vals)
n.samp = sum(keep)

rand.mu = mean(rand[keep,3]) / P.samp
rand.se = mean(rand[keep,4]) / P.samp

dens = density(vals)
# get the middle 95%
q = c(vals[ as.integer(n.samp * 0.025) ] , vals[ as.integer(n.samp * 0.975) ])

plot( dens , lwd=2 , xlab="per-SNP heritability" , ylab="density" , main="500 random h2 estimates from 5k SNPs" , xlim=xlim , bty="n"  )
# shade the middle 95%
x1 = min(which(dens$x >= q[1]))
x2 = max(which(dens$x <  q[2]))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=clr[2], border=NA))
# replot the main curve
lines( dens , lwd=2 )
# true results
abline( v = truth.mu / P , lwd=4 , col=clr[3] )
abline( v = c(truth.mu - 1.96*truth.se,truth.mu + 1.96*truth.se) / P , lty=2 , lwd=2 , col=clr[3] )
# down-sampled results
abline( v = c(rand.mu - 1.96*rand.se,rand.mu + 1.96*rand.se) , lty=3 , lwd=2 )

legend("topright",pt.cex=2,lwd=2,legend=c("Estimated h2 from 50k SNPs","Estimated 95CI from 50k SNPs","Estimated 95CI from 5k SNPs","Middle 95% of estimates from 5k SNPs") , lty=c(1:3,0) , pch=c(NA,NA,NA,15) , col=c(clr[3],clr[3],1,clr[2]), cex=0.75 , bg="white" , box.lwd=0 )
#box()

dev.off()
