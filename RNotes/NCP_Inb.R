RP=1000
n1=500
n2=500
N=n1+n2
f1=0.4
f2=0.6
paraA=matrix(0, RP, 1)
for(i in 1:RP) {
  g1=rbinom(n1, 1, f1)*2
  g2=rbinom(n2, 1, f2)*2
  y=c(rep(1,n1), rep(0, n2))
  ys=scale(y)
  G=c(g1, g2)
  Ga=c(g1, g2)

  modA=lm(ys~Ga)
  paraA[i,1]=summary(modA)$coefficients[2,3]^2
}

ncpA=N*n1/N*n2/N*(mean(g1)/2-mean(g2)/2)^2/(mean(Ga)/2*(1-mean(Ga)/2))
qqplot(main="", rchisq(RP,1, ncp = ncpA), paraA[,1], pch=16, cex=0.5, bty='n',
       xlab=expression(paste("Theoretical ", chi[1]^2)), ylab=expression(paste("Observed ",chi[1]^2)))
abline(a=0, b=1, lty=2, col="grey")
