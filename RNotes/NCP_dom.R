RP=500
n1=300
n2=700
N=n1+n2
f1=0.4
f2=0.6
para=matrix(0, RP, 6)
paraA=matrix(0, RP, 6)
for(i in 1:RP) {
  g1=rbinom(n1, 2, f1)
  g2=rbinom(n2, 2, f2)
  y=c(rep(1,n1), rep(0, n2))
  ys=scale(y)
  G=c(g1, g2)
  Gd=ifelse(G==1, 1, 0)
  Ga=c(g1, g2)
#  Gd=scale(Gd)
  mod=lm(ys~Gd)
  Ecov=sqrt(n1/N*n2/N)*(length(which(g1==1))/n1-length(which(g2==1))/n2)
  EV=mean(Gd)*(1-mean(Gd))
  Eb=Ecov/EV
  b=mod$coefficients[2]
  para[i,1]=Eb
  para[i,2]=b
  para[i,3]=sqrt(1/(N*var(Gd)))
  para[i,4]=summary(mod)$coefficients[2,2]
  para[i,5]=summary(mod)$coefficients[2,3]^2
  para[i,6]=summary(mod)$coefficients[2,4]
  
  modA=lm(ys~Ga)
  paraA[i,5]=summary(modA)$coefficients[2,3]^2
}
#layout(matrix(1:2, 1, 2))
#vF2=f1*(1-f1)*n1/N+f2*(1-f2)*n2/N
#ncpD=n1*n2/N * (2*f1*(1-f1)-2*f2*(1-f2))^2/vF2
#qqplot(main="Dom", rchisq(RP,1, ncp = ncpD), para[,5], pch=16, cex=0.5, bty='n',
#       xlab=expression(paste("Theoretical ", chi[1]^2)), ylab=expression(paste("Obs ",chi[1]^2)))
#abline(a=0, b=1)

ncpA=4*N*n1/N*n2/N*(mean(g1)/2-mean(g2)/2)^2/(2*mean(Ga)/2*(1-mean(Ga)/2))
qqplot(main="", rchisq(RP,1, ncp = ncpA), paraA[,5], pch=16, cex=0.5, bty='n',
       xlab=expression(paste("Theoretical ", chi[1]^2)), ylab=expression(paste("Observed ",chi[1]^2)))
abline(a=0, b=1, lty=2, col="grey")
