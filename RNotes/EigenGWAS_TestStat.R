rep=1
EV=matrix(0, rep, 10)
typeI=matrix(0, rep, 4)
fst=c(0.002, 0.01)
N1=c(100, 500, 1000)
N2=1000
Ml=c(8000, 2000)
M=sum(Ml)
ME=matrix(0, length(N1), 3)
for(s in 1:length(N1)) {
  for(r in 1:rep) {
    P=runif(M, 0.2, 0.8)
    p1=runif(Ml[1], 0.2, 0.8)
    p2=runif(Ml[2], 0.2, 0.8)
    P=c(p1, p2)
    P1=c(rbeta(Ml[1], p1*(1-fst[1])/fst[1], (1-p1)*(1-fst[1])/fst[1]),
         rbeta(Ml[2], p2*(1-fst[2])/fst[2], (1-p2)*(1-fst[2])/fst[2]))
    P2=c(rbeta(Ml[1], p1*(1-fst[1])/fst[1], (1-p1)*(1-fst[1])/fst[1]),
         rbeta(Ml[2], p2*(1-fst[2])/fst[2], (1-p2)*(1-fst[2])/fst[2]))

    Gn=matrix(0, nrow=N1[s]+N2, ncol=M)
    for(i in 1:N1[s]) {
      Gn[i,] = rbinom(M, 2, P1)
    }

    for(i in (N1[s]+1):(N2+N1[s])) {
      Gn[i,] = rbinom(M, 2, P2)
    }
    Frq1=apply(Gn[1:N1[s],], 2, mean)/2
    Frq2=apply(Gn[(N1[s]+1):(N1[s]+N2),], 2, mean)/2
    FrqM=apply(Gn, 2, mean)/2
    Fst=2*(N1[s]/(N1[s]+N2)*(Frq1-FrqM)^2 + N2/(N1[s]+N2) * (Frq2-FrqM)^2)/(FrqM*(1-FrqM))
    FstN=(Frq1-Frq2)^2/(2*FrqM*(1-FrqM))
    plot(Fst, FstN)
  }

  GnS=apply(Gn, 2, scale)
  G=GnS %*% t(GnS)/M
  EigenGN=eigen(G)
  
  RegB=matrix(0, M, 4)
  for(i in 1:M) {
    mod=lm(EigenGN$vectors[,1]~Gn[,i])
    RegB[i,1]=summary(mod)$coefficients[2,1]
    RegB[i,2]=summary(mod)$coefficients[2,2]
  }
  RegB[,3]=RegB[,1]^2/RegB[,2]^2
  qqplot(rchisq(M, 1), RegB[,3], bty='n', pch=16)
  abline(a=0, b=1, col="red")
  gc=median(RegB[,3])/0.455
  abline(a=0, b=gc, col="blue")
  
  ME[s,1]=mean(Fst)*(N1[s]+N2)
  ME[s,2]=gc
  ME[s,3]=EigenGN$values[1]
}
rownames(ME)=N1
barplot(t(ME), beside = T, border = F)
legend("topleft", legend = c("Fst", "GC", "Eigenvalue"), pch=15, col=c("black", "grey50", "grey"), bty='n')



EG=matrix(0, 50, 1)
for(rp in 1:50) {
  N=500
  M=5000
  p=runif(M, 0.2, 0.8)
  X=matrix(0, N, M)
  for(i in 1:M) {
    X[,i] = rbinom(N, 2, p[i])
  }
  mu=apply(X, 2, mean)
  sig=apply(X, 2, var)
  Xs=apply(X, 2, scale)
  G=Xs %*% t(Xs) / M
  rc=rchisq(nrow(G), nrow(G)-1)/nrow(G)
  qqplot(rc, diag(G), xlab=expression(paste(chi[N]^2, "/N")), ylab="Diag (G)", bty='n', pch=16)
  abline(a=0, b=1, col="red")
  mod=lm(sort(diag(G))~sort(rc))
  Eg=eigen(G)
  EG[rp,1]=Eg$values[1]
}
Emu=(sqrt(M-1)+sqrt(N))^2/M
Esd=sqrt(Emu/M)*(1/sqrt(M-1/2)+1/sqrt(N-1/2))^(1/3)

hist((EG/(sum(diag(G))/N)-Emu)/Esd)
x=matrix(rnorm(N*M), N, M)
