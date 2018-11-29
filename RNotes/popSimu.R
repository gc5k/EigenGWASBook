N=c(100, 100, 100, 100, 100) #N
COL=rep(c(1,2,3,4,5), N)

M=5000 #M
fst=c(0.05, 0.02) #Fst

FP=matrix(0, length(N), M)

P0=runif(M, 0.1, 0.9) #ancestry p0
FP[1,]=rbeta(M, P0*(1-fst[1])/fst[1], (1-P0)*(1-fst[1])/fst[1])
FP[2,]=rbeta(M, P0*(1-fst[1])/fst[1], (1-P0)*(1-fst[1])/fst[1])

P1=(FP[1,]+FP[2,])*0.5
FP[3,]=P1
FP[4,]=rbeta(M, P1*(1-fst[2])/fst[2], (1-P1)*(1-fst[2])/fst[2])
FP[5,]=rbeta(M, P1*(1-fst[2])/fst[2], (1-P1)*(1-fst[2])/fst[2])

G=matrix(0, sum(N), M)
cnt=1
for(i in 1:length(N)) {
  for(j in 1:N[i]) {
    G[cnt,]=rbinom(M, 2, FP[i,])
    cnt=cnt+1
  }
}
sprintf("Fst1=%f, Fst2=%f", fst[1], fst[2])
sprintf("%d ind, %d marker",dim(G)[1], dim(G)[2])
Gs=apply(G, 2, scale)
GG=Gs %*% t(Gs)/M
EigenG=eigen(GG)
layout(matrix(1:2, 1, 2))
barplot(EigenG$values[1:10], border = F)
plot(EigenG$vectors[,1], EigenG$vectors[,2], bty='n', col=COL, pch=16, cex=0.5, xlab="E1", ylab="E2")
