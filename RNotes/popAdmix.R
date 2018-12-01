N=c(100, 100, 100) #N
COL=rep(rep(1:length(N)), N)
M=10000 #M
fst=c(0.05, 0.02) #Fst

FP=matrix(0, sum(N), M)

P0=runif(M, 0.1, 0.9) #ancestry p0
z1_A=rbeta(M, P0*(1-fst[1])/fst[1], (1-P0)*(1-fst[1])/fst[1])
z1_B=rbeta(M, P0*(1-fst[1])/fst[1], (1-P0)*(1-fst[1])/fst[1])
Z1=rbind(z1_A, z1_B)

P1=(z1_A+z1_B)*0.5
z2_A=rbeta(M, P1*(1-fst[2])/fst[2], (1-P1)*(1-fst[2])/fst[2]) - P1
z2_B=rbeta(M, P1*(1-fst[2])/fst[2], (1-P1)*(1-fst[2])/fst[2]) - P1
Z2=rbind(z2_A, z2_B)

for(i in 1:N[1]) {
  FP[i,]=z1_A
}

for(i in (N[1]+1):(N[1]+N[2])) {
  w1 = rbeta(1, 5, 2)
#  w1=0.5
  w2 = 1 - w1
  FP[i,] = w1*z1_A + w2*z1_B 
  if(rbinom(1, 1, 0.5)>0) {
    FP[i,] = FP[i,] + w1*z2_A 
  } else {
    FP[i,] = FP[i,] + w1*z2_B
  }
}

for(i in (N[1]+N[2]+1):sum(N)) {
  FP[i,]=z1_B
}

for(i in 1:nrow(FP)) {
  if (length(which(FP[i,]<0)>0)) {
    FP[i, FP[i,]<0]=abs(FP[i,FP[i,]<0])
  }
  if (length(which(FP[i,]>1)>0)) {
    FP[i, FP[i,]>1]=2-(FP[i,FP[i,]>1])
  }
}

G=matrix(0, sum(N), M)
cnt=1
for(i in 1:sum(N)) {
    G[i,]=rbinom(M, 2, FP[i,])
}
sprintf("Fst1=%f, Fst2=%f", fst[1], fst[2])
sprintf("%d ind, %d marker",dim(G)[1], dim(G)[2])
Gs=apply(G, 2, scale)
GG=Gs %*% t(Gs)/M
EigenG=eigen(GG)
layout(matrix(1:2, 1, 2))
barplot(EigenG$values[1:10], border = F)
plot(EigenG$vectors[,1], EigenG$vectors[,2], bty='n', col=COL, pch=16, cex=0.5, xlab="E1", ylab="E2")

