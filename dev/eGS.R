#setup
n=100
m=10000
freq=runif(m, 0.1, 0.9)
#generate sample
X=matrix(0, n, m)
for(i in 1:n) {
  X[i,]=rbinom(m, 2, freq)
  mL=ceiling(runif(runif(1, 1, m*0.01), 1, m))
  X[i,mL]=NA
}

#step 0 freq
Fq=colMeans(X, na.rm = T)/2 #remove missing and estimate frequency

#step 1 make G
G=matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:i) {
    cnt=0
    for (l in 1:m) {
      if (!is.na(X[i,l]) && !is.na(X[j,l])) {
        G[i,j]=G[i,j]+(X[i,l]-2*Fq[l])*(X[j,l]-2*Fq[l])/(2*Fq[l]*(1-Fq[l]))
        cnt = cnt+1
      }
    }
    if (cnt > 0) {
      G[i,j] = G[j,i] = G[i,j]/cnt
    }
  }
}
layout(matrix(1:2, 1, 2))
hist(G[row(G)<col(G)])
hist(diag(G))
ne=-1/mean(G[col(G)<row(G)]) #G_o
me=1/var(G[col(G)<row(G)]) #G_o

#Step 2
Eg=eigen(G)
layout(matrix(1:2, 1, 2))
barplot(Eg$values)
plot(Eg$vectors[,1], Eg$vectors[,2], xlab="PC 1", ylab="PC 2", bty='n')

#linear model
pc=5
eRes=array(NA, dim=c(pc, m, 6))
GC=array(NA, dim=c(2,pc))
for(k in 1:pc) {
  for(l in 1:m) {
    md=lm(Eg$vectors[,k]~X[,l])
    eRes[k, l, 1]=summary(md)$coefficients[2, 1]
    eRes[k, l, 2]=summary(md)$coefficients[2, 2]
    eRes[k, l, 3]=(eRes[k, l, 1]/eRes[k, l, 2])^2
    eRes[k, l, 4]=pchisq(eRes[k, l, 3], 1, lower.tail = F)
  }
  GC[1,k]=median(eRes[k, ,3])/qchisq(0.5, 1, lower.tail = TRUE)
  eRes[k,,5]=eRes[k,,3]/GC[1,k]
  eRes[k,,6]=pchisq(eRes[k,,5], 1, lower.tail = TRUE)
}
GC[2,]=Eg$values[1:pc]

layout(matrix(1:4, 2, 2))
barplot(main="GC", GC, beside = T, bty='n', border = F)
plot(main="PC", Eg$vectors[,1], Eg$vectors[,2], bty="n", pch=16)
hist(eRes[1,,4], main="p-value")
chiSeq=rchisq(m, 1)
qqplot(chiSeq, eRes[1,,3], col='red', pch=16, bty='n', xlab="Theory chisq", ylab="Obs chisq")
points(sort(chiSeq), sort(eRes[1,,5]), col="blue", pch=16)
abline(a=0, b=1)
