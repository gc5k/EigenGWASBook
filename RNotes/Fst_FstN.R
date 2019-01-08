rep=1
EV=matrix(0, rep, 10)
typeI=matrix(0, rep, 4)
fst=0.05
N1=100
N2=1000
M=50000
for(r in 1:rep)
{
  P=runif(M, 0.05, 0.95)
  P1=rbeta(M, P*(1-fst)/fst, (1-P)*(1-fst)/fst)
  P2=rbeta(M, P*(1-fst)/fst, (1-P)*(1-fst)/fst)
  
  Gn=matrix(0, nrow=N1+N2, ncol=M)
  for(i in 1:N1)
  {
    Gn[i,] = rbinom(M, 2, P1)
  }
  
  for(i in (N1+1):(N2+N1))
  {
    Gn[i,] = rbinom(M, 2, P2)
  }
  Frq1=apply(Gn[1:N1,], 2, mean)/2
  Frq2=apply(Gn[(N1+1):(N1+N2),], 2, mean)/2
  FrqM=apply(Gn, 2, mean)/2
  Fst=2*(N1/(N1+N2)*(Frq1-FrqM)^2 + N2/(N1+N2) * (Frq2-FrqM)^2)/(FrqM*(1-FrqM))
  FstN=(Frq1-Frq2)^2/(2*FrqM*(1-FrqM))
  
}
