####mapfile
mapInfoCur=mapInfo=read.table(file = "lr_10000.plk.map", as.is = T)
colnames(mapInfoCur)=c("Chr", "SNP", "Gdis", "BP")

####pedfile
conn <- gzfile("lr_10000.plk.ped.gz", "rt") #repalce the filename with your own ped
#conn <- gzfile("LRWithPosition_2691.plk.ped.gz", "rt") #repalce the filename with your own ped

ped=read.table(conn, as.is = T, sep = "\t")#remove the first 6 cols
pedInfoCur=pedInfo=ped[,c(1:6)]
ped=ped[,-c(1:6)]
hped=ped[,seq(1, ncol(ped), 2)] #inbred, one haploid is enough


##QC stats
freq=colMeans(hped, na.rm = T) #freq
sMiss=array(0, dim=nrow(hped)) #sample-level missing
lMiss=array(0, dim=ncol(hped)) #locus-leve missing
hped_T=t(hped)
for(i in 1:length(sMiss)) {
  sMiss[i]=length(which(is.na(hped_T[,i])))
}
rm(hped_T)

for(i in 1:length(lMiss)) {
  lMiss[i]=length(which(is.na(hped[,i])))
}

layout(matrix(1:3, ncol=3))
plot(main="Frequency", freq, pch=16, cex=0.5)
plot(main="Individual missing rate", sMiss/ncol(hped), pch=16, cex=0.5)
plot(main="Locus missing rate", lMiss/nrow(hped), pch=16, cex=0.5)

#########QC steps
#QC 1 for MAF
QCcut=list("fqMin"=0.03, "fqMax"=0.97, "IndMiss"=0.4)
fQC=which(freq > QCcut$fqMin & freq < QCcut$fqMax) #remove freq <0.03 and > 0.97
#QC 2 for ind miss
iQC=which(sMiss/ncol(hped) < QCcut$IndMiss) #individual genotyping missing rate

hpedQC=hped[iQC,fQC]
hpedQC_T=t(hpedQC)
freqQC=colMeans(hpedQC, na.rm = T)
vQC=freqQC*(1-freqQC)

mapInfoCur=mapInfoCur[fQC,]
pedInfoCur=pedInfo[iQC,]
##test sample size, the maxima of testInd=nrow(hpedQC)
#testInd=100
testInd=nrow(hpedQC)

#########make grm
G=matrix(0, nrow = testInd, ncol = testInd)
for(i in 1:testInd) {
  print(paste0(i, "/", testInd))

##slow  x1=ifelse(!is.na(as.numeric(hpedQC[i,])), as.numeric(hpedQC[i,]), freqQC) 
  ##fast
  x1=ifelse(!is.na(as.numeric(hpedQC_T[,i])), as.numeric(hpedQC_T[,i]), freqQC)
  s1=(x1-freqQC)/sqrt(vQC)
  for(j in i:testInd) {
##slow    x2=ifelse(!is.na(as.numeric(hpedQC[j,])), as.numeric(hpedQC[j,]), freqQC)
    ##fast
    x2=ifelse(!is.na(as.numeric(hpedQC_T[,j])), as.numeric(hpedQC_T[,j]), freqQC)
    
    s2=(x2-freqQC)/sqrt(vQC)
    effM=length(freqQC)-length(which(x1==0 | x2==0))
    G[j,i]=G[i,j]=sum(s1*s2)/effM
  }
}
rm(hpedQC_T)

##save GRM
write.table(G, "G.txt", row.names = F, col.names = F, quote = F)

##eigen
g1=read.table("G.txt", as.is = T)
eG=eigen(g1)
write.table(eG$values, "Gvalue.txt", row.names = F, col.names = F, quote = F)
write.table(eG$vectors, "Gvec.txt", row.names = F, col.names = F, quote = F)

eVec=read.table("Gvec.txt", as.is = T)
eVal=read.table("Gvalue.txt", as.is = T)
layout(matrix(1:2, 1, 2))
barplot(eVal[1:5,1], main="Eigenvalues")
plot(eVec$V1, eVec$V2, pch=16, cex=0.5)

pcIdx=1 #if you want to try another eigenvector 'x' try eVec[,x]
pcY=eVec[,pcIdx]   

fqAll=colMeans(hpedQC, na.rm = T)
MOD=matrix(2, nrow=ncol(hpedQC), 8)
for(i in 1:nrow(MOD)) {
  mod=lm(pcY~hpedQC[,i]) 
  sm=summary(mod)
  if(nrow(sm$coefficients)>1) {
    MOD[i,1:4]=sm$coefficients[2,]

    pcPos=which(pcY > 0 & !is.na(hpedQC[,i]))
    pcNeg=which(pcY <= 0 & !is.na(hpedQC[,i]))
    fqP=mean(hpedQC[pcPos,i])
    fqN=mean(hpedQC[pcNeg,i])
    MOD[i,8]=(length(pcPos)/(length(pcPos)+length(pcNeg))*(fqP-fqAll[i])^2
              +length(pcNeg)/(length(pcPos)+length(pcNeg))*(fqN-fqAll[i])^2)/(fqAll[i]*(1-fqAll[i]))
  }
}

goodSNP=which(MOD[,4]!=2)
MOD1=MOD[goodSNP,]

colnames(MOD1)=c(colnames(sm$coefficients), "Chisq", "PGC", "ChisqGC", "Fst")
gc=qchisq(median(MOD1[,4]), 1, lower.tail = F)/0.455
#using approximating that t^2=chi when df is large
#MOD1[,5]=qchisq(MOD1[,4], 1, lower.tail = F)
MOD1[,5]=MOD1[,3]^2
#.Machine$double.xmin
MOD1[,6]=-1*pchisq(MOD1[,5]/gc, 1, lower.tail = F, log.p = T)/log(10) #-log10(p)
MOD1[,7]=MOD1[,5]/gc

Res=cbind(mapInfoCur[goodSNP,], MOD1)
Res=Res[order(Res$Chr, Res$BP),] #assuming V1=chr, V4=pos in map file

layout(matrix(1:3, 1, 3))
qqplot(-log10(runif(nrow(Res))), -log10(Res$`Pr(>|t|)`), pch=16, cex=0.5,
       main="p-value", bty='n', xlab="E(-log10(p))", ylab="-log10(p)",
       bty='n')
abline(a=0, b=1, col="red", lty=2)
qqplot(-log10(runif(nrow(Res))), Res$PGC, pch=16, cex=0.5,
       main="p-value after GC", xlab="E(-log10(p))", ylab="-log10(p)",
       bty='n')
abline(a=0, b=1, col="red", lty=2)
plot(Res$Fst, Res$ChisqGC, xlab=expression(paste(F[st])), ylab=expression(paste(chi["1.GC"]^2)), 
     bty='l', pch=16, cex=0.5)

print("Summary:")
print(QCcut)
print(paste0(nrow(pedInfoCur), " samples and ", nrow(Res), " SNPs included in the analysis."))
print(paste0("Removed ", nrow(pedInfo)-nrow(pedInfoCur), " samples"))
print(paste0("Removed ", nrow(mapInfo)-nrow(mapInfoCur), " snps"))
print(paste("GC=", gc, ", eigenvalue is ", eVal[pcIdx,1]))
