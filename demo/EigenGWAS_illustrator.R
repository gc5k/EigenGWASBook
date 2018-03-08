install.packages("~/git/EigenGWASFriends_0.1.0.tar.gz", repos = NULL, type = "source")
library(EigenGWASFriends)

#demo
FN="arab"
PC=5

RunEigenGWAS(FN, PC, T, "~/Documents/workspace/FromSVN/GEAR/")

AIM(FN, 5)

####GRM stats
grmStats(FN)

############PC plot
pcMatPlot(FN, c(1,2,3,4,5), ma=0.3)

####EigenValue plot
EigenQQPlot(FN, 1)
EigenValuePlot(FN, 5)
DeepEigenValuePlot(FN, 1, c(0.5,0.1, 0.05, 0.001, 0.0005, 0.0001, 0.00001))

###EigenGWAS plot
#EigenGWASPlot(FN, 1)

EigenGWASPlot(FN, 1)
SWEigenGWASPlot(FN, 1, 10)

####pheno eigen
#SWPhenoEigenGWASPlot("Arab452SetMAFclean_Naive_1.assoc.linear", FN, 10)
