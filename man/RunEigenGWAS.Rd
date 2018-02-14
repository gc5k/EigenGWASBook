\name{RunEigenGWAS}
\alias{RunEigenGWAS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Mapping genes under selection/domestication with EigenGWAS algorithm
}
\description{
 It mappes genes with EigenGWAS algorithm, which could provide the eigenvalues, eigenvectors, lambdaGC and genes related to selection.
}
\usage{
 RunEigenGWAS(dataName, PC, inbred, gearPath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataName}{
  File name which should be inputted into EigenGWAS
}
  \item{PC}{
  The number of PC to be used.
}
  \item{inbred}{
  The type of population. By default inbred = T.
}
  \item{gearPath}{
  The path where gear.jar locates. By default gearPath=".", which means that gear.jar locates in the path where data file is.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{*.No.egwas }{Gives P valus and other statistics for all markers related to PC NO.}
  \item{*.eigenval }{Gives all eigenvalues}
  \item{*.eigenvec }{Gives all eigenvectors}
  \item{*.grm.gz }{Gives the genetic relatedness matrix}
  \item{*.grm.id }{Gives ID of the genetic relatedness matrix}
%% ...
}
\references{
Chen, G.B., et al, EigenGWAS: finding loci under selection through genome-wide association studies of eigenvectors in structured populations, Heredity, 2016, 117:51-61
}
\author{
Guo-Bo Chen, chenguobo@gmail.com
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
RunEigenGWAS("arab", 2, inbred = T, gearPath="your path")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RunEigenGWAS }% use one of  RShowDoc("KEYWORDS")

