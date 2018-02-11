\name{DeepEigenValuePlot}
\alias{DeepEigenValuePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deep evaluation of lambdaGC and comparison with eigenvalue
}
\description{
The default lambdaGC is calculated on the median p-value, but in this function, the p-value can be set at any quantile.
}
\usage{
DeepEigenValuePlot(root, PC, cutV, pdf = "n")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root}{
File name, the working file name will be "root.eigenval" for reading eigenvalues, and "root.PC.egwas" for reading EigenGWAS result files. A string.
}
  \item{PC}{
  The eigenvector for deep evaluation. An integer.
}
  \item{qcut}{
  An array for quantile level. An array, by default qcut = c(0.25, 0.5, 0.75).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
DeepEigenValuePlot("arab", 1, c(0.5, 0.05, 0.005, 0.001, 0.0005, 0.0001))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EigenGWAS }% use one of  RShowDoc("KEYWORDS")
