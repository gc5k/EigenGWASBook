\name{EigenQQPlot}
\alias{EigenQQPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This is the EigenGWAS interface for illustrating QQ plot.
}
\description{
QQ plot for chisq tests with and without correction for lambdaGC.
}
\usage{
EigenQQPlot(root, pc)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root}{
The file name of the eigengwas result file. A string.
  }
  \item{pc}{
The eigenvector to be illustrated. An integer.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A plot is generated.
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
EigenQQPlot("arab", 1) #QQ plot for eigengwas for the 1st eigenvector
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EigenGWAS }% use one of  RShowDoc("KEYWORDS")
