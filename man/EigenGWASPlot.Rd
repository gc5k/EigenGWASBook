\name{EigenGWASPlot}
\alias{EigenGWASPlot}
\title{
EigenGWAS major plot
}
\description{
It gives a three row plot.
}
\usage{
EigenGWASPlot(root, pc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root}{
The file name of the eigengwas result file. A string.
}
  \item{pc}{
The index of the eigenvector to demonstrate.
}
}
\details{
Three plots will be illustrated. The first one is the plot for EigenGWAS with corrected p-value, the second one is the plot for Fst, and the last one is the contrast between p-value and Fst.
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
EigenGWAS("arab", 2) #plot for the second eigenvector
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EigenGWAS }% use one of  RShowDoc("KEYWORDS")
