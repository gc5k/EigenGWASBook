\name{SWEigenGWASPlot}
\alias{SWEigenGWASPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sliding window method for EigenGWAS plot
}
\description{
When sliding window is used, the regional signals may be shown more clearly, especially for polygenic traits.
}
\usage{
SWEigenGWASPlot(root, pc, kb = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root}{
The file name of the eigengwas result file. A string.
}
  \item{pc}{
The eigenvector to be illustrated. An integer.
}
  \item{kb}{
Window size. Default value is 5 for "5 kilo basepair".
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
SWEigenGWASPlot(root, 2, kb = 5) #5 kb window size,
#and the second eigenvector

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EigenGWAS }% __ONLY ONE__ keyword per line
