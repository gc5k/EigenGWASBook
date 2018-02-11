\name{EigenValuePlot}
\alias{EigenValuePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparison between eigenvalue and lambdaGC
}
\description{
It compares the eigenvalue and lambdaGC. Often, the eigenvalue will be very close to its corresponing lambdaGC.
}
\usage{
EigenValuePlot(root, PC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root}{
File name, the working file name will be "root.eigenval" for reading eigenvalues, and "root.*.egwas" for reading EigenGWAS result files. A string.
}
  \item{PC}{
  The number of eigenvalues to be used. By default PC = 1.
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
Assuming the median of the p-value is 0.35, the corresponding lambdagc is qchisq(0.35, 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F) = 1.92.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
EigenValuePlot("arab", 2) #plot for the first two eigenvalues and lambdaGC
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EigenGWAS }% use one of  RShowDoc("KEYWORDS")
