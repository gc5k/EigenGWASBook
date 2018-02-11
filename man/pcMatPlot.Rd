\name{pcMatPlot}
\alias{pcMatPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrixized plot for eigenvectors
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pcMatPlot(root, PC = c(1,2), COL = "grey", ma = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root}{
File name, the working file name will be "root.eigenvec". A string.
}
  \item{PC}{
Indeces of eigenvectors to be used. An array, by default pc=c(1,2).
}
  \item{COL}{
Color to be used. By default, it is "grey". A vector can be specified to give the color for each sample.
}
  \item{ma}{
  Margin of the plot. By default, ma=0.3
}
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
pcMatPlot("arab", c(1:5), col="green", ma=0.2) #using eigenvector 1~5, and the color is green.

pcMatPlot("arab", c(1,2,4,5), col=c(rep("green", 100), rep("red", 195)), ma=0.2) #using eigenvector 1, 2, 4, and 5,
#and the color is green for the first 100 samples and red for the rest 195 samples.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{EigenGWAS}% use one of  RShowDoc("KEYWORDS")
