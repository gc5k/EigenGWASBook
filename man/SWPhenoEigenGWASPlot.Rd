\name{SWPhenoEigenGWASPlot}
\alias{SWPhenoEigenGWASPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EigenGWAS plot agaist conventional GWAS plot.
}
\description{
It is a bottom-up search for loci under selection. The overlapping peaks between the two plots indicate selection of the traits.
}
\usage{
SWPhenoEigenGWASPlot(phe, root, pc, kb = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phe}{
The phenotype GWAS result file. A string.
}
  \item{root}{
The file name of the eigengwas result file. A string.
}
  \item{pc}{
The index of the eigenvector to demonstrate.
}
  \item{kb}{
Window size. Default value is 5 for "5 kilo basepair".
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
SWPhenoEigenGWASPlot("arab.linear", "arab", 1, 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EigenGWAS }% use one of  RShowDoc("KEYWORDS")
