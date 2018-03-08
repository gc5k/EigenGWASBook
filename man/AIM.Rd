\name{AIM}
\alias{AIM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Picking up ancestry informative markers for the samples
}
\description{
It harvests all signals at a region.
}
\usage{
AIM(root, PC, "bonferroni", T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root}{
The file name of the eigengwas result file. A string.
}
  \item{PC}{
The eigenvector to be illustrated. An integer.
}
  \item{pcutoff}{
P-value cutoff for signals. By default, it is 0.05/(number of markers) or bonferroni correction. The user can specify it by giving a cutoff pvalue.
}
  \item{GCcorrection}{
GC correction for p-values. By default, it uses corrected p-values.
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
AIM("arab", 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EigenGWAS }% use one of  RShowDoc("KEYWORDS")
