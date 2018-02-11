\name{grmStats}
\alias{grmStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary statistics of sample
}
\description{
The additive genetic relationship matrix gives the effective number of sample size [ne] and effective number of markers [me].

If the samples are unrelated, the effective number of sample size will be almost identical to the real sample size. The effective number of markers is ofter far smaller than the real marker numers, indicating the linkage disequilibrium betweem markers.
}

\usage{
grmStats(root)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root}{
  The file name of the additive genetic relationship matrix. The workding file name will be "root.grm.gz", and ".grm.gz" is added by default.
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
grmStats("arab")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EigenGWAS }% use one of  RShowDoc("KEYWORDS")
