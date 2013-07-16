\name{read.luminex.xls}
\alias{read.luminex.xls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read a Luminex File
}
\description{
Read a Luminex raw output .xls file 
}
\usage{
read.luminex.xls (file,verbose = FALSE,sheets = NULL, assay_id=NULL, 
    na.strings = c("NA", "#DIV/0!"), ..., perl = "perl") 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
string. The file name.
}
  \item{verbose}{
Boolean. If TRUE, some debug messages are printed.
}
  \item{sheets}{
vector of strings. The names of the sheets to be read. If NULL, all sheets are read.
}
  \item{assay_id}{
string. Used to name the assay that generated the data file. If not provided, a random number is generated as the name.
}
  \item{na.strings}{
vector of strings. Strings to be mapped to NA.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{perl}{
string. Command used to execute perl. If perl.exe is not in the path, the full path can be specified, e.g. "C:/Perl64/bin/perl.exe"
}
}

\examples{

#begin=Sys.time()

# example from https://www.labkey.org/wiki/home/Documentation/page.view?name=luminexFileFormats
dat = read.luminex.xls(paste(system.file(package="nCal")[1],
    '/misc/02-14A22-IgA-Biotin-tiny.xls', sep=""), verbose=TRUE)
out = ncal(log(fi)~expected_conc, dat, return.fits = TRUE, plot.se.profile=FALSE)
out

library(RUnit)
checkEqualsNumeric(unlist(out[1:3,c("fi")]), c(15,45,19.33908), tolerance=1e-3)
checkEqualsNumeric(unlist(out[12,c("fi")]), c(183.73622), tolerance=1e-3)

#end=Sys.time();print(end-begin)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
