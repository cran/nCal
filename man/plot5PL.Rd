\name{plot5PL}
\alias{plot5PL}
\title{ Plot a 5PL function}
\description{ 
  \code{plot5PL} plots a 5PL function
}

\usage{
plot5PL (param, xlim, ylim=NULL, col=NULL, lty=NULL, lwd=1, plot.legend=FALSE, add=FALSE, 
    legend=NULL, main=NULL, x.exp=FALSE, xlab=NULL, ylab=NULL, xaxt="s") 
}

\arguments{

  \item{param}{Can be a vector or a matrix. If matrix, each row represents one curve. First check for classical parameterization, then gh, then ED50 parameterization. e/b/f can be provided in log form.}
  \item{xlim}{range of t, i.e. log (concentration)}
  \item{ylim}{ylim}
  \item{col}{col}
  \item{lty}{lty}
  \item{lwd}{lwd}
  \item{plot.legend}{Boolean. Indicates whether to plot legend.}
  \item{add}{Boolean. Indicates whether to add to an existing plot or create a new plot.}
  \item{legend}{vector. Legend.}
  \item{main}{string. Title.}
  \item{x.exp}{Boolean. If TRUE, x axis in on the concentration scale, not the log scale.}
  \item{xlab}{xlab}
  \item{ylab}{ylab}
  \item{xaxt}{xaxt}

}

\author{ Youyi Fong \email{yfong@fhcrc.org} }

\keyword{ distribution }
