\name{plot.bcrm}
\alias{plot.bcrm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot function for the bcrm object
}
\description{
Plot function for the bcrm object
}
\usage{
\method{plot}{bcrm} (x, assay_id = NULL, add = FALSE, lcol = 1, main = NULL, fit.2 = NULL, 
    lwd = 0.1, points.only = FALSE, all.lines.only = FALSE, t = NULL, ylim = NULL, 
    same.ylim = FALSE, lty3 = NULL, fit.3 = NULL, lcol2 = NULL, lcol3 = NULL, xlab = NULL,
    ylab = NULL, col.outliers = TRUE, lty = 1, cex = 1, log = "x", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a bcrm object.
}
  \item{assay_id}{
a string. Assay identifier.}
  \item{add}{
Boolean. If TRUE, adding to an existing plot.}
  \item{lcol}{
integer. Line color.}
  \item{main}{
string. Plot title.}
  \item{fit.2}{
a bcrm object. A second fit object to be plotted together with x.}
  \item{lwd}{
numeric. Line width.}
  \item{points.only}{
Boolean. If TRUE, only the data points are plotted and not the fitted curves
}
  \item{all.lines.only}{
Boolean. If TRUE, only the fitted curves are plotted.
}
  \item{t}{
a numeric vector. The log concentrations.}
  \item{ylim}{
a numeric vector. The usual ylim of plot function.}
  \item{same.ylim}{
Boolean. If TRUE, all fitted curves are plotted with the same ylim.}
  \item{lty3}{
integer. lty for plotting fit.3.
}
  \item{fit.3}{
bcrm object. A third optional bcrm object to be plotted.}
  \item{lcol2}{
integer. Line color for plotting fit.2.}
  \item{lcol3}{
integer. Line color for plotting fit.3.}
  \item{xlab}{
string. The usual xlab argument of plot function.}
  \item{ylab}{
string. The usual ylab argument of plot function.}
  \item{col.outliers}{
Boolean. If TRUE, outliers are colored differently.}
  \item{lty}{
integer. Line type for plotting x.}
  \item{cex}{
numeric. The usual cex argument of plot function.}
  \item{log}{
string. If it is "x", then the x axis is labeled on the scale of concentration; otherwise, it is labeled on the scale of log concentration.}
  \item{...}{ Not passed to other functions, but needed to keep some similarity with plot.drm. }
}
