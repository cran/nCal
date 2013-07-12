\name{drm.fit}
\alias{drm.fit}
\alias{getVarComponent.drc}
\title{Fit drm}
\description{ 
  \code{drm.fit} fit concentration-response curves using drm function from drc package. 

}
    

\usage{

drm.fit (formula, data, robust="mean", fit.4pl=FALSE, force.fit=TRUE, 
    weighting=FALSE, coln.weight=NULL, pow.weight=1, verbose=FALSE) 

getVarComponent.drc (fit) 

}

\arguments{

  \item{formula}{
a formula object.}
  \item{data}{
a data frame object.}
  \item{robust}{
a string. Passed to drm. See ?drm for more details.}
  \item{fit.4pl}{
boolean. If TRUE, 4PL model is fitted. If FALSE, 5PL model is fitted.}
  \item{force.fit}{
boolean. If TRUE, a fit is returned even if it is determined that the goodness of fit is poor.
}
  \item{weighting}{
boolean. If TRUE, weighted least square optimization is carried out. The weights are in the form of coln.weight^pow.weight. }
  \item{coln.weight}{
string. The column used in the weights. If NULL, default to the response variable in the formula.
}
  \item{pow.weight}{
numeric. The power used in the weights.
}
  \item{verbose}{
Boolean. If TRUE, verbose messages are printed.
}
  \item{fit}{
a drm object.
}

}

\details{
drm.fit differs from drc::drm in several aspects. 

(1) It tries several self start functions in order to get better fits. 

(2) It uses gof.threshold to report lack of fit.

(3) It tried to determine whether the standard deviation of the parameter estimates can be estimated.

(4) It implements weighting through coln.weight and pow.weight.
}

\value{
An object of type drm. 
}

\examples{

# simulate a dataset
set.seed(1)
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
fi=simulate1curve (p.eotaxin[1,], rep(log.conc,each=n.replicate), sd.e=0.2)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="Test", 
    assay_id="Run 1", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate),
    replicate=rep(1:n.replicate, 10))

fit = drm.fit(fi ~ expected_conc, data = dat.std, force.fit=TRUE)
plot(fit, log="xy")
fit

fit.2 = drm.fit(fi ~ expected_conc, data = dat.std, force.fit=TRUE, weighting=TRUE, pow.weight=-1)
plot(fit.2, log="xy")
fit.2



}

\keyword{ distribution }
