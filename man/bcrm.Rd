\name{bcrm}
\alias{bcrm}
\alias{getVarComponent.bcrm}
\alias{vcov.bcrm}
\alias{coef.bcrm}
\alias{print.bcrm}
\alias{get.single.fit}
\title{ Bayesian Concentration-Response Model}
\description{ 
  \code{bcrm} fit concentration-response curves with a Bayesian random effects model using JAGS
}
    

\usage{

bcrm (formula, data, error.model, 
    informative.prior, prior.sensitivity=NULL, 
    n.iter=1e5, jags.seed=1, n.thin=NULL, T.init=NULL, keep.data=TRUE, 
    keep.jags.samples=FALSE, standards.only=TRUE, verbose=FALSE,
    t.unk.truth=NULL, params.true=NULL)

\method{print}{bcrm} (x, ...)
\method{coef}{bcrm} (object, type="gh", ...)
\method{vcov}{bcrm} (object, type="gh", ...)
\method{getVarComponent}{bcrm} (fit) 
get.single.fit (fit, assay_id) 

}

\arguments{

  \item{formula}{formula. Gives the response column and concentration column.}
  \item{data}{ a data frame. Each row represents the measurement from one well/bead_type. See details }
  \item{error.model}{ a string specifying the error.model to fit. See details.}
  \item{informative.prior}{boolean. If TRUE, use informative priors. See details.}
  \item{n.iter}{ a number indicating the number of iterations to run.}
  \item{jags.seed}{ a number to seed the random number generator used within jags.}
  \item{n.thin}{ a number specifying the thinning factor applied to the jags samples.}
  \item{keep.data}{ a boolean. If TRUE, the fit object being returned has an element named "data", which contains ordered input data.}
  \item{keep.jags.samples}{boolean. If TRUE, the fit object being returned has an element named "jags.samples". coef samples will always be saved in "coef.samples".}
  \item{t.unk.truth}{True unknown concentrations, for simulation study use only.}
  \item{params.true}{True curve parameters, for simulation study use only.}
  \item{T.init}{a integer vector. Initial values for mixture indicators.}
  \item{prior.sensitivity}{integer. A number between 1 and 4. Change priors.}
  \item{standards.only}{boolean. If TRUE, data is subset to standard samples only.}
  \item{verbose}{boolean. If TRUE, debug messages are printed.}
  
  \item{x}{bcrm fit object.}
  \item{object}{bcrm fit object.}
  \item{type}{string. 5PL parameterization.}
  \item{fit}{bcrm fit object.}
  \item{...}{...}
  \item{assay_id}{string. Label for the assay run.}

}

\details{

data is expected to contain one to many plates with the same analyte. 
    \itemize{
        \item well_role     Defines the role of a well. This should be from {Standard, Unknown, ...}. Standard wells are used to generate standard curves, and concentrations of the substance in the Unknown well will be estimated
        \item assay_id      Identifies an assay, which is defined to be a collection of Standard and non-Standard wells. Measured fi from the Standard wells are used to create a set of standard curves, one of each bead type.        
                            Based on the standard curves and the measured fi from the non-Standard wells, concentrations of the substance in the non-Standard wells will be estimated. 
                            An assay can be a plate, if every plate has Standard wells; or it can be multiple plates run by one technician at one time, if there are only one set of Standard wells on these plates
        \item dilution      Standard samples are often prepared by starting with one sample and doing serial dilutions. 
                            Unknown samples may be measured at several dilutions so that one of the dilutions may fall into the more reliable region of the standard curve
        \item replicate     Index of technical replicates for a sample. Typical values are 1 or 2. May be used in plotting. Optional
        \item expected_conc Standard sapmles have expected concentrations. If this column is present, the dilution and starting_conc are optional and will not be used. This column does not apply to non-Standard samples
    }
    
Main error.model supported: drc, classical_norm, classical_t4, classical_mixnorm, classical_lar1, gh_norm, gh_mixnorm, gh_lar1 
Also support: classical_replicate_re, gh_replicate_re, gh_tar1

Only two replicates are supported for now for the correlated noise models.

Sometimes jags.model fails with one .RNG.seed. The function will increase it by 1 and try again. The function tries three times before giving up.


}

\value{
An object of type bcrm. 
}

\author{ Youyi Fong \email{yfong@fhcrc.org} }
\references{
Fong, Y., Wakefield, J., DeRosa, S., Frahm, N. (2012) A robust Bayesian random effects model for nonlinear calibration problems, Biometrics, 68:1103-1112.
}
\examples{

set.seed(1)
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
fi=simulate1curve (p.eotaxin[1,], rep(log.conc,each=n.replicate), sd.e=0.3)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="test", 
    assay_id="assay1", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate))
# add unknown
dat.unk=rbind(
      data.frame(fi=exp(6.75), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=1, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(6.70), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=2, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(3), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=3, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(4.4), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=4, well_role="Unknown", dilution=10)
)
dat=rbind(dat.std, dat.unk)
# second plate
fi=simulate1curve (p.eotaxin[2,], rep(log.conc,each=n.replicate), sd.e=0.3)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="test", 
    assay_id="assay2", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate))
# add unknown
dat.unk=rbind(
      data.frame(fi=exp(6.75), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=1, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(6.70), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=2, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(3), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=3, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(4.4), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=4, well_role="Unknown", dilution=10)
)
dat=rbind(dat, dat.std, dat.unk)

fits = bcrm(log(fi)~expected_conc, dat, error.model="gh_norm", informative.prior=TRUE, n.iter=1e4)
par(mfrow=c(1,2))
plot(fits)


\dontrun{
# this block of code takes a few minutes to run

# Example from Fong et al. (2012) 
fits.t4 = bcrm (log(fi)~expected_conc, dat.QIL3, error.model="gh_t4", informative.prior=TRUE)
par(mfrow=c(2,3))
plot(fits.t4)

fits.norm = bcrm (log(fi)~expected_conc, dat.QIL3, error.model="gh_norm", informative.prior=TRUE)
par(mfrow=c(2,3))
plot(fits.norm)


# Example from Fong et al (2013)
fit.bcrm=bcrm(log(fi)~expected_conc, subset(hier.model.ex.2,assay_id %in% c("Run 1","Run 2","Run 3",
    "Run 4")), error.model="gh_t4", informative.prior=T)

for (p in c("Run 1","Run 2","Run 3","Run 4")) {
    fit.drm = drm(log(fi)~expected_conc, data=subset(hier.model.ex.2, assay_id==p), 
        fct=LL.5(), robust="median")
    plot(fit.drm, type="all", col="black", main=p, lty=2)
    plot(get.single.fit(fit.bcrm, p), add=T, log="x", col=1)
    mylegend(x=9,legend=c("drm, median","bcrm, t4"),lty=c(2,1))
}
    
}

}
\keyword{ distribution }
