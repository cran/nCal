### --- Test setup ---

if(FALSE) {
  library("RUnit")
  library("nCal")
}

test.bcrm <- function() {

tolerance.jags=1e-1 # JAGS is not yet reproducible, see http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/6c8c3e6a/

# more stringent tolerance for one system to ensure algorithm accuracy
if (R.Version()$system %in% c("x86_64, mingw32")) {
    tolerance.jags=1e-6
}
 
RNGkind("Mersenne-Twister", "Inversion")

# decreasing curves
set.seed(1)
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
p.1=p.eotaxin[1,]
p.1["b"]= -p.1["b"]
fi=simulate1curve (p.1, rep(log.conc,each=n.replicate), sd.e=0.3)
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
p.2=p.eotaxin[2,]
p.2["b"]= -p.2["b"]
fi=simulate1curve (p.2, rep(log.conc,each=n.replicate), sd.e=0.3)
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

checkException(
    bcrm(log(fi)~expected_conc, dat, error.model="classical_norm", informative.prior=TRUE, n.iter=1e1)
)

fits = bcrm(log(fi)~expected_conc, dat, error.model="gh_norm", informative.prior=TRUE, n.iter=1e1, n.adapt=0)

checkEqualsNumeric(
    mean(fits$median.coef["assay1",])
    , 5.471837
, tolerance=tolerance.jags)




}
