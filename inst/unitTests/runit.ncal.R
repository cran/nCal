### --- Test setup ---

if(FALSE) {
  library("RUnit")
  library("nCal")
}

test.ncal <- function() {

RNGkind("Mersenne-Twister", "Inversion")
#RNGkind("Marsaglia-Multicarry", "Kinderman-Ramage") 


# a dataset without analyte and assay_id column
set.seed(1)
#    print(runif(1))
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
fi=simulate1curve (p.eotaxin[1,], rep(log.conc,each=n.replicate), sd.e=0.2)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)))

checkException (
    ncal(log(fi)~expected_conc, dat.std, return.fits = TRUE)
, msg="check no analyte or assay_id column")

# add analyte and assay_id column
dat.std$analyte="Test"
dat.std$assay_id="Run 1"

checkEquals (
    nrow(ncal(log(fi)~expected_conc, dat.std, return.fits = TRUE))
, 0, msg="check no analyte or assay_id column")

# add unknown
dat.unk=rbind(
  data.frame(fi=exp(6.75), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=1)
, data.frame(fi=exp(11), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=2)
, data.frame(fi=exp(3),    expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=3)
, data.frame(fi=exp(4.4),  expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=4)
)
dat.std=cbind(dat.std, sample_id=NA)
dat=rbind(dat.std, dat.unk)

checkException (
    ncal(log(fi)~expected_conc, dat, return.fits = TRUE)
, msg="check no well_role")

# add well_role
dat.std=cbind(dat.std, well_role="Standard")
dat.unk=cbind(dat.unk, well_role="Unknown")
dat=rbind(dat.std, dat.unk)

out=ncal(log(fi)~expected_conc, dat, return.fits = TRUE, additional.plot.func=function() abline(v=10))
checkEqualsNumeric(unlist(out[1,c("est.log.conc","se")]), c(3.940014, 0.1622546), tolerance=1e-6)

out.norm=ncal(log(fi)~expected_conc, dat, return.fits = TRUE, bcrm.fit=TRUE, bcrm.model="norm", control.jags=list(n.iter=10, n.adapt=0))

checkEqualsNumeric(
    unlist(out.norm[1,c("est.log.conc","se")])
    , 
    c(3.9400136, 0.1253771)
, tolerance=1e-6)

# weighting
out.w = ncal(fi~expected_conc, dat, return.fits = TRUE, plot.se.profile=TRUE, weighting=TRUE, pow.weight=-1)
fit.w=attr(out.w, "fits")[[1]]

checkEqualsNumeric(
    coef(fit.w)
    , 
    c(-1.213899,    76.742779, 31382.356851,   764.912405,     1.129558 )
, tolerance=1e-6)



# test decreasing curve
dat.2 = dat
dat.2$expected_conc=1/dat.2$expected_conc

out.2=ncal(log(fi)~expected_conc, dat.2, return.fits = TRUE)

checkTrue(
    coef(attr(out.2, "fits")[[1]])["b"]>0
)

checkEqualsNumeric(
    unlist(out.2[1:3,c("est.log.conc","se")])
    , 
    c(-3.939381, -9.9034876, 0.6771702, 0.162688, Inf, Inf)
, tolerance=1e-6)


out.2.norm=ncal(log(fi)~expected_conc, dat.2, return.fits = TRUE, bcrm.fit=TRUE, control.jags=list(n.iter=1e1, n.adapt=0), bcrm.model="norm", verbose=FALSE)

checkEqualsNumeric(
    coef(attr(out.2.norm, "fits"))
    ,
    c(4.386389, 10.511270, -4.145641, -1.322411,  2.547013)
, tolerance=1e-6)

checkEqualsNumeric(
    unlist(out.2.norm[1:3,c("est.log.conc","se")])
    , 
    c(-3.9400136,    -9.9034876,     0.6771702,     0.1253771, Inf, Inf)
, tolerance=1e-6)


# test 4PL
out.4pl=ncal(log(fi)~expected_conc, dat, return.fits = TRUE, fit.4pl=TRUE)
    
checkTrue(
    is.na(coef(attr(out.4pl, "fits")[[1]])["f"])
)

checkEqualsNumeric(
    unlist(out.4pl[1:3,c("est.log.conc","se")])
    , 
    c(3.9829620, 9.2103404, -1.3703174, 0.1629386, Inf, Inf)
, tolerance=1e-6)


out.4pl.norm=ncal(log(fi)~expected_conc, dat, return.fits = TRUE, fit.4pl=TRUE, bcrm.fit=TRUE, bcrm.model="norm", control.jags=list(n.iter=10, n.adapt=0))

checkEqualsNumeric(
    coef(attr(out.4pl.norm, "fits"))
    , 
    c(4.386389, 10.511270,  4.145641,  1.322411)
, tolerance=1e-6)

out.2.4pl.norm=ncal(log(fi)~expected_conc, dat.2, return.fits = TRUE, fit.4pl=TRUE, bcrm.fit=TRUE, bcrm.model="norm", control.jags=list(n.iter=10, n.adapt=0))

checkEqualsNumeric(
    coef(attr(out.2.4pl.norm, "fits"))
    , 
    c(4.386389, 10.511270, -4.145641, -1.322411)
, tolerance=1e-6)


out.2.4pl=ncal(log(fi)~expected_conc, dat.2, return.fits = TRUE, fit.4pl=TRUE)

checkEqualsNumeric(
    coef(attr(out.2.4pl, "fits")[[1]])
    , 
    c(0.86770147,  4.25401765, 10.38578190,  0.01207711)
, tolerance=1e-6)


out.4pl.t4=ncal(log(fi)~expected_conc, dat, return.fits = TRUE, fit.4pl=TRUE, bcrm.fit=TRUE, bcrm.model="t4", control.jags=list(n.iter=10, n.adapt=0))

checkEqualsNumeric(
    coef(attr(out.4pl.t4, "fits"))
    , 
    c(4.386389, 10.511270,  4.145641,  1.322411)
, tolerance=1e-6)



}
