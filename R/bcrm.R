# n.iter=1e5; jags.seed=1; n.thin=NULL; keep.data=TRUE; keep.jags.samples=FALSE; t.unk.truth=NULL; params.true=NULL; T.init=NULL; prior.sensitivity=NULL; verbose=F #default
bcrm = function (formula, data, error.model, informative.prior, 
    prior.sensitivity=NULL, n.iter=1e5, jags.seed=1, n.thin=NULL, T.init=NULL, keep.data=TRUE, keep.jags.samples=FALSE, standards.only=TRUE, verbose=FALSE,
    # for simulation study only
    t.unk.truth=NULL, params.true=NULL)
{
    
    load.rjags=try(require(rjags))
    if (!load.rjags) stop("rjags does not load successfully. Is JAGS installed?")
    
    if (standards.only) data=data[data$well_role=="Standard",]
    
    dat.model.frame=model.frame(formula, data)
    
    if (verbose) print(error.model)
    
    if (is.null(data$replicate)) {
        if (error.model %in% c("gh_mixnorm","classical_mixnorm")) {
            stop("data is missing a replicate column")
        } else data$replicate=1 
    } 
    n.replicate=max(data$replicate)
    if (n.replicate>2) stop("Only two replicates are supported for now")
        
    # ordering the data helps later to tell which is which in the indicators
    data=data[order(data$replicate),]        
    data=data[order(dat.model.frame[[2]]),]
    data=data[order(data$assay_id),]
    assay_names=unique(data$assay_id)
    data=rbind(data[data$well_role=="Standard",], data[data$well_role=="Unknown",])
    
    # create seqno column, assuming all curves have the same dilution!
    data$seqno=as.numeric(as.factor(data$expected_conc))
    nDil=max(data$seqno,na.rm=TRUE)
    
    # create i.curve column to index curves, used by the error.model file
    data$i.curve=as.numeric(as.factor(data$assay_id))
    n.curve=max(data$i.curve)
    
    if(is.null(T.init)) T.init=array(0,dim=c(n.curve,n.replicate,nDil)) # T is now a 3-dimensional vector
    
    # make sure sample_id is NA for Standard and is 1:S for Unknown
    data$sample_id[data$well_role=="Standard"]=NA
    data$sample_id[!is.na(data$sample_id)]=as.numeric(as.factor(data$sample_id[!is.na(data$sample_id)]))
            
    # do drm fits, fitted parameter will be used to get initial values for samplers and set the hyperparameter for mean
    fits = attr(ncal(formula, data, return.fits=TRUE, force.fit=TRUE, plot.se.profile=FALSE, auto.layout=FALSE, plot=FALSE), "fits")
    params.drm = mysapply(fits, coef)
    if (startsWith(error.model,"gh_")) {
        theta.init=cla2gh(params.drm)
        theta.init[,4]=log(theta.init[,4])
        theta.init[,5]=log(theta.init[,5])
        var.matrix=var.matrix.gh
        R=R.gh
        mean.distr=mean.distr.gh
        #mean.distr["mean",]=colMeans(theta.init)
        #mean.distr["prec","logf"]=2
    } else if (startsWith(error.model,"classical_")) {
        colnames(params.drm) = substr(colnames(params.drm), 1, 1)
        theta.init=cbind(params.drm[,"c"], params.drm[,"d"], log(params.drm[,"e"]), log(-params.drm[,"b"]), log(params.drm[,"f"]))
        var.matrix=var.matrix.classical
        R=R.classical
        mean.distr=mean.distr.classical
    } else stop("error.model not right: "%+%error.model)
    sd.resid=(sqrt(sapply(fits, function (fit) summary(fit)$resVar)))
    if (verbose) print("bcrm 100")
    if (verbose) print(theta.init)
    
    if (!is.null(prior.sensitivity)) {
        if (prior.sensitivity==1){
            scale.f=2
            R["g","g"]=R["g","g"]*scale.f**2
            R["logh","logh"]=R["logh","logh"]*scale.f**2
        } else if (prior.sensitivity==2){
            scale.f=4
            R["g","g"]=R["g","g"]*scale.f**2
            R["logh","logh"]=R["logh","logh"]*scale.f**2
        } else if (prior.sensitivity==3){
            scale.f=2
            R["c","c"]=R["c","c"]*scale.f**2
            R["d","d"]=R["d","d"]*scale.f**2
            R["logf","logf"]=R["logf","logf"]*scale.f**2
        } else if (prior.sensitivity==4){
            scale.f=4
            R["c","c"]=R["c","c"]*scale.f**2
            R["d","d"]=R["d","d"]*scale.f**2
            R["logf","logf"]=R["logf","logf"]*scale.f**2
        } else {
            stop("prior.sensitivity not supported: "%+%prior.sensitivity)
        }
    }
    
    
    # the meaning of tau in t4 and norm is different    
    if (endsWith(error.model,"t4")) {
        #var.comp=c(2,0.002) # used in the paper, a flat prior on tau, also favors large tao, small variance
        var.comp=c(2,0.02) # "calibrated" to give similar results as drm estimate of the variance component
    } else {
        #var.comp = c(2, 0.02) # used in the paper
        var.comp = c(2, 0.06) # "calibrated" to give similar results as drm estimate of the variance component
    }
    var.comp.2 = c(2, 0.2) # second component 
    
#    nn is for normal mixture and n is for normal. the following two used to use the same error.model file as error.model_classical_t4.txt
#    } else if (error.model=="nn") {
#        mixp=0.05
#        var.comp=c(2,0.02) # empirical Bayes estimate
#        var.comp.2=c(2,0.2) # the second component is expected to have a much smaller preicsion and larger variance
#        dof=c(Inf,Inf)
#    } else if (error.model=="n") {
#        mixp=0
#        dof=c(Inf,Inf)
#        var.comp=c(2,0.02) # empirical Bayes estimate
#        var.comp.2=var.comp # doesn't matter
    
    
#    wd=3 # weak
#    wd=10 # strong
    if(informative.prior) {
        dof.wish=8; R=R
        v=mean.distr[1,]; m=diag(mean.distr[2,])
    } else {
        dof.wish=5; R=diag(1,5)
        v=rep(0,5); m=diag(1e-4, 5)
    }
    
    
    ###############
    # jags.data
    
    if (verbose) print("bcrm 200")
    jags.data = list(
        "t"=log(dat.model.frame[[2]]), "y"=dat.model.frame[[1]], "i.curve"=data$i.curve, "I"=n.curve, "K"=sum(data$well_role=="Standard"), "nDil"=nDil, "seqno"=data$seqno, "var.comp"=var.comp,
        "dof.wish"=dof.wish, "R"= R,         
        "v"=v, "m"=m, 
#        dof.wish.1=wd + 4, "R1"= diag(diag(var.matrix))*(wd-2), 
        "precision.matrix"=diag(1/diag(var.matrix)),
        "var.comp.2"=var.comp.2,        
        "t.min"=min(log(data$expected_conc), na.rm=TRUE), "t.max"=max(log(data$expected_conc), na.rm=TRUE),
        "replicate"=data$replicate
    )
    if (error.model=="classical_replicate_re" | error.model=="gh_replicate_re") {
        jags.data=c(jags.data, list("dof.wish.0"=5, "R0"= diag(1,5)))
    }
    #if (!standards.only) 
    jags.data=c(jags.data, list("J"=sum(data$well_role=="Unknown"), "sample_id"=data$sample_id, "S"=length(setdiff(unique(data$sample_id),NA))))
    #print(jags.data)
    
    ###############
    # jags.init    
    # Note that after JAGS 3.3, if a node is a constant node, it cannot have an init value
    
    jags.inits = list( 
        "tau1"=100, 
        "r"=8,
        #"t.unk.sample"=mean(min(log(data$expected_conc), na.rm=TRUE)),
        "theta"=theta.init, "theta.0"=colMeans(theta.init), "TAO"=diag(1/diag(var.matrix)),
        .RNG.name="base::Mersenne-Twister", .RNG.seed=jags.seed +11
    ) 
    if (!error.model %in% c("classical_norm","classical_t4","classical_replicate_re","classical_tar1","gh_norm","gh_replicate_re","gh_tar1","gh_t4")) {
        jags.inits=c(jags.inits, list("T"=T.init)) # c(jags.init, "T"=T.init) won't work
    }
    if (!error.model %in% c("classical_t4","gh_t4","classical_mixnorm","gh_mixnorm","gh_replicate_re")) {
        jags.inits=c(jags.inits, "alpha"=0.5)    
    }
    if (error.model=="classical_replicate_re" | error.model=="gh_replicate_re"){
        theta.r.init=array(0,dim=c(n.curve,2,5))
        jags.inits=c(jags.inits, list("theta.r"=theta.r.init, "W"=diag(10/diag(var.matrix)) ))
    }   
    if (!error.model %in% c("classical_t4","gh_t4","classical_replicate_re","gh_replicate_re")) jags.inits=c(jags.inits, "mixp"=.05) 
    if (verbose) print(names(jags.inits))
    
    ################
    # trace.var
    
    pname=c("c","d","f","g","h") # these have to be alphabetically ordered, so as to match the order of columns in jags samples
    if (endsWith(error.model,"replicate_re")) {
        trace.var="p"%+% pname
    } else trace.var=pname
    trace.var =c(trace.var,"sigma")
    if (error.model=="classical_t4" | error.model=="gh_t4") {
        # do nothing
    } else if(!contain(error.model,"tar1")){
        trace.var =c(trace.var,"T","mixp","alpha")
    } else {
        trace.var =c(trace.var,"alpha")
    }        
    if (!standards.only & length(setdiff(unique(data$sample_id),NA))>0) trace.var=c(trace.var, "t.unk.sample", "T.unk")
    
    
    ################
    # run jags
    
    cat("  Running jags\n")
    # try three times
    jags.success=F
    for (i.try in 1:3) {
        jags.model.1 = try(suppressWarnings(
            jags.model(file=system.file(package="nCal")[1]%+%"/jags_script/model_"%+%error.model%+%".txt", data=jags.data, inits=jags.inits, n.chains = 1, n.adapt=1e3, quiet=TRUE)
         ), silent=FALSE)
        if (inherits(jags.model.1, "try-error")) {
            print("jags.model fails, try with a different seed")            
            jags.inits$.RNG.seed = jags.inits$.RNG.seed+1
        } else {
            jags.success=T
            break
        }
    }
    if (!jags.success) return (NULL)
    
    if (is.null(n.thin)) {
        # choose a n.thin so that  samples are saved
        n.target = 5e3
        if (n.iter<n.target) n.thin=1 else n.thin = floor(n.iter/n.target)
    }
    n.burnin=n.iter/n.thin/5 # 20% burnin
    #n.burnin=n.iter/n.thin/2 # 50% burnin
    samples = coda.samples(jags.model.1, trace.var, n.iter=n.iter, thin = n.thin)[[1]][-(1:n.burnin),,drop=FALSE]
    
    fit=list()
    attr(fit, "class")="bcrm"         
    if (keep.jags.samples) fit$jags.samples=samples
    
    if("alpha" %in% colnames(samples)) fit$alpha=mean(samples[,"alpha"])
    
    # summarize posterior distributions for parameters: median, sd, 95% CI
    if (!endsWith(error.model,"replicate_re")) {
        samples.2=samples[, regexpr("^[cdghf]\\[.+\\]", colnames(samples))!=-1, drop=F ] # select columns
        if (ncol(samples.2)==0) samples.2=samples[, regexpr("^[cdghf]", colnames(samples))!=-1,drop=F ] # if there is only one curve, the column names are b,c,d,e,f
    } else {
        samples.2=samples[, regexpr("^p[cdghf]\\[.+\\]", colnames(samples))!=-1, drop=F ] # select columns
        if (ncol(samples.2)==0) samples.2=samples[, regexpr("^p[cdghf]", colnames(samples))!=-1,drop=F ] # if there is only one curve, the column names are b,c,d,e,f
        colnames(samples.2)=substr(colnames(samples.2),2,1000)
    }
#    print(str(samples.2))
#    print(pname)
    fit$median.coef =matrix(apply(samples.2, 2, median),                        nrow=n.curve, dimnames=list(assay_names, pname)) # median is better than mean here for abc
    fit$mean.coef=matrix(apply(samples.2, 2, mean),                        nrow=n.curve, dimnames=list(assay_names, pname)) 
    fit$mode.coef=matrix(apply(samples.2, 2, function(x) {den<-density(x); den$x[which(den$y==max(den$y))]}),                        nrow=n.curve, dimnames=list(assay_names, pname)) 
    fit$sd.coef=     matrix(apply(samples.2, 2, sd),                            nrow=n.curve, dimnames=list(assay_names, pname))
    fit$low.coef=    matrix(apply(samples.2, 2, function(x) quantile(x,0.025)), nrow=n.curve, dimnames=list(assay_names, pname))
    fit$high.coef=   matrix(apply(samples.2, 2, function(x) quantile(x,0.975)), nrow=n.curve, dimnames=list(assay_names, pname))
    fit$coefficients=fit$median.coef
    fit$coef.samples=samples.2
    
#    if (!is.null(params.true)) {
#        for (k in 1:n.curve) {
#            f0=FivePL.t.func(params.true[k,])
#            apply(samples.2, 1, function(x) {
#                f1=FivePL.t.func( x )
#                integrate( function(t) abs(f1(t)-f0(t)) , lower=min(log.conc), upper=max(log.conc), subdivisions=1000 )$value
#            })
#        }
#        
#    }
#    fit$abc=
    
    fit$fitted = sapply(1:nrow(data[data$well_role=="Standard",]), function(i.row) {
        FivePL.x(data$expected_conc[i.row], fit$coefficients[data[i.row,"assay_id"],])
    })
    fit$resid=log(data$fi[data$well_role=="Standard"])-fit$fitted
    
    # summarize posterior distributions for the variance components, there may be one or two
    samples.4=samples[,startsWith(colnames(samples),"sigma"),drop=FALSE ] 
    fit$varcomp=apply(samples.4, 2, function(x) median(x))
    
    if(!contain(error.model,"tar1") & error.model!="classical_t4" & error.model!="gh_t4"){
        fit$mixp=mean(samples[,"mixp"])
    
        # summarize mixture indicators
        samples.3=samples[,startsWith(colnames(samples),"T[") ] # select columns
        tmp=apply(samples.3, 2, mean) # cannot do median here, because it is Bernoulli variable
        #the following line makes the assumption that there are the same number of points each curve, this is not true sometimes
        #fit$mixture.indicators=matrix(tmp[1:sum(data$well_role=="Standard")], ncol=n.curve, dimnames=list(1:(sum(data$well_role=="Standard")/n.curve), assay_names))        
        # simple give it to a vector, should have the same length as the number of rows in fit$data
        #if (substr(error.model,1,2)=="nn") 
        fit$mixture.indicators=tmp["T["%+%data$i.curve[data$well_role=="Standard"]%+%","%+%data$replicate[data$well_role=="Standard"]%+%","%+%data$seqno[data$well_role=="Standard"]%+%"]"]
    }
        
    # summarize Unknown concentrations
    samples.6=samples[,startsWith(colnames(samples),"T.unk[") ] # select columns
    tmp=apply(samples.6, 2, mean) # cannot do median here, because it is Bernoulli variable
    unk = data.frame(data[data$well_role=="Unknown",], "mix.ind"=tmp)
    samples.5 = samples[,startsWith(colnames(samples),"t.unk.sample"),drop=FALSE ] 
    if(ncol(samples.5)>0){
        fit$t.unk.mean=  apply(samples.5, 2, function(x) mean(x)) # mean is better than median for unk mse
        fit$t.unk.median=apply(samples.5, 2, function(x) median(x)) 
        fit$t.unk.mode= apply(samples.5, 2, function(x) {den<-density(x); den$x[which(den$y==max(den$y))]} )
        if (!is.null(t.unk.truth)) {    
            fit$t.unk.cp = sapply(1:ncol(samples.5), function(i) {
                x=samples.5[,i]        
                quantile(x,0.025) < t.unk.truth[i] & t.unk.truth[i] < quantile(x,0.975)
            })
            names(fit$t.unk.cp)="sample"%+%1:ncol(samples.5)
            fit$t.unk.mse = sapply(1:ncol(samples.5), function(i) {
                x=samples.5[,i]        
                mean ((x-t.unk.truth[i])**2)
            })
            names(fit$t.unk.mse)="sample"%+%1:ncol(samples.5)
            fit$t.unk.perc.bias = sapply(1:ncol(samples.5), function(i) {
                x=samples.5[,i]        
                mean ((x-t.unk.truth[i])/t.unk.truth[i]*100)
            })
            names(fit$t.unk.mse)="sample"%+%1:ncol(samples.5)
            fit$t.unk.var = sapply(1:ncol(samples.5), function(i) {
                x=samples.5[,i]        
                var (x)
            })
            names(fit$t.unk.mse)="sample"%+%1:ncol(samples.5)
       }
    }
        
    if (keep.data) {
        if(!contain(error.model,"tar1") & error.model!="classical_t4" & error.model!="gh_t4"){
            fit$data=cbind(data, fitted=fit$fitted, resid=fit$resid, mixture.indicators=fit$mixture.indicators)
        } else {
            fit$data=cbind(data, fitted=fit$fitted, resid=fit$resid)
        }
    }
    fit$error.model=error.model
    fit$bad.se=FALSE
    return (fit)        
}



# if assay_id is not null, then only one curve is plotted, e.g. "LMX004-L-RV144"
# when x2 is not null, then the second x is used to plot another line
# only works with two replicates for now
# log="x" means concentration is plotted, otherwise, log concentration is plotted
# assay_id=NULL; add=F; lcol=1; main=NULL; x2=NULL; lwd=1; points.only=FALSE; all.lines.only=FALSE; t=NULL; ylim=NULL
plot.bcrm=function(x, assay_id=NULL, add=FALSE, lcol=1, main=NULL, x2=NULL, lwd=.1, points.only=FALSE, all.lines.only=FALSE, t=NULL, ylim=NULL, same.ylim=FALSE, 
    lty3=NULL,x3=NULL,lcol2=NULL,lcol3=NULL,xlab=NULL,ylab=NULL,col.outliers=TRUE,lty=1,cex=1,log="x",...){
    
    dat=x$data
    dat=dat[dat$well_role=="Standard",]
    if (!is.null(assay_id)) assay_names=assay_id else assay_names=unique(dat$assay_id)
    
    # add mixture.indicators to the data frame
    if (!is.null(x$mixture.indicators) ) {
        if (!any(is.na(x$mixture.indicators))) {
            dat$mixture.indicators=x$mixture.indicators
        }
    } 
    
    if (all.lines.only) {
        t=log(dat$expected_conc)[dat$assay_id==assay_names[1]]
        y=log(dat$fi)[dat$assay_id==assay_names[1]]
        if(log!="x") {
            plot(t,y,main=ifelse(is.null(main),"",main),type="n",ylim=ylim,xlab=xlab,ylab=ylab)
        } else {
            plot(exp(t),y,main=ifelse(is.null(main),"",main),type="n",ylim=ylim,xlab=xlab,ylab=ylab,log="x")
        }
    }
    
    if (same.ylim) ylim=range(log(dat$fi))
    for(a in assay_names) {
        dat.a=dat[dat$assay_id==a,]
        # order points
        dat.a=dat.a[order(dat.a$expected_conc),]
        
        t=log(dat.a$expected_conc)
        
        # plot points
        if (!add & !all.lines.only) {
            y=log(dat.a$fi)
            col=1
            if(col.outliers)
                if (!is.null(x$mixture.indicators) ) {
                    if (!any(is.na(x$mixture.indicators))) {
                        col=ifelse(dat.a$mixture.indicators>.5,2,1)
                    }
                } 
            if (!is.null(dat.a$replicate)) pch=ifelse(dat.a$replicate==1, 1, 19) else pch=1
            if(log!="x") {
                plot(t,y,main=ifelse(is.null(main),a,main),col=col, cex=cex, pch=pch,ylim=ylim,xlab=xlab,ylab=ylab)
            } else {
                plot(exp(t),y,main=ifelse(is.null(main),a,main),col=col, cex=cex, pch=pch,ylim=ylim,xlab=xlab,ylab=ylab,log="x")
            }
        }
                
        # plot lines
        if (!points.only & !is.null(x$coefficients)) {
            t.1=seq(min(t), max(t), length=100)
            x.1=t.1
            if(log=="x") x.1=exp(t.1)
            lines(x.1, FivePL.t(t.1, x$coefficients[a,]), lty=lty, col=lcol, lwd=lwd)
            if (!is.null(x2)) {
                lines(x.1, FivePL.t(t.1, x2$coefficients[a,]), lty=lty, col=ifelse(is.null(lcol2),2,lcol2), lwd=lwd)
            }
            if (!is.null(x3)) {
                lines(x.1, FivePL.t(t.1, x3$coefficients[a,]), lty=ifelse(is.null(lty3),lty,lty3), col=ifelse(is.null(lcol2),3,lcol3), lwd=lwd)
            }
        }
    }
}

print.bcrm=function (x, ...) {
    print(x[c("median.coef")])
    cat("Complete list of fields: \n")
    print(names(x))
}

coef.bcrm=function(object, type="gh", ...) {
    if(type=="gh") {
        object$median.coef[1,c("c","d","g","h","f")]
    } else if (type=="classical") {
        gh2cla(object$median.coef[1,])
    } else stop("type not supported")    
}

# note that vcov for object is not very meaningful b/c the distribution is far from multivariate normal
vcov.bcrm=function(object, type="gh", ...) {
    if (is.null(object$vcov)) {
        if(type=="gh") {
            object$vcov=cov(object$coef.samples[,c("c","d","g","h","f")])
        } else if (type=="classical") {
            object$vcov=cov(gh2cla(object$coef.samples))
        } else stop("type not supported")
    }
    object$vcov
}

getVarComponent.bcrm=function (fit) {
    s=fit$varcomp[1]
    if (endsWith(fit$error.model, "t4")) {
        s^2 * 2 # 2 is the var of standard Student's t
    } else {
        s^2
    }
}

# y is the left hand side of the formula
predict.bcrm=function (fit, y){
    xx=apply(fit$coef.samples,1, function (param) FivePL.x.inv(y, param, incr=TRUE) )
    median(xx)    
}

get.single.fit=function(fit, assay_id) {
    assay_names=rownames(fit$median.coef)
    single.fit=fit
    single.fit$median.coef=fit$median.coef[assay_id,,drop=F]    
    id=match(assay_id,assay_names)
    single.fit$coef.samples=fit$coef.samples[,id+length(assay_names)*0:4]
    dimnames(single.fit$coef.samples)[[2]]=substr(dimnames(single.fit$coef.samples)[[2]],1,1)
    #single.fit$varcomp=fit$varcomp[id+length(assay_names)*0:1]
    single.fit$data=single.fit$data[single.fit$data$assay_id==assay_id,]
    single.fit$mixture.indicators=NULL
    
    single.fit
} 
