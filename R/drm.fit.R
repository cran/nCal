if(getRversion() >= "2.15.1")  utils::globalVariables(c("drm.weights"))


# better ssfct to get better fits
# Use gof.threshold to report lack of fit
# this is only working for log transformed
# this is not working for weighting for two reasons 1) if weights are present, fails 2) gof needs to be computed differently
# robust="mean"; fit.4pl=FALSE; force.fit=F; weighting=F # default
drm.fit=function (formula, data, robust="mean", fit.4pl=FALSE, force.fit=FALSE, weighting=FALSE, pow.weight=1, verbose=FALSE) {
    
    require(drc)
    
    # drm 5pl is sensitive to the ordering of rows
    data=data[order(data$expected_conc),]
    
    gof.threshold=.2 # .2 is empirically determined
    control=drmc(maxIt=5000, method="BFGS", relTol=1e-7, trace=FALSE)
    
    if (weighting) {
        outcome.coln=all.vars(formula)[1]
        eval(eval(substitute(expression( drm.weights <<- data[,outcome.coln%+%".avg"]^pow.weight ))))     
    } else {
        eval(eval(substitute(expression( drm.weights <<- rep(1,nrow(data)) ))))     
    }

    gofs=rep(Inf, 3)
    fits=list()
    
    if (!fit.4pl){
        fit1= try(drm(formula=formula, data=data, robust=robust, fct = LL.5(ssfct=ss.fct.via.LL4), control=control, weights=drm.weights), silent=FALSE)    
        
        if (!inherits(fit1, "try-error")) {
            vcov.=vcov(fit1)
            if(is.null(vcov.)) {
                bad.se = T 
            } else if(all(is.na(diag(vcov.)))) {
                bad.se=T
            } else if (any(diag(vcov.)<0)) {
                bad.se=T
            } else {
                bad.se=F
            }            
            if (bad.se) gof1=Inf
            else gof1 = gof.cost( resid(fit1) ) 
        } else {
            gof1=Inf
        }
        gofs[1]=gof1
        fits[[1]]=fit1
        
        if (gof1>gof.threshold) {
            
            if (verbose) print("ss.fct.via.LL4 fails, try ssfct.drc.1.5.2")
            fit2= try(drm(formula=formula, data=data, robust=robust, fct = LL.5(ssfct=ssfct.drc.1.5.2), control=control, weights=drm.weights), silent=FALSE)
            
            if (!inherits(fit2, "try-error")) {
                vcov.=vcov(fit2)
                if(is.null(vcov.)) {
                    bad.se = T 
                } else if(all(is.na(diag(vcov.)))) {
                    bad.se=T
                } else if (any(diag(vcov.)<0)) {
                    bad.se=T
                } else {
                    bad.se=F
                }            
                if (bad.se) gof2=Inf
                else gof2 = gof.cost( resid(fit2) ) 
            } else gof2=Inf
            gofs[2]=gof2
            fits[[2]]=fit2
            
            if (gof2>gof.threshold) {
                
                if(verbose) print("ssfct.drc.1.5.2 fails, try default ssfct")
                fit3 = try(drm(formula=formula, data=data, robust=robust, fct = LL.5(), control=control, weights=drm.weights), silent=FALSE)
                
                if (!inherits(fit3, "try-error")) {
                    vcov.=vcov(fit3)
                    if(is.null(vcov.)) {
                        bad.se = T 
                    } else if(all(is.na(diag(vcov.)))) {
                        bad.se=T
                    } else if (any(diag(vcov.)<0)) {
                        bad.se=T
                    } else {
                        bad.se=F
                    }            
                    if (bad.se) gof3=Inf
                    else gof3 = gof.cost( resid(fit3) ) 
                } else gof3=Inf
                gofs[3]=gof3
                fits[[3]]=fit3
                if (gof3>gof.threshold) print ("residuals larger than ususal in "%+%data[1,"assay_id"]%+%", "%+%data[1,"analyte"])
                else {
                    
                }
                
            } else {
                
            }
        
        } 
        
        fit=fits[[which.min(gofs)]]
        
    } else {
        fit= try(drm(formula=formula, data=data, robust=robust, fct = LL.4(), control=control, weights=drm.weights), silent=FALSE)  
        if (!inherits(fit, "try-error")) {
            vcov.=vcov(fit)
            if(is.null(vcov.)) {
                bad.se = T 
            } else if(all(is.na(diag(vcov.)))) {
                bad.se=T
            } else if (any(diag(vcov.)<0)) {
                bad.se=T
            } else {
                bad.se=F
            }            
            if (bad.se) gof3=Inf
            else gof3 = gof.cost( resid(fit) ) 
        } else gof3=Inf
        gofs[1]=gof3
        
    }
    
    fit$bad.se=bad.se
    if (force.fit) {
        if (class(fit)=="try-error") return (NULL) 
        else {
            if (verbose) print("return forced fit")
            return (fit)
        }
    } else {
        if (min(gofs)>gof.threshold) {
            print ("residuals larger than ususal")
            return (NULL)
        }
        else return (fit)
    }
}
#fit.drc(fi ~ expected_conc, data = dat)
#fit.drc(MFI ~ Expected.Conc, data = tt, weights= tt$MFI.avgY^-1)


# self start functions
ssfct.drc.1.5.2 <- function (dataFra) {
    x <- dataFra[, 1]
    y <- dataFra[, 2]
    startVal <- rep(0, 5)
    startVal[3] <- max(y) + 0.001
    startVal[2] <- min(y) - 0.001
    startVal[5] <- 1
    indexT2 <- (x > 0)
    x2 <- x[indexT2]
    y2 <- y[indexT2]
    startVal[c(1, 4)] <- find.be2(x2, y2, startVal[2] - 0.001, 
        startVal[3])
#    print (startVal)
    return(startVal)
}
find.be2 <- function(x, y, c, d)
{
#    myprint (x," ")
#    myprint (y)
#    myprint (c)
#    myprint (d)
    logitTrans <- log((d - y)/(y - c))  
    
    lmFit <- lm(logitTrans ~ log(x))
#        eVal <- exp((-coef(logitFit)[1]/coef(logitFit)[2]))
#        bVal <- coef(logitFit)[2]
    
    coefVec <- coef(lmFit)
    bVal <- coefVec[2]        
    eVal <- exp(-coefVec[1]/bVal)    
    
    return(as.vector(c(bVal, eVal)))
}
ss.fct.via.LL4= function (dataFra) {
    x <- dataFra[, 1]
    y <- dataFra[, 2]
    fit0= drm(y ~ x, fct = LL.4(), weights= y^-.5)
    start=c(coef(fit0), 1)
    return (start)
}

# sometimes when y is MFI, the lower asymptote can be below 0, after taking log, we get NaN
gof.cost=function (x) {
    x=ifelse (is.na(x), 10, x)
    mean(abs(x))
} 

getVarComponent.drc=function (fit) {
    summary(fit)$resVar
}
