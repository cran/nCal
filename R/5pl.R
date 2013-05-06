# FIX incr


get.curve.param.list=function(param){
    tmp=substr(names(param),1,1)
    names(param)=tmp
    if(!"b" %in% tmp & !"e" %in% tmp) {
        # not classical parameterization
        if("g" %in% tmp) param=gh2cla(param) else stop("not gh not cla")
    } 
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]; f=param["f"]    
    names(b)=NULL; names(c)=NULL; names(d)=NULL; names(e)=NULL; names(f)=NULL; # so that the outcome is not named incorrectly    
    list(b=b, c=c, d=d, e=e, f=f)
}


############################################################################
# classical parameterization

# following functions are vectorized for t, but not for param
FivePL.t=function (t,param) {
    tmp=sapply(names(param),function(x) strsplit(x,"\\[")[[1]][1] )
    names(param)=tmp
    names(param)=substr(names(param),1,1)
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]; f=param["f"]; g=param["g"]; h=param["h"]; logtao=param["logtao"]
    
    if(!"b" %in% tmp & "logmb"%in%tmp) b=unname(-exp(param["logmb"]))
    if(!"e" %in% tmp & "loge"%in%tmp) e=unname(exp(param["loge"]))
    if(!"f" %in% tmp & "logf"%in%tmp) f=unname(exp(param["logf"]))
    if(!"h" %in% tmp & "logh"%in%tmp) h=unname(exp(param["logh"]))
        
    if(is.na(b) | is.na(e)) {
        if(!is.na(g) & !is.na(h)) {
            pp=gh2cla(c(c,d,g,h=h,f=f)) 
            b=pp["b"]
            e=pp["e"]
        } else if (!is.na(logtao)) {
            pp=ed502cla(c(c,d,b=b,logtao=unname(logtao),f=f))            
            e=pp["e"]
        } else  {
            stop("parameterization not recognized")
        }
    } 
    
    out=(d-c)/{1+exp(b*t-b*log(e))}^f+c
    names(out)=rep("y", length(out))
    out
}

FivePL.t.func=function (param) {
    tmp=substr(names(param),1,1)
    names(param)=tmp
    if(!"b" %in% tmp & !"e" %in% tmp) {
        # not classical parameterization
        if("g" %in% tmp) param=gh2cla(param) else stop("not gh not cla")
    } 
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]; f=param["f"]    
    names(b)=NULL; names(c)=NULL; names(d)=NULL; names(e)=NULL; names(f)=NULL; # so that the outcome is not named incorrectly
    return (function (t) (d-c)/{1+exp(b*t-b*log(e))}^f+c)
}

FivePL.x=function (x,param) {
    tmp=substr(names(param),1,1)
    names(param)=tmp
    if(!"b" %in% tmp & !"e" %in% tmp) {
        # not classical parameterization
        if("g" %in% tmp) param=gh2cla(param) else stop("not gh not cla")
    } 
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]; f=param["f"]    
    out=(d-c)/{1+(x/e)^b}^f+c
    names(out)=rep("y", length(out))
    out
}

FivePL.x.inv = function (y,param, incr) {
    if (!is.matrix(param)) param=matrix(param, nrow=1, dimnames=list(NULL, substr(names(param),1,1)))
    flag=ifelse(incr,1,-1)
    if(!"b" %in% colnames(param) & !"e" %in% colnames(param)) {
        if("g" %in% colnames(param)) param=gh2cla(param) else stop("not gh not cla")
    } 
    b=param[,"b"]; c=param[,"c"]; d=param[,"d"]; e=param[,"e"]; f=param[,"f"]    
    out=(((d-c)/(y-c))^(1/f)-1)^(1/b)*e
    out = ifelse(y>c & y<d, out, {
        ifelse(y<c, -Inf, Inf)*flag
    })
    names(out)=rep("x", length(out))
    out        
}

FivePL.x.inv.func=function (param) {
    tmp=substr(names(param),1,1)
    names(param)=tmp
    if(!"b" %in% tmp & !"e" %in% tmp) {
        # not classical parameterization
        if("g" %in% tmp) param=gh2cla(param) else stop("not gh not cla")
    } 
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]; f=param["f"]    
    names(b)=NULL; names(c)=NULL; names(d)=NULL; names(e)=NULL; names(f)=NULL; # so that the outcome is not named incorrectly
    return (function (y)  (((d-c)/(y-c))^(1/f)-1)^(1/b)*e  )
}

FivePL.t.inv = function (y,param, incr) {
    log(FivePL.x.inv(y, param, incr))        
}

FivePL.t.inv.func=function (param) {
    param.list=get.curve.param.list(param)
    with(param.list, function (y) log((((d-c)/(y-c))^(1/f)-1)^(1/b)*e))
}



FourPL.x.inv <- function (y, param, incr=TRUE)  ###CHANGE 9 6/16/2011
{
    flag=ifelse(incr,1,-1)
    tmp = substr(names(param), 1, 1)
    names(param) = tmp
    if (!"b" %in% tmp & !"e" %in% tmp) {
        if ("g" %in% tmp)
            param = gh2cla(param)
        else stop("not gh not cla")
    }
    b = param["b"]
    c = param["c"]
    d = param["d"]
    e = param["e"]
    out = (((d - c)/(y - c)) - 1)^(1/b) * e
    out = ifelse(y > c & y < d, out, {
        ifelse(y < c, -Inf, Inf) * flag
    })
    names(out) = rep("x", length(out))
    out
}

FourPL.x <- function (x, param) ###Change 10 6/16/2011
{
    tmp = substr(names(param), 1, 1)
    names(param) = tmp
    if (!"b" %in% tmp & !"e" %in% tmp) {
        if ("g" %in% tmp)
            param = gh2cla(param)
        else stop("not gh not cla")
    }
    b = param["b"]
    c = param["c"]
    d = param["d"]
    e = param["e"]
    out = (d - c)/{
        1 + (x/e)^b
    } + c
    names(out) = rep("y", length(out))
    out
}

FourPL.t.func=function (param) {
    tmp=substr(names(param),1,1)
    names(param)=tmp
    if(!"b" %in% tmp & !"e" %in% tmp) {
        # not classical parameterization
        if("g" %in% tmp) param=gh2cla(param) else stop("not gh not cla")
    } 
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]    
    names(b)=NULL; names(c)=NULL; names(d)=NULL; names(e)=NULL; # so that the outcome is not named incorrectly
    return (function (t) (d-c)/{1+exp(b*t-b*log(e))}+c)
}







# simulate one curve, return FI, a vectorized function
simulate1curve=function(param, t, sd.e=0.1) {
    .mean=FivePL.t(t, param)
    y = rnorm (n=length(.mean), mean=.mean, sd=sd.e)
    exp(y)
}

# x can be a single number or a vector
# return a matrix, where each row corresponds to an x
vpl1.deriv = function(x,param){
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]; f=param["f"]    
    u=(x/e)^{b}
    y=c+((d-c)/((1+u)^{f}))
    res=matrix(
        c(  "c" = 1-(1/((1+u)^{f})), 
            "d" = (1/((1+u)^{f})), 
            "e" = (y-c)*((b*f)/e)*(u/(1+u)),
            "loge" = (y-c)*(b*f)*(u/(1+u)),
            "f" = -(y-c)*log(1+u),
            "b" = -(y-c)*(f/b)*(u/(1+u))*log(u)
        ), 
        nrow=length(x)
    )
    colnames(res)=c("c","d","e","loge","f","b")
    res
}
## test
#vpl1.deriv(c(1,100), coef(fit))

# return a list of derivative functions
vpl1.deriv.func = function(param){
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]; f=param["f"]    
    f.c=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname(1-(1/((1+u)^{f})))
    }
    f.d=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname((1/((1+u)^{f})))
    }
    f.e=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname((y-c)*((b*f)/e)*(u/(1+u)))
    }
    f.loge=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname((y-c)*b*f*(u/(1+u)))
    }
    f.f=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname(-(y-c)*log(1+u))
    }
    f.b=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname(-(y-c)*(f/b)*(u/(1+u))*log(u))
    }
    list("c"=f.c, "d"=f.d, "e"=f.e, "loge"=f.loge, "f"=f.f, "b"=f.b)
}


############################################################################
# parameterization 2: b,c,d,e,f -> b,c,d,tao,f, where tao is the ED50

# x can be a single number or a vector
# return a matrix, where each row corresponds to an x
vpl2.deriv = function(x,param){
    c=param["c"]; d=param["d"]; logtao=param["logtao"]; b=param["b"]; f=param["f"]    
    t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
    matrix(
        c(  "c" = 1-(1/((1+u)^{f})), 
            "d" = (1/((1+u)^{f})), 
            "logtao" = b*f*(y-c)*(u/(1+u)),
            "b" = -f*(t-logtao)*(y-c)*(u/(1+u)),
            "f" = ((log(2))/f)*((2^{1/f})/(2^{1/f}-1))*(y-c)*(u/(1+u))
        ), 
        nrow=length(x)
    )
}
## test
#vpl2.deriv(c(1,100), theta)

# return a list of derivative functions
vpl2.deriv.func = function(param){
    c=param["c"]; d=param["d"]; logtao=param["logtao"]; b=param["b"]; f=param["f"]    
    f.c=function(x) {
        t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
        unname(1-(1/((1+u)^{f})))
    }
    f.d=function(x) {
        t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
        unname((1/((1+u)^{f})))
    }
    f.logtao=function(x) {
        t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
        unname(b*f*(y-c)*(u/(1+u)))
    }
    f.b=function(x) {
        t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
        unname(-f*(t-logtao)*(y-c)*(u/(1+u)))
    }
    f.f=function(x) {
        t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
        unname(((log(2))/f)*((2^{1/f})/(2^{1/f}-1))*(y-c)*(u/(1+u)))
    }
    list("c"=f.c, "d"=f.d, "logtao"=f.logtao, "b"=f.b, "f"=f.f)
}


############################################################################
# parameterization 3: b,c,d,e,f -> b,c,d,g,h, where g is the inflection point, and h is the hill slope at the inflection point

# vectorized functions for converting parameters between parameterizations
# param can be a vector or a matrix where each row is a parameter value
cla2gh=function(param){
    is.v=FALSE
    if(is.vector(param)) {
        is.v=TRUE
        tmp=substr(names(param),1,1)
        param=matrix(param, nrow=1)
        colnames(param)=tmp
    } else {
        colnames(param)=substr(colnames(param),1,1)
    }
    
    if(ncol(param)==4) param=cbind(param, "f"=rep(1,nrow(param)))
    
    if(!"b" %in% colnames(param) & "logmb"%in%colnames(param)) b=unname(-exp(param[,"logmb"])) else b=param[,"b"]
    if(!"e" %in% colnames(param) & "loge"%in%colnames(param)) e=unname(exp(param[,"loge"])) else e=param[,"e"]
    if(!"f" %in% colnames(param) & "logf"%in%colnames(param)) f=unname(exp(param[,"logf"])) else f=param[,"f"]    
    c=param[,"c"]; d=param[,"d"]; 
    g=log(e)-(1/b)*log(f)
    h=-b*(d-c)/(1+(1/f))^(f+1)
    if(is.v) c(c,d,g=unname(g),h=unname(h),f) else cbind(c,d,g=unname(g),h=unname(h),f)
}

gh2cla=function(param){
    is.v=FALSE
    if(is.vector(param)) {
        is.v=TRUE
        tmp=substr(names(param),1,1)
        param=matrix(param, nrow=1)
        colnames(param)=tmp
    }
    
    if(ncol(param)==4) param=cbind(param, "f"=rep(1,nrow(param)))
    
    c=param[,"c"]; d=param[,"d"]; g=param[,"g"]
    if(!"f" %in% colnames(param) & "logf"%in%colnames(param)) f=unname(exp(param[,"logf"])) else f=param[,"f"]    
    if(!"h" %in% colnames(param) & "logh"%in%colnames(param)) h=unname(exp(param[,"logh"])) else h=param[,"h"]    
    
    b=-(h/(d-c))*(1+(1/f))^{f+1}
    e=exp(g+(1/b)*log(f))
    if(is.v) c(b=unname(b),c,d,e=unname(e),f) else cbind(b=unname(b),c,d,e=unname(e),f)
}

cla2ed50=function(param){
    is.v=FALSE
    if(is.vector(param)) {
        is.v=TRUE
        tmp=substr(names(param),1,1)
        param=matrix(param, nrow=1)
        colnames(param)=tmp
    }else colnames(param)=substr(colnames(param),1,1)
    
    if(ncol(param)==4) param=cbind(param, "f"=rep(1,nrow(param)))
    
    b=param[,"b"]; c=param[,"c"]; d=param[,"d"]; e=param[,"e"]; f=param[,"f"]    
    tao=e*(2^{1/f}-1)^{1/b}
    if(is.v) c(b,c,d,logtao=unname(log(tao)),f) else cbind(b,c,d,logtao=unname(log(tao)),f)
}

ed502cla=function(param){
    is.v=FALSE
    if(is.vector(param)) {
        is.v=TRUE
#        tmp=substr(names(param),1,1)
        tmp=names(param)
        param=matrix(param, nrow=1)
        colnames(param)=tmp
    }
    
    if(ncol(param)==4) param=cbind(param, "f"=rep(1,nrow(param)))
    
    c=param[,"c"]; d=param[,"d"]; b=param[,"b"]; f=param[,"f"]; logtao=param[,"logtao"]
    if(!"f" %in% colnames(param) & "logf"%in%colnames(param)) f=unname(exp(param[,"logf"])) else f=param[,"f"]    
    if(!"b" %in% colnames(param) & "logmb"%in%colnames(param)) b=unname(-exp(param[,"logmb"])) else b=param[,"b"]
    
    loge=logtao-(1/b)*log(2^{1/f}-1)
    e=exp(loge)
    if(is.v) c(b=unname(b),c,d,e=unname(e),f) else cbind(b=unname(b),c,d,e=unname(e),f)
}


# x can be a single number or a vector
# return a matrix, where each row corresponds to an x
vpl3.deriv = function(x,param){
    c=param["c"]; d=param["d"]; g=param["g"]; h=param["h"]; f=param["f"]    
    t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
    matrix(
        c(  "c" = 1-(1/((1+u)^{f})), 
            "d" = (1/((1+u)^{f})), 
            "g" = -h*((y-c)/(d-c))*f*(1+(1/f))^{f+1}*(u/(1+u)),
            "h" = (t-g)*((y-c)/(d-c))*f*(1+(1/f))^{f+1}*(u/(1+u)),
            "f" = -(y-c)*{log(1+u)-(u/(1+u))*(log(u)+log(f)+1)+f*(u/(1+u))*(log(u)+log(f))*log(1+(1/f))}
        ), 
        nrow=length(x)
    )
}
## test
#vpl3.deriv(c(1,100), theta)

# return a list of derivative functions
vpl3.deriv.func = function(param){
    c=param["c"]; d=param["d"]; g=param["g"]; h=param["h"]; f=param["f"]    
    f.c=function(x) {
        t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
        unname(1-(1/((1+u)^{f})))
    }
    f.d=function(x) {
        t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
        unname((1/((1+u)^{f})))
    }
    f.g=function(x) {
        t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
        unname(-h*((y-c)/(d-c))*f*(1+(1/f))^{f+1}*(u/(1+u)))
    }
    f.h=function(x) {
        t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
        unname((t-g)*((y-c)/(d-c))*f*(1+(1/f))^{f+1}*(u/(1+u)))
    }
    f.f=function(x) {
        t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
        unname(-(y-c)*{log(1+u)-(u/(1+u))*(log(u)+log(f)+1)+f*(u/(1+u))*(log(u)+log(f))*log(1+(1/f))})
    }
    list("c"=f.c, "d"=f.d, "g"=f.g, "h"=f.h, "f"=f.f)
}

# ylim=NULL; col=NULL; lty=NULL; lwd=1; plot.legend=FALSE; add=FALSE; legend=NULL; main=NULL # default 
plot5PL=function(param, xlim, ylim=NULL, col=NULL, lty=NULL, lwd=1, plot.legend=FALSE, add=FALSE, legend=NULL, main=NULL, x.exp=FALSE, xlab=NULL, ylab=NULL, xaxt="s") {
    
    if (!is.matrix(param) & !is.vector(param)) stop("param has to be vector or matrix")
    if(is.vector(param)) {
        tmp=names(param)
        param=matrix(param, nrow=1)
        dimnames(param)[[2]]=tmp
    }
    
    t.1=seq(xlim[1], xlim[2], length=1000)
    if (is.null(ylim)) ylim=range(apply(param,1,function(x) FivePL.t(t.1, x)))
    if (!add) {
        if (!x.exp) {
            plot(1,1,xlim=xlim, ylim=ylim, type="n", xlab=ifelse(is.null(xlab),"t",xlab), ylab=ifelse(is.null(ylab),"y",ylab), main=main, xaxt=xaxt)
        } else {
            plot(1,1,xlim=xlim, ylim=ylim, type="n", xlab=ifelse(is.null(xlab),"x",xlab), ylab=ifelse(is.null(ylab),"y",ylab), main=main, xaxt="n")
            axis(side=1, at=seq(xlim[1], xlim[2], length=10), labels=round(exp(seq(xlim[1], xlim[2], length=10)),1))
        }
    }
    if (is.null(col)) col=1:nrow(param) else if(length(col)==1) col=rep(col, nrow(param))
    if (is.null(lty)) lty=rep(1, nrow(param))
    for (i in 1:nrow(param)) {
        lines(t.1, FivePL.t(t.1, param[i,]), col=col[i], lty=lty[i], lwd=lwd)
    }
    
    if(is.null(legend)) legend=1:nrow(param)
    if (!add & plot.legend) mylegend(x=9, legend=legend, col=col, lty=lty)
}

#Treat out of bound concentration estimates
#y: a number. The readout
#p: a vector of number. Parameters for a 5pl/4pl curve.
#t.range: a vector of two numbers. The range of log standard samples concentrations.
#If y is less than lower asymptote, return t.range[1]+log(1/2), i.e. log of half of smallest standard samples concentration.
#If y is higher than upper asymptote, return t.range[2], i.e. log of largest standard samples concentration
treat.out.of.bound=function(y, p, t.range){
    if (y<get.curve.param.list(p)$c) {
        t.0=t.range[1]+log(1/2) # half of the smallest standard concentration
    } else if (y>get.curve.param.list(p)$d) {
        t.0=t.range[2] # the largest standard concentration
    } else stop("treat.out.of.bound: this cannot be right.")
    t.0
}

get.abs.dev = function(p1, p2, t.range, y.range, vpl=TRUE) {
    if (vpl) {
        f.0  =FivePL.t.inv.func(p1)
        f.hat=FivePL.t.inv.func(p2)
    } else  {
        stop ("not supported yet")
    }
    integrate( function(y) {
        t.0 =   f.0  (y)
        t.0 = sapply (1:length(y), function (i) { # as y is a vector and treat.out.of.bound is not vectorized, we need to loop through the lenght of y
            if (is.nan(t.0[i])) treat.out.of.bound (y[i], p1, t.range) else t.0[i]
        })
        t.hat = f.hat(y)
        t.hat = sapply (1:length(y), function (i) {
            if (is.nan(t.hat[i])) treat.out.of.bound (y[i], p2, t.range) else t.hat[i]
        })
        abs(t.hat - t.0)
    }, lower=y.range[1], upper=y.range[2], subdivisions=1000 )$value/(y.range[2]-y.range[1])
}

# get area between two curves betwee two curves
get.abc = function(p1, p2, t.range, vpl=TRUE) {
    if (vpl) {
        f1=FivePL.t.func( p1 )
        f0=FivePL.t.func( p2 )
    } else  {
        f1=FourPL.t.func( p1 )
        f0=FourPL.t.func( p2 )
    }
    integrate( function(t) abs(f1(t)-f0(t)) , lower=t.range[1], upper=t.range[2], subdivisions=1000 )$value/(t.range[2]-t.range[1])
}

# get S1 betwee two curves
get.S1 = function(p1, p2, t.range, vpl=TRUE) {
    if (vpl) {
        f1=FivePL.t.func( p1 )
        f0=FivePL.t.func( p2 )
    } else  {
        f1=FourPL.t.func( p1 )
        f0=FourPL.t.func( p2 )
    }
    integrate( function(t) (f1(t)-f0(t))^2 , lower=t.range[1], upper=t.range[2], subdivisions=1000 )$value/(t.range[2]-t.range[1])
}

# get S2 betwee two curves, percent bias
get.S2 = function(p1, p2, t.range, vpl=TRUE) {
    if (vpl) {
        f0=FivePL.t.func( p1 )
        f1.inv=FivePL.x.inv.func(p2)
    } else  {
        stop ("not supported yet")
    }
    integrate( function(t) {
        x0= exp(t)
        x1= f1.inv(f0(t))
        abs(x1-x0)/x0 * 100
    }, lower=t.range[1], upper=t.range[2], subdivisions=1000 )$value/(t.range[2]-t.range[1])
}

ED5PL = function (param, tao) {
    
    names(param)=substr(names(param),1,1)
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]; f=param["f"]; g=param["g"]; h=param["h"]; logtao=param["logtao"]
    if(is.na(b) | is.na(e)) {
        if(!is.na(g) & !is.na(h)) {
            pp=gh2cla(c(c,d,g,h=h,f=f)) 
            b=pp["b"]
            e=pp["e"]
        } else if (!is.na(logtao)) {
            pp=ed502cla(c(c,d,b=b,logtao=unname(logtao),f=f))            
            e=pp["e"]
        } else  {
            stop("parameterization not recognized")
        }
    } 
    
    unname(e*(tao^{-1/f}-1)^{1/b})
}
