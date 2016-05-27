
##' @rdname estimates-methods
##' @aliases estimates
setMethod("estimates","Bacon", function(object){
    if(nrow(object@traces) != object@niter) {
        return(NULL)
        message("No traces to calculate estimates!")
    } else
        return(object@estimates)
})

##' @rdname inflation-methods
##' @aliases inflation
setMethod("inflation","Bacon", function(object){ return(object@estimates[,7]) })

##' @rdname bias-methods
##' @aliases bias
setMethod("bias","Bacon", function(object){ return(object@estimates[,4]) })

##' @rdname tstat-methods
##' @aliases tstat
setMethod("tstat","Bacon", function(object, corrected){
    if(!corrected | any(is.na(bias(object))) | any(is.na(inflation(object))))
        teststatistics <- object@teststatistics
    else {
        teststatistics <- t(t(object@teststatistics) - bias(object))
        teststatistics <- t(t(teststatistics)/inflation(object))
    }
    return(teststatistics)
})

##' @rdname pval-methods
##' @aliases pval
setMethod("pval","Bacon", function(object, corrected){
    if(!corrected | any(is.na(bias(object))) | any(is.na(inflation(object))))
        pvalues <- 2*pnorm(-abs(object@teststatistics))
    else {
        teststatistics <- t(t(object@teststatistics) - bias(object))
        teststatistics <- t(t(teststatistics)/inflation(object))
        pvalues <- 2*pnorm(-abs(teststatistics), mean=0, sd=1, lower.tail=TRUE, log.p=FALSE)
    }
    return(pvalues)
})

##' @rdname es-methods
##' @aliases es
setMethod("es","Bacon", function(object, corrected){
    if(nrow(object@effectsizes) == 1)
        stop("Effect-sizes not provided!")
    if(!corrected | any(is.na(bias(object))) | any(is.na(inflation(object))))
        effectsizes <- object@effectsizes
    else
        effectsizes <- object@effectsizes - t(t(object@standarderrors)*bias(object))
    return(effectsizes)
})

##' @rdname se-methods
##' @aliases se
setMethod("se","Bacon", function(object, corrected){
    if(nrow(object@standarderrors) == 1)
        stop("Standard errors not provided!")

    if(!corrected |  any(is.na(bias(object))) | any(is.na(inflation(object))))
        standarderrors <- object@standarderrors
    else
        standarderrors <- t(t(object@standarderrors)*inflation(object))
    return(standarderrors)
})

##' @rdname traces-methods
##' @aliases traces
setMethod("traces", "Bacon", function(object, burnin=TRUE, index=1){

    if(burnin)
        gstraces <- object@traces[,,index]
    else
        gstraces <- object@traces[-c(1:object@nburnin),,index]

    thetahat <- estimates(object)[index,]

    op <- par(mfcol=c(3, 3), mar=c(2,4,2,2))
    for(i in 1:9) {
        if(burnin)
            plot(1:object@niter, gstraces[,i],
                 ylab=colnames(gstraces)[i], type="l", xlab="", main = "", lwd=0.3)
        else
            plot((object@nburnin+1):object@niter, gstraces[,i],
                 ylab=colnames(gstraces)[i], type="l", xlab="", main = "", lwd=0.3)
        abline(h=thetahat[i], col=2, lwd=2)
    }
    par(op)
})

##' @rdname posteriors-methods
##' @aliases posteriors
setMethod("posteriors", "Bacon", function(object, thetas, index, alphas, xlab, ylab, ...){

    if(any(!(thetas %in% colnames(object@traces[,, index]))))
        stop("'thetas' should be two of: ", paste(colnames(object@traces[,, index]), collapse=", "), "!")

    gstraces <- object@traces[-c(1:object@nburnin), thetas, index]

    if(xlab=="") xlab <- thetas[1]
    if(ylab=="") ylab <- thetas[2]

    plot(gstraces, pch=20, xlab=xlab, ylab=ylab, bty='n',
         main=c("median at:", round(estimates(object)[thetas], 3)))
    points(estimates(object)[index, thetas], col=3, pch=17, cex=2)

    for(alpha in alphas)
        lines(ellipse(cov(gstraces), centre=colMeans(gstraces), level=alpha), col="blue", ...)
})

##' @rdname fit-methods
##' @aliases fit
setMethod("fit", "Bacon", function(object, index,...){
    plotnormmix(tstat(object, corrected=FALSE)[, index], estimates(object)[index, ], ...)
})

##' @rdname meta-methods
##' @aliases meta
setMethod("meta", "Bacon",function(object, corrected=TRUE, ...){
    ES <- es(object, corrected=corrected)
    SE <- se(object, corrected=corrected)
    W <- 1/SE^2
    V <- 1/rowSums(W)
    TS <- rowSums(ES*W)*V
    Z <- TS/sqrt(V)

    ##TODO maybe add more from:
    ##http://www.netstorm.be/home/meta_analysis#metaAnalysisU
    ##Q <- rowSums((ES - TS)^2*W)
    ##Qc <- qchisq(0.95, df= ncol(ES)-1)
    ##Qp <- 1 - pchisq(Q, df= ncol(ES)-1)
    ##P <- 2*pnorm(-abs(Z))

    object@teststatistics <- cbind(object@teststatistics, meta=Z)
    object@effectsizes <- cbind(object@effectsizes, meta=TS)
    object@standarderrors <- cbind(object@standarderrors, meta=sqrt(V))
    object@estimates <- rbind(object@estimates, meta=c(NA,NA,NA,0, NA, NA, 1, NA, NA))
    invisible(object)
})

##' @rdname topTable-methods
##' @aliases topTable
setMethod("topTable", "Bacon", function(object,
                                        number=10,
                                        adjust.method="bonf",
                                        sort.by=c("pval", "eff.size")){

    sort.by <- match.arg(sort.by)

    ##only return top table after meta-analysis has been performed
    if(!grepl("meta", colnames(pv)))
        stop("First run fixed-effect meta-analysis using the 'meta'-function!")


    pv <- pval(object)
    n <- ncol(pv)
    padj <- p.adjust(pv[,n], method=adjust.method)
    tst <- tstat(object)

    eff <- es(object)
    std <- se(object)

    data <- cbind(eff[, -n], std[, -n], pv[, -n], tst[,-n])
    colnames(data) <- paste(rep(c("eff.size", "std.err", "pval", "tstat"), each=n-1),
                            colnames(tst)[-n], sep=".")

    ##reorder columns cohorts together
    data <- data[, order(gsub(".*\\.", "", colnames(data)))]

    meta <- cbind(eff[,n], std[,n], padj, pv[,n], tst[,n])
    colnames(meta) <- paste(c("eff.size", "std.err", "pval.adj", "pval.org", "tstat"),
                            "meta", sep=".")

    if(number==-1)
        number <- nrow(pv)

    if(sort.by=="pval")
        topId <- order(pv[,n])[1:number]
    else if(sort.by=="eff.size")
        topId <- order(eff[,n])[1:number]

    tt <- cbind(meta[topId,], data[topId,])

    tt <- as.matrix(tt)
    colnames(tt) <- c(colnames(meta), colnames(data))


    rownames(tt) <- rownames(pv)[topId]
    invisible(tt)
})
