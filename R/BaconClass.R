##' An S4 class container for storing Gibbs Sampler input and output
##'
##' @slot teststatistics numeric vector or matrix of test-statistics
##' @slot effectsizes numeric vector or matrix of effect-sizes
##' @slot standarderrors numeric vector or matrix of standard errors
##' @slot traces array of Gibbs Sampler traces
##' @slot estimates vector or matrix of parameter estimates
##' @slot priors list of parameters of for the prior distributions
##' @slot niter number of iterations
##' @slot nburnin length of the burnin period
##' @import methods
setClass("Bacon",
         representation(teststatistics = "matrix",
                        effectsizes    = "matrix",
                        standarderrors = "matrix",
                        traces         = "array",
                        estimates      = "matrix",
                        priors         = "list",
                        niter          = "integer",
                        nburnin        = "integer",
                        na.exclude     = "logical"),
         prototype(teststatistics = matrix(1),
                   effectsizes    = matrix(1),
                   standarderrors = matrix(1),
                   traces         = array(1),
                   estimates      = matrix(1),
                   priors         = list(),
                   niter          = integer(1),
                   nburnin        = integer(1),
                   na.exclude     = logical(1)),
         validity = function(object){
             if(niter <= nburnin)
                 return("niter should be > nburnin!")
             else
                 return(TRUE)}
         )

##' Method to extract the estimated parameters from the 'bacon'-object
##' @name estimates
##' @param object 'bacon'-object
##' @return vector or matrix of estimates
##' @rdname estimates-methods
##' @seealso \code{\link{bacon}}
##' @exportMethod estimates
##' @examples
##' y <- rnormmix(2000, c(0.9, 0, 1, 0, 4, 1))
##' ##nbins = 100 to speed up the calculations
##' bc <- bacon(y, nbins=100)
##' estimates(bc)
setGeneric("estimates", function(object){ standardGeneric("estimates") })

##' Method to extract the estimated inflation from the 'bacon'-object
##' @name inflation
##' @param object 'bacon'-object
##' @return vector or matrix of inflation
##' @rdname inflation-methods
##' @seealso \code{\link{bacon}}
##' @exportMethod inflation
##' @examples
##' y <- rnormmix(2000, c(0.9, 0, 1, 0, 4, 1))
##' ##nbins = 100 to speed up the calculations
##' bc <- bacon(y, nbins=100)
##' inflation(bc)
setGeneric("inflation", function(object){ standardGeneric("inflation") })

##' Method to extract the estimated bias from the 'bacon'-object
##' @name bias
##' @param object 'bacon'-object
##' @return vector or matrix of inflation
##' @rdname bias-methods
##' @seealso \code{\link{bacon}}
##' @exportMethod bias
##' @examples
##' y <- rnormmix(2000, c(0.9, 0, 1, 0, 4, 1))
##' ##nbins = 100 to speed up the calculations
##' bc <- bacon(y, nbins=100)
##' bias(bc)
setGeneric("bias", function(object){ standardGeneric("bias") })

##' Method to extract inflation- and bias-corrected test-statistics
##' @name tstat
##' @param object 'bacon'-object
##' @param corrected optional return uncorrected
##' @return vector or matrix of test-statistics
##' @rdname tstat-methods
##' @seealso \code{\link{bacon}}
##' @exportMethod tstat
##' @examples
##' y <- rnormmix(2000, c(0.9, 0, 1, 0, 4, 1))
##' ##nbins = 100 to speed up the calculations
##' bc <- bacon(y, nbins=100)
##' head(tstat(bc))
setGeneric("tstat", function(object, corrected=TRUE){ standardGeneric("tstat") })

##' Method to extract inflation- and bias-corrected P-values
##' @name pval
##' @rdname pval-methods
##' @param object 'bacon'-object
##' @param corrected optional return uncorrected
##' @return vector or matrix of P-values
##' @seealso \code{\link{bacon}}
##' @exportMethod pval
##' @examples
##' y <- rnormmix(2000, c(0.9, 0, 1, 0, 4, 1))
##' bc <- bacon(y, nbins=100) #nbins = 100 to speed up the calculations
##' head(pval(bc))
setGeneric("pval", function(object, corrected=TRUE){ standardGeneric("pval") })

##' Method to extract inflation- and bias-corrected effect-sizes
##' @name es
##' @rdname es-methods
##' @param object 'bacon'-object
##' @param corrected optional return uncorrected
##' @return vector or matrix of effect-sizes
##' @seealso \code{\link{bacon}}
##' @exportMethod es
##' @examples
##' es <- replicate(6, rnormmix(2000, c(0.9, 0, 1, 0, 4, 1)))
##' se <- replicate(6, 0.8*sqrt(4/rchisq(2000,df=4)))
##' bc <- bacon(NULL, es, se)
##' head(es(bc))
setGeneric("es", function(object, corrected=TRUE){ standardGeneric("es") })

##' Method to extract inflation- and bias-corrected standard errors
##' @name se
##' @rdname se-methods
##' @param object 'bacon'-object
##' @param corrected optional return uncorrected
##' @return vector or matrix of standard-errors
##' @seealso \code{\link{bacon}}
##' @exportMethod se
##' @examples
##' es <- replicate(6, rnormmix(2000, c(0.9, 0, 1, 0, 4, 1)))
##' se <- replicate(6, 0.8*sqrt(4/rchisq(2000,df=4)))
##' bc <- bacon(NULL, es, se)
##' head(se(bc))
setGeneric("se", function(object, corrected=TRUE){ standardGeneric("se") })

##' Method to plot Gibbs sampling traces
##' @name traces
##' @rdname traces-methods
##' @param object 'bacon'-object
##' @param burnin include burnin period default true
##' @param index if multiple sets of test-statsistics where provided
##' @return plot of the Gibbs Sampler traces
##' @seealso \code{\link{bacon}}
##' @exportMethod traces
##' @importFrom graphics abline curve hist lines par points
##' @examples
##' y <- rnormmix(2000, c(0.9, 0, 1, 0, 4, 1))
##' ##nbins = 100 to speed up the calculations
##' bc <- bacon(y, nbins=100)
##' traces(bc)
setGeneric("traces", function(object, burnin=TRUE, index=1){ standardGeneric("traces") })

##' Method to plot posterior distribution
##' @name posteriors
##' @rdname posteriors-methods
##' @param object 'bacon'-object
##' @param thetas which thetas to plot
##' @param index if multiple sets of test-statsistics where provided
##' @param alphas significance level confidence ellipses
##' @param xlab optional xlab
##' @param ylab  optional ylab
##' @param ... additional plotting parameters
##' @return plot of the Gibbs Sampler posterior probabilities
##' @import ellipse
##' @seealso \code{\link{bacon}}
##' @exportMethod posteriors
##' @examples
##' y <- rnormmix(2000, c(0.9, 0, 1, 0, 4, 1))
##' ##nbins = 100 to speed up the calculations
##' bc <- bacon(y, nbins=100)
##' posteriors(bc)
setGeneric("posteriors", function(object,
                                  thetas = c("sigma.0", "p.0"), index = 1,
                                  alphas=c(0.95, 0.9, 0.75), xlab="", ylab="", ...){
    standardGeneric("posteriors")
})

##' Method to plot mixture fit
##' @name fit
##' @rdname fit-methods
##' @param object 'bacon'-object
##' @param index if multiple sets of test-statsistics where provided
##' @param col line color default 'grey75'
##' @param border  border color 'grey75'
##' @param ... additional plotting parameters
##' @return plot of the Gibbs Sampler mixture fit
##' @seealso \code{\link{bacon}}
##' @exportMethod fit
##' @examples
##' y <- rnormmix(2000, c(0.9, 0, 1, 0, 4, 1))
##' ##nbins = 100 to speed up the calculations
##' bc <- bacon(y, nbins=100)
##' fit(bc)
setGeneric("fit", function(object, index=1, ...){ standardGeneric("fit") })

##' Perform fixed meta-analysis using inflation and bias corrected
##' effect-sizes and standard errors
##'
##' TODO maybe add idea's from http://www.netstorm.be/home/meta_analysis#metaAnalysisU
##' @name meta
##' @rdname meta-methods
##' @title fixed meta-analysis
##' @param object 'bacon'-object
##' @param corrected optional return uncorrected
##' @param ... additional arguments
##' @return object of class 'bacon' with added fixed-effect meta-analysis
##' test-statistics, effect-sizes and standard-errors
##' @import stats
##' @seealso \code{\link{bacon}}
##' @exportMethod meta
##' @examples
##' es <- replicate(6, rnormmix(2000, c(0.9, 0, 1, 0, 4, 1)))
##' se <- replicate(6, 0.8*sqrt(4/rchisq(2000,df=4)))
##' bc <- bacon(NULL, es, se)
##' mbc <- meta(bc)
setGeneric("meta", function(object, corrected=TRUE, ...){ standardGeneric("meta") })

##' Extract top features after meta analysis
##'
##' @name topTable
##' @rdname topTable-methods
##' @title topTable
##' @param object 'bacon'-object
##' @param number return specified number of top features, n=-1 return all features
##' @param adjust.method P-value multiple testing adjustment method default bonferroni
##' @param sort.by order results by pval or eff.size
##' @return table with top features
##' @seealso \code{\link{bacon}}
##' @exportMethod topTable
##' @examples
##' es <- replicate(6, rnormmix(2000, c(0.9, 0, 1, 0, 4, 1)))
##' se <- replicate(6, 0.8*sqrt(4/rchisq(2000,df=4)))
##' bc <- bacon(NULL, es, se)
##' mbc <- meta(bc)
##' topTable(mbc)
setGeneric("topTable", function(object, number=10, adjust.method="bonf", sort.by=c("pval", "eff.size")){ standardGeneric("topTable") })

setMethod("initialize", "Bacon",
          function(.Object, teststatistics, effectsizes, standarderrors, 
                            niter, nburnin, priors, na.exclude) {

              if(!is.null(teststatistics)){
                  .Object@teststatistics <- as.matrix(teststatistics)
              }
              else if(is.null(effectsizes) & is.null(standarderrors))
                  stop("Need to provide test-statistics or both effect-sizes and standard errors!")
              else {
                  .Object@effectsizes <- as.matrix(effectsizes)
                  .Object@standarderrors <- as.matrix(standarderrors)
                  .Object@teststatistics <- as.matrix(effectsizes/standarderrors)
              }

              if(!all(is.finite(.Object@teststatistics)) & !na.exclude)
                  stop("Non finite value(s) in test statistics and na.exclude = FALSE!")
              
              .Object@estimates <- matrix(nrow=ncol(.Object@teststatistics), ncol=9,
                                          dimnames=list(colnames(.Object@teststatistics),
                                                        paste0(rep(c("p.", "mu.", "sigma."), each=3), 0:2)))
              .Object@traces <- array(dim=c(niter, 9, ncol(.Object@teststatistics)),
                                      dimnames=list(NULL,
                                                    paste0(rep(c("p.", "mu.", "sigma."), each=3), 0:2),
                                                    colnames(.Object@teststatistics)))

              .Object@niter <- niter
              .Object@nburnin <- nburnin
              .Object@priors <- priors
              .Object@na.exclude <- na.exclude
              .Object
          })

setMethod("show", "Bacon", function (object) {
    cat(sprintf("%s-object containing %s set(s) of %s test-statistics.\n",
                class(object),
                ncol(tstat(object)),
                nrow(tstat(object))))

    cat(sprintf("...estimated bias: %s.\n",
                paste(signif(bias(object), 2), collapse=",")))

    cat(sprintf("...estimated inflation: %s.\n\n",
                paste(signif(inflation(object), 2), collapse=",")))

    cat(sprintf("Empirical null estimates are based on %s iterations with a burnin-period of %s.\n",
                object@niter,
                object@nburnin))

    ##maybe add cat(sprintf("...prior parameters...\n"))
})


n2mfcol <- function(n){
    if(n<1)
        stop("Not a (postive) integer!")
    switch(as.character(n),
           "1" = c(1, 1),
           "2" = c(1, 2),
           "3" = c(2, 2),
           "4" = c(2, 2),
           "5" = c(2, 3),
           "6" = c(2, 3),
           "7" = c(3, 3),
           "8" = c(3, 3),
           stop("Too many sets of statistics to visualize nicely!"))
}

.hist <- function(object, trim, ...) {

    mu <- bias(object)
    sigma <- inflation(object)

    tstats <- tstat(object, corrected=FALSE)
    
    if(is.null(colnames(tstats))) colnames(tstats) <- LETTERS[1:ncol(tstats)]
    
    nas <- is.na(tstats)
    tstats <- na.omit(tstats)
    column <- rep(colnames(tstats), each=nrow(tstats))
    
    stdnorm <- apply(tstats, 2, dnorm, mean=0, sd=1)
    empnull <- 0*stdnorm
    
    for(i in 1:ncol(tstats)) empnull[,i] <- dnorm(tstats[,i], mean=mu[i], sd=sigma[i])
    
    data <- data.frame(tstats = as.vector(tstats),
                       column = column[as.vector(!nas)],
                       stdnorm = as.vector(stdnorm),
                       empnull = as.vector(empnull))

    nrow <- n2mfcol(nlevels(factor(data$column)))[1]
    ncol <- n2mfcol(nlevels(factor(data$column)))[2]
   
    qnts <- quantile(data$tstats, prob=c(trim, 1-trim))
    data <- data[data$tstat > qnts[1] & data$tstat < qnts[2],]
    
    gp <- ggplot(data, aes(x = tstats, ...))
    gp <- gp + geom_histogram(aes(y = ..density..), colour = "grey", fill = "grey")
    gp <- gp + facet_wrap(~column, nrow=nrow, ncol=ncol)
    gp <- gp + xlab("test-statistics")
    gp <- gp + geom_line(aes(y = stdnorm), col=1) ##std. norm. black
    gp <- gp + geom_line(aes(y = empnull), col=2) ##emp. null red
    gp
}

.qq <- function(object, ...){
    
    pvalues  <- pval(object, corrected=FALSE)
    
    if(is.null(colnames(pvalues))) colnames(pvalues) <- LETTERS[1:ncol(pvalues)]
    
    nas <- is.na(pvalues)
    pvalues <- na.omit(pvalues)
    column <- rep(colnames(pvalues), each=nrow(pvalues))
  
    d1 <- data.frame(pvalues = as.vector(pvalues),
                     column = column[as.vector(!nas)],
                     bacon = "uncorrected")

    pvalues <- pval(object, corrected=TRUE)
    
    if(is.null(colnames(pvalues))) colnames(pvalues) <- LETTERS[1:ncol(pvalues)]
    
    nas <- is.na(pvalues)
    pvalues <- na.omit(pvalues)
    column <- rep(colnames(pvalues), each=nrow(pvalues))
    
    d2 <- data.frame(pvalues = as.vector(pvalues),
                     column = column[as.vector(!nas)],
                     bacon = "corrected")

    data <- rbind(d1, d2)

    gp <- ggplot(data, aes(sample=-log10(pvalues), colour=column, ...))
    gp <- gp + stat_qq(distribution=stats::qexp, dparams=list(rate=1/log10(exp(1))))
    gp <- gp + xlab(expression(paste("Expected -log"[10], plain(P))))
    gp <- gp + ylab(expression(paste("Observed -log"[10], plain(P))))
    gp <- gp + geom_abline(slope=1, intercept=0)
    gp <- gp + facet_wrap(~bacon)
    gp
}

##' simple ggplot2 plotting function for 'bacon'-object
##' @title plot hist or qq
##' @param x 'bacon'-object
##' @param y NULL
##' @param type hist or qq
##' @return either qq-plot of P-values or histogram of Test-statistics
##' @export
##' @import ggplot2
setMethod("plot", "Bacon", function(x, y, type=c("hist", "qq")) {
    type <- match.arg(type)
    switch(type,
           hist = .hist(x, trim=0.01),
           qq = .qq(x))
})
