.bacon <- function(i, object, niter, nbins, trim, level, verbose, priors){

    ##TODO add some kind of trimming?
    
    tstats <- tstat(object)[,i]
    medy <- median(tstats)
    mady <- mad(tstats)
    binned <- !is.null(nbins)

    if(verbose & binned)
        message("Use multinomial weighted sampling...")
    else if(verbose & !binned)
        message("Use fast random weighted sampling...")
    
    if(!is.null(nbins)){
        ##q <- range(tstats) ##sensitive to outlying values
        q <- quantile(tstats, prob=c(1 - trim, trim))
        tstats <- tstats[tstats > q[1] & tstats < q[2]]
        breaks <- seq(q[1], q[2], length = nbins + 1)
        x <- breaks[-c(nbins+1)] + 0.5*(q[2] - q[1])/(nbins) #identical to h$mids
        h <- hist(tstats, breaks = breaks, plot=FALSE)
        w <- h$counts
        tstats <- h$mids
    }
    else
        w <- 1.0+0*tstats

    output <- .C("bacon",
                 y = as.vector(tstats, mode="double"),
                 w = as.vector(w, mode="double"),
                 medy = as.double(medy),
                 mady = as.double(mady),
                 n = as.integer(length(tstats)),
                 niter = as.integer(niter),
                 level = as.double(level),
                 binned = as.integer(binned),
                 verbose = as.integer(verbose),
                 gibbsmu = vector(3*niter, mode="double"),
                 gibbssig = vector(3*niter, mode="double"),
                 gibbsp = vector(3*niter, mode="double"),
                 alpha = as.double(priors$sigma$alpha),
                 beta = as.double(priors$sigma$beta),
                 lambda = as.vector(priors$mu$lambda, mode="double"),
                 tau = as.vector(priors$mu$tau, mode="double"),
                 gamma = as.vector(priors$epsilon$gamma, mode="double"),
                 PACKAGE="bacon")

    gibbsp <- matrix(output$gibbsp, ncol=3, byrow=TRUE,
                     dimnames=list(NULL, paste0("p.", 0:2)))
    gibbsmu <- matrix(output$gibbsmu, ncol=3, byrow=TRUE,
                      dimnames=list(NULL, paste0("mu.", 0:2)))
    gibbssig <- matrix(output$gibbssig, ncol=3, byrow=TRUE,
                       dimnames=list(NULL, paste0("sigma.", 0:2)))

    return(cbind(gibbsp, gibbsmu, gibbssig))
}

##' Gibbs Sampler Algorithm to fit a three component normal mixture to
##' z-scores
##'
##'
##' @title Gibbs sampler
##' @param teststatistics numeric vector or matrix of test-statistics
##' @param effectsizes numeric vector or matrix of effect-sizes
##' @param standarderrors numeric vector or matrix of standard errors
##' @param niter number of iterations
##' @param nburnin length of the burnin period
##' @param nbins default 1000 else bin test-statistics
##' @param trim default 0.999 trimming test-statistics 
##' @param level significance leve used to determine prop. null for
##'     starting values
##' @param verbose default FALSE
##' @param priors list of parameters of for the prior distributions
##' @return object of class-Bacon
##' @examples
##' ##simulate some test-statistic from a normal mixture
##' ##and run bacon
##' y <- rnormmix(2000, c(0.9, 0, 1, 0, 4, 1))
##' bc <- bacon(y)
##' ##extract all estimated mixture parameters
##' estimates(bc)
##' ##extract inflation
##' inflation(bc)
##' ##extract bias
##' bias(bc)
##'
##' ##extract bias and inflation corrected test-statistics
##' head(tstat(bc))
##'
##' ##inspect the Gibbs Sampling output
##' traces(bc)
##' posteriors(bc)
##' fit(bc)
##'
##' ##simulate multiple sets of test-statistic from a normal mixture
##' ##and run bacon
##' y <- matrix(rnormmix(10*2000, c(0.9, 0, 1, 0, 4, 1)), ncol=10)
##' bc <- bacon(y)
##' ##extract all estimated mixture parameters
##' estimates(bc)
##' ##extract only the inflation
##' inflation(bc)
##' ##extract only the bias
##' bias(bc)
##' ##extract bias and inflation corrected P-values
##' head(pval(bc))
##' ##extract bias and inflation corrected test-statistics
##' head(tstat(bc))
##' @author mvaniterson
##' @references Implementation is based on a version from Zhihui Liu
##'     \url{https://macsphere.mcmaster.ca/handle/11375/9368}
##' @export
##' @importFrom BiocParallel bplapply bpworkers bpparam
##' @useDynLib bacon
bacon <- function(teststatistics=NULL, effectsizes=NULL, standarderrors=NULL,
                  niter=5000L, nburnin = 2000L, nbins=1000, trim =0.999, level=0.05, verbose=FALSE,
                  priors = list(sigma = list(alpha = 1.28,
                                             beta = 0.36), ##original uses 0.36*mad(teststatistics)
                                mu = list(lambda = c(0.0, 3.0, -3.0),
                                          tau = c(1000.0, 100.0, 100.0)),
                                epsilon = list(gamma = c(90.0, 5.0, 5.0)))){

    ##create new Bacon-object
    object <- new("Bacon", teststatistics = teststatistics,
                  effectsizes = effectsizes,
                  standarderrors = standarderrors,
                  niter = niter,
                  nburnin = nburnin,
                  priors = priors)

    ##run the Gibbs Sampler
    if(ncol(tstat(object)) > 1){
        nworkers <- bpworkers(bpparam())
        if(nworkers <= 1) {
            message("Did you registered a biocparallel back-end?\n Continuing serial!")
            for(i in 1:ncol(tstat(object)))
                object@traces[,,i] <- .bacon(i, object, niter, nbins, trim, level, verbose, priors)
        }
        else{
            message("Detected ", nworkers, " workers!\n Running in parallel!")
            ret <- bplapply(1:ncol(tstat(object)), .bacon, object=object,
                            niter=niter, nbins=nbins, trim=trim, level=level, verbose=verbose, priors=priors)
            object@traces <- simplify2array(ret)
        }
    } else
        object@traces[,,1] <- .bacon(1, object, niter, nbins, trim=trim, level, verbose, priors)

    ##summarize traces
    for(i in 1:ncol(tstat(object)))
        object@estimates[i,] <- apply(object@traces[-c(1:object@nburnin),,i], 2, mean)

    return(object)
}
