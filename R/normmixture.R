##' sample from a normal mixture
##'
##' details follow
##' @title sample from a normal mixture
##' @param n size
##' @param theta parameters
##' @param shuffle shuffle return vectors or keep nulls and alternative ordered (null, alts)
##' @return n samples from a normal mixture with parameters theta
##' @author mvaniterson
##' @examples
##' n <- 2000
##' theta <- c(0.8, 0, 1, 0, 4, 1)
##' x <- rnormmix(n, theta)
##' @export
rnormmix <- function(n, theta, shuffle=TRUE){
    epsilon1 <- theta[1]
    mu0 <- theta[2]
    sigma0 <- theta[3]
    mu1 <- theta[4]
    sigma1 <- theta[5]
    sigma2 <- theta[6]

    U <- runif(n, min=0, max=1)
    n1 <- sum(U < epsilon1)
    n2 <- sum(U >=  epsilon1)
    mu <- rnorm(n2, mean=mu1, sd=sigma1)

    if(shuffle)
        return(sample(c(rnorm(n1, mean=mu0, sd=sigma0), mu + rnorm(n2, mean=0, sd=sigma2))))
    else
        return(c(rnorm(n1, mean=mu0, sd=sigma0), mu + rnorm(n2, mean=0, sd=sigma2)))
}

##' density of a k-component normal mixture
##'
##' details follow
##' @title density of a k-component normal mixture
##' @param x x like dnorm(x, ...
##' @param theta parameters of the mixture proportion, mean and sd
##' @return density of a k-component normal mixture
##' @author mvaniterson
##' @examples
##' n <- 2000
##' theta <- c(0.8, 0, 1, 0, 4, 1)
##' x <- rnormmix(n, theta)
##' hist(x, freq=FALSE, n=100)
##' curve(dnormmix(x, theta), add=TRUE, lwd=2)
##' @export
dnormmix <- function(x, theta){
    if(length(theta) %% 3 != 0)
        stop("Length of theta should be a multiple of three!")
    ncomp <- length(theta)/3
    y <- 0*x
    for(k in 1:ncomp)
        y <- y + theta[k]*dnorm(x, mean=theta[k+3], sd=theta[k+6], log=FALSE)
    y
}

##' plot normal mixtures
##'
##' details follow
##' @title plot normal mixtures
##' @param x vector of test statistics
##' @param theta parameters describing the mixture components
##' @param ... arguments passed to hist
##' @return return plot with histogram of the data and mixture and individual components
##' @author mvaniterson
##' @examples
##' n <- 2000
##' theta <- c(0.8, 0, 1, 0, 4, 1)
##' x <- rnormmix(n, theta)
##' plotnormmix(x, theta)
##' @export
plotnormmix <- function(x, theta, ...) {
    if(length(theta) %% 3 != 0)
        stop("Length of theta should be a multiple of three!")
    x <- data.frame(x = x)
    theta <- data.frame(y = theta)
    fit <- ggplot(x, aes(x=x , y = after_stat(density))) +
    geom_histogram(fill = "grey", color="black", binwidth = 1.5) +
    geom_line(aes(x=x, y =dnorm(x, mean(x), sd(x))), lwd=1) +
    geom_line(aes(x=x, 
              y=theta["p.0",]*dnorm(x, theta["mu.0",], theta["sigma.0",])), 
              color="red",
              lwd = 1.5) +
    geom_line(aes(x=x, 
              y=theta["p.1",]*dnorm(x, theta["mu.1",], theta["sigma.1",])), 
              color="green",
              lwd=1.5) +
    geom_line(aes(x=x, 
              y=theta["p.2",]*dnorm(x, theta["mu.2",], theta["sigma.2",])), 
              color="blue",
              lwd=1.5) +
    theme_cowplot(font_size = 12) +
    xlab("Test Statistics") +
    ylab("Density")
  
    return(fit)
}
