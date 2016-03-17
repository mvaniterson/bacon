Controlling Inflation and Bias in Epigenome- and Transcriptome-wide Association Studies using Empirical Calibration
================
Maarten van Iterson, Erik van Zwet, Eline Slagboom and Bas Heijmans Department of Molecular Epidemiology, Leiden University Medical Center, Leiden, The Netherlands
2 March 2016

Introduction
============

*[bacon](http://bioconductor.org/packages/release/bioc/html/bacon.html)* can be used to remove inflation and bias often observed in epigenome- and transcriptome-wide association studies(Iterson et al. 2016).

To this end *[bacon](http://bioconductor.org/packages/release/bioc/html/bacon.html)* constructs an empirical null distribution using a Gibbs Sampling algorithm by fitting a three-component normal mixture on z-scores. One component is forced, using prior knowledge, to represent the null distribution with mean and standard deviation representing the bias and inflation. The other two components are necessary to capture the amount of true associations present in the data, which we assume unknown but small.

*[bacon](http://bioconductor.org/packages/release/bioc/html/bacon.html)* provides functionality to inspect the output of the Gibbs Sampling algorithm, i.e., traces, posterior distributions and the mixture fit, are provided. Furthermore, inflation- and bias-corrected test-statistics or P-values are extracted easily. In addition, functionality for performing fixed-effect meta-analysis are provided.

The function `bacon` requires a vector or a matrix of z-scores, e.g., those extracted from association analyses using a linear regression approach. For fixed-effect meta-analysis a matrix of effect-sizes and standard errors is required.

The vignette illustrates the use of *[bacon](http://bioconductor.org/packages/release/bioc/html/bacon.html)* using simulated z-scores, effect-sizes and standard errors to avoid long runtimes. If multiple sets of test-statisics or effect-sizes and standard errors are provided, the Gibbs Sampler Algorithm can be executed on multiple nodes to reduce computation time using functionality provide by *[BiocParallel](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html)*-package.

One set of test-statistics
==========================

A vector containing 2000 z-scores is generated from a normal mixture distribution 90% of the z-scores were drawn from \(\mathcal{N}(0,1)\) and the remaining z-scores from \(\mathcal{N}(\mu,1)\), where \(\mu \sim \mathcal{N}(4,1)\).

The function `bacon` executes the Gibbs Sampler Algorithm and stores all input and output in a `Bacon`-object. Several accessor-functions are available to access data contained in the `Bacon`-object.

``` r
library(bacon)
set.seed(12345)
y <- rnormmix(2000, c(0.9, 0, 1, 0, 4, 1))
bc <- bacon(y)
bc
```

    ## Bacon-object containing 1 set(s) of 2000 test-statistics.
    ## ...estimated bias: -0.0092.
    ## ...estimated inflation: 0.97.
    ## 
    ## Emprical null estimates are based on 5000 iterations with a burnin-period of 2000.

``` r
estimates(bc)
```

    ##            p.0        p.1        p.2         mu.0     mu.1      mu.2
    ## [1,] 0.9065106 0.04109918 0.05239019 -0.009192046 2.977157 -3.017523
    ##        sigma.0  sigma.1  sigma.2
    ## [1,] 0.9739519 2.795116 3.133882

``` r
inflation(bc)
```

    ##   sigma.0 
    ## 0.9739519

``` r
bias(bc)
```

    ##         mu.0 
    ## -0.009192046

Several methods are provided to inspect the output of the Gibbs Sampler Algorithm, such as traces of all estimates, the posterior distributions, provide as a scatter plot between two parameters, the actual fit of the three component mixture to the histogram of z-scores.

``` r
traces(bc, burnin=FALSE)
```

![](bacon_files/figure-markdown_github/gsoutput-1.png)<!-- -->

``` r
posteriors(bc)
```

![](bacon_files/figure-markdown_github/gsoutput-2.png)<!-- -->

``` r
fit(bc, n=100)
```

![](bacon_files/figure-markdown_github/gsoutput-3.png)<!-- -->

There is also a generic plot function that can generate two types of plots. A histogram of the z-scores with on top the standard normal distribution and the bacon estimated empirical null distribution. Or a quantile-quantile plot of the \(-log_{10}\) transformed P-values.

``` r
plot(bc, type="hist", n=100, xlab="", ylab="", main="")
```

![](bacon_files/figure-markdown_github/plothist1-1.png)<!-- -->

``` r
plot(bc, type="qq", pch=1, xlab="expected", ylab="observed", main="")
```

![](bacon_files/figure-markdown_github/plothist1-2.png)<!-- -->

Multiple sets of test-statistics
================================

Matrices containing \(2000\times6\) effect-sizes and standard errors are generated to simulated data for a fixed-effect meta-analyses.

By default the function `bacon` detects the number of cores/nodes registered, as described in the *[BiocParallel](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html)*, to perform bacon in parallel. To run the vignette in general we set it here for convenience to 1 node.

``` r
library(bacon)
set.seed(12345)
es <- replicate(6, rnormmix(2000, c(0.9, 0, 1, 0, 4, 1)))
se <- replicate(6, 0.8*sqrt(4/rchisq(2000,df=4)))
colnames(es) <- colnames(se) <- LETTERS[1:ncol(se)]
rownames(es) <- rownames(se) <- 1:2000
head(rownames(es))
```

    ## [1] "1" "2" "3" "4" "5" "6"

``` r
head(colnames(es))
```

    ## [1] "A" "B" "C" "D" "E" "F"

``` r
library(BiocParallel)
register(MulticoreParam(1, log=TRUE))
bc <- bacon(NULL, es, se)
```

    ## Did you registered a biocparallel back-end?
    ##  Continuing serial!

``` r
bc
```

    ## Bacon-object containing 6 set(s) of 2000 test-statistics.
    ## ...estimated bias: -0.014,-0.023,-0.00037,0.028,0.00098,-0.0089.
    ## ...estimated inflation: 1.1,1.1,1.1,1.1,1,1.
    ## 
    ## Emprical null estimates are based on 5000 iterations with a burnin-period of 2000.

``` r
estimates(bc)
```

    ##         p.0        p.1        p.2          mu.0     mu.1      mu.2
    ## A 0.8830234 0.03506520 0.08191139 -0.0143163201 2.868507 -2.297921
    ## B 0.8792335 0.03045660 0.09030987 -0.0229185423 2.963691 -2.232250
    ## C 0.8802886 0.05452202 0.06518936 -0.0003660301 2.798501 -2.832849
    ## D 0.8923815 0.07744259 0.03017588  0.0281181199 1.792805 -2.998751
    ## E 0.8670583 0.09495013 0.03799162  0.0009763801 1.628463 -2.913473
    ## F 0.8614458 0.08444307 0.05411111 -0.0088889235 2.402517 -2.727277
    ##    sigma.0  sigma.1   sigma.2
    ## A 1.112229 2.460015 5.1908447
    ## B 1.093458 1.422252 5.2182474
    ## C 1.082837 3.941034 3.6474522
    ## D 1.133225 5.431706 0.6814031
    ## E 1.048074 5.110939 1.5047098
    ## F 1.048279 4.614447 3.2496472

``` r
inflation(bc)
```

    ##        A        B        C        D        E        F 
    ## 1.112229 1.093458 1.082837 1.133225 1.048074 1.048279

``` r
bias(bc)
```

    ##             A             B             C             D             E 
    ## -0.0143163201 -0.0229185423 -0.0003660301  0.0281181199  0.0009763801 
    ##             F 
    ## -0.0088889235

``` r
head(tstat(bc))
```

    ##            A          B           C           D          E           F
    ## 1  1.6729718 -1.0820172 -1.25647093 -3.51574318  0.2136967 -0.04993011
    ## 2 -0.9705573 -1.3746489  0.26593945  0.24936918  0.3173021  0.34490647
    ## 3 -2.3941099 -3.5741371  0.06156933 -0.87574774 -0.1632296  0.05507560
    ## 4 -0.6631808 -0.5541443  0.53077341  0.03065278 -0.4603291  1.01139602
    ## 5  0.2382283  1.2571103 -0.89428739  1.00604737 -0.0535097 -0.08733374
    ## 6  0.9410900 -2.1828069  0.97804957  1.05891670  0.3546004  1.83774199

``` r
head(pval(bc))
```

    ##            A            B         C            D         E          F
    ## 1 0.09433285 0.2792448686 0.2089453 0.0004385249 0.8307836 0.96017809
    ## 2 0.33176878 0.1692403135 0.7902858 0.8030752236 0.7510144 0.73016469
    ## 3 0.01666076 0.0003513848 0.9509058 0.3811671883 0.8703377 0.95607823
    ## 4 0.50721480 0.5794801073 0.5955758 0.9755464456 0.6452800 0.31182692
    ## 5 0.81170406 0.2087137054 0.3711681 0.3143927840 0.9573258 0.93040624
    ## 6 0.34665876 0.0290500302 0.3280498 0.2896377181 0.7228890 0.06610043

``` r
head(se(bc))
```

    ##           A         B         C         D         E         F
    ## 1 0.9399676 1.3921374 1.2634769 0.8932216 0.9446594 0.5745967
    ## 2 0.5686220 0.6185539 0.9605522 2.2718792 1.3006636 1.3560164
    ## 3 2.3006394 0.5083181 0.9391743 1.2413721 0.7490925 0.8564455
    ## 4 0.6100582 0.7596804 1.2295672 1.5054981 1.2763231 0.5644831
    ## 5 0.7086247 0.6847576 0.6606352 1.3753840 0.9339212 0.7392906
    ## 6 0.6982982 0.6258524 0.8701762 1.9655763 1.0804439 0.6274508

``` r
head(es(bc))
```

    ##            A          B           C           D           E           F
    ## 1  1.5725393 -1.5063166 -1.58752197 -3.14033762  0.20187056 -0.02868968
    ## 2 -0.5518803 -0.8502944  0.25544872  0.56653664  0.41270326  0.46769884
    ## 3 -5.5079835 -1.8167985  0.05782433 -1.08712878 -0.12227405  0.04716925
    ## 4 -0.4045789 -0.4209725  0.65262155  0.04614771 -0.58752872  0.57091599
    ## 5  0.1688144  0.8608158 -0.59079776  1.38370144 -0.04997384 -0.06456501
    ## 6  0.6571614 -1.3661150  0.85107548  2.08138152  0.38312579  1.15309274

The accessor-function return as expected matrices of estimates. For the plotting functions an additional index of the ith study or z-score is required.

``` r
traces(bc, burnin=FALSE, index=3)
```

![](bacon_files/figure-markdown_github/traces2-1.png)<!-- -->

``` r
posteriors(bc, index=3)
```

![](bacon_files/figure-markdown_github/traces2-2.png)<!-- -->

``` r
fit(bc, index=3, n=100)
```

![](bacon_files/figure-markdown_github/traces2-3.png)<!-- -->

``` r
plot(bc, type="hist", n=100, xlab="", ylab="")
```

![](bacon_files/figure-markdown_github/plothist2-1.png)<!-- -->

``` r
plot(bc, type="qq", xlab="expected", ylab="observed")
```

![](bacon_files/figure-markdown_github/qqplot2-1.png)<!-- -->

Fixed-effect meta-analysis
==========================

The following code chunk show how the neta-analysis P-values from bias and inflation corrected effect-sizes and standard errors can be add to the qqplot.

``` r
fm <- meta(bc)
hist(fm$pval)
```

![](bacon_files/figure-markdown_github/meta-1.png)<!-- -->

``` r
pl <- plot(bc, type="qq", pch=1, xlab="expected", ylab="observed", main="")
points(pl$x, sort(-log10(fm$pval)), col=1, pch=16)
```

![](bacon_files/figure-markdown_github/meta-2.png)<!-- -->

Session Info
============

    ## R Under development (unstable) (2016-02-24 r70217)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 14.04.3 LTS
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] bacon_0.99.0       ellipse_0.3-8      BiocParallel_1.4.3
    ## [4] rmarkdown_0.9.5    BiocStyle_1.8.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.9         futile.options_1.0.0 formatR_1.2.1       
    ##  [4] magrittr_1.5         evaluate_0.8         stringi_1.0-1       
    ##  [7] futile.logger_1.4.1  lambda.r_1.1.7       tools_3.3.0         
    ## [10] stringr_1.0.0        yaml_2.1.13          parallel_3.3.0      
    ## [13] htmltools_0.3        knitr_1.12.3

References
==========

Iterson, M van, E van Zwet, P Slagboom, and B.T Heijmans. 2016. “Controlling Inflation and Bias in Epigenome- and Transcriptome-Wide Association Studies Using Empirical Calibration.” *Manuscript in Preparation*.
