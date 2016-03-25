Controlling bias and inflation in association studies using the empirical null distribution
================

    ## 
    ## Attaching package: 'BiocStyle'

    ## The following objects are masked from 'package:rmarkdown':
    ## 
    ##     html_document, md_document, output_format, pdf_document

<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").className = "title";
});
</script>
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>
**Package**: *[bacon](http://bioconductor.org/packages/bacon)*<br /> **Authors**: Maarten van Iterson \[aut, cre\], Erik van Zwet \[ctb\]<br /> **Modified**: Fri Mar 25 11:54:05 2016<br /> **Compiled**: Fri Mar 25 12:53:56 2016

Introduction
============

*[bacon](http://bioconductor.org/packages/bacon)* can be used to remove inflation and bias often observed in epigenome- and transcriptome-wide association studies(Iterson et al. 2016).

To this end *[bacon](http://bioconductor.org/packages/bacon)* constructs an empirical null distribution using a Gibbs Sampling algorithm by fitting a three-component normal mixture on z-scores. One component is forced, using prior knowledge, to represent the null distribution with mean and standard deviation representing the bias and inflation. The other two components are necessary to capture the amount of true associations present in the data, which we assume unknown but small.

*[bacon](http://bioconductor.org/packages/bacon)* provides functionality to inspect the output of the Gibbs Sampling algorithm, i.e., traces, posterior distributions and the mixture fit, are provided. Furthermore, inflation- and bias-corrected test-statistics or P-values are extracted easily. In addition, functionality for performing fixed-effect meta-analysis are provided.

The function `bacon` requires a vector or a matrix of z-scores, e.g., those extracted from association analyses using a linear regression approach. For fixed-effect meta-analysis a matrix of effect-sizes and standard errors is required.

The vignette illustrates the use of *[bacon](http://bioconductor.org/packages/bacon)* using simulated z-scores, effect-sizes and standard errors to avoid long runtimes. If multiple sets of test-statisics or effect-sizes and standard errors are provided, the Gibbs Sampler Algorithm can be executed on multiple nodes to reduce computation time using functionality provide by *[BiocParallel](http://bioconductor.org/packages/BiocParallel)*-package.

One set of test-statistics
==========================

A vector containing 2000 z-scores is generated from a normal mixture distribution 90% of the z-scores were drawn from \(\mathcal{N}(0,1)\) and the remaining z-scores from \(\mathcal{N}(\mu,1)\), where \(\mu \sim \mathcal{N}(4,1)\).

The function `bacon` executes the Gibbs Sampler Algorithm and stores all input and output in a `Bacon`-object. Several accessor-functions are available to access data contained in the `Bacon`-object.

``` r
set.seed(12345)
y <- rnormmix(2000, c(0.9, 0, 1, 0, 4, 1))
bc <- bacon(y)
bc
```

    ## Bacon-object containing 1 set(s) of 2000 test-statistics.
    ## ...estimated bias: -0.0082.
    ## ...estimated inflation: 0.97.
    ## 
    ## Emprical null estimates are based on 5000 iterations with a burnin-period of 2000.

``` r
estimates(bc)
```

    ##            p.0        p.1        p.2         mu.0    mu.1      mu.2
    ## [1,] 0.9063948 0.04170951 0.05189565 -0.008240754 2.96292 -3.043919
    ##        sigma.0  sigma.1  sigma.2
    ## [1,] 0.9740429 2.810962 3.105406

``` r
inflation(bc)
```

    ##   sigma.0 
    ## 0.9740429

``` r
bias(bc)
```

    ##         mu.0 
    ## -0.008240754

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
plot(bc, type="hist")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 40 rows containing non-finite values (stat_bin).

    ## Warning: Removed 1 rows containing missing values (geom_bar).

    ## Warning: Removed 40 rows containing missing values (geom_path).

    ## Warning: Removed 40 rows containing missing values (geom_path).

![](bacon_files/figure-markdown_github/plothist1-1.png)<!-- -->

``` r
plot(bc, type="qq")
```

![](bacon_files/figure-markdown_github/plothist1-2.png)<!-- -->

Multiple sets of test-statistics
================================

Matrices containing \(2000\times6\) effect-sizes and standard errors are generated to simulated data for a fixed-effect meta-analyses.

By default the function `bacon` detects the number of cores/nodes registered, as described in the *[BiocParallel](http://bioconductor.org/packages/BiocParallel)*, to perform bacon in parallel. To run the vignette in general we set it here for convenience to 1 node.

``` r
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
    ## ...estimated bias: -0.022,-0.025,-0.0014,0.028,0.0038,-0.0099.
    ## ...estimated inflation: 1.1,1.1,1.1,1.1,1,1.1.
    ## 
    ## Emprical null estimates are based on 5000 iterations with a burnin-period of 2000.

``` r
estimates(bc)
```

    ##         p.0        p.1        p.2         mu.0     mu.1      mu.2  sigma.0
    ## A 0.8874249 0.02731197 0.08526312 -0.021866990 2.958835 -2.077888 1.114897
    ## B 0.8792706 0.02904402 0.09168536 -0.024518736 3.000491 -2.184858 1.091183
    ## C 0.8810674 0.05296603 0.06596655 -0.001382182 2.812989 -2.789718 1.085869
    ## D 0.8931000 0.07631524 0.03058475  0.028493259 1.797394 -2.996150 1.133673
    ## E 0.8668525 0.09937511 0.03377243  0.003803909 1.461691 -2.949467 1.043746
    ## F 0.8631620 0.07754931 0.05928872 -0.009896142 2.611823 -2.679531 1.054271
    ##    sigma.1  sigma.2
    ## A 1.384168 5.973949
    ## B 1.201077 5.274064
    ## C 3.913841 3.890048
    ## D 5.442560 0.680674
    ## E 5.294036 1.078277
    ## F 4.503480 3.759275

``` r
inflation(bc)
```

    ##        A        B        C        D        E        F 
    ## 1.114897 1.091183 1.085869 1.133673 1.043746 1.054271

``` r
bias(bc)
```

    ##            A            B            C            D            E 
    ## -0.021866990 -0.024518736 -0.001382182  0.028493259  0.003803909 
    ##            F 
    ## -0.009896142

``` r
head(tstat(bc))
```

    ##            A          B           C           D          E           F
    ## 1  1.6757419 -1.0828059 -1.25202621 -3.51468705  0.2118737 -0.04869098
    ## 2 -0.9614628 -1.3760475  0.26613256  0.24893989  0.3159087  0.34390169
    ## 3 -2.3816096 -3.5801198  0.06233318 -0.87573314 -0.1666154  0.05571797
    ## 4 -0.6548216 -0.5538327  0.53022692  0.03030979 -0.4649469  1.00660348
    ## 5  0.2444309  1.2611968 -0.89085413  1.00531956 -0.0564406 -0.08588204
    ## 6  0.9456110 -2.1858898  0.97625398  1.05816803  0.3533617  1.82825320

``` r
head(pval(bc))
```

    ##            A            B         C            D         E          F
    ## 1 0.09378876 0.2788945974 0.2105603 0.0004402725 0.8322056 0.96116557
    ## 2 0.33631952 0.1688069416 0.7901371 0.8034072773 0.7520718 0.73092022
    ## 3 0.01723716 0.0003434367 0.9502975 0.3811751256 0.8676727 0.95556648
    ## 4 0.51258261 0.5796933203 0.5959546 0.9758199920 0.6419695 0.31412536
    ## 5 0.80689715 0.2072379541 0.3730074 0.3147430018 0.9549908 0.93156019
    ## 6 0.34434703 0.0288236674 0.3289386 0.2899788451 0.7238173 0.06751156

``` r
head(se(bc))
```

    ##           A         B         C         D         E         F
    ## 1 0.9422218 1.3892419 1.2670153 0.8935741 0.9407586 0.5778809
    ## 2 0.5699857 0.6172674 0.9632422 2.2727758 1.2952927 1.3637669
    ## 3 2.3061567 0.5072609 0.9418045 1.2418620 0.7459992 0.8613406
    ## 4 0.6115212 0.7581003 1.2330106 1.5060923 1.2710527 0.5677095
    ## 5 0.7103241 0.6833334 0.6624854 1.3759268 0.9300647 0.7435161
    ## 6 0.6999728 0.6245508 0.8726132 1.9663520 1.0759823 0.6310371

``` r
head(es(bc))
```

    ##            A          B           C           D           E           F
    ## 1  1.5789205 -1.5042793 -1.58633631 -3.14063331  0.19932203 -0.02813758
    ## 2 -0.5480200 -0.8493892  0.25635011  0.56578457  0.40919429  0.46900175
    ## 3 -5.4923650 -1.8160547  0.05870567 -1.08753972 -0.12429497  0.04799215
    ## 4 -0.4004373 -0.4198608  0.65377540  0.04564933 -0.59097202  0.57145836
    ## 5  0.1736251  0.8618178 -0.59017781  1.38324614 -0.05249341 -0.06385468
    ## 6  0.6619020 -1.3651992  0.85189207  2.08073085  0.38021093  1.15369562

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
plot(bc, type="hist")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 240 rows containing non-finite values (stat_bin).

    ## Warning: Removed 6 rows containing missing values (geom_bar).

    ## Warning: Removed 43 rows containing missing values (geom_path).

    ## Warning: Removed 43 rows containing missing values (geom_path).

![](bacon_files/figure-markdown_github/plothist2-1.png)<!-- -->

``` r
plot(bc, type="qq")
```

![](bacon_files/figure-markdown_github/qqplot2-1.png)<!-- -->

Fixed-effect meta-analysis
==========================

The following code chunk shows how to perform fixed-effect meta-analysis and the inspection of results.

``` r
bcm <- meta(bc)
head(pval(bcm))
```

    ##             A            B         C            D         E          F
    ## 1 0.064832387 2.277949e-01 0.1735394 7.621167e-05 0.8220210 0.95117635
    ## 2 0.274043274 1.270002e-01 0.7736507 7.560214e-01 0.7387323 0.72433638
    ## 3 0.007425897 8.456294e-05 0.9471362 3.348949e-01 0.8649312 0.96104228
    ## 4 0.452095884 5.294460e-01 0.5657142 9.498823e-01 0.6301734 0.29310415
    ## 5 0.802086138 1.764782e-01 0.3326783 2.427275e-01 0.9560542 0.91999577
    ## 6 0.301888697 1.596452e-02 0.2897356 2.194059e-01 0.7094284 0.05516461
    ##          meta
    ## 1 0.171408974
    ## 2 0.286601859
    ## 3 0.003041214
    ## 4 0.985273782
    ## 5 0.656391353
    ## 6 0.376440192

``` r
plot(bcm, type="qq")
```

![](bacon_files/figure-markdown_github/meta-1.png)<!-- -->

``` r
(topTable(bcm))
```

    ##      eff.size.meta std.err.meta pval.adj.meta pval.org.meta tstat.meta
    ## 348      -5.374985    0.2963590  3.269017e-70  1.634508e-73 -18.136737
    ## 416      -4.172012    0.3205998  2.060282e-35  1.030141e-38 -13.013145
    ## 846       3.450712    0.2717292  1.196514e-33  5.982569e-37  12.699084
    ## 1512     -3.250383    0.2714998  9.968262e-30  4.984131e-33 -11.971953
    ## 251      -3.986363    0.3632577  1.020179e-24  5.100893e-28 -10.973925
    ## 298      -3.059337    0.2906307  1.303643e-22  6.518217e-26 -10.526546
    ## 1471      2.950354    0.3013521  2.476434e-19  1.238217e-22   9.790386
    ## 1691     -3.180895    0.3264487  3.916377e-19  1.958188e-22  -9.743936
    ## 1890     -3.110623    0.3302108  9.012161e-18  4.506081e-21  -9.420112
    ## 1350      3.040821    0.3255549  1.918992e-17  9.594958e-21   9.340423
    ##       eff.size.A std.err.A        pval.A      tstat.A   eff.size.B
    ## 348  -11.0847513 0.3839902 3.076025e-183 -28.86727342 -0.002347721
    ## 416   -1.4946442 0.7346860  4.191140e-02  -2.03439882 -0.713180529
    ## 846    5.0336849 0.3291048  8.245950e-53  15.29507995 -0.291598513
    ## 1512  -7.2369898 0.3770892  4.340610e-82 -19.19171675  0.305888081
    ## 251   -0.4157304 0.8258287  6.146761e-01  -0.50340991  1.394683042
    ## 298    0.2092830 0.7198036  7.712424e-01   0.29075021 -6.847817543
    ## 1471   0.9123215 1.0002536  3.617212e-01   0.91209023  1.394658968
    ## 1691  -9.5521842 0.5086909  1.143625e-78 -18.77797337 -0.685003121
    ## 1890   0.5313915 0.5741804  3.547171e-01   0.92547839 -2.639799193
    ## 1350  -0.0253185 1.1099642  9.818017e-01  -0.02281019 -0.674084505
    ##      std.err.B       pval.B      tstat.B eff.size.C std.err.C       pval.C
    ## 348  0.6217137 9.969870e-01  -0.00377621 -0.0924156 1.0514122 9.299588e-01
    ## 416  0.6071969 2.401765e-01  -1.17454569 -1.4933336 1.0967759 1.733346e-01
    ## 846  0.6321306 6.445872e-01  -0.46129470  7.8275355 0.8878722 1.185568e-18
    ## 1512 1.1346987 7.874861e-01   0.26957648 -0.4332764 1.6394504 7.915630e-01
    ## 251  1.0739642 1.940706e-01   1.29863084 -8.6737136 0.4865015 4.227732e-71
    ## 298  0.4595883 3.301036e-50 -14.89989384 -1.2122887 1.3571105 3.717037e-01
    ## 1471 0.6142030 2.316631e-02   2.27068073 -0.1177449 0.8191734 8.857088e-01
    ## 1691 0.9903737 4.891501e-01  -0.69166126  0.3424911 0.9464873 7.174604e-01
    ## 1890 1.3692436 5.386444e-02  -1.92792514 -0.6460778 0.7607509 3.957348e-01
    ## 1350 0.5802341 2.453388e-01  -1.16174577  4.1729192 1.1888716 4.481354e-04
    ##           tstat.C  eff.size.D std.err.D    pval.D     tstat.D eff.size.E
    ## 348   -0.08789664  0.58542791 0.7697885 0.4469529  0.76050484 -1.6929514
    ## 416   -1.36156684 -1.32193935 0.9029136 0.1431716 -1.46408182 -9.1117195
    ## 846    8.81606106  0.03076455 0.8011290 0.9693676  0.03840150  2.1437613
    ## 1512  -0.26428149 -0.42489336 0.4689393 0.3648971 -0.90607327  0.3541021
    ## 251  -17.82874930  0.34576929 0.8080852 0.6687333  0.42788715  0.7894370
    ## 298   -0.89328665  0.54502896 0.8285481 0.5106589  0.65781209 -2.2480981
    ## 1471  -0.14373623 -0.20322159 0.7001583 0.7716243 -0.29025090  7.0540388
    ## 1691   0.36185498  0.42698157 0.9247673 0.6442837  0.46171785 -0.5142322
    ## 1890  -0.84926327 -0.17489012 0.8228693 0.8316882 -0.21253693 -1.6622526
    ## 1350   3.50998295  0.08360362 0.8805529 0.9243589  0.09494446  7.0010324
    ##      std.err.E       pval.E     tstat.E  eff.size.F std.err.F       pval.F
    ## 348  0.8750208 5.302028e-02  -1.9347557  1.18337261 1.3723951 3.885400e-01
    ## 416  0.4919889 1.419758e-76 -18.5201746  0.19057584 1.3295419 8.860221e-01
    ## 846  1.2190856 7.866257e-02   1.7584994  0.09806573 0.9011307 9.133411e-01
    ## 1512 0.9157899 6.990057e-01   0.3866630 -0.92131675 0.6243603 1.400467e-01
    ## 251  1.2763918 5.362516e-01   0.6184911 -1.59444026 1.6455646 3.325791e-01
    ## 298  1.1915960 5.921040e-02  -1.8866278 -1.96129684 0.4701802 3.027700e-05
    ## 1471 0.4907888 7.659986e-47  14.3728614  1.67861541 0.8257847 4.207759e-02
    ## 1691 0.5621529 3.603202e-01  -0.9147551  1.80041635 1.1701208 1.238877e-01
    ## 1890 1.2262967 1.752564e-01  -1.3555061 -8.67048813 0.5424235 1.632961e-57
    ## 1350 0.4756816 4.948163e-49  14.7178961 -1.00026586 1.1975375 4.035663e-01
    ##          tstat.F
    ## 348    0.8622682
    ## 416    0.1433395
    ## 846    0.1088252
    ## 1512  -1.4756172
    ## 251   -0.9689320
    ## 298   -4.1713728
    ## 1471   2.0327520
    ## 1691   1.5386585
    ## 1890 -15.9847211
    ## 1350  -0.8352689

Session Info
============

    ## R Under development (unstable) (2016-03-21 r70361)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 14.04.4 LTS
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
    ## [1] BiocStyle_1.9.5     bacon_0.99.5        ellipse_0.3-8      
    ## [4] BiocParallel_1.5.20 ggplot2_2.1.0       rmarkdown_0.9.5    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.3      digest_0.6.9     plyr_1.8.3       grid_3.4.0      
    ##  [5] gtable_0.2.0     formatR_1.3      magrittr_1.5     evaluate_0.8.3  
    ##  [9] scales_0.4.0     stringi_1.0-1    labeling_0.3     tools_3.4.0     
    ## [13] stringr_1.0.0    munsell_0.4.3    parallel_3.4.0   yaml_2.1.13     
    ## [17] colorspace_1.2-6 htmltools_0.3.5  knitr_1.12.3

References
==========

Iterson, M van, E van Zwet, P Slagboom, and B.T Heijmans. 2016. “Controlling Bias and Inflation in Association Studies Using the Empirical Null Distribution.” *Manuscript in Preparation*.
