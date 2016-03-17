#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <R_ext/Print.h>

void wsmultinom(double *y, double *w, int n, double *mu, double *sig, double *p,
                double *nj, double *sj, double *sj2){

  double prob[3] = {0.0, 0.0, 0.0};
  double prob_tot = 0.;
  int rN[3];

  for(int j = 0; j < 3; j++)
    nj[j] = sj[j] = sj2[j] = 0.0;

  for(int i = 0; i < n; i++){
    
    prob_tot = 0.;
    for(int j = 0; j < 3; j++) {
      prob[j] = p[j]*dnorm(y[i], mu[j], sig[j], 0);
      prob_tot += prob[j];
    }
	
	//normalized probabilities
    for(int j = 0; j < 3; j++)
      prob[j] = prob[j]/prob_tot;

    rmultinom(w[i], prob, 3, rN);

    for(int j = 0; j < 3; j++) {
      nj[j] += (rN[j] == 0 ? 0.0 : rN[j]);
      sj[j] += (rN[j] == 0 ? 0.0 : rN[j]*y[i]);
      sj2[j] += (rN[j] == 0 ? 0.0 : rN[j]*y[i]*y[i]);
    }
  }
}

void fwrs(double *y, int n, double *mu, double *sig, double *p,
          double *nj, double *sj, double *sj2){

  //  Fast weighted random sampling according to: Efraimidisa and Spirakisb
  //  http://www.sciencedirect.com/science/article/pii/S002001900500298X#

  int l = 0;
  double s = 0.0;
  double maxval = 0.0;
  double val = 0.0;
  int z[n];
  double iweights[3] = {0.0, 0.0, 0.0};

  for(int i = 0; i < n; i++) {
    s = 0.0;
    for(int j = 0; j < 3; j++) {
      iweights[j] = p[j]*dnorm(y[i], mu[j], sig[j], 0);
      s += iweights[j];
    }

    for(int j = 0; j < 3; j++)
      iweights[j] = s/(iweights[j]+DBL_EPSILON);

    l = 0;
    maxval = pow(runif(0.0, 1.0), iweights[0]);
    for(int j = 1; j < 3; j++) {
      val = pow(runif(0.0, 1.0), iweights[j]);
      if(maxval < val) {
        maxval = val;
        l = j;
      }
    }

    z[i] = l;
  }

  for(int j = 0; j < 3; j++) {
    nj[j] = sj[j] = sj2[j] = 0.0;
    for(int i = 0; i < n; i++) {
      nj[j] += (z[i] == j ? 1.0 : 0.0);
      sj[j] += (z[i] == j ? y[i] : 0.0);
      sj2[j] += (z[i] == j ? y[i]*y[i] : 0.0);
    }
  }

}

void bacon(double *y, double *w, double *medyR, double *madyR, int *nR, int *niterR, double *levelR, int *binnedR, int *verboseR,
           double *gibbsmu, double *gibbssig, double *gibbsp,
           double *alpha, double *beta, double *lambda, double *tau, double *gamma){

  /*
    Fit three component normal mixture using Gibbs sampler of niter
    starting values use median and mad of y
    output values are full traces of the mixture parameters: mu, sigma, and proportion
  */

  int n = *nR;
  int niter = *niterR;
  int binned = *binnedR;
  int verbose = *verboseR;
  double medy = *medyR;
  double mady = *madyR;
  double level = *levelR;

  double s = 0.0;
  for(int i = 0; i < n; i++)
    s += w[i];

  double threshold = qnorm(1.0 -level/(2.0*s), 0, 1, 0, 0);

  /* starting values */

  double p[3] = {0.0, 0.0, 0.0};
  double mu[3] = {medy, medy-threshold*mady, medy+threshold*mady};
  double sig[3] = {mady, mady, mady};

  for(int i = 0; i < n; i++) {
    p[1] += (y[i] > medy - threshold*mady) ? w[i]/s : 0.0;
    p[2] += (y[i] < medy + threshold*mady) ? w[i]/s : 0.0;
  }
  p[0] = 1.0 - p[1] - p[2];

  /* working variables */

  double nj[3] = {0.0, 0.0, 0.0};
  double sj[3] = {0.0, 0.0, 0.0};
  double sj2[3] = {0.0, 0.0, 0.0};

  if(verbose) {
    Rprintf("threshold = %1.4f\n", threshold);
    Rprintf("Starting values:\np0 = %1.4f, p1 = %1.4f, p2 = %1.4f\n", p[0], p[1], p[2]);
    Rprintf("mu0 = %1.4f, mu1 = %1.4f, mu2 = %1.4f\n", mu[0], mu[1], mu[2]);
    Rprintf("sigma0 = %1.4f, sigma1 = %1.4f, sigma2 = %1.4f\n", sig[0], sig[1], sig[2]);
  }

  for(int k = 0; k < niter; k++) {

    if(binned)
      wsmultinom(y, w, n, mu, sig, p, nj, sj, sj2);
    else
      fwrs(y, n, mu, sig, p, nj, sj, sj2);

    for(int j = 0; j < 3; j++) {
      gibbsmu[3*k+j] = mu[j] = rnorm((lambda[j]*tau[j] + sj[j])/(nj[j] + tau[j]), sqrt(sig[j]*sig[j]/(tau[j] + nj[j])));
      sj2[j] = sj2[j] + mu[j]*(nj[j]*mu[j] - 2*sj[j]);
    }

    s = 0.0;
    for(int j = 0; j < 3; j++) {
      gibbssig[3*k+j] = sig[j] = sqrt(1.0/rgamma(*alpha + 0.5*(nj[j] + 1.0), 1.0/(*beta + 0.5*tau[j]*(mu[j] - lambda[j])*(mu[j] - lambda[j]) + 0.5*sj2[j])));
      gibbsp[3*k+j] = p[j] = rgamma(nj[j] + gamma[j], 1);
      s += p[j];
    }

    for(int j = 0; j < 3; j++)
      gibbsp[3*k+j] = p[j] = p[j]/s;
  }
}
