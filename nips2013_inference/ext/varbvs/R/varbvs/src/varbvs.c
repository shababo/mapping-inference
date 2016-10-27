#include "varbvs.h"
#include "sigmoid.h"
#include "vectorops.h"
#include <math.h>

// Function definitions.
// -----------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the
// variational lower bound for Bayesian variable selection in linear
// regression.
void varbvsupdate (const double* x, double xy, double d, double sigma, 
		   double sa, double logodds, double* alpha, double* mu, 
		   double* Xr, Size n) {

  // Compute the variational estimate of the posterior variance.
  double s = sa*sigma/(sa*d + 1);
  
  // Update the variational estimate of the posterior mean.
  double r = (*alpha) * (*mu);
  *mu = s/sigma * (xy + d*r - dot(x,Xr,n));
  
  // Update the variational estimate of the posterior inclusion
  // probability.
  double SSR = (*mu) * (*mu) / s;
  *alpha = sigmoid(logodds + (log(s/(sa*sigma)) + SSR)/2);
  
  // Update Xr = X*r.
  double rnew = (*alpha) * (*mu);
  add(Xr,rnew - r,x,n);
}

    
// Execute a single coordinate ascent update to maximize the
// variational lower bound for Bayesian variable selection in linear
// regression but this does it without a matrix multiply that is linear in T
void varbvsupdate_coord_fast (const double* XtXk, double xyk, double dk, double sigma_n_sq, 
		   double sigma_s_sq, double logoddsk, double* alpha, double* mu, 
		   Size p, Index k) {

  // Compute the variational estimate of the posterior variance.
  double s_sq = sigma_s_sq*sigma_n_sq/(sigma_s_sq*dk + 1);
  
  // Update the variational estimate of the posterior mean.
  
  double sum_x_alpha_mu = 0;
  
  for (Index i = 0; i < p; i++)
      sum_x_alpha_mu += XtXk[i] * alpha[i] * mu[i]; 
  
  mu[k] = s_sq/sigma_n_sq * (xyk + dk*alpha[k]*mu[k] - sum_x_alpha_mu);
  
  // Update the variational estimate of the posterior inclusion
  // probability.
  double SSR = mu[k] * mu[k] / s_sq;
  alpha[k] = sigmoid(logoddsk + (log(s_sq/(sigma_s_sq*sigma_n_sq)) + SSR)/2);
  
}

// Execute a single coordinate ascent update to maximize the
// variational lower bound for Bayesian variable selection in linear
// regression with nonzero prior slab mean.
void varbvsupdate_general (const double* x, double xy, double d, double sigma, 
		   double sa, double logodds, double* alpha, double* mu, 
		   double* Xr, Size n, double eta) {

  // Compute the variational estimate of the posterior variance.
  double s = sa*sigma/(sa*d + sigma);
  
  // Update the variational estimate of the posterior mean.
  double r = (*alpha) * (*mu);
  *mu = s*(eta/sa + 1/sigma * (xy + d*r - dot(x,Xr,n)));
  
  
  // Update the variational estimate of the posterior inclusion
  // probability.
  double SSR = (*mu) * (*mu) / s - eta*eta/sa;
  *alpha = sigmoid(logodds + (log(s/sa) + SSR)/2);
  
  // Update Xr = X*r.
  double rnew = (*alpha) * (*mu);
  add(Xr,rnew - r,x,n);
}