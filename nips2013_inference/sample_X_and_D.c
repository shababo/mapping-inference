#include "mex.h"
#include <math.h>

/* sample uniform random variable */
const double rand_unif() {
    const int numerator = rand();
    return ((double)numerator)/((double)RAND_MAX);
}

/* sample gaussian random variable */
bool new_gaussian_sample_needed = true;
double next_gaussian_sample = 0;
double rand_gaussian() {
    if (!new_gaussian_sample_needed) {
        new_gaussian_sample_needed = true;
        return next_gaussian_sample;
    }
    
    double x, y, rds, c;

    /* rejection sample from a unit disk */
    do {
        x = rand_unif()*2-1;
        y = rand_unif()*2-1;
        rds = x*x + y*y;
    } while (rds == 0 || rds > 1);

    /* Box-Muller transform constant */
    c = sqrt(-2*log(rds)/rds);

    next_gaussian_sample = y*c;
    new_gaussian_sample_needed = false;
    
    return x*c;
}

/* compute alpha waveform */
const double alpha_waveform(const double t, const double offset, const double tau, const double gmax) {
    return (offset <= t) ? gmax*(t - offset) / tau * exp(-(t - offset - tau)/tau) : 0;
}

/* sample new X and D */
void sample_X_D(
    double *X_s, double *D_s, 
    const int N, const int K, const int T,
    double *Y, double *w_s, double *pi_nk, double *d_mean_nk, double *d_sigma_nk, double *t, double tau, double gmax, double sigma_n
) {

    int n,t_ix,k;
    double Y_n[T], Y_n_xk0[T], Y_n_xk1[T], Y_n_star[T];
    double alpha_tmp;

    double logR_num, logR_denom, odds_ratio;
    double d_nk_star, log_accept_ratio_num, log_accept_ratio_denom, accept_ratio;

    for (n=0; n<N; n++) {
        
        /* calculate baseline Y_n */
        for (t_ix=0; t_ix<T; t_ix++) {
            Y_n[t_ix] = 0;
            for (k=0; k<K; k++) {
                Y_n[t_ix] = Y_n[t_ix] + w_s[k]*X_s[n+N*k]*alpha_waveform(t[t_ix],D_s[n+N*k],tau,-gmax);
            }
        }
        
        for (k=0; k<K; k++) {

            logR_num = 0;
            logR_denom = 0; 
            for (t_ix=0; t_ix<T; t_ix++) {
                alpha_tmp = alpha_waveform(t[t_ix],D_s[n+N*k],tau,-gmax);
                Y_n[t_ix] = Y_n[t_ix] - w_s[k]*X_s[n+N*k]*alpha_tmp;
                Y_n_xk0[t_ix] = Y_n[t_ix];
                Y_n_xk1[t_ix] = Y_n[t_ix] + w_s[k]*alpha_tmp;
                logR_num += pow(Y_n_xk1[t_ix] - Y[n+N*t_ix], 2);
                logR_denom += pow(Y_n_xk0[t_ix] - Y[n+N*t_ix], 2);
            }
            logR_num = log(pi_nk[n+N*k]) - 0.5*logR_num/pow(sigma_n,2);
            logR_denom = log(1 - pi_nk[n+N*k]) - 0.5*logR_denom/pow(sigma_n,2);
            odds_ratio = exp(logR_num - logR_denom);

            if (odds_ratio == INFINITY || rand_unif() < odds_ratio/(1+odds_ratio)) {
                X_s[n+N*k] = 1;
                for (t_ix=0; t_ix<T; t_ix++) {
                    Y_n[t_ix] += w_s[k]*X_s[n+N*k]*alpha_waveform(t[t_ix],D_s[n+N*k],tau,-gmax);
                }
            } else {
                X_s[n+N*k] = 0;
            }
            
            /* get sample for D_s(n,k)  */
            d_nk_star = D_s[n+N*k] + rand_gaussian()*5;
            if (X_s[n+N*k] == 1 && w_s[k] != 0) {
                
                log_accept_ratio_num = 0; 
                log_accept_ratio_denom = 0; 
                for (t_ix=0; t_ix<T; t_ix++) {
                    Y_n_star[t_ix] = Y_n[t_ix] - w_s[k]*alpha_waveform(t[t_ix],D_s[n+N*k],tau,-gmax) + w_s[k]*alpha_waveform(t[t_ix],d_nk_star,tau,-gmax);
                    log_accept_ratio_num += pow(Y_n_star[t_ix]-Y[n+N*t_ix],2);
                    log_accept_ratio_denom += pow(Y_n[t_ix]-Y[n+N*t_ix],2);
                }
                log_accept_ratio_num = -.5*log_accept_ratio_num/pow(sigma_n,2) - .5*pow(d_nk_star - d_mean_nk[n+N*k],2)/pow(d_sigma_nk[n+N*k],2);
                log_accept_ratio_denom = -.5*log_accept_ratio_denom/pow(sigma_n,2) - .5*pow(D_s[n+N*k] - d_mean_nk[n+N*k],2)/pow(d_sigma_nk[n+N*k],2);
                accept_ratio = exp(log_accept_ratio_num - log_accept_ratio_denom);
            } else {
                accept_ratio = exp(-0.5/pow(d_sigma_nk[n+N*k],2)*(pow(d_nk_star - d_mean_nk[n+N*k],2) - pow(D_s[n+N*k] - d_mean_nk[n+N*k],2)));
            }
            
            if(rand_unif() < accept_ratio) {
                D_s[n+N*k] = d_nk_star;
            }
        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *Y,*w_s,*pi_nk,*d_mean_nk,*d_sigma_nk,*t,tau,gmax,sigma_n;

    double *X_s, *X_s_prev;
    double *D_s, *D_s_prev;

    /* Allocate memory for return values */
    mwSize ndim = mxGetNumberOfDimensions(prhs[1]);
    const mwSize *dims = mxGetDimensions(prhs[1]);
    plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);

    /* Extract matrix dimensions */
    const int N = (int)dims[0];
    const int K = (int)dims[1];
    const int T = (int)mxGetDimensions(prhs[0])[1];

    /* Assign pointers to each input and output */
    Y = mxGetPr(prhs[0]);
    X_s_prev = mxGetPr(prhs[1]);
    D_s_prev = mxGetPr(prhs[2]);
    w_s = mxGetPr(prhs[3]);
    pi_nk = mxGetPr(prhs[4]);
    d_mean_nk = mxGetPr(prhs[5]);
    d_sigma_nk = mxGetPr(prhs[6]);
    t = mxGetPr(prhs[7]);
    tau = mxGetScalar(prhs[8]);
    gmax = mxGetScalar(prhs[9]);
    sigma_n = mxGetScalar(prhs[10]);

    /* Create pointers to return arguments, and copy memory */
    X_s = mxGetPr(plhs[0]);
    D_s = mxGetPr(plhs[1]);
    int i;
    for (i=0; i<K*N; i++) {
        X_s[i] = (double)X_s_prev[i];
        D_s[i] = D_s_prev[i];
    }

    /* Call the subroutine. */
    sample_X_D(X_s,D_s,N,K,T,Y,w_s,pi_nk,d_mean_nk,d_sigma_nk,t,tau,gmax,sigma_n);
}
