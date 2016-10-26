#include "mex.h"
#include <math.h>

/* compute alpha waveform */
const double alpha_waveform(const double t, const double offset, const double tau, const double gmax) {
    return (offset <= t) ? gmax*(t - offset) / tau * exp(-(t - offset - tau)/tau) : 0;
}

/* compute M and r */
void compute_M_and_r(
        double *M, double *r, 
        const int N, const int K, const int T,
        const double * const Y, const bool * const X_s, const double * const D_s, const double * const t, 
        const double tau, const double gmax, const double sigma_n, const double * const sign_w) {
    int n,j,k,t_ix;
    double alpha_k_Y_n, alpha_kk, alpha_jk, alpha_k_tmp;
    
    for (n=0; n<N; n++) {
        for (k=0; k<K; k++) {
            if (X_s[n+N*k] == 1) {
                for (j=k; j<K; j++) {

                    alpha_k_Y_n = 0;
                    alpha_jk = 0;

                    for (t_ix=0; t_ix<T; t_ix++) {
                        alpha_k_tmp = alpha_waveform(t[t_ix],D_s[n+N*k],tau,-gmax*sign_w[k]);
                        alpha_jk += alpha_k_tmp*alpha_waveform(t[t_ix], D_s[n+N*j], tau, -gmax*sign_w[k]);
                        if (j==k) {
                            alpha_k_Y_n += alpha_k_tmp*Y[n+N*t_ix];
                        }
                    }
                    
                    if (j==k) {
                        r[k] += X_s[n+N*k]*alpha_k_Y_n/pow(sigma_n,2);
                    }
                    
                    M[k+K*j] += X_s[n+N*k]*X_s[n+N*j]*alpha_jk/pow(sigma_n, 2);
                    M[j+K*k] = M[k+K*j];
                }
            }
        }
    }

    /* alternate loop formulation: not any better. */
    /*
        for (n=0; n<N; n++) {
            for (k=0; k<K; k++) {
                alpha_kk = 0;
                alpha_k_Y_n = 0;
                for (t_ix=0; t_ix<T; t_ix++) {
                    alpha_k_tmp = alpha_waveform(t[t_ix],D_s[n+N*k],tau,-gmax);
                    alpha_k_Y_n += alpha_k_tmp*Y[n+N*t_ix];
                    alpha_kk += pow(alpha_k_tmp,2);
                }
                r[k] += X_s[n+N*k]*alpha_k_Y_n/pow(sigma_n,2);
                M[k+K*k] += pow(X_s[n+N*k],2)*alpha_kk/pow(sigma_n, 2);
                 
                for (j=k+1; j<K; j++) {
                    alpha_jk = 0;
                    for (t_ix=0; t_ix<T; t_ix++) {
                        alpha_jk += alpha_waveform(t[t_ix], D_s[n+N*k], tau, -gmax)*alpha_waveform(t[t_ix], D_s[n+N*j], tau, -gmax);
                    }
                    M[k+K*j] += X_s[n+N*k]*X_s[n+N*j]*alpha_jk/pow(sigma_n, 2);
                    M[j+K*k] = M[k+K*j];
                }
            }
        }
    */
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    /* Extract matrix dimensions */
    const mwSize *dims = mxGetDimensions(prhs[1]);
    const int N = (int)dims[0];
    const int K = (int)dims[1];
    const int T = (int)mxGetDimensions(prhs[0])[1];

    double *M,*r,*Y,*D_s,*t,tau,gmax,sigma_n,*sign_w;
    bool X_s[N*K];
    mxLogical *X_s_logical;

    /* Assign pointers to each input and output */
    Y = mxGetPr(prhs[0]);
    X_s_logical = (mxLogical*)mxGetData(prhs[1]);
    D_s = mxGetPr(prhs[2]);
    t = mxGetPr(prhs[3]);
    tau = mxGetScalar(prhs[4]);
    gmax = mxGetScalar(prhs[5]);
    sigma_n = mxGetScalar(prhs[6]);
    sign_w = mxGetPr(prhs[7]);
    
    /* Allocate memory for return values */
    int dims_KK[2] = { K, K };
    int dims_K1[2] = { K, 1 };
    plhs[0] = mxCreateNumericArray(2, (mwSize*)dims_KK, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, (mwSize*)dims_K1, mxDOUBLE_CLASS, mxREAL);

    /* Create pointers for M and r, and zero them out */
    M = mxGetPr(plhs[0]);
    r = mxGetPr(plhs[1]);
    int i;
    for (i=0; i<K; i++) { r[i] = 0; }
    for (i=0; i<(K*K); i++) { M[i] = 0; }

    /* Cast X_s */
    for (i=0; i<N*K; i++) { X_s[i] = X_s_logical[i] > 0; }

    /* Do something */
    compute_M_and_r(M,r,N,K,T,Y,X_s,D_s,t,tau,gmax,sigma_n,sign_w);
}
