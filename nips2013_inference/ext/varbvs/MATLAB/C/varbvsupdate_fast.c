// For a description of this C code, see varbvsupdate.m.
#include "types.h"
#include "vectorops.h"
#include "doublevectormatlab.h"
#include "singlematrixmatlab.h"
#include "varbvs.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// MEX-file gateway routine. Note that varbvsupdate.m checks the
// inputs, so we do not have to do it here.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {
    
//  MATLAB call:
//     [alpha mu Xr] = varbvsupdate_fast(XtX,double(sigma_n_sq),...
//              double(sigma_s_sq),double(Xty),double(d),double(alpha),...
//              double(mu),double(Xr),double(I-1),double(logodds),X);
        
        

  // GET INPUTS.
  const SingleMatrix XtX        = getSingleMatrix(prhs[0]);
  const double       sigma_n_sq = *mxGetPr(prhs[1]);
  const double       sigma_s_sq = *mxGetPr(prhs[2]);
  const DoubleVector Xty        = getDoubleVector(prhs[3]);
  const DoubleVector d          = getDoubleVector(prhs[4]);
  const DoubleVector alpha0     = getDoubleVector(prhs[5]);
  const DoubleVector mu0        = getDoubleVector(prhs[6]);
  const DoubleVector I          = getDoubleVector(prhs[7]);
  const DoubleVector logodds    = getDoubleVector(prhs[8]);


  // Get the number of samples (n), the number of variables (p), and
  // the number of coordinate ascent updates (m).
  const Size p = XtX.nc;
  const Size m = I.n;

  // INITIALIZE OUTPUTS.
  DoubleVector alpha = createMatlabVector(p,&plhs[0]);
  DoubleVector mu    = createMatlabVector(p,&plhs[1]);

  copyDoubleVector(alpha0,alpha);
  copyDoubleVector(mu0,mu);

  // This is storage for a column of matrix X.
  double* XtXk = malloc(sizeof(double)*p);

  // RUN COORDINATE ASCENT UPDATES.
  // Repeat for each coordinate ascent update.
  for (Index j = 0; j < m; j++) {
    Index k = (Index) I.elems[j];

    // Copy the kth column of matrix X.
    copyColumn(XtX.elems,XtXk,k,p);

    // Perform the update.
    varbvsupdate_coord_fast(XtXk,Xty.elems[k],d.elems[k],sigma_n_sq,sigma_s_sq,logodds.elems[k],
		 alpha.elems,mu.elems,p,k);
  }

  // Free the dynamically allocated memory.
  free(XtXk);
}
