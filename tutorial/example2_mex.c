#include "mex.h"
#include "math.h"
#include "matrix.h"

/*
 * example2_mex.c
 *
 *EXAMPLE2 Compute log odds for Aspen Yoo's word recognition memory model.
 *
 * ================ INPUT VARIABLES ====================
 * M: number of features. [scalar] (integer)
 * SIGMA: memory noise. [scalar] (double)
 * NS: number of samples of SNew and SOld. [scalar] (integer)
 * NNEW: number of new words. [scalar] (integer)
 * NOLD: number of old words. [scalar] (integer)
 * SNEW: new words across S simulations. [Nnew*nS,M] (double)
 * X: noisy memories. [Nold,M] (double)
 * 
 * ================ OUTPUT VARIABLES ==================
 * D_NEW: log odds of new trials. [Nnew*nS,1] (double)
 *
 * This is a MEX-file for MATLAB.
 */

void example2( double *d, int M, double sigma, int nS, int Nold, int N, double *S, int Srows, double* X)
{
	
    mwSize i,j,k;
    double J,Jss,tmp;
    double *S0,*X0,*d0,sum,SUM;
        
    /* store initial position */
    S0 = S;     
    d0 = d;
        
    /* Compute d = M/2*log(1+1/sigma^2) + 0.5*sum(SNew.^2,2) */

    J = 1/(sigma*sigma) + 1;
    
    /* this summation is awkward because MATLAB is column-major */
    for (i=0; i<Srows; i++) {
        sum = 0.;
        *d = 0.5 * M * log(J);
        S = S0 + i;   /* initialize pointer at current row */
        for (j=0; j<M; j++) {
            sum += (*S)*(*S);
            S += Srows; /* move pointer to next column */
        }
        *(d++) += 0.5 * sum;
    }
        
/* Compute d = d + log(squeeze(mean(exp( ...
 *        -0.5*J*sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) ...
 *                    - repmat(X/J/sigma^2,[1,1,Nnew*nS])).^2,2) ...
 *        )))) */
    
    Jss = (1/J)/(sigma*sigma);
    
    X0 = X;
    d = d0;
    
    for (k=0; k < Srows; k++) {
        S = S0 + k;
        SUM = 0.;
        
        for (i=0; i < Nold; i++) {
            X = X0 + i;
            
            sum = 0.;            
            for (j=0; j < M; j++) {
                tmp = S[j*Srows] - (X[j*Nold] * Jss);
                sum += tmp * tmp;
            }
            SUM += exp(-0.5 * J * sum);
        }    
        *(d++) += log( SUM / (double) Nold );
        
    }
        
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *d_new, sigma, *SNew, *X;
	int M, nS, Nnew, Nold;

	/* Get first input (M, scalar int) */
	M = (int) mxGetScalar(prhs[0]);

	/* Get second input (SIGMA, scalar double) */
	sigma = (double) mxGetScalar(prhs[1]);

	/* Get third input (NS, scalar int) */
	nS = (int) mxGetScalar(prhs[2]);

	/* Get fourth input (NNEW, scalar int) */
	Nnew = (int) mxGetScalar(prhs[3]);

	/* Get fifth input (NOLD, scalar int) */
	Nold = (int) mxGetScalar(prhs[4]);

	/* Get sixth input (SNEW, Nnew*nS-by-M double) */
	SNew = (double*) mxGetPr(prhs[5]);

	/* Get seventh input (X, Nold-by-M double) */
	X = (double*) mxGetPr(prhs[6]);

	/* Pointer to first output (D_NEW, Nnew*nS-by-1 double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (Nnew*nS), (mwSize) (1), mxREAL);
	d_new = mxGetPr(plhs[0]);

	/* Call the C subroutine */
    example2(d_new, M, sigma, nS, Nnew, Nold, SNew, Nnew*nS, X);

}
