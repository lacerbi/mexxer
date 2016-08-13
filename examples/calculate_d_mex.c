#include "mex.h"
#include "math.h"
#include "matrix.h"

/*
 * calculate_d_mex.c
 *
 *% calculate log odds
 * function used as a basis for Luigi to code up C code!
 *
 * ================ INPUT VARIABLES ====================
 * M: number of features. [scalar] (integer)
 * SIGMA: memory noise. [scalar] (double)
 * NS: number of samples of SNew and SOld. [scalar] (integer)
 * NNEW: number of new words. [scalar] (integer)
 * NOLD: number of old words. [scalar] (integer)
 * SNEW: new words across S simulations. [Nnew*nS,M] (double)
 * SOLD: old words across S simulations. [Nold*nS,M] (double)
 * X: noisy memories. [Nold,M] (double)
 * 
 * ================ OUTPUT VARIABLES ==================
 * D_NEW: log odds of new trials. [Nnew*nS,1] (double)
 * D_OLD: log odds of old trials. [Nold*nS,1] (double)
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 13-Aug-2016 with MEXXER v0.1 
 * (https://github.com/lacerbi/mexxer).
*/

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

void calculate_d( double *d_new, double *d_old, int M, double sigma, int nS, int Nnew, int Nold, double *SNew, double *SOld, double *X )
{
	
	/* Write your main calculations here... */
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *d_new, *d_old, sigma, *SNew, *SOld, *X;
	int M, nS, Nnew, Nold;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<8 || nrhs>8 )
		mexErrMsgIdAndTxt( "MATLAB:calculate_d:invalidNumInputs",
			"Eight inputs required.");
	if ( nlhs<2 || nlhs>2 )
		mexErrMsgIdAndTxt( "MATLAB:calculate_d:invalidNumOutputs",
			"Two outputs required.");

	/* Get 1st input (M, scalar int) */
	M = (int) mxGetScalar(prhs[0]);

	/* Get 2nd input (SIGMA, scalar double) */
	sigma = (double) mxGetScalar(prhs[1]);

	/* Get 3rd input (NS, scalar int) */
	nS = (int) mxGetScalar(prhs[2]);

	/* Get 4th input (NNEW, scalar int) */
	Nnew = (int) mxGetScalar(prhs[3]);

	/* Get 5th input (NOLD, scalar int) */
	Nold = (int) mxGetScalar(prhs[4]);

	/* Get 6th input (SNEW, Nnew*nS-by-M double) */
	SNew = (double*) mxGetPr(prhs[5]);

	/* Get 7th input (SOLD, Nold*nS-by-M double) */
	SOld = (double*) mxGetPr(prhs[6]);

	/* Get 8th input (X, Nold-by-M double) */
	X = (double*) mxGetPr(prhs[7]);

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip this part) */
	if ( ARGSCHECK ) {
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || (mxGetN(prhs[0])*mxGetM(prhs[0])!=1) )
			mexErrMsgIdAndTxt("MATLAB:calculate_d:MNotScalar", "Input M must be a scalar.");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetN(prhs[1])*mxGetM(prhs[1])!=1) )
			mexErrMsgIdAndTxt("MATLAB:calculate_d:sigmaNotScalar", "Input SIGMA must be a scalar.");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || (mxGetN(prhs[2])*mxGetM(prhs[2])!=1) )
			mexErrMsgIdAndTxt("MATLAB:calculate_d:nSNotScalar", "Input NS must be a scalar.");

		if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || (mxGetN(prhs[3])*mxGetM(prhs[3])!=1) )
			mexErrMsgIdAndTxt("MATLAB:calculate_d:NnewNotScalar", "Input NNEW must be a scalar.");

		if ( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || (mxGetN(prhs[4])*mxGetM(prhs[4])!=1) )
			mexErrMsgIdAndTxt("MATLAB:calculate_d:NoldNotScalar", "Input NOLD must be a scalar.");
	}

	/* Pointer to 1st output (D_NEW, Nnew*nS-by-1 double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (Nnew*nS), (mwSize) (1), mxREAL);
	d_new = mxGetPr(plhs[0]);

	/* Pointer to 2nd output (D_OLD, Nold*nS-by-1 double) */
	plhs[1] = mxCreateDoubleMatrix((mwSize) (Nold*nS), (mwSize) (1), mxREAL);
	d_old = mxGetPr(plhs[1]);

	/* Call the C subroutine */
	calculate_d(d_new, d_old, M, sigma, nS, Nnew, Nold, SNew, SOld, X);

}
