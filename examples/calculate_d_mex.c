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

	/* Get seventh input (SOLD, Nold*nS-by-M double) */
	SOld = (double*) mxGetPr(prhs[6]);

	/* Get eighth input (X, Nold-by-M double) */
	X = (double*) mxGetPr(prhs[7]);

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
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

		dims_SNew = (mwSize*) mxGetDimensions(prhs[5]);
		if ( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) )
				mexErrMsgIdAndTxt("MATLAB:calculate_d:SNewNotReal", "Input SNEW must be real.");
		if ( dims_SNew[0] != ((mwSize) (Nnew*nS)) )
			mexErrMsgIdAndTxt("MATLAB:calculate_d:SNewWrongSize", "The first dimension of input SNEW has the wrong size (should be Nnew*nS).");
		if ( dims_SNew[1] != ((mwSize) (M)) )
			mexErrMsgIdAndTxt("MATLAB:calculate_d:SNewWrongSize", "The second dimension of input SNEW has the wrong size (should be M).");

		dims_SOld = (mwSize*) mxGetDimensions(prhs[6]);
		if ( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) )
				mexErrMsgIdAndTxt("MATLAB:calculate_d:SOldNotReal", "Input SOLD must be real.");
		if ( dims_SOld[0] != ((mwSize) (Nold*nS)) )
			mexErrMsgIdAndTxt("MATLAB:calculate_d:SOldWrongSize", "The first dimension of input SOLD has the wrong size (should be Nold*nS).");
		if ( dims_SOld[1] != ((mwSize) (M)) )
			mexErrMsgIdAndTxt("MATLAB:calculate_d:SOldWrongSize", "The second dimension of input SOLD has the wrong size (should be M).");

		dims_X = (mwSize*) mxGetDimensions(prhs[7]);
		if ( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) )
				mexErrMsgIdAndTxt("MATLAB:calculate_d:XNotReal", "Input X must be real.");
		if ( dims_X[0] != ((mwSize) (Nold)) )
			mexErrMsgIdAndTxt("MATLAB:calculate_d:XWrongSize", "The first dimension of input X has the wrong size (should be Nold).");
		if ( dims_X[1] != ((mwSize) (M)) )
			mexErrMsgIdAndTxt("MATLAB:calculate_d:XWrongSize", "The second dimension of input X has the wrong size (should be M).");
	}

	/* Pointer to first output (D_NEW, Nnew*nS-by-1 double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (Nnew*nS), (mwSize) (1), mxREAL);
	d_new = mxGetPr(plhs[0]);

	/* Pointer to second output (D_OLD, Nold*nS-by-1 double) */
	plhs[1] = mxCreateDoubleMatrix((mwSize) (Nold*nS), (mwSize) (1), mxREAL);
	d_old = mxGetPr(plhs[1]);

	/* Call the C subroutine */
	calculate_d(d_new, d_old, M, sigma, nS, Nnew, Nold, SNew, SOld, X);

}
