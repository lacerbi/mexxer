#include "mex.h"
#include "math.h"
#include "matrix.h"

/*
 * example0_mex.c
 *
 *EXAMPLE0 Return square of an array.
 *  Y = EXAMPLE0(X) takes one array X and computes its square.
 *
 * ================ INPUT VARIABLES ====================
 * X: array to be squared.    [M-by-N] (double)
 * 
 * ================ OUTPUT VARIABLES ==================
 * Y: squared array.          [M-by-N] (double)
 *
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 15-Aug-2016 with MEXXER v0.1 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

void example0( double *y, double *x, mwSize M, mwSize N )
{
	
	/* Write your main calculations here... */
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *y, *x;
	mwSize *dims_x;
	mwSize M, N;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<1 || nrhs>1 )
		mexErrMsgIdAndTxt( "MATLAB:example0:invalidNumInputs",
			"One inputs required.");
	if ( nlhs<1 || nlhs>1 )
		mexErrMsgIdAndTxt( "MATLAB:example0:invalidNumOutputs",
			"One outputs required.");

	/* Get first input (X, M-by-N double) */
	x = (double*) mxGetPr(prhs[0]);
	dims_x = (mwSize*) mxGetDimensions(prhs[0]);
	M = dims_x[0];
	N = dims_x[1];

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
	if ( ARGSCHECK ) {
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
				mexErrMsgIdAndTxt("MATLAB:example0:xNotReal", "Input X must be real.");
		if ( dims_x[0] != ((mwSize) (M)) )
			mexErrMsgIdAndTxt("MATLAB:example0:xWrongSize", "The first dimension of input X has the wrong size (should be M).");
		if ( dims_x[1] != ((mwSize) (N)) )
			mexErrMsgIdAndTxt("MATLAB:example0:xWrongSize", "The second dimension of input X has the wrong size (should be N).");
	}

	/* Pointer to first output (Y, M-by-N double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (M), (mwSize) (N), mxREAL);
	y = mxGetPr(plhs[0]);

	/* Call the C subroutine */
	example0(y, x, M, N);

}
