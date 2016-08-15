#include "mex.h"
#include "math.h"
#include "matrix.h"

/*
 * example1_mex.c
 *
 *EXAMPLE1 Outer product of two vectors and the sum along the 1st dimension.
 *  Z = EXAMPLE1(X,Y) takes two vectors X (column vector) and Y (row vector), 
 *  computes the outer product, that is Z_ij = X_i * Y_j, and returns the
 *  sum along the first dimension (a row vector Z_j).
 *
 *  We explicitly define input and outputs for the sake of MEXXER 
 *  (you can download MEXXER from here: https://github.com/lacerbi).
 *
 * ================ INPUT VARIABLES ====================
 * X: column vector.    [M-by-1] (double)
 * Y: row vector. [1-by-N] (double)
 * 
 * ================ OUTPUT VARIABLES ==================
 * Z: summed vector.    [1-by-N] (double)
 *
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 13-Aug-2016 with MEXXER v0.1 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

void example1( double *z, double *x, double *y, mwSize M, mwSize N )
{
	mwSize i,j;
    double sum, *x0;
    
    x0 = x;     /* store initial pointer position */
    
    for (i=0; i < M; i++) {
        sum = 0.;
        for (j=0; j < N; j++) {
            sum += x[i] * y[j];
        }
        z[j] = sum;
    }
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *z, *x, *y;
	mwSize *dims_x, *dims_y;
	mwSize M, N;

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<2 || nrhs>2 )
		mexErrMsgIdAndTxt( "MATLAB:example1:invalidNumInputs",
			"Two inputs required.");
	if ( nlhs<1 || nlhs>1 )
		mexErrMsgIdAndTxt( "MATLAB:example1:invalidNumOutputs",
			"One outputs required.");

	/* Get first input (X, M-by-1 double) */
	x = (double*) mxGetPr(prhs[0]);
	dims_x = (mwSize*) mxGetDimensions(prhs[0]);
	M = dims_x[0];

	/* Get second input (Y, 1-by-N double) */
	y = (double*) mxGetPr(prhs[1]);
	dims_y = (mwSize*) mxGetDimensions(prhs[1]);
	N = dims_y[1];

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
	if ( ARGSCHECK ) {
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
				mexErrMsgIdAndTxt("MATLAB:example1:xNotReal", "Input X must be real.");
		if ( dims_x[0] != ((mwSize) (M)) )
			mexErrMsgIdAndTxt("MATLAB:example1:xWrongSize", "The first dimension of input X has the wrong size (should be M).");
		if ( dims_x[1] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:example1:xWrongSize", "The second dimension of input X has the wrong size (should be 1).");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
				mexErrMsgIdAndTxt("MATLAB:example1:yNotReal", "Input Y must be real.");
		if ( dims_y[0] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:example1:yWrongSize", "The first dimension of input Y has the wrong size (should be 1).");
		if ( dims_y[1] != ((mwSize) (N)) )
			mexErrMsgIdAndTxt("MATLAB:example1:yWrongSize", "The second dimension of input Y has the wrong size (should be N).");
	}

	/* Pointer to first output (Z, 1-by-N double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (1), (mwSize) (N), mxREAL);
	z = mxGetPr(plhs[0]);

	/* Call the C subroutine */
	example1(z, x, y, M, N);

}
