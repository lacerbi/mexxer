#include "mex.h"

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
 */

void example0( double *y, double *x, mwSize M, mwSize N )
{
	mwSize i,j;
    
    for ( i=0; i < M; i++ ) {
       for ( j=0; j < N; j++ ) {
            y[i + j*M] = x[i + j*M]*x[i + j*M];
        }
    }	
	    
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *y, *x;
	mwSize *dims_x;
	mwSize M, N;

	/* Get input (X, M-by-N double) */
	x = (double*) mxGetPr(prhs[0]);
	dims_x = (mwSize*) mxGetDimensions(prhs[0]);
	M = dims_x[0];
	N = dims_x[1];
    
	/* Pointer to output (Y, M-by-N double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (M), (mwSize) (N), mxREAL);
	y = mxGetPr(plhs[0]);

	/* Call the C subroutine */
	example0(y, x, M, N);

}