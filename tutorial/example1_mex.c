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
 */

void example1( double *z, double *x, double *y, mwSize M, mwSize N )
{
    mwSize i,j;
    
    for (j=0; j<N; j++) {
        z[j] = 0.0;
        for (i=0; i<M; i++) {
            z[j] += x[i]*y[j];
        }
    }
    
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *z, *x, *y;
	mwSize *dims_x, *dims_y;
	mwSize M, N;
    
	/* Get first input (X, M-by-1 double) */
	x = (double*) mxGetPr(prhs[0]);
	dims_x = (mwSize*) mxGetDimensions(prhs[0]);
	M = dims_x[0];

	/* Get second input (Y, 1-by-N double) */
	y = (double*) mxGetPr(prhs[1]);
	dims_y = (mwSize*) mxGetDimensions(prhs[1]);
	N = dims_y[1];

	/* Pointer to first output (Z, 1-by-N double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (1), (mwSize) (N), mxREAL);
	z = mxGetPr(plhs[0]);

	/* Call the C subroutine */
	example1(z, x, y, M, N);

}
