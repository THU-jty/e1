const char* dgemm_desc = "Naive, three-loop dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */    
void square_dgemm (int n, double* A, double* B, double* C)
{

    for (int j = 0; j < n; ++j) 
    {
		for( int k = 0; k < n; k++ ){
			double tmp = B[k+j*n];
			for (int i = 0; i < n; ++i)	
				C[i+j*n] += A[i+k*n] * tmp;
		}
    }
}

// void square_dgemm (int n, double* A, double* B, double* C)
// {
  // for( int i = 0; i < n; i ++ ){
	  // for( int j = 0; j < i; j ++ ){
		  // double tmp = B[i+j*n];
		  // B[i+j*n] = B[j+i*n];
		  // B[j+i*n] = tmp;
	  // }
  // }
  // /* For each row i of A */
  // for (int i = 0; i < n; ++i)
    // /* For each column j of B */
    // for (int j = 0; j < n; ++j) 
    // {
      // /* Compute C(i,j) */
	  // double tmp = C[i+j*n];
      // for( int k = 0; k < n; k++ )
	// tmp += A[i+k*n] * B[k+j*n];
		// C[i+j*n] = tmp;
    // }
// }
