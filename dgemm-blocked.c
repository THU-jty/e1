const char* dgemm_desc = "kji";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 65
#endif

#include<stdlib.h>
#define min(a,b) (((a)<(b))?(a):(b))
#define rep(j) for(int j = 0; j < lda; j += BLOCK_SIZE)
#define repj for (int j = 0; j < N; ++j)
#define repi for (int i = 0; i < M; ++i)
#define repk for (int k = 0; k < K; ++k)
	
/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int M, int N, int K, double* A, double* B, double* C)
{
 repk
	repj{
		double tmp = B[k+j*K];
		repi{
			C[i+j*lda] += A[i+k*lda]*tmp;
		}
	}
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
// void square_dgemm (int lda, double* A, double* B, double* C)
// {
  // /* For each block-row of A */
  // for (int j = 0; j < lda; j += BLOCK_SIZE)
    // /* For each block-column of B */
    // for (int k = 0; k < lda; k += BLOCK_SIZE)
      // /* Accumulate block dgemms into block of C */
      // for (int i = 0; i < lda; i += BLOCK_SIZE)
      // {
	// /* Correct block dimensions if block "goes off edge of" the matrix */
	// int M = min (BLOCK_SIZE, lda-i);
	// int N = min (BLOCK_SIZE, lda-j);
	// int K = min (BLOCK_SIZE, lda-k);

	// /* Perform individual block dgemm */
	// do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
      // }
// }

void square_dgemm (int lda, double* A, double* B, double* C)
{
  /* For each block-row of A */
  double *buf = (double *)malloc( sizeof(double)*BLOCK_SIZE*BLOCK_SIZE );
  int N, M, K;
  rep(k){
    /* For each block-column of B */
	K = min (BLOCK_SIZE, lda-k);
	rep(j){
		int N = min( BLOCK_SIZE, lda-j );
		int cnt = 0;
		for( int p = 0; p < N; p ++ ){
			for( int s = 0; s < K; s ++ ){
				buf[ cnt++ ] = B[ (j+p)*lda+k+s ];
			}
		}
		rep(i){
			M = min (BLOCK_SIZE, lda-i);
			do_block(lda, M, N, K, A + i + k*lda, buf, C + i + j*lda);
		}
	}
  }
  free(buf);
 }
 
 /*
 C[i+j*lda] += A[i+k*lda] * B[k+j*K];
 
 i j k
 repi
	repj{
	 double tmp = C[i+j*lda];
	 repk
	   tmp += A[i+k*lda] * B[k+j*K];
	 C[i+j*lda] = tmp;  
	}
 
 i k j
 repi
	repk{
		double tmp = A[i+k*lda];
		repj{
			C[i+j*lda] += tmp * B[k+j*K];
		}
	}
 
 j i k
 repj
	repi{
	 double tmp = C[i+j*lda];
	 repk
	   tmp += A[i+k*lda] * B[k+j*K];
	 C[i+j*lda] = tmp; 
	}
 
 j k i
   repj
	repk{
		double tmp = B[k+j*K];
		repi{
			C[i+j*lda] += A[i+k*lda]*tmp;
		}
	}
 
 k i j
 repk
    repi{
		double tmp = A[i+k*lda];
		repj{
			 C[i+j*lda] += tmp * B[k+j*K];
		}
	}
 
 k j i
 repk
	repj{
		double tmp = B[k+j*K];
		repi{
			C[i+j*lda] += A[i+k*lda]*tmp;
		}
	}
 
 */
