const char* dgemm_desc = "kji";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE_I 2
#define BLOCK_SIZE_J 65
#define BLOCK_SIZE_K 65

#define INNER_I 2
#define INNER_J 2
#define INNER_K 2
#endif

#include<stdlib.h>
#define min(a,b) (((a)<(b))?(a):(b))
#define perj for(int j = 0; j < lda; j += BLOCK_SIZE_J)
#define peri for(int i = 0; i < lda; i += BLOCK_SIZE_I)
#define perk for(int k = 0; k < lda; k += BLOCK_SIZE_K)
#define repj for (int j = 0; j < J; ++j)
#define repi for (int i = 0; i < I; ++i)
#define repk for (int k = 0; k < K; ++k)
	
/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int I, int J, int K, double* A, double* B, double* C)
{
 repi
	repj{
		double tmp = C[i*J+j];
		repk{
			tmp += A[i*K+k]*B[k+j*K];
		}
		C[i*J+j] = tmp;
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
  double *buf = (double *)malloc( sizeof(double)*BLOCK_SIZE_J*BLOCK_SIZE_K );
  double *cuf = (double *)malloc( sizeof(double)*BLOCK_SIZE_K*lda );
  double *auf = (double *)malloc( sizeof(double)*BLOCK_SIZE_J*lda );
  int I, J, K;
  perk{
    /* For each block-column of B */
	K = min (BLOCK_SIZE_K, lda-k);
	int a_cnt = 0;
	for( int i = 0; i < lda; i ++ ){
		for( int kk = 0; kk < K; kk ++ )
			auf[a_cnt++] = A[ (k+kk)*lda+i ];
	}
	int t = 0;
	perj{
		J = min (BLOCK_SIZE_J, lda-j);
		int cnt = 0;
		for( int p = 0; p < J; p ++ ){
			for( int s = 0; s < K; s ++ ){
				buf[ cnt++ ] = B[ (j+p)*lda+k+s ];
			}
		}
		peri{
			I = min( BLOCK_SIZE_I, lda-i );
			do_block(lda, I, J, K, auf + K*i, buf, cuf + i*J );
		}
		cnt = 0;
		for( int ii = 0; ii < lda; ii ++ )
			for( int jj = 0; jj < J; jj ++ ){
				double tmp = cuf[ cnt ];
				cuf[ cnt++ ] = 0.0;
				C[ (jj+j)*lda + ii ] += tmp;
			}
		
	}
	
  }
  free(buf);
  free(auf);
  free(cuf);
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