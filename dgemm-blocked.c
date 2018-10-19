const char* dgemm_desc = "kji";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE_I 8
#define BLOCK_SIZE_J 64
#define BLOCK_SIZE_K 128

#define MR 8

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
		double tmp = 0.0;
		repk{
			tmp += A[i*K+k]*B[k+j*K];
		}
		C[i+j*I] = tmp;
	}
}

void square_dgemm (int lda, double* A, double* B, double* C)
{
  /* For each block-row of A */
  // double *buf = (double *)_mm_malloc( sizeof(double)*BLOCK_SIZE_J*BLOCK_SIZE_K, 64);
  // double *cuf = (double *)_mm_malloc( sizeof(double)*BLOCK_SIZE_J*BLOCK_SIZE_I, 64 );
  // double *auf = (double *)_mm_malloc( sizeof(double)*BLOCK_SIZE_K*lda, 64 );
  double *buf = (double *)malloc( sizeof(double)*BLOCK_SIZE_J*BLOCK_SIZE_K);
  double *cuf = (double *)malloc( sizeof(double)*BLOCK_SIZE_J*BLOCK_SIZE_I);
  double *auf = (double *)malloc( sizeof(double)*BLOCK_SIZE_K*lda);
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
			for( int p = 0; p < J; p += MR ){
				int P = min( MR, J-p );
				do_block(lda, I, P, K, auf + K*i, buf + K*p, cuf + p*I );
			}
			int cnt = 0;
            for( int jj = 0; jj < J; jj ++ ){
                for( int ii = 0; ii < I; ii ++ ){
					C[ (jj+j)*lda + ii+i ] += cuf[ cnt++ ];
				}
			}
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
