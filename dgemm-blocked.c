const char* dgemm_desc = "kji";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE_I 8
#define BLOCK_SIZE_J 48
#define BLOCK_SIZE_K 128

#define MR 6

#define DI 4
#define DJ 3
#define DK 4

#endif

#include<stdlib.h>
#include <immintrin.h>
#define min(a,b) (((a)<(b))?(a):(b))
#define perj for(int j = 0; j < lda; j += BLOCK_SIZE_J)
#define peri for(int i = 0; i < lda; i += BLOCK_SIZE_I)
#define perk for(int k = 0; k < lda; k += BLOCK_SIZE_K)
#define repj for (int j = 0; j < J; ++j)
#define repi for (int i = 0; i < I; ++i)
#define repk for (int k = 0; k < K; ++k)
#define MM
	
/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
 
static void do_block (int lda, int I, int J, int K, double* A, double* B, double* C)
{
	int II = I/DI*DI;
	int JJ = J/DJ*DJ;
	int KK = K/DK*DK;

	__m256d a, c[DI][DJ], b[DJ];
	for( int i = 0; i < II; i += DI ){
		for( int j = 0; j < JJ; j += DJ ){
			for( int jj = 0; jj < DJ; jj ++ ){
				for( int ii = 0; ii < DI; ii ++ ){
					c[ii][jj] = _mm256_setzero_pd();
				}
			}
			
			for( int k = 0; k < KK; k += DK ){
				for( int jj = 0; jj < DJ; jj ++ ){
					b[jj] = _mm256_load_pd( B+(j+jj)*BLOCK_SIZE_K+k );
				}
				for( int ii = 0; ii < DI; ii ++ ){
					a = _mm256_load_pd( A+(i+ii)*BLOCK_SIZE_K+k );
					for( int jj = 0; jj < DJ; jj ++ ){
						c[ii][jj] = _mm256_fmadd_pd( b[jj], a, c[ii][jj] );
					}
				}
			}
			for( int jj = 0; jj < DJ; jj ++ )
				for( int ii = 0; ii < DI; ii ++ ){
					double d[4];
					_mm256_store_pd( d, c[ii][jj] );
					C[BLOCK_SIZE_I*(j+jj)+i+ii] = d[0]+d[1]+d[2]+d[3];
				}
		}
	}
	
	//(1,1)*(1,1)
	for( int j = 0; j < JJ; j ++ ){
		for( int i = 0; i < II; i ++ ){
			double tmp = C[i+j*BLOCK_SIZE_I];
			for( int k = KK; k < K; k ++ ){
				tmp += A[i*BLOCK_SIZE_K+k]*B[k+j*BLOCK_SIZE_K];
			}
			C[i+j*BLOCK_SIZE_I] = tmp;
		}
	}
	
	//(1,2)*(1,2)
	for( int j = 0; j < JJ; j ++ ){
		for( int i = II; i < I; i ++ ){
			double tmp = 0.0;
			for( int k = 0; k < K; k ++ ){
				tmp += A[i*BLOCK_SIZE_K+k]*B[k+j*BLOCK_SIZE_K];
			}
			C[i+j*BLOCK_SIZE_I] = tmp;
		}
	}
	
	//(1,2)*(2,2)
	for( int j = JJ; j < J; j ++ ){
		for( int i = 0; i < I; i ++ ){
			double tmp = 0.0;
			for( int k = 0; k < K; k ++ ){
				tmp += A[i*BLOCK_SIZE_K+k]*B[k+j*BLOCK_SIZE_K];
			}
			C[i+j*BLOCK_SIZE_I] = tmp;
		}
	}
	
}

// static void do_block (int lda, int I, int J, int K, double* A, double* B, double* C)
// {
   // repj
	// repi{
	 // double tmp = 0.0;
	 // repk
	   // tmp += A[i*BLOCK_SIZE_K+k]*B[k+j*BLOCK_SIZE_K];
	 // C[i+j*BLOCK_SIZE_I] = tmp;
	// }
// }

void square_dgemm (int lda, double* A, double* B, double* C)
{
  /* For each block-row of A */
#ifdef MM
  double *buf = (double *)_mm_malloc( sizeof(double)*BLOCK_SIZE_J*BLOCK_SIZE_K, 4096);
  double *cuf = (double *)_mm_malloc( sizeof(double)*BLOCK_SIZE_J*BLOCK_SIZE_I, 4096 );
  double *auf = (double *)_mm_malloc( sizeof(double)*BLOCK_SIZE_K*lda, 4096 );
  //double *d = (double *)_mm_malloc( sizeof(double)*BLOCK_SIZE_K, 64 );
#else
  
  // double *auf = (double *)malloc( sizeof(double)*BLOCK_SIZE_K*lda);
  // double *buf = (double *)malloc( sizeof(double)*BLOCK_SIZE_J*BLOCK_SIZE_K);
  // double *cuf = (double *)malloc( sizeof(double)*BLOCK_SIZE_J*BLOCK_SIZE_I);
  
  double *buf = (double *)malloc( sizeof(double)*BLOCK_SIZE_J*BLOCK_SIZE_K);
  double *cuf = (double *)malloc( sizeof(double)*BLOCK_SIZE_J*BLOCK_SIZE_I);
  double *auf = (double *)malloc( sizeof(double)*BLOCK_SIZE_K*lda);
#endif
  
  
  int I, J, K;
  perk{
    /* For each block-column of B */
	K = min (BLOCK_SIZE_K, lda-k);
	int a_cnt = 0;
	for( int i = 0; i < lda; i ++ ){
		for( int kk = 0; kk < K; kk ++ )
			auf[ i*BLOCK_SIZE_K+kk ] = A[ (k+kk)*lda+i ];
	}
	int t = 0;
	perj{
		J = min (BLOCK_SIZE_J, lda-j);
		int cnt = 0;
		for( int p = 0; p < J; p ++ ){
			for( int s = 0; s < K; s ++ ){
				buf[ p*BLOCK_SIZE_K+s ] = B[ (j+p)*lda+k+s ];
			}
		}
		peri{
			I = min( BLOCK_SIZE_I, lda-i );
			for( int p = 0; p < J; p += MR ){
				int P = min( MR, J-p );
				do_block(lda, I, P, K, auf + BLOCK_SIZE_K*i, \
				buf + BLOCK_SIZE_K*p, cuf + p*BLOCK_SIZE_I );
				// do_block(lda, I, P, K, auf + K*i, buf + K*p, cuf );
				// int cnt = 0;
				// for( int jj = 0; jj < P; jj ++ ){
					// for( int ii = 0; ii < I; ii ++ ){
						// double tmp = cuf[ cnt ++ ];
						// C[ (jj+j+p)*lda + ii+i ] += tmp;
					// }
				// }
			}
			int cnt = 0;
            for( int jj = 0; jj < J; jj ++ ){
                for( int ii = 0; ii < I; ii ++ ){
					double tmp = cuf[ jj*BLOCK_SIZE_I + ii ];
					C[ (jj+j)*lda + ii+i ] += tmp;
				}
			}
		}

	}

  }
#ifdef MM
  _mm_free(auf);
  _mm_free(buf);
  _mm_free(cuf);
  //_mm_free(d);
#else  
  free(buf);
  free(auf);
  free(cuf);
#endif  
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
