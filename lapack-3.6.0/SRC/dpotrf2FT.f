*> \brief \b DPOTRF2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE DPOTRF2( UPLO, N, A, LDA, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DPOTRF2 computes the Cholesky factorization of a real symmetric
*> positive definite matrix A using the recursive algorithm.
*>
*> The factorization has the form
*>    A = U**T * U,  if UPLO = 'U', or
*>    A = L  * L**T,  if UPLO = 'L',
*> where U is an upper triangular matrix and L is lower triangular.
*>
*> This is the recursive version of the algorithm. It divides
*> the matrix into four submatrices:
*>
*>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
*>    A = [ -----|----- ]  with n1 = n/2
*>        [  A21 | A22  ]       n2 = n-n1
*>
*> The subroutine calls itself to factor A11. Update and scale A21
*> or A12, update A22 then calls itself to factor A22.
*> 
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          N-by-N upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky
*>          factorization A = U**T*U or A = L*L**T.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the leading minor of order i is not
*>                positive definite, and the factorization could not be
*>                completed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup doublePOcomputational
*
*  =====================================================================
      SUBROUTINE DPOTRF2FT( UPLO, N, A, LDA, INFO,CHKA,LDCA,CHKV,LDCV )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N,LDCA,LDCV,LDAR
      DOUBLE PRECISION   ZERO, ALPHA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ),CHKA(LDCA,*),CHKV(2,16),
      DOUBLE PRECISION CHKAR(16,16)
*     ..
      EXTERNAL DPOTRF2, DGEMM

      INTEGER           J

*  =====================================================================
*
      ZERO=0.0D+0
      ALPHA=-1.0D+0
      LDCR=16
      LDAR=16

      CALL DGEMM('No transpose','No transpose',2,N,N,ALPHA,CHKV,LDCV,
     &            A,LDA,ZERO,CHKAR,LDAR)

      CALL DPOTRF2( UPLO, N, A, LDA, INFO)
      DO J=1, N
         CHKA(1,J)=CHKA(1,J)/A(J,J)
         CHKA(2,J)=CHKA(2,J)/A(J,J)
         CHKA(1,J+1:N)=CHKA(1,J+1:N)-CHKA(1,J)*A(J+1:N,J)
         CHKA(2,J+1:N)=CHKA(2,J+1:N)-CHKA(2,J)*A(J+1:N,J)

      END DO

      PRINT *, "POTRF2 UPDATED MATRIX"

      DO I=1, N
         Print 100, ( A(I,J), J=1,N )
      end do

      PRINT *, "POTRF2 UPDATED CHKSUM"
      DO I=1, 2
         Print 100, ( CHKA(I,J), J=1,N )
      end do
      
 100  format (1x, 16(1x,f7.1))

      RETURN

*     End of DPOTRF2FT .
      END
