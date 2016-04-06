*> \brief \b DGEMM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,M,N
*       CHARACTER TRANSA,TRANSB
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEMM  performs one of the matrix-matrix operations
*>
*>    C := alpha*op( A )*op( B ) + beta*C,
*>
*> where  op( X ) is one of
*>
*>    op( X ) = X   or   op( X ) = X**T,
*>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )
*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n',  op( A ) = A.
*>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c',  op( A ) = A**T.
*> \endverbatim
*>
*> \param[in] TRANSB
*> \verbatim
*>          TRANSB is CHARACTER*1
*>           On entry, TRANSB specifies the form of op( B ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSB = 'N' or 'n',  op( B ) = B.
*>
*>              TRANSB = 'T' or 't',  op( B ) = B**T.
*>
*>              TRANSB = 'C' or 'c',  op( B ) = B**T.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry,  M  specifies  the number  of rows  of the  matrix
*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N  specifies the number  of columns of the matrix
*>           op( B ) and the number of columns of the matrix C. N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry,  K  specifies  the number of columns of the matrix
*>           op( A ) and the number of rows of the matrix op( B ). K must
*>           be at least  zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by m  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*>           least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*>           part of the array  B  must contain the matrix  B,  otherwise
*>           the leading  n by k  part of the array  B  must contain  the
*>           matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*>           least  max( 1, n ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*>           supplied as zero then C need not be set on input.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*>           Before entry, the leading  m by n  part of the array  C must
*>           contain the matrix  C,  except when  beta  is zero, in which
*>           case C need not be set on entry.
*>           On exit, the array  C  is overwritten by the  m by n  matrix
*>           ( alpha*op( A )*op( B ) + beta*C ).
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, m ).
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
*> \ingroup double_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGEMMFT(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,
     &                   LDC, CHKA,LDCA,CHKB,LDCB,CHKC,LDCC,CHKV,LDCV)
*
*  -- Reference BLAS level3 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA,ZERO
      INTEGER K,LDA,LDB,LDC,M,N,LDCA,LDCB,LDCC,LDCV,LDAR,LDVR,LDCR
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*),CHKA(LDCA,*)
      DOUBLE PRECISION CHKB(LDCB,*),CHKC(LDCC,*),CHKV(2,16)
      DOUBLE PRECISION CHKAR(16,16),CHKBR(16,16),CHKCR(16,16)
*     ..
      EXTERNAL DGEMM,checkFT


*  =====================================================================
*     .. Local Scalars ..
      INTEGER J

      LDAR=16
      LDBR=16
      LDCR=16
      ZERO=0.0D+0
     

      CALL DGEMM('No transpose','No transpose',2,K,N,BETA,CHKV,LDCV,
     &           B,LDB,ZERO,CHKBR,LDBR)
      
      DO 10 J=1, M, N
         CALL DGEMM('No transpose','No transpose',2,K,N,BETA,CHKV,LDCV,
     &              A(J,1),LDA,ZERO,CHKAR((J/N)*2+1,1),LDAR)
 10   CONTINUE

      DO 20 J=1, M, N
         CALL DGEMM('No transpose','No transpose',2,N,N,BETA,CHKV,LDCV,
     &              C(J,1),LDC,ZERO,CHKCR((J/N)*2+1,1),LDCR)
 20   CONTINUE   

      PRINT *, "GEMM NEW CHECKSUM OF A"
      DO I=1,(M/N)*2
         PRINT 100, (CHKAR(I,J),J=1,K)
      END DO

      PRINT *,"GEMM OLD CHECKSUM OF A"
      DO I=1,(M/N)*2
         PRINT 100, (CHKA(I,J),J=1,K)
      END DO

      PRINT *, "GEMM NEW CHECKSUM OF B"
      DO I=1,2
         PRINT 100, ( CHKBR(I,J),J=1,K)
      END DO

      PRINT *, "GEMM OLD CHECKSUM OF B"
      DO I=1,2
         PRINT 100, (CHKB(I,J), J=1,K)
      END DO

      PRINT *, "GEMM NEW CHECKSUM OF C"
      DO I=1,(M/N)*2
         PRINT 100, (CHKCR(I,J), J=1,N)
      END DO

      PRINT *, "GEMM OLD CHECKSUM OF C"
      DO I=1,(M/N)*2
         PRINT 100, (CHKC(I,J), J=1,N)
      END DO

      checkFT(A,LDA,M,N,K,CHKA,LDCA,CHKAR,LDAR)
      checkFT(B,LDB,N,N,K,CHKB,LDCB,CHKBR,LDBR)
      checkFT(C,LDC,M,N,N,CHKC,LDCC,CHKCR,LDCR)

      CALL DGEMM(TRANSA, TRANSB, M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      CALL DGEMM(TRANSA,TRANSB,(M/N)*2,N,K,ALPHA,CHKA,LDCA,B,LDB,
     &            BETA,CHKC,LDCC )

*      PRINT *, "GEMM UPDATED MATRIX"

*      DO I=1, M
*         Print 100, ( C(I,J), J=1,N )
*      end do

*      PRINT *, "GEMM UPDATED CHKSUM"
*      DO I=1, (M/N)*2
*         Print 100, ( CHKC(I,J), J=1,N )
*      end do
      
 100  format (1x, 16(1x,f7.1))


      RETURN
*
*     End of DGEMMFT .
*
      END
