*  =====================================================================
      SUBROUTINE GENERATE_CHECKSUM(N, B, A, LDA, CHECK_VECTOR, LDV, CHECK_MATRIX, LDM)
*
*  -- Reference BLAS level3 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER N,B, LDA, LDV, LDM
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),CHECK_VECTOR(LDV,*),CHECK_MATRIX(LDM,*)
*     ..
*
*  =====================================================================
*

*     .. External Subroutines ..
      EXTERNAL DGEMM

*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
*     ..
*     .. Parameters ..
      INTEGER NUMBER_CHECKSUM


      PARAMETER (NUMBER_CHECKSUM=2)
      DO 10 I = 1,N/B
          CALL DGEMM（'N', 'N', NUMBER_CHECKSUM, I*B, B, 1.0, CHECK_VECTOR, LDV, A(I+(I-1)*B,1), LDA, 0.0, CHECK_MATRIX, LDM）
10    CONTINUE



*
      RETURN
*
*     End of  GENERATE_CHECKSUM.
*
      END
