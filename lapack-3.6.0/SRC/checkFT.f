      SUBROUTINE checkFT(A,LDA,M,B,N, CHKA, LDCA,CHKAR,LDAR)

      INTEGER           LDA, M, B, N, LDCA, LDAR,X,Y,I,J
      DOUBLE PRECISION  C1, C2, R1, R2, d1, d2,e, d1abs

*     ..Array Arguments ..
      DOUBLE PRECISION  A( LDA,*),CHKA(LDCA,*),CHKAR(LDAR,*)

*  ==============================================================
      e = 10.0D-10

      DO  I=1,M/B
         DO  J=1,N
            C1 = CHKA(I*2-1,J)
            C2 = CHKA(I*2,J)
            
            R1 = CHKAR(I*2-1,J)
            R2 = CHKAR(I*2,J)
            
            d1 = C1-R1
            d2 = C2-R2
            d1abs = d1

            IF ( d1 .LT. 0.0D+0)
               d1abs = -1 * d1
            END IF

*            IF (d1abs.GT. e)
*          print *,"maybe error"
*          x = NINT(d2/d1)
*               A((I-1)*B+X,J) += d1
*            end if
         end do 
      end do
      print *, "done"
      
      return
      end
