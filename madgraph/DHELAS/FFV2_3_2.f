C     This File is Automatically generated by ALOHA 
C     
      SUBROUTINE FFV2_3_2(F1, V3, COUP1,COUP2, M2, W2, F2)
      IMPLICIT NONE
      DOUBLE COMPLEX F1(*)
      DOUBLE COMPLEX F2(*)
      DOUBLE COMPLEX V3(*)
      DOUBLE COMPLEX COUP1,COUP2
      DOUBLE COMPLEX DENOM
      DOUBLE PRECISION M2, W2
      DOUBLE PRECISION P2(0:3)
      DOUBLE COMPLEX TMP(6)
      INTEGER I

      CALL FFV2_2(F1, V3, COUP1, M2, W2, F2)
      CALL FFV3_2(F1, V3, COUP2, M2, W2, TMP)
      DO I=1,4
        F2(I) = F2(I) + TMP(I)
      ENDDO
      END


