C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,2)*P(2,1) - P(-1,1)*P(-1,2)*Metric(1,2)
C     
      SUBROUTINE VVS3_3(V1, V2, COUP, M3, W3, S3)
      IMPLICIT NONE
      DOUBLE COMPLEX V1(*)
      DOUBLE COMPLEX V2(*)
      DOUBLE COMPLEX S3(*)
      DOUBLE COMPLEX COUP
      DOUBLE COMPLEX DENOM
      DOUBLE PRECISION M3, W3
      DOUBLE PRECISION P2(0:3),P3(0:3),P1(0:3)

      S3(2)= V1(5)+V2(5)
      S3(3)= V1(6)+V2(6)
      P2(0) =  DBLE(V2(5))
      P2(1) =  DBLE(V2(6))
      P2(2) =  DIMAG(V2(6))
      P2(3) =  DIMAG(V2(5))
      P3(0) = - DBLE(S3(2))
      P3(1) = - DBLE(S3(3))
      P3(2) = - DIMAG(S3(3))
      P3(3) = - DIMAG(S3(2))
      P1(0) =  DBLE(V1(5))
      P1(1) =  DBLE(V1(6))
      P1(2) =  DIMAG(V1(6))
      P1(3) =  DIMAG(V1(5))

      DENOM =1D0/(( (M3*( -M3+(0, 1)*W3))+( (P3(0)**2)-(P3(1)**2)
     $ -(P3(2)**2)-(P3(3)**2))))
      S3(1)= COUP*DENOM*( (V2(1)*( (V1(1)*( (0, 1)*(P2(1)*P1(1))
     $ +(0, 1)*(P2(2)*P1(2))+(0, 1)*(P2(3)*P1(3))))+(P1(0)*( (0, 
     $ -1)*(V1(2)*P2(1))+(0, -1)*(V1(3)*P2(2))+(0, -1)*(V1(4)
     $ *P2(3))))))+( (V2(2)*( (V1(2)*( (0, 1)*(P2(0)*P1(0))+(0, 
     $ -1)*(P2(2)*P1(2))+(0, -1)*(P2(3)*P1(3))))+(P1(1)*( (0, 
     $ -1)*(V1(1)*P2(0))+(0, 1)*(V1(3)*P2(2))+(0, 1)*(V1(4)*P2(3))))))
     $ +( (V2(3)*( (V1(3)*( (0, 1)*(P2(0)*P1(0))+(0, -1)*(P2(1)*P1(1))
     $ +(0, -1)*(P2(3)*P1(3))))+(P1(2)*( (0, -1)*(V1(1)*P2(0))
     $ +(0, 1)*(V1(2)*P2(1))+(0, 1)*(V1(4)*P2(3))))))+(V2(4)*( (V1(4)
     $ *( (0, 1)*(P2(0)*P1(0))+(0, -1)*(P2(1)*P1(1))+(0, -1)*(P2(2)
     $ *P1(2))))+(P1(3)*( (0, -1)*(V1(1)*P2(0))+(0, 1)*(V1(2)*P2(1))
     $ +(0, 1)*(V1(3)*P2(2)))))))))
      END


