      subroutine hvvhxx(v1,v2,gc,smass,swidth , jsvv)
c
c- by RF - Mar. 2006
c
c
c This subroutine computes an off-shell (pseudo-)scalar current from
c two incoming vectors using the scalar effective coupling.
c
c input:
c       complex v1(6)          : first  vector
c       complex v2(6)          : second vector
c       complex gc(2)          : coupling constant: gc(1) scalar
c                                                   gc(2) pseudo-scalar
c       real smass             : mass of the outgoing (pseudo-)scalar
c       real swidth            : width of the outgoing (pseudo-)scalar
c
c output:
c       complex jsvv(3)        : (pseudo-)scalar current
c     
c
      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex v1(DIM),v2(DIM),jsvv(DIM),vertex1,vertex2,dj
      double complex v12,p2v1,p1v2	
      double complex v13,v14,v23,v24,v34
      double precision p12,p13,p14,p23,p24,p34
      double precision p1(0:3),p2(0:3),q(4),q2
      double precision smass,swidth
      double complex gc(2)

      p1(0) = dble( v1(5))
      p1(1) = dble( v1(6))
      p1(2) = dimag(v1(6))
      p1(3) = dimag(v1(5))

      p2(0) = dble( v2(5))
      p2(1) = dble( v2(6))
      p2(2) = dimag(v2(6))
      p2(3) = dimag(v2(5))

      vertex1 = (0D0,0D0)
      vertex2 = (0D0,0D0)
      
      jsvv(2) = v1(5) + v2(5)
      jsvv(3) = v1(6) + v2(6)

      q(1) = -dble( jsvv(2))
      q(2) = -dble( jsvv(3))
      q(3) = -dimag(jsvv(3))
      q(4) = -dimag(jsvv(2))

      q2 = q(1)**2 - q(2)**2 - q(3)**2 - q(4)**2

      dj = dcmplx( q2-smass**2, smass*swidth )


      if (gc(1).NE.(0D0,0D0)) then

         v12  = v1(1)*v2(1) - v1(2)*v2(2) - v1(3)*v2(3) - v1(4)*v2(4)
         p12  = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
         p2v1 = v1(1)*p2(0) - v1(2)*p2(1) - v1(3)*p2(2) - v1(4)*p2(3)
         p1v2 = p1(0)*v2(1) - p1(1)*v2(2) - p1(2)*v2(3) - p1(3)*v2(4)	

         vertex1 = gc(1)*(v12*p12 - p2v1*p1v2)
      endif

      if (gc(2).NE.(0D0,0D0)) then
          p12 = p1(0)*p2(1) - p1(1)*p2(0)
          p13 = p1(0)*p2(2) - p1(2)*p2(0)
          p14 = p1(0)*p2(3) - p1(3)*p2(0)
          p23 = p1(1)*p2(2) - p1(2)*p2(1)
          p24 = p1(1)*p2(3) - p1(3)*p2(1)
          p34 = p1(2)*p2(3) - p1(3)*p2(2)

          v12 = v1(1)*v2(2) - v1(2)*v2(1)
          v13 = v1(1)*v2(3) - v1(3)*v2(1)
          v14 = v1(1)*v2(4) - v1(4)*v2(1)
          v23 = v1(2)*v2(3) - v1(3)*v2(2)
          v24 = v1(2)*v2(4) - v1(4)*v2(2)
          v34 = v1(3)*v2(4) - v1(4)*v2(3)

          vertex2 = - gc(2)*( v12*p34 - v13*p24 + v14*p23
     &                       +v23*p14 - v24*p13 + v34*p12 )
      endif
       
      jsvv(1) = (vertex1 + vertex2) /dj


      return
      end
