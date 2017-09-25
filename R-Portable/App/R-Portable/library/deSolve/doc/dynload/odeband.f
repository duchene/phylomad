c ==========================================================================
c Example 1 of help file of lsode:
c a simple function with banded jacobian - upper and lower band = 1
c note that number of rows of PD  = nupper + 2*nlower + 1
c ==========================================================================


c Rate of change  
      subroutine derivsband (neq, t, y, ydot,out,IP)
      integer          neq, IP(*)
      DOUBLE PRECISION  T, Y(5), YDOT(5), out(*)
        ydot(1) = 0.1*y(1) -0.2*y(2)
        ydot(2) = -0.3*y(1) +0.1*y(2) -0.2*y(3)
        ydot(3) =           -0.3*y(2) +0.1*y(3) -0.2*y(4)
        ydot(4) =                     -0.3*y(3) +0.1*y(4) -0.2*y(5)
        ydot(5) =                               -0.3*y(4) +0.1*y(5)
      RETURN
      END

c The jacobian matrix 
      subroutine jacband (neq, t, y, ml, mu, pd, nrowpd,RP,IP)
      INTEGER  NEQ, ML, MU, nrowpd, ip(*)
      DOUBLE PRECISION  T, Y(5), PD(nrowpd,5), rp(*)

        PD(:,:) = 0.D0

        PD(1,1) =  0.D0
        PD(1,2) = -.02D0
        PD(1,3) = -.02D0
        PD(1,4) = -.02D0
        PD(1,5) = -.02D0

        PD(2,:) = 0.1D0

        PD(3,1) = -0.3D0
        PD(3,2) = -0.3D0
        PD(3,3) = -0.3D0
        PD(3,4) = -0.3D0
        PD(3,5) = 0.D0

      RETURN
      END
