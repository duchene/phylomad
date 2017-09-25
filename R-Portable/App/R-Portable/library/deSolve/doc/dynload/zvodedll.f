
C The program below uses ZVODE to solve the following system of 2 ODEs:
C dz/dt = i*z; dw/dt = -i*w*w*z,z(0) = 1; w(0) = 1/2.1,  t = 0 to 2*pi.
C Solution: w = 1/(z + 1.1), z = exp(it).  As z traces the unit circle,
C w traces a circle of radius 10/2.1 with center at 11/2.1.


      SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
      INTEGER NEQ, IPAR(*)
      DOUBLE COMPLEX Y(NEQ), YDOT(NEQ), RPAR(*), CMP
      DOUBLE PRECISION T
      character(len=100) msg

c the imaginary unit i
      CMP = DCMPLX(0.0D0,1.0D0)

      YDOT(1) = CMP*Y(1)
      YDOT(2) = -CMP*Y(2)*Y(2)*Y(1)

      RETURN
      END

      SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      INTEGER NEQ, ML, MU, NRPD, IPAR(*)
      DOUBLE COMPLEX Y(NEQ), PD(NRPD,NEQ), RPAR(*), CMP
      DOUBLE PRECISION T
c the imaginary unit i
      CMP = DCMPLX(0.0D0,1.0D0)

      PD(2,3) = -2.0D0*CMP*Y(1)*Y(2)
      PD(2,1) = -CMP*Y(2)*Y(2)
      PD(1,1) = CMP
      RETURN
      END
