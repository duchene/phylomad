C  file dede_lfF.f
C  Initializer for parameter common block
      subroutine initmod(odeparms)
      external odeparms
      double precision parms(5)
      common /myparms/parms

      call odeparms(5, parms)
      return
      end

C  Derivatives and one output variable
      subroutine derivs(neq, t, y, ydot, yout, ip)
      integer neq, ip(*)
      double precision t, y(neq), ydot(neq), yout(*)

      double precision N, P, ytau(2), tlag
      integer nr(2)

      double precision f, g, e, m, tau
      common /myparms/f, g, e, m, tau

      if (ip(1) < 2) call rexit("nout should be at least 2")

      N = y(1)
      P = y(2)
      nr(1) = 0
      nr(2) = 1
      ytau(1) = 1.0
      ytau(2) = 1.0
 
      tlag = t - tau
      if (tlag .GT. 0.0) call lagvalue(tlag, nr, 2, ytau)

      ydot(1) = f * N - g * N * P
      ydot(2) = e * g * ytau(1) * ytau(2) - m * P

      yout(1) = ytau(1)
      yout(2) = ytau(2)

      return
      end

