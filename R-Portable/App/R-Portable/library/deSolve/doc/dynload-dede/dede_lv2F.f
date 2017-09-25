C  file dede_lf2F.f
C  Initializer for parameter common block
      subroutine initmod(odeparms)
      external odeparms
      double precision parms(6)
      common /myparms/parms

      call odeparms(6, parms)
      return
      end

C  Derivatives and one output variable
      subroutine derivs(neq, t, y, ydot, yout, ip)
      integer neq, ip(*)
      double precision t, y(neq), ydot(neq), yout(*)

      double precision N, P, ytau(2), tlag1, tlag2
      integer nr(2)

      double precision f, g, e, m, tau1, tau2
      common /myparms/f, g, e, m, tau1, tau2

      if (ip(1) < 2) call rexit("nout should be at least 2")

      N = y(1)
      P = y(2)
      nr(1) = 0
      nr(2) = 1
      ytau(1) = 1.0
      ytau(2) = 1.0
 
      tlag1 = t - tau1
      tlag2 = t - tau2
      if (min(tlag1, tlag2) .GE. 0.0) then
        call lagvalue(tlag1, nr(1), 1, ytau(1))
        call lagvalue(tlag2, nr(2), 1, ytau(2))
      endif

      ydot(1) = f * N - g * N * P
      ydot(2) = e * g * ytau(1) * ytau(2) - m * P

      yout(1) = ytau(1)
      yout(2) = ytau(2)

      return
      end

      double precision function getlag(t0, t, tau, ydef, nr)
      double precision t0, t, tau, ydef
      integer nr
      double precision tlag, y
      tlag = t - tau
      y = ydef
      if (tlag .GE. t0) call lagvalue(tlag, nr, 1, y)
      getlag = y
      return
      end


      subroutine derivs2(neq, t, y, ydot, yout, ip)
      integer neq, ip(*)
      double precision t, y(neq), ydot(neq), yout(*)

      double precision N, P, ytau(2), getlag

      double precision f, g, e, m, tau1, tau2
      common /myparms/f, g, e, m, tau1, tau2

      if (ip(1) < 2) call rexit("nout should be at least 2")

      N = y(1)
      P = y(2)
      ytau(1) = 1.0
      ytau(2) = 1.0
      ytau(1) = getlag(0.0, t, tau1, ytau(1), 0)
      ytau(2) = getlag(0.0, t, tau2, ytau(2), 1)
 
      ydot(1) = f * N - g * N * P
      ydot(2) = e * g * ytau(1) * ytau(2) - m * P

      yout(1) = ytau(1)
      yout(2) = ytau(2)

      return
      end

