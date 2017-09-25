C  file dedesimpleF.f
C  Initializer for parameter common block
      subroutine initmod(odeparms)
      external odeparms
      double precision parms(2)
      common /myparms/parms

      call odeparms(2, parms)
      return
      end

C  Derivatives and one output variable
      subroutine derivs(neq, t, y, ydot, yout, ip)
      integer neq, ip(*)
      double precision t, y(neq), ydot(neq), yout(*)
      double precision tau, k, ytau(1), tlag
      integer nr(1)
      common /myparms/tau, k

      if (ip(1) < 1) call rexit("nout should be at least 1")
      nr(1) = 0
      ytau(1) = 1.0
      tlag = t - tau
      if (tlag .GT. 0.0) call lagvalue(tlag, nr, 1, ytau)
      yout(1) = ytau(1)
      ydot(1) = k * ytau(1)
      return
      end

