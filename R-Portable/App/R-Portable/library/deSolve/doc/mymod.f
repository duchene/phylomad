c file mymodf.f
       subroutine initmod(odeparms)
        external odeparms
        double precision parms(3)
        common /myparms/parms

         call odeparms(3, parms)
        return
       end

       subroutine derivs (neq, t, y, ydot, yout, ip)
        double precision t, y, ydot, k1, k2, k3
        integer neq, ip(*)
        dimension y(3), ydot(3), yout(*)
        common /myparms/k1,k2,k3

          if(ip(1) < 1) call rexit("nout should be at least 1")

          ydot(1) = -k1*y(1) + k2*y(2)*y(3)
          ydot(3) = k3*y(2)*y(2)
          ydot(2) = -ydot(1) - ydot(3)

          yout(1) = y(1) + y(2) + y(3)
        return
       end

       subroutine jac (neq, t, y, ml, mu, pd, nrowpd, yout, ip)
        integer neq, ml, mu, nrowpd, ip
        double precision y(*), pd(nrowpd,*), yout(*), t, k1, k2, k3
        common /myparms/k1, k2, k3

          pd(1,1) = -k1
          pd(2,1) = k1
          pd(3,1) = 0.0
          pd(1,2) = k2*y(3)
          pd(2,2) = -k2*y(3) - 2*k3*y(2)
          pd(3,2) = 2*k3*y(2)
          pd(1,3) = k2*y(2)
          pd(2,3) = -k2*y(2)
          pd(3,3) = 0.0
        return
       end
c end of file mymodf.f
