c -------- odefor2.f -> odefor2.dll ------
c compile in R with: system("g77 -shared -o odefor2.dll odefor2.f")
c or with system("R CMD SHLIB odefor2.f")
c fortran source without initialiser
c Rate of change and 3 output variables

      subroutine derivsfor2 (neq, t, y, ydot,out,IP)
      integer          neq, IP(*)
      double precision t, y(neq), ydot(neq), out(*), k1, k2, k3

      k1 = 0.04
      k2 = 1e4
      k3 = 3e7
      if(IP(1) < 1) call rexit("nout should be at least 1")      

      ydot(1) = -k1*y(1) + k2*y(2)*y(3)
      ydot(3) = k3*y(2)*y(2)
      ydot(2) = -ydot(1) - ydot(3)
      out(1)=y(1)+y(2)+y(3)
      out(2)=y(1)*2
      out(3)=k3
      return
      end

c The jacobian matrix 
      subroutine jacfor2 (neq, t, y, ml, mu, pd, nrowpd,RP,IP)
      integer neq, ml, mu, nrowpd ,IP(*)
      double precision y(*), pd(nrowpd,*), t, RP(*), k1, k2, k3
      k1 = 0.04
      k2 = 1e4
      k3 = 3e7

      pd(1,1) = -k1
      pd(2,1) =  k1
      pd(3,1) = 0.0
      pd(1,2) = k2*y(3)
      pd(2,2) = -k2*y(3) - 2*k3*y(2)
      pd(3,2) = 2*k3*y(2)
      pd(1,3) = k2*y(2)
      pd(2,3) = -k2*y(2)
      pd(3,3) = 0.0
                         

      return
      end

