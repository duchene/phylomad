c -------- odefor.f -> odefor.dll ------
c compile in R with: system("g77 -shared -o odefor.dll odefor.f")
c or with system("R CMD SHLIB odefor.f")
    
c Initialiser for parameter common block
      subroutine odefor(odeparms)

      external odeparms
      integer N
      double precision parms(3)
      common /myparms/parms

      N = 3
      call odeparms(N, parms)
      return
      end

c Rate of change and 3 output variables
      subroutine derivsfor (neq, t, y, ydot,out,IP)
      integer          neq, IP(*)
      double precision t, y(neq), ydot(neq), out(*), k1, k2, k3
      common /myparms/k1,k2,k3
      if(IP(1) < 3) call rexit("nout should be at least 3")      
      ydot(1) = -k1*y(1) + k2*y(2)*y(3)
      ydot(3) = k3*y(2)*y(2)
      ydot(2) = -ydot(1) - ydot(3)
      out(1)= y(1)+y(2)+y(3)
      out(2)= y(1)*2
      out(3)= IP(1)
      return
      end

c The jacobian matrix 
      subroutine jacfor (neq, t, y, ml, mu, pd, nrowpd,RP,IP)
      integer neq, ml, mu, nrowpd ,IP(*)
      double precision y(*), pd(nrowpd,*), t, RP(*), k1, k2, k3
      common /myparms/k1,k2,k3
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

