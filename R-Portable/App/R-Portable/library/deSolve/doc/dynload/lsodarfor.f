c----------------------------------------------------------------
c----------------------------------------------------------------
c--- The root model example of lsodar
c----------------------------------------------------------------
c----------------------------------------------------------------


c -------- lsodarfor.f -> lsodarfor.dll ------
c compile in R with: system("g77 -shared -o lsodarfor.dll lsodarfor.f")
c or with system("R CMD SHLIB lsodarfor.f")
    
c----------------------------------------------------------------
c Initialiser for parameter common block
c----------------------------------------------------------------
      subroutine lsodarfor(odeparms)

      external odeparms
      integer, parameter :: N = 3
      double precision parms(N)
      common /myparms/parms

      call odeparms(N, parms)
      return
      end

c----------------------------------------------------------------
c rate of change and 1 output variable
c----------------------------------------------------------------
      subroutine modfor(neq, t, y, ydot,out,IP)
      integer          neq, IP(*)
      double precision t, y(neq), ydot(neq), out(*), aa, bb, cc
      common /myparms/aa,bb,cc
 
      if(IP(1) < 1) call rexit("nout should be at least 1") 
      ydot(1) = aa*y(1) + bb*y(2)*y(3)
      ydot(3) = cc*y(2)*y(2)
      ydot(2) = -ydot(1) - ydot(3)
      out(1)=y(1)+y(2)+y(3)

      return
      end

c----------------------------------------------------------------
c The root function
c----------------------------------------------------------------
      subroutine myroot(neq, t, y, ng, gout)
      integer :: neq, ng
      double precision :: t, y(neq), gout(ng)

      gout(1) = y(1) - 1.e-4
      gout(2) = y(3) - 1e-2
      
      return
      end

