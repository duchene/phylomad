c -------- scoc.f -> scoc.dll ------
c compile in R with: system("g77 -shared -o scoc.dll SCOC.f")
c or with system("R CMD SHLIB scoc.f")
    
c Initialiser for parameter common block
      subroutine scocpar(odeparms)

      external odeparms
      integer N
      double precision parms(1)
      common /myparms/parms

      N = 1
      call odeparms(N, parms)
      return
      end

c Initialiser for forcing common block
      subroutine scocforc(odeforcs)

      external odeforcs
      integer N
      double precision forcs(1)
      common /myforcs/forcs

      N = 1
      call odeforcs(N, forcs)
      return
      end

c Rate of change and output variables
      subroutine scocder (neq, t, y, ydot,out,IP)
      integer          neq, IP(*)
      double precision t, y(neq), ydot(neq), out(*), k, depo

      common /myparms/k
      common /myforcs/depo
      
      if(IP(1) < 2) call rexit("nout should be at least 2")

      ydot(1) = -k*y(1) + depo

      out(1)= k*y(1)
      out(2)= depo

      return
      end
                      

