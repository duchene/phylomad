c----------------------------------------------------------------
c----------------------------------------------------------------
c--- The chemical model example of daspk
c----------------------------------------------------------------
c----------------------------------------------------------------


c -------- daspkdll.f -> daspkdll.dll ------
c compile in R with: system("g77 -shared -o daspkfor.dll daspkfor.f")
c or with system("R CMD SHLIB daspkfor.f")
    
c----------------------------------------------------------------
c Initialiser for parameter common block
c----------------------------------------------------------------
      subroutine daspkfor(daspkparms)

      external daspkparms
      integer, parameter :: N = 4
      double precision parms(N)
      common /myparms/parms

      call daspkparms(N, parms)
      return
      end

c----------------------------------------------------------------
c residual of rate of change and 1 output variable
c----------------------------------------------------------------
      subroutine resfor (t, y, ydot, cj, delta, ires, out, ipar)
      integer            :: ires, ipar(*)
      integer, parameter :: neq = 3
      double precision   :: t, y(neq), ydot(neq), delta(neq), out(*)
      double precision   :: K, ka, r, prod, ra, rb
      common /myparms/K,ka,r,prod

      if(IPar(1) < 1) call rexit("nout should be at least 1")
      ra  = ka* y(3)            ! forward rate
      rb  = ka/K *y(1) * y(2)   ! backward rate

      ! residuals of rates of changes
      delta(3) = -ydot(3)  -  ra + rb + prod
      delta(1) = -ydot(1)  +  ra - rb
      delta(2) = -ydot(2)  +  ra - rb - r*y(2)
      out(1)   = y(1) + y(2) + y(3)
      return
      end

c----------------------------------------------------------------
c The jacobian matrix 
c----------------------------------------------------------------
      subroutine resjacfor (t, y, dy, pd, cj, out, ipar)
      integer, parameter :: neq = 3
      integer            :: ipar(*)
      double precision   :: K, ka, r, prod
      double precision   :: pd(neq,neq),y(neq),dy(neq),out(*)
      common /myparms/K,ka,r,prod

      ! residuals of rates of changes
      !res1 = -dD  - ka*D + ka/K *A*B + prod
      PD(1,1) = ka/K *y(2)
      PD(1,2) = ka/K *y(1)
      PD(1,3) = -ka -cj
      !res2 = -dA  + ka*D - ka/K *A*B
      PD(2,1) = -ka/K *y(2) -cj
      PD(2,2) = -ka/K *y(2)
      PD(2,3) = ka
      !res3 = -dB  + ka*D - ka/K *A*B - r*B
      PD(3,1) = -ka/K *y(2)
      PD(3,2) = -ka/K *y(2) -r -cj
      PD(3,3) = ka

      return
      end

