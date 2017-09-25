c----------------------------------------------------------------
c The chemical model example of daspk but with the
c production rate a forcing function rather than 
c a parameter...
c----------------------------------------------------------------


c -------- ChemicalDAE.f -> ChemicalDAE.dll ------
c compile in R with: system("g77 -shared -o ChemicalDAE.dll ChemicalDAE.f")
c or with system("R CMD SHLIB ChemicalDAE.f")
    
c----------------------------------------------------------------
c Initialiser for parameter common block
c----------------------------------------------------------------
      subroutine initparms(daspkparms)

      external daspkparms
      double precision parms(3)
      common /myparms/parms

      call daspkparms(3, parms)
      return
      end

c----------------------------------------------------------------
c Initialiser for forcing common block
c----------------------------------------------------------------
      subroutine initforcs(daspkforcs)

      external daspkforcs

      double precision forcs(1)
      common /myforcs/forcs

      call daspkforcs(1, forcs)
      return
      end


c----------------------------------------------------------------
c residual of rate of change and 1 output variable
c----------------------------------------------------------------
      subroutine chemres (t, y, ydot, cj, delta, ires, out, ipar)
      integer            :: ires, ipar(*)
      integer, parameter :: neq = 3
      double precision   :: t, y(neq), ydot(neq), delta(neq), out(*)
      double precision   :: K, ka, r, prod, ra, rb

      common / myparms / K, ka, r
      common / myforcs / prod

      if(IPar(1) < 2) call rexit("nout should be at least 2")
      
      ra  = ka* y(3)            ! forward rate
      rb  = ka/K *y(1) * y(2)   ! backward rate

      ! residuals of rates of changes
      delta(3) = -ydot(3)  -  ra + rb + prod
      delta(1) = -ydot(1)  +  ra - rb
      delta(2) = -ydot(2)  +  ra - rb - r*y(2)
      out(1)   = y(1) + y(2) + y(3)
      out(2)   = prod

      return
      end
