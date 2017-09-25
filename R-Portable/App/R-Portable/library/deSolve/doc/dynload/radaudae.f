c----------------------------------------------------------------
c----------------------------------------------------------------
c--- The car axis problem of radau
c----------------------------------------------------------------
c----------------------------------------------------------------


c -------- radaudae.f -> radaudae.dll ------
c compile in R with: system("g77 -shared -o radaudae.dll radaudae.f")
c or with system("R CMD SHLIB radaudae.f")
    
c----------------------------------------------------------------
c Initialiser for parameter common block
c----------------------------------------------------------------
      subroutine initcaraxis(daeparms)

      external daeparms
      integer, parameter :: N = 8
      double precision parms(N)
      common /myparms/parms

      call daeparms(N, parms)
      return
      end

c----------------------------------------------------------------
c rate of change
c----------------------------------------------------------------
      subroutine caraxis(neq, t, y, ydot, out, ip)
      implicit none
      integer          neq, IP(*)
      double precision t, y(neq), ydot(neq), out(*)
      double precision eps, M, k, L, L0, r, w, g
      common /myparms/ eps, M, k, L, L0, r, w, g

      double precision xl, yl, xr, yr, ul, vl, ur, vr, lam1, lam2
      double precision yb, xb, Ll, Lr, dxl, dyl, dxr, dyr
      double precision dul, dvl, dur, dvr, c1, c2

c expand state variables 
       xl = y(1)
       yl = y(2)
       xr = y(3) 
       yr = y(4) 
       ul = y(5) 
       vl = y(6) 
       ur = y(7) 
       vr = y(8) 
       lam1 = y(9) 
       lam2 = y(10)

       yb = r * sin(w * t)
       xb = sqrt(L * L - yb * yb)
       Ll = sqrt(xl**2 + yl**2)
       Lr = sqrt((xr - xb)**2 + (yr - yb)**2)
        
       dxl = ul
       dyl = vl
       dxr = ur
       dyr = vr
        
       dul = (L0-Ll) * xl/Ll      + 2 * lam2 * (xl-xr) + lam1*xb
       dvl = (L0-Ll) * yl/Ll      + 2 * lam2 * (yl-yr) + lam1*yb - k*g
       dur = (L0-Lr) * (xr-xb)/Lr - 2 * lam2 * (xl-xr)
       dvr = (L0-Lr) * (yr-yb)/Lr - 2 * lam2 * (yl-yr) - k*g
        
       c1  = xb * xl + yb * yl
       c2  = (xl - xr)**2 + (yl - yr)**2 - L * L
        
c function values in ydot
       ydot(1)  = dxl
       ydot(2)  = dyl
       ydot(3)  = dxr
       ydot(4)  = dyr
       ydot(5)  = dul
       ydot(6)  = dvl
       ydot(7)  = dur
       ydot(8)  = dvr
       ydot(9)  = c1
       ydot(10) = c2
      return
      end

