C  file satres.f
C  Initializer for parameter common block
      subroutine initmod(odeparms)
      external odeparms
      double precision parms(16)
      common /myparms/parms

      call odeparms(16, parms)
      return
      end
C Initializer for forcing common block
      subroutine initforc(odeforcs)
      external odeforcs
      double precision forcs(1)
      common /myforcs/forcs

      call odeforcs(1, forcs)
      return
      end

C  Compartments are:
C  y(1)  central compartment
C  y(2)  second compartment
C  y(3)  filtrate compartment
C  y(4)  'Gut'
C  y(5)  Total eliminated
C  y(6)  AUC central compartment

C  Derivatives and one output variable
      subroutine derivs(neq, t, y, ydot, out, ip)
      integer neq, IP(*)
      double precision t, y(neq), ydot(neq), out(*)
      double precision Vc, Vt, kd, ka, Tm, KT, Kfil, Vfil, free, BW,
     $     Dose, DoseInt, Qd, Qfil, MaxTime, TDose, TDoseRt
      common /myparms/Vc, Vt, kd, ka, Tm, KT, Kfil, Vfil, free, BW,
     $     Dose, DoseInt, Qd, Qfil, MaxTime, TDose
      common /myforcs/TDoseRt

      if (ip(1) < 1) call rexit("nout should be at least 1")

      ydot(1) = (ka * y(4) - Qd * free * y(1) + Qd * y(2) -
     $     Qfil * y(1) * free) / Vc + Tm * y(3) / (KT + y(3))
      ydot(2) = (free * Qd * y(1) - Qd * y(2)) / Vt
      ydot(3) = (Vc * kfil * y(1) * free - Vc * Tm * y(3) / (KT + y(3))-
     $     Vc * kfil * y(3)) / Vfil
      ydot(4) = -ka * y(4) + TDoseRt
      ydot(5) = Vc * kfil * y(3)
      ydot(6) = y(1)

      out(1) = y(1) * Vc + y(2) * Vt + y(3) * Vfil + y(4) + y(5)
      return
      end
