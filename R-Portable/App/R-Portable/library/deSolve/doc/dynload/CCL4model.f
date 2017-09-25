
c the CCl4 inhalation model
c based on the demo in odesolve


c -------- ccl4model.f -> ccl4model.dll ------
c compile in R with: system("g77 -shared -o ccl4model.dll ccl4model.f")
c or with system("R CMD SHLIB ccl4model.f")
    

c=======================================================================
c=======================================================================
c                      Model initialisation
c=======================================================================
c=======================================================================

c=======================================================================
c Initialise primary parameter common block
c=======================================================================

      subroutine initccl4(odeparms)

      external odeparms
      integer N
c parameters are divided into primary and derived parameters
      double precision pars(21), derivedpars(15)
      common /myparms/pars,derivedpars

       N = 21
       call odeparms(N, pars)

       call derived()           
      return
      end

c=======================================================================
c In this "event", state variable 1 is increased with 1. DOES NOT WORK...
c=======================================================================

      subroutine eventfun(n, t, y)
      integer n
      double precision t, y(7)

      y(1) = y(1) + 1
      
      end subroutine
      
      
c=======================================================================
c Calculate derived parameters from primary parameters
c=======================================================================

      subroutine derived  
      implicit none
      double precision  BW,QP,QC,VFC,VLC,VMC,QFC,QLC,QMC,PLA,PFA,         &
     &                  PMA,PTA,PB,MW,VMAX,KM,CONC,KL,RATS,VCHC,          &
     &                  VCH,VM,VT,VF,VL,PM,PT,PF,PL,AI0,VTC,QM,QT,QF,QL
      common /myparms/  BW,QP,QC,VFC,VLC,VMC,QFC,QLC,QMC,PLA,PFA,         &
     &                  PMA,PTA,PB,MW,VMAX,KM,CONC,KL,RATS,VCHC,
     &                  VCH,VM,VT,VF,VL,PM,PT,PF,PL,AI0,VTC,QM,QT,QF,QL
                        

c     Fraction viscera (kg/(kg BW))
      VTC = 0.91 - (VLC+VFC+VMC)

c     Net chamber volume
      VCH = VCHC - RATS*BW 
      VM  = VMC*BW
      VT  = VTC*BW
      VF  = VFC*BW
      VL  = VLC*BW

c     Initial amt. in chamber (mg)         
      AI0 = CONC*VCH*MW/24450. 
      PL  = PLA/PB
      PF  = PFA/PB
      PT  = PTA/PB
      PM  = PMA/PB
      
      QF = QFC*QC
      QL = QLC*QC
      QM = QMC*QC
    
      QT = QC - (QF+QL+QM)

      return 
      end subroutine derived

c=======================================================================
c The dynamic model
c=======================================================================

      subroutine derivsccl4 (neq, t, y, ydot,out,IP)
      implicit none
      integer          neq, IP(*), i
      double precision t, y(neq), ydot(neq), out(*)

      double precision  BW,QP,QC,VFC,VLC,VMC,QFC,QLC,QMC,PLA,PFA,         &
     &                  PMA,PTA,PB,MW,VMAX,KM,CONC,KL,RATS,VCHC,          &
     &                  V(5), P(4),AI0,VTC,Q(4)

c here we lump parameters Vx, Qx and Px into vectors    
      common /myparms/  BW,QP,QC,VFC,VLC,VMC,QFC,QLC,QMC,PLA,PFA,         &
     &                  PMA,PTA,PB,MW,VMAX,KM,CONC,KL,RATS,VCHC,          &
     &                  V,    P,   AI0,VTC,  Q 


      double precision tconc(5), vconc(5), dose, mass, cp, ca, cx, RAM



      
c check if provision has been made for at least 3 output variables 
      if (IP(1) < 3) call rexit("nout should be at least 3") 
      
c y = AI, AAM, AT, AF, AL CLT, AM
c where clt = the area under the concentration-time curve in the liver 
c       AM  = total amount metabolised

c     concentrations
      do i =1,5
        tconc(i) = y(i)/v(i)
      enddo
   
c     vconc(1) is conc in mixed venous blood
      vconc(1) = 0.d0
      do i = 2,5
       vconc(i) = tconc(i)/P(i-1)
       vconc(1) = vconc(1) + vconc(i)*Q(i-1)/QC
      enddo
      
c     CA is conc in arterial blood
      CA = (QC * Vconc(1) + QP * tconc(1))/ (QC + QP/PB)

c     Exhaled chemical
      CX = CA/PB

c     metabolisation rate 
      RAM = VMAX*Vconc(5)/(KM + Vconc(5))
  
c     the rate of change
      ydot(1) = RATS*QP*(CX - tconc(1)) - KL*y(1)
      do i = 2,5
        ydot(i) = Q(i-1)*(CA-vconc(i))
      enddo
      ydot(5) = ydot(5) - RAM
      ydot(6) = tconc(5)
      ydot(7) = RAM
      
c the mass balance (MASS=AAM+AT+AF+AL+AM), should be constant
      DOSE =  AI0 - y(1)
      MASS = (y(2)+y(3)+y(4)+y(5)+y(7))*RATS
      CP   = tconc(1)*24450.0/MW
      
      out(1) = DOSE
      out(2) = MASS
      out(3) = CP
        
      return
      end
