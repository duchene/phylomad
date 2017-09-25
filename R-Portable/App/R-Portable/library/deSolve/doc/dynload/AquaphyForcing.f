

c the Aquaphy algal model with forcing function light intensity

c -------- Aquaphy2.f -> Aquaphy2.dll ------
c compile in R with: system("g77 -shared -o Aquaphy Aquaphy")
c or with system("R CMD SHLIB AquaphyForcing")


c=======================================================================
c=======================================================================
c                      Model initialisation
c=======================================================================
c=======================================================================

c=======================================================================
c Initialise parameter common block
c=======================================================================

      subroutine initaqparms(odeparms)

      external odeparms
      double precision pars(16)
      common /myparms/pars

       call odeparms(16, pars)

      return
      end


      subroutine initaqforc(odeforc)

      external odeparms
      double precision forcs(2)
      common /myforcs/forcs

       call odeforc(2, forcs)

      return
      end

c=======================================================================
c In this "event", state variable 1 is increased with 1. DOES NOT WORK...
c=======================================================================

      subroutine eventfun(n, t, y)
      integer n
      double precision t, y(n)

      y(1) = y(1) + 1
      
      end subroutine
      

c=======================================================================
c Algal dynamics
c=======================================================================

      subroutine  aquaphy2 (neq, t, y, ydot,out,IP)
      implicit none
      
      integer           neq, ip(*)
      double precision  t, y(*), ydot(*), out(*)
      
c parameters      
      double precision  maxPhotoSynt,rMortPHY,alpha,pExudation,          &
     &     maxProteinSynt,ksDIN,minpLMW,maxpLMW,minQuotum,               &
     &     maxStorage,respirationRate,pResp,catabolismRate,              &
     &     rNCProtein,inputDIN,rChlN

      common /myparms/  maxPhotoSynt,rMortPHY,alpha,pExudation,          &
     &     maxProteinSynt,ksDIN,minpLMW,maxpLMW,minQuotum,               &
     &     maxStorage,respirationRate,pResp,catabolismRate,              &
     &     rNCProtein,inputDIN,rChlN
	double precision PAR, dilutionRate
      common /myforcs/PAR, dilutionRate
		  
c variables
      double precision ::                                                &                      
     &  DIN,PROTEIN,RESERVE,LMW,dLMW,dRESERVE,dPROTEIN,dDIN,             &
     &  PhytoC,PhytoN,NCratio,Chlorophyll,TotalN,ChlCratio,PartLMW,      &      
     &  hourofday, Limfac,PhotoSynthesis,Exudation,MonodQuotum,          &
     &  ProteinSynthesis,Storage,Respiration,Catabolism

c ------------------------------------------------------------------------
      if(ip(1) < 6) call rexit("nout should at least be 6")
      
      DIN      = y(1)
      PROTEIN  = y(2)
      RESERVE  = y(3)
      LMW      = y(4)

c the output variables
      PhytoC           = PROTEIN + RESERVE + LMW       ! all components contain carbon
      PhytoN           = PROTEIN * rNCProtein          ! only proteins contain nitrogen
      NCratio          = PhytoN / PhytoC
      Chlorophyll      = PhytoN * rChlN
      TotalN           = PhytoN + DIN
      ChlCratio        = Chlorophyll / PhytoC

c the rates, in mmol/hr
      PartLMW         = LMW / PhytoC
      Limfac          = min(1.d0,(maxpLMW -PartLMW)/(maxpLMW-minpLMW))
      Limfac          = max(0.d0,Limfac)
      PhotoSynthesis  = maxPhotoSynt*Limfac *                             &
     &                 (1.d0-exp(alpha*PAR/maxPhotoSynt)) * PROTEIN
      Exudation       = pExudation * PhotoSynthesis
      MonodQuotum     = max(0.d0,LMW / PROTEIN - minQuotum)
      ProteinSynthesis= maxProteinSynt*MonodQuotum                        &
     &                  * DIN / (DIN+ksDIN)            * PROTEIN
      Storage         = maxStorage    *MonodQuotum     * PROTEIN
      Respiration     = respirationRate * LMW                             &
     &                + pResp * ProteinSynthesis
      Catabolism      = catabolismRate  * RESERVE

c the rates of change of state variables; includes dilution effects (last term)
      dLMW     = PhotoSynthesis + Catabolism                              &
     &         - Exudation - Storage  - Respiration - ProteinSynthesis    &
     &         - dilutionRate * LMW

      dRESERVE = Storage - Catabolism    - dilutionRate * RESERVE

      dPROTEIN = ProteinSynthesis        - dilutionRate * PROTEIN

      dDIN     = -ProteinSynthesis * rNCProtein                            &
     &          - dilutionRate * (DIN - inputDIN)

c the vector with rate of changes
      ydot(1) = dDIN
      ydot(2) = dPROTEIN
      ydot(3) = dRESERVE
      ydot(4) = dLMW
      
c the ordinary variables
      out(1) = PAR
      out(2) = TotalN
      out(3) = PhotoSynthesis
      out(4) = NCratio
      out(5) = ChlCratio
      out(6) = Chlorophyll

      return
      end subroutine Aquaphy2



      
