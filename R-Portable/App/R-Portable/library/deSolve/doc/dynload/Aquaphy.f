

c the Aquaphy algal model

c -------- Aquaphy.f -> Aquaphy.dll ------
c compile in R with: system("g77 -shared -o Aquaphy Aquaphy")
c or with system("R CMD SHLIB Aquaphy")


c=======================================================================
c=======================================================================
c                      Model initialisation
c=======================================================================
c=======================================================================

c=======================================================================
c Initialise parameter common block
c=======================================================================

      subroutine initaquaphy(odeparms)

      external odeparms
      double precision pars(19)
      common /myparms/pars

       call odeparms(19, pars)

      return
      end

c=======================================================================
c Algal dynamics
c=======================================================================

      subroutine  aquaphy (neq, t, y, ydot,out,IP)
      implicit none
      
      integer           neq, ip(*)
      double precision  t, y(*), ydot(*), out(*)
      
c parameters      
      double precision  maxPhotoSynt,rMortPHY,alpha,pExudation,          &
     &     maxProteinSynt,ksDIN,minpLMW,maxpLMW,minQuotum,               &
     &     maxStorage,respirationRate,pResp,catabolismRate,              &
     &     dilutionRate,rNCProtein,inputDIN,rChlN,parMean,               &
     &     dayLength
      common /myparms/  maxPhotoSynt,rMortPHY,alpha,pExudation,          &
     &     maxProteinSynt,ksDIN,minpLMW,maxpLMW,minQuotum,               &
     &     maxStorage,respirationRate,pResp,catabolismRate,              &
     &     dilutionRate,rNCProtein,inputDIN,rChlN,parMean,               &
     &     dayLength
c variables
      double precision ::                                                &                      
     &  DIN,PROTEIN,RESERVE,LMW,dLMW,dRESERVE,dPROTEIN,dDIN,PAR,         &
     &  PhytoC,PhytoN,NCratio,Chlorophyll,TotalN,ChlCratio,PartLMW,      &      
     &  hourofday, Limfac,PhotoSynthesis,Exudation,MonodQuotum,          &
     &  ProteinSynthesis,Storage,Respiration,Catabolism

c ------------------------------------------------------------------------
      if(ip(1) < 6) call rexit("nout should at least be 6")
      
      DIN      = y(1)
      PROTEIN  = y(2)
      RESERVE  = y(3)
      LMW      = y(4)

c PAR, on-off function depending on the hour within a day
      hourofday       = mod(t,24.d0)
      if (hourofday  < dayLength) THEN
       PAR = parMean
      else
       PAR = 0.d0
      endif

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
      end subroutine Aquaphy



      
