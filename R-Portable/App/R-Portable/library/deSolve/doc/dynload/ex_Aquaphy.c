/* file ex_aquaphy.c  
 The Aquaphy algal model

 -------- ex_Aquaphy.c -> ex_Aquaphy.dll ------
 compile in R with: system("gcc -shared -o Aquaphy Aquaphy")
 or with system("R CMD SHLIB ex_Aquaphy")
*/

#include <R.h>
static double parms[19];
#define maxPhotoSynt     parms[0]
#define rMortPHY         parms[1]
#define alpha            parms[2]
#define pExudation       parms[3]
#define maxProteinSynt   parms[4]
#define ksDIN            parms[5] 
#define minpLMW          parms[6]
#define maxpLMW          parms[7]
#define minQuotum        parms[8]
#define maxStorage       parms[9]
#define respirationRate  parms[10]
#define pResp            parms[11]
#define catabolismRate   parms[12]
#define dilutionRate     parms[13]
#define rNCProtein       parms[14]
#define inputDIN         parms[15]
#define rChlN            parms[16]
#define parMean          parms[17]
#define dayLength        parms[18]

static double forcs[1];
#define Light            forcs[0]

#define DIN              y[0]
#define PROTEIN          y[1]
#define RESERVE          y[2]
#define LMW              y[3]

#define dDIN             ydot[0]
#define dPROTEIN         ydot[1]
#define dRESERVE         ydot[2]
#define dLMW             ydot[3]

#define PAR              out[0]
#define TotalN           out[1]
#define PhotoSynthesis   out[2]
#define NCratio          out[3]
#define ChlCratio        out[4]
#define Chlorophyll      out[5]

/*=======================================================================
c=======================================================================
c                      Model initialisation
c=======================================================================
c=======================================================================

c=======================================================================
c Initialise parameter common block
c=======================================================================
*/
void iniaqua(void (* odeparms)(int *, double *))
{
    int N=19;
    odeparms(&N, parms);
}

/*
c=======================================================================
c Initialise forcing function common block
c=======================================================================
*/
void initaqforc(void (* odeforc)(int *, double *))
{
    int N=1;
    odeforc(&N, forcs);
}

/*
c=======================================================================
c Algal dynamics - light an on-off function
c=======================================================================
*/
void aquaphy (int *neq, double *t, double *y, double *ydot, double *out, int *ip)

{     
double PhytoC,PhytoN,PartLMW,Limfac,Exudation,MonodQuotum,hourofday,                     
       ProteinSynthesis,Storage,Respiration,Catabolism;
        
      if(ip[0] < 6) error("nout should at least be 6");
      

/* PAR, on-off function depending on the hour within a day*/
      hourofday       = fmod(*t,24.0);
      if (hourofday  < dayLength)  
       PAR = parMean;
      else
       PAR = 0.0;
       

/*  the output variables - all components contain carbon
  only proteins contain nitrogen   */

      PhytoC           = PROTEIN + RESERVE + LMW;        
      PhytoN           = PROTEIN * rNCProtein;          
      NCratio          = PhytoN / PhytoC;
      Chlorophyll      = PhytoN * rChlN;
      TotalN           = PhytoN + DIN;
      ChlCratio        = Chlorophyll / PhytoC;

/* the rates, in mmol/hr */
      PartLMW         = LMW / PhytoC;
      Limfac          = fmin(1.0,(maxpLMW -PartLMW)/(maxpLMW-minpLMW));
      Limfac          = fmax(0.0,Limfac);
      PhotoSynthesis  = maxPhotoSynt*Limfac *                             
                       (1.0-exp(alpha*PAR/maxPhotoSynt)) * PROTEIN;
      Exudation       = pExudation * PhotoSynthesis;
      MonodQuotum     = fmax(0.0,LMW / PROTEIN - minQuotum);
      ProteinSynthesis= maxProteinSynt*MonodQuotum                        
                       * DIN / (DIN+ksDIN)            * PROTEIN;
      Storage         = maxStorage    *MonodQuotum     * PROTEIN;
      Respiration     = respirationRate * LMW        
                     + pResp * ProteinSynthesis;
      Catabolism      = catabolismRate  * RESERVE;

/* the rates of change of state variables; includes dilution effects (last term) */
      dLMW     = PhotoSynthesis + Catabolism                              
              - Exudation - Storage  - Respiration - ProteinSynthesis     
              - dilutionRate * LMW;

      dRESERVE = Storage - Catabolism    - dilutionRate * RESERVE;

      dPROTEIN = ProteinSynthesis        - dilutionRate * PROTEIN;

      dDIN     = -ProteinSynthesis * rNCProtein                            
               - dilutionRate * (DIN - inputDIN);

}
/* Algal dynamics with forcings
c======================================================================= */

void aquaphyforc (int *neq, double *t, double *y, double *ydot, double *out, int *ip)
{             

double PhytoC,PhytoN,PartLMW,Limfac,Exudation,MonodQuotum,                      
       ProteinSynthesis,Storage,Respiration,Catabolism;

      if(ip[0] < 6) error("nout should at least be 6");

      PAR = Light;
      
/*  the output variables - all components contain carbon
  only proteins contain nitrogen   */

      PhytoC           = PROTEIN + RESERVE + LMW;       
      PhytoN           = PROTEIN * rNCProtein;          
      NCratio          = PhytoN / PhytoC;
      Chlorophyll      = PhytoN * rChlN;
      TotalN           = PhytoN + DIN;
      ChlCratio        = Chlorophyll / PhytoC;

/* the rates, in mmol/hr */
      PartLMW         = LMW / PhytoC;
      Limfac          = fmin(1.0,(maxpLMW -PartLMW)/(maxpLMW-minpLMW));
      Limfac          = fmax(0.0,Limfac);
      PhotoSynthesis  = maxPhotoSynt*Limfac *                              
                       (1.0-exp(alpha*PAR/maxPhotoSynt)) * PROTEIN;
      Exudation       = pExudation * PhotoSynthesis;
      MonodQuotum     = fmax(0.0,LMW / PROTEIN - minQuotum);
      ProteinSynthesis= maxProteinSynt*MonodQuotum                         
                        * DIN / (DIN+ksDIN)            * PROTEIN;
      Storage         = maxStorage    *MonodQuotum     * PROTEIN;
      Respiration     = respirationRate * LMW                              
                      + pResp * ProteinSynthesis;
      Catabolism      = catabolismRate  * RESERVE;

/* the rates of change of state variables; includes dilution effects (last term) */
      dLMW     = PhotoSynthesis + Catabolism                              
               - Exudation - Storage  - Respiration - ProteinSynthesis    
               - dilutionRate * LMW;

      dRESERVE = Storage - Catabolism    - dilutionRate * RESERVE;

      dPROTEIN = ProteinSynthesis        - dilutionRate * PROTEIN;

      dDIN     = -ProteinSynthesis * rNCProtein                            
                - dilutionRate * (DIN - inputDIN);

}

