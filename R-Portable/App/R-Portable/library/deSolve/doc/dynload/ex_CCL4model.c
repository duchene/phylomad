/* c the CCl4 inhalation model
 -------- ex_ccl4model.c -> ex_ccl4model.dll ------
 compile in R with: system("gcc -shared -o ex_ccl4model.dll ex_ccl4model.c")
 or with system("R CMD SHLIB ex_ccl4model.c")
*/

#include <R.h>
static double parms[21];
#define BW     parms[0]
#define QP     parms[1]
#define QC     parms[2]
#define VFC    parms[3]
#define VLC    parms[4]
#define VMC    parms[5]
#define QFC    parms[6]
#define QLC    parms[7]
#define QMC    parms[8]
#define PLA    parms[9]
#define PFA    parms[10]
#define PMA    parms[11]
#define PTA    parms[12]
#define PB     parms[13]
#define MW     parms[14]
#define VMAX   parms[15]
#define KM     parms[16]
#define CONC   parms[17]
#define KL     parms[18]
#define RATS   parms[19]
#define VCHC   parms[20]

double V[5], P[4], AI0, VTC, Q[4];

#define DOSE   out[0]
#define MASS   out[1]
#define CP     out[2]

/* 
c=======================================================================
c=======================================================================
c                      Model initialisation
c=======================================================================
c=======================================================================

c=======================================================================
2c Initialise primary parameter common block
c=======================================================================
 */

void initccl4(void (* odeparms)(int *, double *))
{
    void derived();
    
    int N=21;
    odeparms(&N, parms);
    derived();           

}

/*=======================================================================
 In this "event", state variable 1 is increased with 1. DOES NOT WORK...
=======================================================================
*/

void eventfun(int *n, double *t, double *y)
{
      y[0] = y[0] + 1;
}
      
/*=======================================================================
c Calculate derived parameters from primary parameters
c=======================================================================
*/

void derived ()
{
//     Fraction viscera (kg/(kg BW))
      VTC = 0.91 - (VLC+VFC+VMC);

//     Net chamber volume

      V[0] = VCHC - RATS*BW; 
      V[1] = VMC*BW;
      V[2] = VTC*BW;
      V[3] = VFC*BW;
      V[4] = VLC*BW;

//     Initial amt. in chamber (mg)         
      AI0 = CONC*V[0]*MW/24450.; 
      P[0]  = PMA/PB;
      P[1]  = PTA/PB;
      P[2]  = PFA/PB;
      P[3]  = PLA/PB;
      
      Q[2] = QFC*QC;
      Q[3] = QLC*QC;
      Q[0] = QMC*QC;
    
      Q[1] = QC - (Q[0]+Q[3]+Q[2]);

}

/*=======================================================================
c The dynamic model
c=======================================================================
*/
void derivsccl4 (int *neq, double *t, double *y, double *ydot,
             double *out, int *ip)
{             
      double vconc[5], tconc[5], CA, CX, RAM;
      int i;
      
      if (ip[0] < 3) error("nout should be at least 3"); 
      
/*c y = AI, AAM, AT, AF, AL CLT, AM
 where clt = the area under the concentration-time curve in the liver 
       AM  = total amount metabolised

     concentrations */
      for (i =0; i<5; i++) {
        tconc[i] = y[i]/V[i];
      }
   
/*     vconc(1) is conc in mixed venous blood */
      vconc[0] = 0.0;
      for (i = 1; i<5; i++){
       vconc[i] = tconc[i]/P[i-1];
       vconc[0] = vconc[0] + vconc[i]*Q[i-1]/QC ;
      }
      
/*     CA is conc in arterial blood  */
      CA = (QC * vconc[0] + QP * tconc[0])/ (QC + QP/PB);

/*     Exhaled chemical   */
      CX = CA/PB;

/*     metabolisation rate   */
      RAM = VMAX*vconc[4]/(KM + vconc[4]);
  
/*      the rate of change    */
      ydot[0] = RATS*QP*(CX - tconc[0]) - KL*y[0];
      for ( i = 1; i<5; i++)
        ydot[i] = Q[i-1]*(CA-vconc[i]);      

      ydot[4] = ydot[4] - RAM;
      ydot[5] = tconc[4];
      ydot[6] = RAM;
      
/*  the mass balance (MASS=AAM+AT+AF+AL+AM), should be constant */
      DOSE =  AI0 - y[0];
      MASS = (y[1]+y[2]+y[3]+y[4]+y[6])*RATS;
      CP   = tconc[0]*24450.0/MW;
}
