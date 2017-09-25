#include <R.h>
static double parms[16];
static double forc[1];

#define Vc      parms[0]  /* Volume of central compartment (L)                */
#define Vt      parms[1]  /* Volume of second compartment  (L)                */
#define kd      parms[2]  /* 1st order rate constant central <-> second       */
                          /*   cmpt (1/hr)                                    */
#define ka      parms[3]  /* absorption 1st order rate constant (1/hr)        */
#define Tm      parms[4]  /* 0 order resorption rate in the limit of          */
                          /*  increasing filtrate PFOA concentrations         */
                          /*  (mg/L/hr)                                       */
#define KT      parms[5]  /* Filtrate cmpt concentration at which             */
                          /*  resorption rate is half maximal                 */
                          /*  (mg/L)                                          */
#define kfil    parms[6]  /* 1st order rate constant central -> filtrate      */
                          /* cmpartment (1/hour)                              */
#define Vfil    parms[7]  /* Volume of filtrate compartment  (L)              */
#define free    parms[8]  /* Free fraction PFOA in central compartment (-)    */
#define BW      parms[9]  /* bodyweight (kg)                                  */
#define Dose    parms[10] /* dose (mg/kg/day)                                 */
#define Doseint parms[11] /* interval between doses (hours)                   */
#define Qd      parms[12] /* Clearance (kd * Vc) central <-> 2nd cmpt (L/hr)  */
#define Qfil    parms[13] /* rate of flow to filtrate compartment             */
#define MaxTime parms[14] /* Duration of simulation                           */
#define TDose   parms[15] /* actual dose (dose * BW) (mg/day)                 */

#define TDoseRt forc[0]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
	int N = 16;
	odeparms(&N, parms);
}

void initforc(void (* odeforcs)(int *, double *)) {
	int N = 1;
	odeforcs(&N, forc);
}

/* Compartments are:
   Cn for the central compartment
   Tc for the second comparment
   Fc for the filtrate compartment
   Gt for the gut
 Elim for total eliminated
  AUC for AUC in the central compartment
*/

#define Cn y[0]
#define Tc y[1]
#define Fc y[2]
#define Gt y[3]
#define Elim y[4]
#define AUC y[5]

#define Cn_dot ydot[0]
#define Tc_dot ydot[1]
#define Fc_dot ydot[2]
#define Gt_dot ydot[3]
#define Elim_dot ydot[4]
#define AUC_dot ydot[5]

#define MassBal yout[0]

/* Derivatives and one output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
	     double *yout, int *ip) {
  if (ip[0] < 1) error("nout should be at least 1");
  Cn_dot = (ka * Gt - Qd * free * Cn + Qd * Tc) / Vc - kfil * Cn * free +
    Tm * Fc/(KT + Fc);
  Tc_dot = (Qd * free * Cn - Qd * Tc) / Vt;
  Fc_dot = (Vc * kfil * Cn * free - Vc * Tm * Fc/(KT + Fc) -
    Vc * kfil * Fc) / Vfil;
  Gt_dot = -ka * Gt + TDoseRt;
  Elim_dot = Vc * kfil * Fc;
  AUC_dot = Cn;
  /* Total amount in all compartments, for mass balance */
  MassBal = Cn * Vc + Tc * Vt + Fc * Vfil + Gt + Elim;
}
