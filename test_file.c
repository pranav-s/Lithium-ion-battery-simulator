#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>
#include <ida/ida_band.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
//#include "half_cell_solver.c"

#define GRID        5
#define N_VAR       5    
#define N           GRID*N_VAR 
#define ZERO        RCONST(0.0)
#define ONE         RCONST(1.0)
#define T_PLUS      RCONST(0.363)
#define R           RCONST(8.314)
#define T           RCONST(298.15)
#define F           RCONST(96487.0)
#define BRUGG       RCONST(4.0)
#define I           RCONST(2.0)
 
typedef struct {
  realtype dx;
  realtype coeff;
  realtype eps;
  realtype sigma;
  realtype diff_coeff;
  realtype radius;
  realtype k;
  realtype c_s_max;
  realtype c_s_0;
  realtype c_0;
  realtype l;
  realtype interfac_area;
  realtype sigma_eff;
  realtype diff_coeff_eff;
} *Material_Data;

realtype ocp_anode(realtype c,realtype c_max);
realtype ocp_cathode(realtype c,realtype c_max);
realtype kappa(realtype c, realtype eps);
realtype Rsinh(realtype x);
realtype Rlog(realtype x);
static void InitAnodeData(Material_Data data_anode);

int main(void){
  realtype b,l_a;
  Material_Data data_anode;
  b = Rlog(I);
  data_anode=NULL;
  data_anode = (Material_Data) malloc(sizeof *data_anode);
  InitAnodeData(data_anode);
  l_a = (data_anode->l)/(data_anode->l + data_anode->radius);
  int i = (int)(l_a*RCONST(GRID));
  printf("%3.8f\n",b);
  printf("%d\n",i);
  printf("%3.8f\n",data_anode->coeff);
  printf("%3.8f\n",T_PLUS);
  printf("%3.8f\n",R);
  printf("%3.8f\n",data_anode->interfac_area);
  return 0;

}

static void InitAnodeData(Material_Data data_anode)
{  
  data_anode->dx = ONE/(GRID - ONE);
  data_anode->coeff = ONE/(data_anode->dx); 
  data_anode->sigma = RCONST(100.0);
  data_anode->eps = RCONST(0.385);
  data_anode->sigma_eff = (data_anode->sigma)*(1-data_anode->eps);
  data_anode->diff_coeff = RCONST(1.0e-14);
  data_anode->diff_coeff_eff = (data_anode->diff_coeff)*SUNRpowerR((data_anode->eps),BRUGG);
  data_anode->k = RCONST(2.334e-11);
  data_anode->c_0 = RCONST(1000.0);
  data_anode->c_s_max = RCONST(51554.0);
  data_anode->c_s_0 = (data_anode->c_s_max)*RCONST(0.4955);
  data_anode->l = RCONST(0.00008);
  data_anode->radius = RCONST(0.000002);
  data_anode->interfac_area = RCONST(3.0)*(1-data_anode->eps)/(data_anode->radius);
}

realtype kappa(realtype c, realtype eps)
{
    realtype exp_term, polynom_term;
    exp_term = SUNRpowerR(eps,BRUGG);
    polynom_term = RCONST(4.1253e-2)+RCONST(5.007e-4)*c-RCONST(4.7212e-7)*SUNRpowerI(c,2)+
                   RCONST(1.5094)*SUNRpowerI(c,3)-RCONST(1.6018)*SUNRpowerI(c,4);
    
    return polynom_term*exp_term;
  
}

realtype Rlog(realtype x)
{
  if (x <= ZERO)
  {
    return(ZERO);
  }
 
 #if defined(SUNDIALS_USE_GENERIC_MATH)
  return((realtype) log((double) x));
 #elif defined(SUNDIALS_DOUBLE_PRECISION)
  return(log(x));
 #elif defined(SUNDIALS_SINGLE_PRECISION)
  return(log(x));
 #elif defined(SUNDIALS_EXTENDED_PRECISION)
  return(log(x));
 #endif
}


realtype ocp_anode(realtype c,realtype c_max)
{
  realtype soc = c/c_max;
  realtype numer,denom;
  numer = -RCONST(4.656)+RCONST(88.669)*SUNRpowerI(soc,2)-RCONST(401.119)*SUNRpowerI(soc,4)+
           RCONST(342.909)*SUNRpowerI(c,6)-RCONST(462.471)*SUNRpowerI(c,8)+RCONST(433.434)*SUNRpowerI(c,10);
  
  denom = -RCONST(1.0)+RCONST(18.933)*SUNRpowerI(soc,2)-RCONST(79.532)*SUNRpowerI(soc,4)+
           RCONST(37.311)*SUNRpowerI(c,6)-RCONST(73.083)*SUNRpowerI(c,8)+RCONST(95.96)*SUNRpowerI(c,10);
           
  return (numer/denom);
}

realtype Rsinh(realtype x)
{
 
 #if defined(SUNDIALS_USE_GENERIC_MATH)
  return((realtype) sinh((double) x));
 #elif defined(SUNDIALS_DOUBLE_PRECISION)
  return(sinh(x));
 #elif defined(SUNDIALS_SINGLE_PRECISION)
  return(sinh(x));
 #elif defined(SUNDIALS_EXTENDED_PRECISION)
  return(sinh(x));
 #endif
}

realtype ocp_cathode(realtype c,realtype c_max)
{
  realtype soc = c/c_max;
  realtype expr;
  expr = RCONST(0.7222)+RCONST(0.1387)*soc+RCONST(0.029)*SUNRpowerR(soc,RCONST(0.5))-
         RCONST(0.0172)/soc+RCONST(0.0019)*SUNRpowerR(soc,RCONST(-1.5))+
         RCONST(0.2808)*SUNRpowerR(RCONST(10.0),RCONST(0.90)-RCONST(15.0)*soc)-
         RCONST(0.7984)*SUNRpowerR(RCONST(10.0),-RCONST(0.4108)+RCONST(0.4465)*soc);
  return expr;
}

