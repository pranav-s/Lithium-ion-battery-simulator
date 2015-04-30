///\file This file contains all the functions used by full_cell_solver_new.c


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>
#include <ida/ida_spgmr.h>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>


#define ZERO        RCONST(0.0)
#define ONE         RCONST(1.0)
#define TWO         RCONST(2.0)
#define T_PLUS      RCONST(0.363)///<Transport number in anode
#define R           RCONST(8.314)///<Universal Gas constant
#define T           RCONST(298.15)///<Temperature
#define F           RCONST(96487.0)///<Faraday's constant
#define BRUGG       RCONST(4.0) ///<Bruggman's coefficient. Used to account for effect of porosity in constants

realtype kappa(realtype c, realtype eps)
{
    realtype exp_term, polynom_term;
    exp_term = SUNRpowerR(eps,BRUGG);
    polynom_term = RCONST(4.1253e-2)+RCONST(5.007e-4)*c-RCONST(4.7212e-7)*SUNRpowerI(c,2)+
                   RCONST(1.5094)*SUNRpowerI(c,3)-RCONST(1.6018)*SUNRpowerI(c,4);
    
    return polynom_term*exp_term;
  
}

///Returns the open circuit potential in the anode at a given point in either electrode domain as a function of the state of charge 

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

///Returns the open circuit potential in the cathode at a given point in either electrode domain as a function of the state of charge 
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
///Returns logarithm of a realtype number
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
///Returns sinh of a realtype number
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

