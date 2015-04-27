#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>
#include <ida/ida_band.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h> 

realtype ocp_anode(realtype c,realtype c_max); //Consider using std library functions
realtype ocp_cathode(realtype c,realtype c_max);

int main(void){
  realtype b;
  b=ocp_cathode(RCONST(5),RCONST(55));
  printf("%3.2f\n",b);
  return 0;

}


realtype ocp_anode(realtype c,realtype c_max) //Consider using std library functions
{
  realtype soc = c/c_max;
  realtype numer,denom;
  numer = -RCONST(4.656)+RCONST(88.669)*SUNRpowerI(soc,2)-RCONST(401.119)*SUNRpowerI(soc,4)+
           RCONST(342.909)*SUNRpowerI(c,6)-RCONST(462.471)*SUNRpowerI(c,8)+RCONST(433.434)*SUNRpowerI(c,10);
  
  denom = -RCONST(1.0)+RCONST(18.933)*SUNRpowerI(soc,2)-RCONST(79.532)*SUNRpowerI(soc,4)+
           RCONST(37.311)*SUNRpowerI(c,6)-RCONST(73.083)*SUNRpowerI(c,8)+RCONST(95.96)*SUNRpowerI(c,10);
           
  return (numer/denom);
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

