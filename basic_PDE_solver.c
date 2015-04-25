#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>
#include <ida/ida_band.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

//Define all the global problem constants. Need to later find a way around global variables. Can use a local struct in which default values have been initialized

#define alpha       RCONST(1.0)
#define N           50
#define delta       RCONST(0.02)
#define RTOL        RCONST(1.e-3)  /* Relative tolerance                */
#define ATOL        RCONST(1.e-3)  /* Absolute tolerance                */
#define ZERO        RCONST(0.0)
#define ONE         RCONST(1.0)
#define TWO         RCONST(2.0)
#define inv         RCONST(1/(delta*delta))

//Declare prototypes of all functions defined in this file


static int check_flag(void *flagvalue, char *funcname, int opt);
static void SetInitialProfiles(N_Vector cc,N_Vector cp,N_Vector id, N_Vector res);
static void PrintOutput(void *mem, N_Vector c, realtype t);
static int res_basic(realtype tt, N_Vector cc, N_Vector cp, N_Vector resval);


int main()
{
//Declare any necessary local constants
N_Vector cc,cp,id,res;
long int mu,ml,iter;
void *mem;
realtype rtol, atol, t0, tout,tret,incr, *ccp;
int iout,retval;

//Initialize all the relevant N-vectors

  cc  = N_VNew_Serial(N);
  if(check_flag((void *)cc, "N_VNew_Serial", 0)) return(1);

  cp  = N_VNew_Serial(N);
  if(check_flag((void *)cp, "N_VNew_Serial", 0)) return(1);
  
  id  = N_VNew_Serial(N);
  if(check_flag((void *)id, "N_VNew_Serial", 0)) return(1);
  
  res = N_VNew_Serial(N);
  if(check_flag((void *)res, "N_VNew_Serial", 0)) return(1);
  
  SetInitialProfiles(cc,cp,id,res);  
  
  
  t0 = ZERO;
  rtol = RTOL; 
  atol = ATOL;
  iter= 10;
  incr= RCONST(0.1);
  tout= RCONST(0.01);

/* Call IDACreate and IDAMalloc to initialize IDA. */
  
  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);

  retval = IDASetId(mem, id);
  if(check_flag(&retval, "IDASetId", 1)) return(1);

  retval = IDAInit(mem, res_basic, t0, cc, cp);
  if(check_flag(&retval, "IDAInit", 1)) return(1);

  retval = IDASStolerances(mem, rtol, atol);
  if(check_flag(&retval, "IDASStolerances", 1)) return(1);

  /* Call IDABand to specify the IDA linear solver. */

  mu = 1;
  ml = 1;
  retval = IDABand(mem, N, mu, ml);
  if(check_flag(&retval, "IDABand", 1)) return(1);

  /* Call IDACalcIC (with default options) to correct the initial values. */
  //printf("%d\n",N);
  //printf("%f\n",inv);
  ccp = NV_DATA_S(cp);
  printf("%f\n",ccp[49]);
  //retval = IDACalcIC(mem, IDA_YA_YDP_INIT, tout);
  //if(check_flag(&retval, "IDACalcIC", 1)) return(1);
  
  for (iout = 1; iout <= iter; iout++) {
    
    retval = IDASolve(mem, tout, &tret, cc, cp, IDA_NORMAL);
    if(check_flag(&retval, "IDASolve", 1)) return(retval);
    
    PrintOutput(mem, cc, tret);
    
    tout=tout+incr;
    
  }
  //Free memory

  IDAFree(&mem);

  N_VDestroy_Serial(cc);
  N_VDestroy_Serial(cp);
  N_VDestroy_Serial(id);
  
  return(0);
}

static void SetInitialProfiles(N_Vector cc,N_Vector cp, N_Vector id, N_Vector res)
{
    realtype *ccv,*cpv,*idv,xfact;
    int i;
    ccv = NV_DATA_S(cc);
    cpv = NV_DATA_S(cp);
    idv = NV_DATA_S(id);
    
    N_VConst(ZERO, cp);
    N_VConst(ONE, id);
    
    for (i=0; i<N; i++){
         xfact = (ONE/(N-ONE)) * i;  
         ccv[i]=RCONST(16.0) * xfact * (ONE - xfact);;
    }
    res_basic(ZERO, cc, cp, res);
    N_VScale(-ONE, res, cp);
    
    cpv[0]=ZERO;
    cpv[N-1]=ZERO;
    idv[0]=ZERO;
    idv[N-1]=ZERO;
    
}


static int res_basic(realtype tt, N_Vector cc, N_Vector cp, N_Vector resval)
{   
    realtype *ccv,*cpv,*res;
    int i;
    N_VScale(ONE, cc, resval);
    ccv = NV_DATA_S(cc);
    cpv = NV_DATA_S(cp);
    res = NV_DATA_S(resval);
    //res[0]=ccv[0]-ZERO;
    //res[N-1]=ccv[N-1]-ZERO;
    
    for (i=1; i<N-1;i++)
     {
        res[i]=cpv[i]-alpha*inv*(ccv[i+1]+ccv[i-1]-RCONST(2.0)*ccv[i]);
        
      }
return 0;

}


static void PrintOutput(void *mem, N_Vector c, realtype t)
{     
      realtype *conc;
      int i;
      conc = NV_DATA_S(c);
      printf("Output at time %3.2f",t);
      printf("\n");
      for (i=0;i<N;i++){
          
          printf("%4.6F",conc[i]);
          printf("\n");  
      }
      printf("\n\n\n\n\n\n");
}


static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  if (opt == 0 && flagvalue == NULL) {
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  } 
  else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
              funcname, *errflag);
      return(1); 
    }
  } 
  else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1); 
  }

  return(0);
}
