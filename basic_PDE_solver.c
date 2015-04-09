#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>
#include <ida/ida_band.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

//Define all the global problem constants. Need to later find a way around global variables. Can use a local struct in which default values have been initialized

#define alpha       RCONST(0.5)
#define N           50
#define delta       RCONST(1/RCONST(N))
#define RTOL        RCONST(1.e-5)  /* Relative tolerance                */
#define ATOL        RCONST(1.e-5)  /* Absolute tolerance                */
#define ZERO        RCONST(0.)

//Declare prototypes of all functions defined in this file


static int check_flag(void *flagvalue, char *funcname, int opt);
static void SetInitialProfiles(N_Vector cc,N_Vector cp);
static void PrintOutput(void *mem, N_Vector c, realtype t);
static int res_basic(realtype tt, N_Vector cc, N_Vector cp, N_Vector resval, void *user_data);


int main(void)
{
//Declare any necessary local constants
N_Vector cc,cp;
long int mu,ml,iter;
void *mem;
realtype rtol, atol, t0, tout,tret,incr;
int iout,retval;

//Initialize all the relevant N-vectors

  cc  = N_VNew_Serial(N);
  if(check_flag((void *)cc, "N_VNew_Serial", 0)) return(1);

  cp  = N_VNew_Serial(N);
  if(check_flag((void *)cp, "N_VNew_Serial", 0)) return(1);
  

  SetInitialProfiles(cc,cp);  
  
  
  t0 = ZERO;
  rtol = RTOL; 
  atol = ATOL;
  iter= 10;
  incr= RCONST(0.01);
  tout=RCONST(0.01);

/* Call IDACreate and IDAMalloc to initialize IDA. */
  
  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);

  retval = IDAInit(mem, res_basic, t0, cc, cp);
  if(check_flag(&retval, "IDAInit", 1)) return(1);

  retval = IDASStolerances(mem, rtol, atol);
  if(check_flag(&retval, "IDASStolerances", 1)) return(1);

  /* Call IDABand to specify the IDA linear solver. */

  mu = 1;
  ml = 1;
  retval = IDABand(mem, N, mu, ml);
  if(check_flag(&retval, "IDABand", 1)) return(1);
  
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
  
  return(0);
}

static void SetInitialProfiles(N_Vector cc,N_Vector cp)
{
    realtype *ccv,*cpv;
    int i;
    ccv = NV_DATA_S(cc);
    cpv = NV_DATA_S(cp);
    for (i=0; i<N; i++){
         ccv[i]=0;
	 cpv[i]=0;
    }


}


static int res_basic(realtype tt, N_Vector cc, N_Vector cp, N_Vector resval, void *user_data)
{   
    realtype *ccv,*cpv,*res;
    int i;
    ccv = NV_DATA_S(cc);
    cpv = NV_DATA_S(cp);
    res = NV_DATA_S(resval);
    res[0]=ccv[0];
    res[N-1]=ccv[N-1]-1;
    for (i=1; i<N-1;i++)
     {
        res[i]=cpv[i]-alpha*(ccv[i+1]+ccv[i-1]-2*ccv[i])/(delta*delta);
        
      }

}



static void PrintOutput(void *mem, N_Vector c, realtype t)
{     
      realtype *conc;
      int i;
      conc = NV_DATA_S(c);
      printf("Output at time %8.2Le",t);
      printf("\n");
      for (i=0;i<N;i++){
          
          printf("%12.4le",conc[i]);
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
