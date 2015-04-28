/*
 * -----------------------------------------------------------------
 * $Date: 2015-04-26$
 * -----------------------------------------------------------------
 * Programmer: Pranav Shetty
 *              
 * -----------------------------------------------------------------
 * Example problem towards building a full fledged LiB solver: A half cell solver.
 *
 * This example solves a set of 5 PDE's to obtain the relevant parameters of a half cell.
 * This version uses the band solver IDABand, and IDACalcIC.
 * y & yp denote the set of variables and its derivaive respectively. A generic label y has been chosen in order to remain
 * consistent with ida's documentation and also because the  encompasses 5 different variables
 *
 *
 * The system is solved with IDA using the banded linear system
 * solver, half-bandwidths equal to 1, and default
 * difference-quotient Jacobian.
 * IDACalcIC is called to compute correct values at the boundary,
 * given incorrect values as input initial guesses. The constraints
 * u >= 0 are posed for some components.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>
#include <ida/ida_band.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

/* Problem Constants */

#define GRID        100
#define N_VAR       5    
#define N           GRID*N_VAR 
#define ZERO        RCONST(0.0)
#define ONE         RCONST(1.0)
#define TWO         RCONST(2.0)
#define T_PLUS      RCONST(0.363)
#define R           RCONST(8.314)
#define T           RCONST(298.15)
#define F           RCONST(96487.0)
#define BRUGG       RCONST(4.0)  //confirm if this is always constant
#define I           RCONST(2.0)

//Problem constants defined in preprocessor in this file. Need to keep all problem constants in a unique location in the next version
//Might have to refactor to relabel anode & cathode

//Look at struct brackets again
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
  realtype diff_coeff_solid;
} *Material_Data;

typedef struct {
  realtype dx_a;
  realtype coeff_a;
  realtype eps_a;
  realtype sigma_a;
  realtype diff_coeff_a;
  realtype radius_a;
  realtype k_a;
  realtype c_s_max_a;
  realtype c_s_0_a;
  realtype c_0_a;
  realtype l_a;
  realtype interfac_area_a;
  realtype sigma_eff_a;
  realtype diff_coeff_eff_a;
  realtype diff_coeff_solid_a;
  
  realtype dx_s;
  realtype coeff_s;
  realtype eps_s;
  realtype diff_coeff_s;
  realtype c_0_s;
  realtype l_s;
  realtype diff_coeff_eff_s;
  
  realtype dx_c;
  realtype coeff_c;
  realtype eps_c;
  realtype sigma_c;
  realtype diff_coeff_c;
  realtype diff_coeff_solid_c;
  realtype radius_c;
  realtype k_c;
  realtype c_s_max_c;
  realtype c_s_0_c;
  realtype c_0_c;
  realtype l_c;
  realtype interfac_area_c;
  realtype sigma_eff_c;
  realtype diff_coeff_eff_c;
  
  int sep_indicator, cath_indicator;
} *Cell_Data;

/* Prototypes of functions called by IDA */

int half_cell_residuals(realtype tres, N_Vector y, N_Vector yp, N_Vector resval, void *user_data);
static void InitAnodeData(Material_Data data_anode);
static void InitSepData(Material_Data data_sep);
static void InitCathodeData(Material_Data data_cathode);
static void InitCellData(Material_Data data_anode,Material_Data data_sep,Material_Data data_cathode,Cell_Data data);
realtype ocp_anode(realtype c,realtype c_max);
realtype ocp_cathode(realtype c,realtype c_max);
realtype SUNRpowerI(realtype a, int b);
realtype SUNRpowerR(realtype a, realtype b);
realtype Rsinh(realtype x);
realtype Rlog(realtype x); 

/* Prototypes of private functions */

//static void PrintHeader(realtype rtol, realtype atol);
//static void PrintOutput(void *mem, realtype t, N_Vector u);
static void SetInitialProfile(Material_Data data_anode,Material_Data data_sep,Material_Data data_cathode, N_Vector y,
                              N_Vector yp, N_Vector id, N_Vector res, N_Vector constraints);

static int check_flag(void *flagvalue, char *funcname, int opt);
realtype kappa(realtype c, realtype eps);
realtype ocp_anode(realtype c_surf, realtype c_max);
realtype ocp_cathode(realtype c_surf, realtype c_max);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(void)
{
  void *mem;
  Material_Data data_anode,data_cathode,data_sep;
  Cell_Data data;
  N_Vector y, yp, constraints, id, res;
  int ier, iout;
  long int mu, ml, netf, ncfn;
  realtype rtol, atol, t0, t1, tout, tret;
  mem = NULL;
  data_anode = data_cathode = data_sep = NULL;
  data = NULL;
  y = yp = constraints = id = res = NULL;
  
  data_anode = (Material_Data) malloc(sizeof *data_anode);
  if(check_flag((void *)data_anode, "malloc", 2)){
    return(1);
  }
  
  data_sep = (Material_Data) malloc(sizeof *data_sep);
  if(check_flag((void *)data_sep, "malloc", 2)){
    return(1);
  }
  
  data_cathode = (Material_Data) malloc(sizeof *data_cathode);
  if(check_flag((void *)data_cathode, "malloc", 2)){
    return(1);
  }
  
  data = (Cell_Data) malloc(sizeof *data);
  if(check_flag((void *)data, "malloc", 2)){
    return(1);
  }
    
  InitAnodeData(data_anode);
  InitSepData(data_sep);
  InitCathodeData(data_cathode);
  InitCellData(data_anode, data_sep,data_cathode,data);

  
  /* Create vectors y, yp, res, constraints, id. */
  y = N_VNew_Serial(N);
  if(check_flag((void *)y, "N_VNew_Serial", 0)){
    return(1);
  }
  
  yp = N_VNew_Serial(N);
  if(check_flag((void *)yp, "N_VNew_Serial", 0)){ 
    return(1);
  }
  
  res = N_VNew_Serial(N);
  if(check_flag((void *)res, "N_VNew_Serial", 0)){
    return(1);
  }
  
  constraints = N_VNew_Serial(N);
  if(check_flag((void *)constraints, "N_VNew_Serial", 0)){
    return(1);
  }
  
  id = N_VNew_Serial(N);
  if(check_flag((void *)id, "N_VNew_Serial", 0)){
    return(1);
  }
  
  
/* Create and load problem data block. */ 
    
  //Need to initialize remaining constants by hardcoding them here. Need to take value from cffi in later iteration  
  //Do it in a separate function
  //Need to deal with units of all constants

  /* Initialize y, yp, id. */
  SetInitialProfile(data_anode,data_sep,data_cathode, y, yp, id, res, constraints);

  /* Set remaining input parameters. */
  t0   = ZERO;
  t1   = RCONST(0.01);
  rtol = RCONST(1.0e-3);
  atol = RCONST(1.0e-3);

  /* Call IDACreate and IDAMalloc to initialize solution */
  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)){
    return(1);
  }
  
  ier = IDASetUserData(mem, data);
  if(check_flag(&ier, "IDASetUserData", 1)){
    return(1);
  }

  ier = IDASetId(mem, id);
  if(check_flag(&ier, "IDASetId", 1)){ 
    return(1);
  }

  ier = IDASetConstraints(mem, constraints);
  if(check_flag(&ier, "IDASetConstraints", 1)){
    return(1);
  }

  ier = IDAInit(mem, half_cell_residuals, t0, y, yp);
  if(check_flag(&ier, "IDAInit", 1)){
    return(1);
  }

  ier = IDASStolerances(mem, rtol, atol);
  if(check_flag(&ier, "IDASStolerances", 1)){
    return(1);
  }
  
  /* Call IDABand to specify the linear solver. */
  mu = 1; ml = 1; 
  ier = IDABand(mem, N, mu, ml);
  if(check_flag(&ier, "IDABand", 1)) return(1);
 
  /* Call IDACalcIC to correct the initial values. */
  
  ier = IDACalcIC(mem, IDA_YA_YDP_INIT, t1);
  if(check_flag(&ier, "IDACalcIC", 1)){ 
    return(1);
  }
  /* Print output heading. */
  //PrintHeader(rtol, atol);  
  
  //PrintOutput(mem, t0, uu); //Change


  /* Loop over output times, call IDASolve, and print results. */
  
  for (tout = t1, iout = 1; iout <= 10; iout++, tout *= RCONST(2.0)) {
    
    ier = IDASolve(mem, tout, &tret, y, yp, IDA_NORMAL);
    if(check_flag(&ier, "IDASolve", 1)){
      return(1);
    }
    //PrintOutput(mem, tret, y);
  
  }
  
  /* Print remaining counters and free memory. 
  ier = IDAGetNumErrTestFails(mem, &netf);
  check_flag(&ier, "IDAGetNumErrTestFails", 1);
  ier = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_flag(&ier, "IDAGetNumNonlinSolvConvFails", 1);
  printf("\n netf = %ld,   ncfn = %ld \n", netf, ncfn);
  */
  
  IDAFree(&mem);
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(id);
  N_VDestroy_Serial(res);
  N_VDestroy_Serial(constraints);
  free(data_anode);
  free(data_sep);
  free(data_cathode);
  free(data);
  
  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * half_cell_residuals: Lithium ion battery system residual function                       
 * 
 */

int half_cell_residuals(realtype tres, N_Vector y, N_Vector yp, N_Vector resval, 
                        void *user_data)
{
  long int j, i;
  int sep_indicator,cath_indicator;
  realtype *yv, *ypv, *resv;
  Cell_Data data;
  
  yv = NV_DATA_S(y); 
  ypv = NV_DATA_S(yp); 
  resv = NV_DATA_S(resval);
  data = (Cell_Data)user_data;
  
  sep_indicator = data->sep_indicator;
  cath_indicator = data->cath_indicator;

  for(i=0; i<N;i++){
      j=i/GRID;  
      switch(j){
          case 0: //c
                 {
                  if(i%GRID==0){
                      resv[i] = yv[i+1]-yv[i];
                  }
                  
                  else if(i%GRID>0 && i%GRID<sep_indicator){
                      resv[i] = (data->eps_a)*ypv[i]-(data->diff_coeff_eff_a)*SUNRpowerI(data->coeff_a,2)*
                                (yv[i+1]+yv[i-1]-TWO*yv[i])-(data->interfac_area_a)*(ONE-T_PLUS)*yv[i+3*GRID];//check here
                      
                  }
                  
                  else if(i%GRID==sep_indicator){
                      resv[i] = (data->diff_coeff_eff_s)*yv[i+1]-(data->diff_coeff_eff_a)*yv[i];
                  }
                  
                  else if(i%GRID>sep_indicator && i%GRID<cath_indicator){
                      resv[i] = (data->eps_s)*ypv[i]-(data->diff_coeff_eff_s)*SUNRpowerI(data->coeff_s,2)*
                                (yv[i+1]+yv[i-1]-TWO*yv[i]);
                      
                  }
                  
                  else if(i%GRID==cath_indicator){
                      resv[i] = (data->diff_coeff_eff_c)*yv[i+1]-(data->diff_coeff_eff_s)*yv[i];
                  }
                  
                  else if(i%GRID>cath_indicator && i%GRID<GRID-1){
                      resv[i] = (data->eps_c)*ypv[i]-(data->diff_coeff_eff_c)*SUNRpowerI(data->coeff_c,2)*
                                (yv[i+1]+yv[i-1]-TWO*yv[i])-(data->interfac_area_c)*(ONE-T_PLUS)*yv[i+3*GRID];;
                     
                  }
                  
                  else if(i%GRID==GRID-1){
                      resv[i] = yv[i]-yv[i-1];
                  }
                  
                  break;
                 }
          
          case 1: //phi1
                                   
                {
                  if(i%GRID==0){
                      resv[i] = (data->coeff_a)*(yv[i+1]-yv[i])+I/(data->sigma_eff_a);
                  }
                  //Might need additional BC
                  else if(i%GRID>0 && i%GRID<=sep_indicator){
                      resv[i] = (data->sigma_eff_a)*SUNRpowerI(data->coeff_a,2)*(yv[i+1]+yv[i-1] - TWO*y[i]) -
                                (data->interfac_area_a)*F*yv[i+2*GRID];//Check
                      
                  }
                  /*
                  else if(i%GRID==sep_indicator){
                      resv[i] = 
                  }
                  */
                  else if(i%GRID>sep_indicator && i%GRID<cath_indicator){
                      resv[i] = yv[i] - ZERO;
                      
                  }
                  
                  else if(i%GRID==cath_indicator){
                      resv[i] = yv[i]-yv[i-1];
                  }
                  
                  else if(i%GRID>cath_indicator && i%GRID<GRID-1){
                      resv[i] = (data->sigma_eff_c)*SUNRpowerI(data->coeff_c,2)*(yv[i+1]+yv[i-1]-TWO*y[i])-
                                (data->interfac_area_c)*F*yv[i+2*GRID];//Check
                     
                  }
                  
                  else if(i%GRID==GRID-1){
                      resv[i] = (data->coeff_c)*(yv[i]-yv[i-1])+I/(data->sigma_eff_c);
                  }
                  
                  break;
                }
                  
          
          case 2://phi2
                 {
                  if(i%GRID==0){
                      resv[i] = kappa(yv[i+1-2*GRID],data->eps_a)*yv[i+1]-kappa(yv[i-2*GRID],data->eps_a)*yv[i];
                  }
                  
                  else if(i%GRID>0 && i%GRID<sep_indicator){
                      resv[i] = I+(data->sigma_eff_a)*(data->coeff_a)*(yv[i+1-GRID]-yv[i-GRID])+
                                kappa(yv[i-2*GRID],data->eps_a)*(data->coeff_a)*(yv[i+1]-yv[i])-
                                TWO*kappa(yv[i-2*GRID],data->eps_a)*R*T*(ONE/F)*(data->coeff_a)*(ONE-T_PLUS)*
                                (Rlog(yv[i+1-2*GRID])-Rlog(yv[i-2*GRID]));
                                
                                                                                      
                  }
                  
                  else if(i%GRID==sep_indicator){
                      resv[i] = kappa(yv[i+1-2*GRID],data->eps_s)*(yv[i+1]-yv[i])-kappa(yv[i-2*GRID],data->eps_a)*
                                (yv[i]-yv[i-1]);
                  }
                  
                  else if(i%GRID>sep_indicator && i%GRID<cath_indicator){
                      resv[i] = I+kappa(yv[i-2*GRID],data->eps_a)*(data->coeff_a)*(yv[i+1]-yv[i])-
                                TWO*kappa(yv[i-2*GRID],data->eps_a)*R*T*(ONE/F)*(data->coeff_a)*(ONE-T_PLUS)*
                                (Rlog(yv[i+1-2*GRID])-Rlog(yv[i-2*GRID]));
                      
                  }
                  
                  else if(i%GRID==cath_indicator){
                      resv[i] = kappa(yv[i+1-2*GRID],data->eps_c)*(yv[i+1]-yv[i])-kappa(yv[i-2*GRID],data->eps_s)*
                                (yv[i]-yv[i-1]);
                  }
                  
                  else if(i%GRID>cath_indicator && i%GRID<GRID-1){
                      resv[i] = I+(data->sigma_eff_c)*(data->coeff_c)*(yv[i+1-GRID]-yv[i-GRID])+
                                kappa(yv[i-2*GRID],data->eps_c)*(data->coeff_c)*(yv[i+1]-yv[i])-
                                TWO*kappa(yv[i-2*GRID],data->eps_c)*R*T*(ONE/F)*(data->coeff_c)*(ONE-T_PLUS)*
                                (Rlog(yv[i+1-2*GRID])-Rlog(yv[i-2*GRID]));
                     
                  }
                  
                  else if(i%GRID==GRID-1){
                      resv[i] = yv[i]-yv[i-1];
                  }
                  
                  break;
                }
          
          case 3://j
                 {
                  if(i%GRID>=0 && i%GRID<=sep_indicator){
                      resv[i] = yv[i] - TWO*(data->k_a)*SUNRpowerR(yv[i-3*GRID],RCONST(0.5))*
                                Rsinh(RCONST(0.5)*(F/(R*T))*(yv[i-2*GRID]-yv[i-GRID]-
                                ocp_anode(yv[i+GRID],data->c_s_max_a)))*          
                                (data->c_s_max_a -(yv[i+GRID]-yv[i]*(data->radius_a)/(RCONST(5.0)*(data
                                ->diff_coeff_solid_a))));
                      
                  }
                  else if(i%GRID>sep_indicator && i%GRID<cath_indicator){
                      resv[i]=yv[i] - ZERO;
                      
                  }
                  else if(i%GRID>=cath_indicator && i%GRID<=GRID-1){
                      resv[i] = yv[i] - TWO*(data->k_s)*SUNRpowerR(yv[i-3*GRID],RCONST(0.5))*
                                Rsinh(RCONST(0.5)*(F/(R*T))*(yv[i-2*GRID]-yv[i-GRID]-
                                ocp_cathode(yv[i+GRID],data->c_s_max_c)))*          
                                (data->c_s_max_c -(yv[i+GRID]-yv[i]*(data->radius_c)/(RCONST(5.0)*(data
                                ->diff_coeff_solid_c))));;
                     
                  }
                  
                  break;
                 }  
           
          case 4:
                 {
                  if(i%GRID>=0 && i%GRID<=sep_indicator){
                      resv[i] = ypv[i]+3*yv[i-GRID]/(data->radius_a);
                      
                  }
                  else if(i%GRID>sep_indicator && i%GRID<=cath_indicator){
                      resv[i]=yv[i] - ZERO;
                      
                  }
                  else if(i%GRID>cath_indicator && i%GRID<=GRID-1){
                      resv[i] = ypv[i]+3*yv[i-GRID]/(data->radius_c);
                     
                  }
 
                  break;
                 }
  
           }
  
  }
  
  return(0);

}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * SetInitialProfile: routine to initialize y, yp,constraints and id vectors       
 */
/*y[0-(GRID-1)]       concentration
 *y[GRID-2*GRID-1]    ph1
 *y[2*GRID-3*GRID-1]  phi2
 *y[3*GRID- 4*GRID-1] j
 *y[4*GRID-5*GRID-1]  c_s
 */
static void SetInitialProfile(Material_Data data_anode,Material_Data data_sep,Material_Data data_cathode,
                             N_Vector y, N_Vector yp, N_Vector res,
                             N_Vector id, N_Vector constraints)
{
  realtype *ydata, *ypdata, *iddata,*constraintdata, c_0,l_a,l_c,phi1_a,phi1_c;
  long int i,j;
  int sep_indicator,cath_indicator;
  ydata = NV_DATA_S(y);
  ypdata = NV_DATA_S(yp);
  iddata = NV_DATA_S(id);
  constraintdata = NV_DATA_S(constraints);
  c_0 = data_anode->c_0;
  
  l_a = (data_anode->l)/(data_anode->l + data_sep->l + data_cathode->l);
  l_s = (data_anode->l+data_sep->l)/(data_anode->l + data_sep->l + data_cathode->l);
  
  sep_indicator = (int)(l_a*RCONST(GRID));
  cath_indicator = (int)(l_s*GRID);

  
  phi1_a = ocp_anode(data_anode->c_s_0,data_anode->c_s_max);
  phi1_c = ocp_anode(data_cathode->c_s_0,data_cathode->c_s_max);
  
  /* Initialize y and yp on all grid points. */ 
  //Future refactoring: Can move all assignments to one block and pass arguments to it
  //j, c_s and phi1 one not used in separator. Initialized to zero for now. May cause issues later.
  //Need to force to zero value in the equation as well
  for(i=0; i<N;i++){
      j=i/GRID;  
      switch(j){
          case 0:
                 {
                  ydata[i]=c_0;
                  ypdata[i]=ZERO;
                  iddata[i]=ONE;   
                  constraintdata[i]=ONE;
                  break;
                 }
          
          case 1:
                 {
                  
                  if(i%GRID>=0 && i%GRID<=sep_indicator){
                      ydata[i]=phi1_a;
                      ypdata[i]=ZERO;
                      iddata[i]=ZERO;   
                      constraintdata[i]=ZERO;
                  }
                  else if(i%GRID>sep_indicator && i%GRID<=cath_indicator){
                      ydata[i]=ZERO;
                      ypdata[i]=ZERO;
                      iddata[i]=ZERO;   
                      constraintdata[i]=ZERO;
                  }
                  else if(i%GRID>cath_indicator && i%GRID<=GRID-1){
                      ydata[i]=phi1_c;
                      ypdata[i]=ZERO;
                      iddata[i]=ZERO;
                      constraintdata[i]=ZERO;
                  }
                  

                  break;
                 }
          
          case 2:
                 {
                  ydata[i]=ZERO;
                  ypdata[i]=ZERO;
                  iddata[i]=ZERO;   
                  constraintdata[i]=ZERO;
                  break;
                 }
          
          case 3:
                 {
                  ydata[i]=ZERO;
                  ypdata[i]=ZERO;
                  iddata[i]=ZERO;   
                  constraintdata[i]=ZERO;
                  break;
                 }  
           
          case 4:
                 {
                  if(i%GRID>=0 && i%GRID<=sep_indicator){
                      ydata[i]=data_anode->c_s_0;
                      ypdata[i]=ZERO;
                      iddata[i]=ONE;   
                      constraintdata[i]=ONE;
                  }
                  else if(i%GRID>sep_indicator && i%GRID<=cath_indicator){
                      ydata[i]=ZERO;
                      ypdata[i]=ZERO;
                      iddata[i]=ONE;   
                      constraintdata[i]=ONE;
                  }
                  else if(i%GRID>cath_indicator && i%GRID<=GRID-1){
                      ydata[i]=data_cathode->c_s_0;
                      ypdata[i]=ZERO;
                      iddata[i]=ONE;   
                      constraintdata[i]=ONE;
                  }
 
                  break;
                 }
  
           }
  
  }
}    

static void InitAnodeData(Material_Data data_anode)
{  
  data_anode->dx = ONE/(RCONST(GRID) - ONE);
  data_anode->coeff = ONE/(data_anode->dx);
  data_anode->sigma = RCONST(100.0);
  data_anode->eps = RCONST(0.385);
  data_anode->sigma_eff = (data_anode->sigma)*(1-data_anode->eps);
  data_anode->diff_coeff = RCONST(1.0e-14);
  data_anode->diff_coeff_eff = (data_anode->diff_coeff)*SUNRpowerR((data_anode->eps),BRUGG);
  data_anode->diff_coeff_solid = RCONST(1.0e-14);
  data_anode->k = RCONST(2.334e-11);
  data_anode->c_0 = RCONST(1000.0);
  data_anode->c_s_max = RCONST(51554.0);
  data_anode->c_s_0 = (data_anode->c_s_max)*RCONST(0.4955);
  data_anode->l = RCONST(80.0e-6);
  data_anode->radius = RCONST(2.0e-6);
  data_anode->interfac_area = RCONST(3.0)*(1-data_anode->eps)/(data_anode->radius);
}

static void InitSepData(Material_Data data_sep)
{  
  data_sep->dx = ONE/(GRID - ONE);
  data_sep->coeff = ONE/(data_sep->dx);
  data_sep->eps = RCONST(0.724);
  data_sep->sigma_eff = (data_sep->sigma)*(ONE-data_sep->eps);
  data_sep->diff_coeff = RCONST(7.5e-10);
  data_sep->diff_coeff_eff = (data_sep->diff_coeff)*SUNRpowerR((data_sep->eps),BRUGG);
  
  data_sep->c_0 = RCONST(1000.0);
  data_sep->l = RCONST(25.0e-6);
}


static void InitCathodeData(Material_Data data_cathode)
{  
  data_cathode->dx = ONE/(RCONST(GRID) - ONE);
  data_cathode->coeff = ONE/(data_cathode->dx);
  data_cathode->sigma = RCONST(100.0);
  data_cathode->eps = RCONST(0.0326);
  data_cathode->sigma_eff = (data_cathode->sigma)*(ONE-data_cathode->eps);
  data_cathode->diff_coeff = RCONST(3.9e-14);
  data_cathode->diff_coeff_eff = (data_cathode->diff_coeff)*SUNRpowerR((data_cathode->eps),BRUGG);
  data_cathode->diff_coeff_solid = RCONST(3.9e-14);
  data_cathode->k = RCONST(5.0307e-11);
  data_cathode->c_0 = RCONST(1000.0);
  data_cathode->c_s_max = RCONST(30555.0);
  data_cathode->c_s_0 = (data_cathode->c_s_max)*RCONST(0.8551);
  data_cathode->l = RCONST(88.0e-6);
  data_cathode->radius = RCONST(2.0e-6);
  data_cathode->interfac_area = RCONST(3.0)*(ONE-data_cathode->eps)/(data_cathode->radius);
}

static void InitCellData(Material_Data data_anode,Material_Data data_sep,Material_Data data_cathode,Cell_Data data)
{
  realtype l_a, l_s;
  
  data->dx_a = data_anode->dx;
  data->coeff = data_anode->coeff;
  data->eps_a = data_anode->eps;
  data->sigma_a = data_anode->sigma;
  data->sigma_eff_a = data_anode->sigma_eff;
  data->diff_coeff_a = data_anode->diff_coeff;
  data->diff_coeff_eff_a = data_anode->diff_coeff_eff;
  data->diff_coeff_solid_a = data_anode->diff_coeff_solid;
  data->k_a = data_anode->k;
  data->c_0_a = data_anode->c_0;
  data->c_s_max_a = data_anode->c_s_max;
  data->c_s_0_a = data_anode->c_s_0;
  data->l_a = data_anode->l;
  data->radius_a = data_anode->radius;
  data->interfac_area_a = data_anode->interfac_area;
  
  data->dx_s = data_sep->dx;
  data->coeff_s = data_sep->coeff;
  data->eps_s = data_sep->eps;
  data->diff_coeff_s = data_sep->diff_coeff;
  data->diff_coeff_eff_s = data_sep->diff_coeff_eff;
  data->c_0_s = data_sep->c_0;
  data->l_s = data_sep->l;
  
  data->dx_c = data_cathode->dx;
  data->coeff_c = data_cathode->coeff;
  data->eps_c = data_cathode->eps;
  data->sigma_c = data_cathode->sigma;
  data->sigma_eff_c = data_cathode->sigma_eff;
  data->diff_coeff_c = data_cathode->diff_coeff;
  data->diff_coeff_eff_c = data_cathode->diff_coeff_eff;
  data->diff_coeff_solid_c = data_cathode->diff_coeff_solid;
  data->k_c = data_cathode->k;
  data->c_0_c = data_cathode->c_0;
  data->c_s_max_c = data_cathode->c_s_max;
  data->c_s_0_c = data_cathode->c_s_0;
  data->l_c = data_cathode->l;
  data->radius_c = data_cathode->radius;
  data->interfac_area_c = data_cathode->interfac_area;
  
  l_a = (data_anode->l)/(data_anode->l + data_sep->l + data_cathode->l);
  l_s = (data_anode->l+data_sep->l)/(data_anode->l + data_sep->l + data_cathode->l);
  
  data->sep_indicator = (int)(l_a*RCONST(GRID));
  data->cath_indicator = (int)(l_s*GRID);
}

realtype kappa(realtype c, realtype eps)
{
    realtype exp_term, polynom_term;
    exp_term = SUNRpowerR(eps,BRUGG);
    polynom_term = RCONST(4.1253e-2)+RCONST(5.007e-4)*c-RCONST(4.7212e-7)*SUNRpowerI(c,2)+
                   RCONST(1.5094)*SUNRpowerI(c,3)-RCONST(1.6018)*SUNRpowerI(c,4);
    
    return polynom_term*exp_term;
  
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


/* May not be required
realtype FLOOR(realtype x)
{
  #if defined(SUNDIALS_USE_GENERIC_MATH)
  return(floor((double) x));
  #elif defined(SUNDIALS_DOUBLE_PRECISION)
  return((floor(x));
  #elif defined(SUNDIALS_SINGLE_PRECISION)
  return(floor(x));
  #elif defined(SUNDIALS_EXTENDED_PRECISION)
  return(floor(x));
  #endif

*/
/* 
 * Print first lines of output (problem description)
 

static void PrintHeader(realtype rtol, realtype atol)
{
  printf("\nidaHeat1D_bnd: Heat equation, serial example problem for IDA\n");
  printf("          Discretized heat equation on 2D unit square.\n");
  printf("          Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("          Mesh dimensions: %d", MGRID);
  printf("        Total system size: %d\n\n", NEQ);
#if defined(SUNDIALS_EXTENDED_PRECISION) 
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION) 
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("Constraints set to force all solution components >= 0. \n");
  printf("Linear solver: IDABAND, banded direct solver \n");
  printf("       difference quotient Jacobian, half-bandwidths = %d \n",MGRID);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("IDACalcIC called with input boundary values = %Lg \n",BVAL);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("IDACalcIC called with input boundary values = %g \n",BVAL);
#else
  printf("IDACalcIC called with input boundary values = %g \n",BVAL);
#endif
  // Print output table heading and initial line of table. 
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time       umax     k  nst  nni  nje   nre   nreLS    h      \n" );
  printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . \n");
}
*/
/*
 * Print Output
 */
/*
static void PrintOutput(void *mem, realtype t, N_Vector uu)
{
  int ier;
  realtype umax, hused;
  long int nst, nni, nje, nre, nreLS;
  int kused;

  umax = N_VMaxNorm(uu);
  
  ier = IDAGetLastOrder(mem, &kused);
  check_flag(&ier, "IDAGetLastOrder", 1);
  ier = IDAGetNumSteps(mem, &nst);
  check_flag(&ier, "IDAGetNumSteps", 1);
  ier = IDAGetNumNonlinSolvIters(mem, &nni);
  check_flag(&ier, "IDAGetNumNonlinSolvIters", 1);
  ier = IDAGetNumResEvals(mem, &nre);
  check_flag(&ier, "IDAGetNumResEvals", 1);
  ier = IDAGetLastStep(mem, &hused);
  check_flag(&ier, "IDAGetLastStep", 1);
  ier = IDADlsGetNumJacEvals(mem, &nje);
  check_flag(&ier, "IDADlsGetNumJacEvals", 1);
  ier = IDADlsGetNumResEvals(mem, &nreLS);
  check_flag(&ier, "IDADlsGetNumResEvals", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION) 
  printf(" %5.2Lf %13.5Le  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2Le \n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION) 
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e \n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused);
#else
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e \n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused);
#endif

}

 */
/* 
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
              funcname, *errflag);
      return(1); 
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  }

  return(0);
}
