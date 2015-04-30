///\file
 /** -----------------------------------------------------------------
 * $Date: 2015-04-26$
 * -----------------------------------------------------------------
 * Programmer: Pranav Shetty
 *              
 * -----------------------------------------------------------------
 * Building a full fledged LiB solver
 *
 * This code solves a set of 5 PDE's to obtain the relevant parameters of a lithium ion cell using a 1D model.
 * This version uses the Krylov solver IDASpgmr.
 * y & yp denote the set of variables and its derivaive respectively. A generic label y has been chosen in order to remain
 * consistent with ida's documentation and also because it encompasses 5 different variables over the entire domain
 *
 * The system is solved with IDA using the banded linear system
 * solver, half-bandwidths equal to 1, and default
 * difference-quotient Jacobian.
 * IDACalcIC is called to compute correct values at the boundary, given incorrect values as input initial guesses 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>
#include <ida/ida_spgmr.h>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include "functions.c"

/* Problem Constants */

#define GRID        50   ///<Number of grid points the domain is split into
#define N_VAR       5    ///<Number of variables in domain
#define N           GRID*N_VAR 
#define ZERO        RCONST(0.0)
#define ONE         RCONST(1.0)
#define TWO         RCONST(2.0)
#define T_PLUS      RCONST(0.363)///<Transport number in anode
#define R           RCONST(8.314)///<Universal Gas constant
#define T           RCONST(298.15)///<Temperature
#define F           RCONST(96487.0)///<Faraday's constant
#define BRUGG       RCONST(4.0) ///<Bruggman's coefficient. Used to account for effect of porosity in constants
#define I           RCONST(2.0)
#define EPS         RCONST(0.01)

///Struct defined in order to store material properties in each region of the cell(anode, cathode and separator)


typedef struct {
  realtype dx; ///< Defines the spacing between grid points
  realtype coeff; ///<Takes the value 1/dx
  realtype eps; ///< Porosity of material. Dimensionless
  realtype sigma;///< Conductivity. Unit: S/m
  realtype diff_coeff;///< Diffusion coefficient. Unit: m^2/s
  realtype radius;///< Average radius of electrode particles: Unit:m
  realtype k;///<Reaction constant for the Butler Volmer equation
  realtype c_s_max;///<Maximum lithium intercalation capacity of the electrode. Unit: mol/m^3
  realtype c_s_0;///<Initial concentration of lithium in the solid phase. Unit: mol/m^3
  realtype c_0;///<Initial concentration of lithium in the electrolyte. Unit: mol/m^3
  realtype l;///<Dimension of the electrode region. Unit: m
  realtype interfac_area;///< Interfacial area of electrode particles. Unit: m^2/m^3
  realtype sigma_eff; ///<Effective conductivity taking into account porosity of medium. Unit: S/m
  realtype diff_coeff_eff; ///<Effective diffusion coefficient taking into account porosity of medium. Unit: m^2/s
  realtype diff_coeff_solid;///<Diffusion coefficient of lithium intercalation in the electrode
} *Material_Data;

///This struct transfers the material data from each of the threee electrode Material data blocks into a single data block which can then be passed to IDA
/**Each variable used in Cell_Data is an analogue of a variable used in Material_Data with a suffix to indicate which region 
   it came from
   
*/
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
  
  int sep_indicator;///< Gives an integer value on the grid between anode and separator
  int cath_indicator;///< Gives an integer value on the grid between cathode and separator
} *Cell_Data;

static void InitAnodeData(Material_Data data_anode);
static void InitSepData(Material_Data data_sep);
static void InitCathodeData(Material_Data data_cathode);
static void InitCellData(Material_Data data_anode,Material_Data data_sep,Material_Data data_cathode,Cell_Data data);
realtype kappa(realtype c, realtype eps);
realtype ocp_anode(realtype c_surf, realtype c_max);
realtype ocp_cathode(realtype c_surf, realtype c_max);

realtype Rsinh(realtype x);
realtype Rlog(realtype x); 

/* Prototypes of private functions */

void SetInitialProfile(Material_Data data_anode,Material_Data data_sep,Material_Data data_cathode, N_Vector y,
                              N_Vector yp, N_Vector id, N_Vector constraints);
                              int half_cell_residuals(realtype tres, N_Vector y, N_Vector yp, N_Vector resval, 
                              void*user_data);
int half_cell_residuals(realtype tres, N_Vector y, N_Vector yp, N_Vector resval, void *user_data);
static int check_flag(void *flagvalue, char *funcname, int opt);

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
  long int mu, ml,i;
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
    
  
  /* Initialize y, yp, id. */
  SetInitialProfile(data_anode,data_sep,data_cathode, y, yp, id, constraints);

  /* Set remaining input parameters. */
  t0   = ZERO;
  t1   = RCONST(0.001);
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
  /*
  ier = IDASetConstraints(mem, constraints);
  if(check_flag(&ier, "IDASetConstraints", 1)){
    return(1);
  }
  */
  ier = IDAInit(mem, half_cell_residuals, t0, y, yp);
  if(check_flag(&ier, "IDAInit", 1)){
    return(1);
  }

  ier = IDASStolerances(mem, rtol, atol);
  if(check_flag(&ier, "IDASStolerances", 1)){
    return(1);
  }
  
  /* Call IDABand to specify the linear solver. */
  //mu = 1; ml = 1; 
  
  ier = IDASpgmr(mem, 0);
  if(check_flag(&ier, "IDASpgmr", 1)) return(1);
  
  /*
  ier = IDABand(mem, N, mu,ml);
  if(check_flag(&ier, "IDABand", 1)) return(1);
  */
  /*
  ier = IDADense(mem, N);
  if(check_flag(&ier, "IDADense", 1)) return(1);
  */
  /* Call IDACalcIC to correct the initial values. */
  
  ier = IDACalcIC(mem, IDA_YA_YDP_INIT, t1);
  if(check_flag(&ier, "IDACalcIC", 1)){ 
    return(1);
  }

  /* Loop over output times, call IDASolve, and print results. */
  realtype *yv;
  yv = NV_DATA_S(y);
  
  for (tout = t1, iout = 1; iout <= 10; iout++, tout *= RCONST(2.0)) {
    
    ier = IDASolve(mem, tout, &tret, y, yp, IDA_NORMAL);
    if(check_flag(&ier, "IDASolve", 1)){
      return(1);
    }
    
  }
  
  for(i=0;i<4*GRID;i++){
      printf("%3.4f\n",yv[i]);
    }
  
 
  
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

/// Lithium ion battery system residual function. This defines the governing PDE equations for this system
 /** The vector y is split into 5 regions, each region representing one of the variables used in the equation \n
 y[0-(GRID-1)]   .......concentration \n
 y[GRID-2*GRID-1]....... ph1 \n
 y[2*GRID-3*GRID-1].....  phi2\n
 y[3*GRID- 4*GRID-1]...... j \n
 y[4*GRID-5*GRID-1]  .....c_s\n
 While looping over y we define equations in each of the 3 regions along with boundary conditions at the boundary of each
 of the three regions
 This function is called by IDA every time it has to solve at a given time t
 */

int half_cell_residuals(realtype tres, N_Vector y, N_Vector yp, N_Vector resval, 
                        void *user_data)
{
  long int j, i, N_t;
  int sep_indicator,cath_indicator;
  realtype *yv, *ypv, *resv;
  Cell_Data data;
  yv = NV_DATA_S(y); 
  ypv = NV_DATA_S(yp); 
  resv = NV_DATA_S(resval);
  data = (Cell_Data)user_data;
  N_t = N;
  sep_indicator = data->sep_indicator;
  cath_indicator = data->cath_indicator;
  
  
  for(i=0; i<2;i++){
      j=i/GRID;  
      switch(j){
          case 0: //c
                 {
                  if(i%GRID==0){
                      resv[i] = yv[i+1]-yv[i];
                      
                  }
                  
                  else if(i%GRID>0 && i%GRID<sep_indicator){
                      resv[i] = (data->eps_a)*ypv[i]-(data->diff_coeff_eff_a)*SUNRpowerI(data->coeff_a,2)*
                                (yv[i+1]+yv[i-1]-TWO*yv[i])-(data->interfac_area_a)*(ONE-T_PLUS)*yv[i+3*GRID];
                      
                  }
                  
                  else if(i%GRID==sep_indicator){
                      resv[i] = (data->diff_coeff_eff_s)*(yv[i+1]-yv[i])-(data->diff_coeff_eff_a)*(yv[i]-yv[i-1]);
                  }
                  
                  else if(i%GRID>sep_indicator && i%GRID<cath_indicator){
                      resv[i] = (data->eps_s)*ypv[i]-(data->diff_coeff_eff_s)*SUNRpowerI(data->coeff_s,2)*
                                (yv[i+1]+yv[i-1]-TWO*yv[i]);
                      
                  }
                  
                  else if(i%GRID==cath_indicator){
                      resv[i] = (data->diff_coeff_eff_c)*(yv[i+1]-yv[i])-(data->diff_coeff_eff_s)*(yv[i]-yv[i-1]);
                  }
                  
                  else if(i%GRID>cath_indicator && i%GRID<GRID-1){
                      resv[i] = (data->eps_c)*ypv[i]-(data->diff_coeff_eff_c)*SUNRpowerI(data->coeff_c,2)*
                                (yv[i+1]+yv[i-1]-TWO*yv[i])-(data->interfac_area_c)*(ONE-T_PLUS)*yv[i+3*GRID];
                     
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
                 
                  else if(i%GRID>0 && i%GRID<sep_indicator){
                      resv[i] = (data->sigma_eff_a)*SUNRpowerI(data->coeff_a,2)*(yv[i+1]+yv[i-1]-TWO*yv[i]) -
                                (data->interfac_area_a)*F*yv[i+2*GRID];
                      
                  }
                  
                  else if(i%GRID==sep_indicator){
                      resv[i] = yv[i]-yv[i-1];
                  }
                  
                  else if(i%GRID>sep_indicator && i%GRID<cath_indicator){
                      resv[i] = yv[i] - ZERO;
                      
                  }
                  
                  else if(i%GRID==cath_indicator){
                      resv[i] = yv[i]-yv[i-1]; 
                  }
                  
                  else if(i%GRID>cath_indicator && i%GRID<GRID-1){
                      resv[i] = (data->sigma_eff_c)*SUNRpowerI(data->coeff_c,2)*(yv[i+1]+yv[i-1]-TWO*yv[i])-
                                (data->interfac_area_c)*F*yv[i+2*GRID];
                     
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
                      resv[i] = RCONST(0.5)*I+(data->sigma_eff_a)*(data->coeff_a)*(yv[i+1-GRID]-yv[i-1-GRID])+
                                kappa(yv[i-2*GRID],data->eps_a)*(data->coeff_a)*(yv[i+1]-yv[i-1])-
                                TWO*kappa(yv[i-2*GRID],data->eps_a)*(R*T/F)*(data->coeff_a)*(ONE-T_PLUS)*
                                ((yv[i+1-2*GRID])-(yv[i-1-2*GRID]))/(yv[i-1-2*GRID]);
                                
                                                                                      
                  }
                  
                  else if(i%GRID==sep_indicator){
                      resv[i] = kappa(yv[i+1-2*GRID],data->eps_s)*(yv[i+1]-yv[i])-kappa(yv[i-2*GRID],data->eps_a)*
                                (yv[i]-yv[i-1]);
                  }
                  
                  else if(i%GRID>sep_indicator && i%GRID<cath_indicator){
                      resv[i] = I+kappa(yv[i-2*GRID],data->eps_a)*(data->coeff_a)*(yv[i+1]-yv[i])-
                                TWO*kappa(yv[i-2*GRID],data->eps_a)*(R*T/F)*(data->coeff_a)*(ONE-T_PLUS)*
                                ((yv[i+1-2*GRID])-(yv[i-2*GRID]))/(yv[i-2*GRID]);
                      
                  }
                  
                  else if(i%GRID==cath_indicator){
                      resv[i] = kappa(yv[i+1-2*GRID],data->eps_c)*(yv[i+1]-yv[i])-kappa(yv[i-2*GRID],data->eps_s)*
                                (yv[i]-yv[i-1]);
                  }
                  
                  else if(i%GRID>cath_indicator && i%GRID<GRID-1){
                      resv[i] = I+(data->sigma_eff_c)*(data->coeff_c)*(yv[i+1-GRID]-yv[i-GRID])+
                                kappa(yv[i-2*GRID],data->eps_c)*(data->coeff_c)*(yv[i+1]-yv[i])-
                                TWO*kappa(yv[i-2*GRID],data->eps_c)*(R*T/F)*(data->coeff_c)*(ONE-T_PLUS)*
                                ((yv[i+1-2*GRID])-(yv[i-2*GRID]))/(yv[i-2*GRID]);
                     
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
                      resv[i] = yv[i] - TWO*(data->k_c)*SUNRpowerR(yv[i-3*GRID],RCONST(0.5))*
                                Rsinh(RCONST(0.5)*(F/(R*T))*(yv[i-2*GRID]-yv[i-GRID]-
                                ocp_cathode(yv[i+GRID],data->c_s_max_c)))*          
                                (data->c_s_max_c -(yv[i+GRID]-yv[i]*(data->radius_c)/(RCONST(5.0)*(data
                                ->diff_coeff_solid_c))));
                     
                  }
                  
                  break;
                 }  
           
          case 4://c_s
                 {
                  if(i%GRID>=0 && i%GRID<=sep_indicator){
                      resv[i] = ypv[i]+RCONST(3.0)*yv[i-GRID]/(data->radius_a);
                      
                  }
                  else if(i%GRID>sep_indicator && i%GRID<=cath_indicator){
                      resv[i]=yv[i] - ZERO;
                      
                  }
                  else if(i%GRID>cath_indicator && i%GRID<=GRID-1){
                      resv[i] = ypv[i]+RCONST(3.0)*yv[i-GRID]/(data->radius_c);
                     
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


 ///Function to initialize y, yp,constraints and id vectors       
 
 /** The vector y is split into 5 regions, each region representing one of the variables used in the equation \n
 y[0-(GRID-1)]   .......concentration.Set to 1000 mol/m^3 \n 
 y[GRID-2*GRID-1]....... ph1. Set to corresponding values of ocp_anode & cathode \n
 y[2*GRID-3*GRID-1].....  phi2. Set to 0\n
 y[3*GRID- 4*GRID-1]...... j. Set to 0 \n
 y[4*GRID-5*GRID-1]  .....c_s. Set to initial values defined by user\n
 While looping over y we define equations in each of the 3 regions along with boundary conditions at the boundary of each
 of the three regions
 This function is called in int main once to set up the values of y and yp before IDA can operate on it
 */
 
void SetInitialProfile(Material_Data data_anode,Material_Data data_sep,Material_Data data_cathode,
                             N_Vector y, N_Vector yp,
                             N_Vector id, N_Vector constraints)
{
  realtype *ydata, *ypdata, *iddata,*constraintdata, c_0,l_a,l_s,phi1_a,phi1_c;
  long int i,j, N_t;
  N_t=N;
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
  phi1_c = ocp_cathode(data_cathode->c_s_0,data_cathode->c_s_max);
  
  /* Initialize y and yp on all grid points. */ 
  //Future refactoring: Can move all assignments to one block and pass arguments to it
  //j, c_s and phi1 one not used in separator. Initialized to zero for now. May cause issues later
  //Need to force to zero value in the equation as well
  for(i=0; i<N_t;i++){
      j=i/GRID; //temp  
      switch(j){
          case 0:
                 {
                  ydata[i]=c_0;
                  ypdata[i]=ZERO;
                  iddata[i]=ONE;   
                  //constraintdata[i]=TWO;
                  
                  break;
                 }
          
          case 1:
                 {
                  
                  if(i%GRID>=0 && i%GRID<=sep_indicator){
                      ydata[i]=phi1_a;
                      ypdata[i]=ZERO;
                      iddata[i]=ZERO;   
                      //constraintdata[i]=ZERO;
                  }
                  else if(i%GRID>sep_indicator && i%GRID<cath_indicator){
                      ydata[i]=ZERO;
                      ypdata[i]=ZERO;
                      iddata[i]=ZERO;   
                      //constraintdata[i]=ZERO;
                  }
                  else if(i%GRID>cath_indicator && i%GRID<=GRID-1){
                      ydata[i]=phi1_c;
                      ypdata[i]=ZERO;
                      iddata[i]=ZERO;
                      //constraintdata[i]=ZERO;
                  }
                  

                  break;
                 }
          
          case 2:
                 {
                  ydata[i]=ZERO;
                  ypdata[i]=ZERO;
                  iddata[i]=ZERO;   
                  //constraintdata[i]=ZERO;
                  break;
                 }
          
          case 3:
                 {
                  ydata[i]=ZERO;
                  ypdata[i]=ZERO;
                  iddata[i]=ZERO;  
                  //constraintdata[i]=ZERO;
                  break;
                 }  
           
          case 4:
                 {
                  if(i%GRID>=0 && i%GRID<=sep_indicator){
                      ydata[i]=data_anode->c_s_0;
                      ypdata[i]=ZERO;
                      iddata[i]=ONE;   
                      //constraintdata[i]=ONE;
                      
                  }
                  else if(i%GRID>sep_indicator && i%GRID<cath_indicator){
                      ydata[i]=ZERO;
                      ypdata[i]=ZERO;
                      iddata[i]=ONE;   
                      //constraintdata[i]=ONE;
                  }
                  else if(i%GRID>cath_indicator && i%GRID<=GRID-1){
                      ydata[i]=data_cathode->c_s_0;
                      ypdata[i]=ZERO;
                      iddata[i]=ONE;   
                      //constraintdata[i]=TWO;
                      
                  }
 
                  break;
                 }
  
           }
  
  }
}    
///Initializes user supplied data into data_anode struct
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

///Initializes user supplied data into data_sep struct
static void InitSepData(Material_Data data_sep)
{  
  data_sep->dx = ONE/(RCONST(GRID) - ONE);
  data_sep->coeff = ONE/(data_sep->dx);
  data_sep->eps = RCONST(0.724);
  data_sep->sigma_eff = (data_sep->sigma)*(ONE-data_sep->eps);
  data_sep->diff_coeff = RCONST(7.5e-10);
  data_sep->diff_coeff_eff = (data_sep->diff_coeff)*SUNRpowerR((data_sep->eps),BRUGG);
  
  data_sep->c_0 = RCONST(1000.0);
  data_sep->l = RCONST(25.0e-6);
}

///Initializes user supplied data into data_cathode struct
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

///Initializes data from the 3 Material_data structs into struct of Cell_Data
static void InitCellData(Material_Data data_anode,Material_Data data_sep,Material_Data data_cathode,Cell_Data data)
{
  realtype l_a, l_s;
  
  data->dx_a = data_anode->dx;
  data->coeff_a = data_anode->coeff;
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

///Returns the electrolyte conductivity at a particular point as a function of electrolyte concentration and the porosity of the medium

 
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
