#include "cvode.h"
#include "nvector_serial.h"
#include "cvode_dense.h"
#include "sundials_dense.h"
#include "main_TF.h"

///For pre-process
#define _STEADY_STATE_ 1e-5 //When the change of dynamics is less than the value,it is assumed as steady-state

static double *X;
static realtype *y_zero;
static CVDlsDenseJacFn dFdy = NULL;

// Input your ordinary differential equations
static void diffeq(int track, realtype t, N_Vector y, N_Vector ydot)
{

    double ERK, IKK, ERK_wt, IKK_wt, ERK_pt, IKK_pt, inTF;

///*
    if(t<_OFFSET_){
	ERK_wt = 0.05;
	IKK_wt = 0.05;
	ERK_pt = 0.05;
	IKK_pt = 0.05;
    } else if(t<_OFFSET_+5 && t>=_OFFSET_){
	ERK_wt = 71.4125*(t-_OFFSET_)/5;
	IKK_wt = 89.5233*(t-_OFFSET_)/5;
	ERK_pt = 0.05;
	IKK_pt = 0.05;
    } else if(t<_OFFSET_+10 && t>=_OFFSET_+5){
	ERK_wt = 71.4125-((71.4125-58.63)*(t-_OFFSET_-5)/5);
	IKK_wt = 89.5233-((89.5233-54.5133)*(t-_OFFSET_-5)/5);
	ERK_pt = 0.05;
	IKK_pt = 0.05;
    } else if(t<_OFFSET_+15 && t>=_OFFSET_+10){
	ERK_wt = 58.63-((58.63-50.1)*(t-_OFFSET_-10)/5);
	IKK_wt = 54.5133-((54.5133-42.9967)*(t-_OFFSET_-10)/5);
	ERK_pt = 0.05;
	IKK_pt = 0.05;
    } else if(t<_OFFSET_+30 && t>=_OFFSET_+15){
	ERK_wt = ((54.3475-50.1)*(t-_OFFSET_-15)/15)+50.1;
	IKK_wt = 42.9967-((42.9967-42.3367)*(t-_OFFSET_-15)/15);
	ERK_pt = 0.05;
	IKK_pt = 0.05;
    } else if(t<_OFFSET_+45 && t>=_OFFSET_+30){
	ERK_wt = ((87.8325-54.3475)*(t-_OFFSET_-30)/15)+54.3475;
	IKK_wt = 42.3367-((42.3367-18.83)*(t-_OFFSET_-30)/15);
	ERK_pt = 0.05;
	IKK_pt = 0.05;
    } else if(t<_OFFSET_+60 && t>=_OFFSET_+45){
	ERK_wt = 87.8325-((87.8325-84.9225)*(t-_OFFSET_-45)/15);
	IKK_wt = 18.83-((18.83-17.0633)*(t-_OFFSET_-45)/15);
	ERK_pt = 0.05;
	IKK_pt = 0.05;
    } else if(t<_OFFSET_+90 && t>=_OFFSET_+60){
	ERK_wt = 84.9225-((84.9225-75.4775)*(t-_OFFSET_-60)/30);
	IKK_wt = 17.0633-((17.0633-4.6033)*(t-_OFFSET_-60)/30);
	ERK_pt = 0.05;
	IKK_pt = 0.05;
    } else if(t<_OFFSET_+120 && t>=_OFFSET_+90){
	ERK_wt = 75.4775-((75.4775-51.2975)*(t-_OFFSET_-90)/30);
	IKK_wt = 4.6033-((4.6033-2.1567)*(t-_OFFSET_-90)/30);
	ERK_pt = 0.05;
	IKK_pt = 0.05;
    } else if(t>=_OFFSET_+120){
	ERK_wt = 51.2975;
	IKK_wt = 2.1567;
	ERK_pt = 0.05;
	IKK_pt = 0.05;
    }// */

    if(track==_No_perturbation_) {
	ERK = ERK_wt;
	IKK = IKK_wt;
    } else if(track==_erk_inhibition_) {
	ERK = ERK_pt;
	IKK = IKK_wt;
    } else if(track==_ikk_inactivation_) {
	ERK = ERK_wt;
	IKK = IKK_pt;
    }

    ERK = ERK/100;
    IKK = IKK/100;

//// Mass conservation equation of TF

    inTF = 1 - Y(aTF);

//// Rate formula

    realtype V[] = {// this variable is each reaction term
	0,										//[0]
	X[k1]*inTF/(X[km1] + inTF),							//[1] :: Activation of TF
	X[k2]*(X[omegaI2]*IKK + X[omegaE2]*ERK)*inTF/(X[km1] + inTF),			//[2] :: Activation of TF
	X[k3]*Y(aTF)/(X[km3] + Y(aTF)),							//[3] :: Activation of TF
	};

//// Differential equation

    dY(aTF) = V[1] + V[2] - V[3];

}

static int f_wt(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
	diffeq(_No_perturbation_, t, y, ydot);
	return 0;
}

static int f_erk(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
	diffeq(_erk_inhibition_, t, y, ydot);
	return 0;
}

static int f_ikk(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
	diffeq(_ikk_inactivation_, t, y, ydot);
	return 0;
}

///////////////////////////////////////////////////////////////////////////////

static N_Vector make_atol(int neq, realtype tol)
{
	N_Vector atol = N_VNew_Serial(neq);
	int i;

	for (i = 0; i < neq; ++i) NV_Ith_S(atol, i) = tol;

	return atol;
}

static N_Vector make_vector(int neq, realtype val[])
{
	N_Vector vec = N_VNew_Serial(neq);
	int i;

	for (i = 0; i < neq; ++i) NV_Ith_S(vec, i) = val[i];

	return vec;
}

static int isproper(int neq, N_Vector y)
{
	for (int i = 0; i < neq; ++i)
    	if (Y(i) < -ATOL)
      		return 0;

	return 1;
}

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}

//track:experiment perturbation, *params:parameter, *tstep:time step, tnum:time step size, **yout:result of simulation
int simulation(int track, double *params, realtype tstep[], int tnum, double **yout)
{

	int i,j;
	int rv = 0;
	realtype reltol = RTOL;
  	void *cvode_mem = NULL;
  	realtype t,comptime;
	double bs[M_n];//values before stimulation for checking steady-state
	int flag0;

	//initialization of CVODE setting
  	y_zero = params;
  	X = params + M_n;

  	N_Vector abstol = make_atol(M_n, ATOL);
  	N_Vector y = make_vector(M_n, y_zero);

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, (char *)"CVodeCreate", 0)) return(1);
	CVodeSetMaxNumSteps(cvode_mem, 327670);//Increase max step

	if (track==_No_perturbation_) {
	    flag0 = CVodeInit(cvode_mem, f_wt, T0, y);
	} else if (track==_erk_inhibition_) {
	    flag0 = CVodeInit(cvode_mem, f_erk, T0, y);
	} else if (track==_ikk_inactivation_) {
	    flag0 = CVodeInit(cvode_mem, f_ikk, T0, y);
	}

	if (check_flag(&flag0, (char *)"CVodeInit", 1)) return(1);
	flag0 = CVodeSVtolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag0, (char *)"CVodeSVtolerances", 1)) return(1);

	/* Call CVDense to specify the CVDENSE dense linear solver */
	flag0 = CVDense(cvode_mem, M_n);
	if (check_flag(&flag0, (char *)"CVDense", 1)) return(1);
	/* Set the Jacobian routine to Jac (user-supplied) */
	flag0 = CVDlsSetDenseJacFn(cvode_mem, dFdy);
	if (check_flag(&flag0, (char *)"CVDlsSetDenseJacFn", 1)) return(1);
	CVodeSetMaxErrTestFails(cvode_mem, 20);

	//get values before stimulation for checking steady-state
    comptime = _OFFSET_ - _OFFSET_/10;

    int flag = CVode(cvode_mem, comptime, y, &t, CV_NORMAL);
    if (flag != CV_SUCCESS) { rv = 1; goto end; }
    if (!isproper(M_n, y)) { rv = 1; goto end; }
	for(i=0;i<M_n;i++) {
	    bs[i]=Y(i);
	}

  	for (i = 0; i < tnum; ++i) {
  		comptime = _OFFSET_ + tstep[i];
    	int flag = CVode(cvode_mem, comptime, y, &t, CV_NORMAL);
    	if (flag != CV_SUCCESS) { rv = 1; goto end; }
    	if (!isproper(M_n, y)) { rv = 1; goto end; }

  		//check whether steady state or not
  		if(i==0){
  			for(j=0;j<M_n;j++){
  				if(fabs(Y(j)-bs[j]) > _STEADY_STATE_){
  					rv = 1; goto end;
  				}
  			}
  		}

		for(j=0;j<M_n;j++) {
		    if (Y(j) < 0) {rv = 1; goto end;}
		}

  		//save simulation data
  		for(j=0;j<M_n;j++)
  			yout[i][j] = Y(j);
  	}


	end:
		N_VDestroy_Serial(y);
		N_VDestroy_Serial(abstol);
		if (cvode_mem != NULL) CVodeFree(&cvode_mem);

	return rv;
}

