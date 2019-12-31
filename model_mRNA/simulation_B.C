#include "cvode.h"
#include "nvector_serial.h"
#include "cvode_dense.h"
#include "sundials_dense.h"
#include "main_mRNA.h"

///For pre-process
#define _STEADY_STATE_ 1e-5 //When the change of dynamics is less than the value,it is assumed as steady-state

static double *X;
static realtype *y_zero;
static CVDlsDenseJacFn dFdy = NULL;

// Input your ordinary differential equations
static void diffeq(int track, realtype t, N_Vector y, N_Vector ydot)
{

    double ERK, IKK, ERK_wt, IKK_wt, ERK_pt, IKK_pt;
    double inTF1, inTF2, Prom;
    double rk5, rk9, rk10;

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

//// Mass conservation equation

    inTF1 = 1 - Y(aTF1);
    inTF2 = 1 - Y(aTF2);

    Prom = 1 - Y(TF1_P) - Y(TF2_P) - Y(TF1_TF2_P);

//// Rate formula

    rk5 = X[k5]*X[k4];
    rk9 = X[k9]*X[k4];
    rk10 = X[k10]*rk9;

    realtype V[] = {// this variable is each reaction term
	0,										//[0]
	X[k1]*inTF1/(X[km1] + inTF1),							//[1] :: Activation of TF1
	X[k2]*(X[omegaI2]*IKK + X[omegaE2]*ERK)*inTF1/(X[km1] + inTF1),			//[2] :: Activation of TF1
	X[k3]*Y(aTF1)/(X[km3] + Y(aTF1)),						//[3] :: Activation of TF1
	X[k4]*pow(Y(aTF1), X[h4])*Prom,							//[4] :: TF1 - promoter intereaction
	rk5*Y(TF1_P),									//[5] :: TF1 - promoter intereaction
	X[k6]*inTF2/(X[km6] + inTF2),							//[6] :: Activation of TF2
	X[k7]*(X[omegaI7]*IKK + X[omegaE7]*ERK)*inTF2/(X[km6] + inTF2),			//[7] :: Activation of TF2
	X[k8]*Y(aTF2)/(X[km8] + Y(aTF2)),						//[8] :: Activation of TF2
	rk9*pow(Y(aTF2), X[h9])*Prom,							//[9] :: TF2 - promoter intereaction
	rk10*Y(TF2_P),									//[10] :: TF2 - promoter intereaction
	X[k11]*pow(Y(aTF1), X[h4])*Y(TF2_P),						//[11] :: TF1 - TF2_promoter intereaction
	X[k12]*pow(Y(aTF2), X[h9])*Y(TF1_P),						//[12] :: TF2 - TF1_promoter intereaction
	X[xi4]*X[k5]*X[k11]*Y(TF1_TF2_P),						//[13] :: TF1 - TF2_promoter intereaction
	X[xi9]*X[k10]*X[k12]*Y(TF1_TF2_P),						//[14] :: TF2 - TF1_promoter intereaction
	X[k19]*Y(TF1_P)/4,								//[15] :: mRNA production
	X[k19]*Y(TF2_P)/4,								//[16] :: mRNA production
	X[k19]*Y(TF1_TF2_P)/2,								//[17] :: mRNA production
	X[k20]*Y(mRNA)									//[18] :: mRNA production
	};

//// Differential equation

    dY(aTF1) = V[1] + V[2] - V[3];
    dY(aTF2) = V[6] + V[7] - V[8];
    dY(TF1_P) = V[4] - V[5] - V[12] + V[14];
    dY(TF2_P) = V[9] - V[10] - V[11] + V[13];
    dY(TF1_TF2_P) = V[11] + V[12] - V[13] - V[14];
    dY(mRNA) = V[15] + V[16] + V[17] - V[18];

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

