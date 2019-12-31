#ifndef __MAIN_TF__
#define __MAIN_TF__


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundials_types.h"

#define RTOL    1e-9
#define ATOL    1e-9
#define Y(i) NV_Ith_S(y, (i))
#define dY(i) NV_Ith_S(ydot, (i))
#define T0 RCONST(0.0)

#define _No_perturbation_  0
#define _erk_inhibition_   1
#define _ikk_inactivation_ 2

#define _OFFSET_ 500.0 // Time until steady-state from initial state

#define STEPSIZE 1 // Step size in time-course simulation
#define ENDTIME 90 // End time of simulation
#define MEMALLOC 91 // = ENDTIME/STEPSIZE + 1

//contant value, search parameter
  enum {
		k1,
		km1,
		k2,
		omegaI2,
		omegaE2,
		k3,
		km3,
		K_n};
//dependent variable
  enum {
		aTF,
  		M_n};			//[0]

#endif

