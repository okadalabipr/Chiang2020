#ifndef __MAIN_RNA__
#define __MAIN_RNA__


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

#define _OFFSET_ 5000.0 // Time until steady-state from initial state

#define STEPSIZE 1 // Step size in time-course simulation
#define ENDTIME 90 // End time of simulation
#define MEMALLOC 91 // = ENDTIME/STEPSIZE + 1

//contant value, search parameter
  enum {
		k1,		//[0]	Natural maximal activation rate of TF1 :: Activation of TF1
		km1,		//[1]	Half-maximal activation concentration of the inactive TF1 :: Activation of TF1
		k2,		//[2]	Coefficient of TF1 relative level to proteins :: Activation of TF1
		omegaI2,	//[3]	Coefficient of the active IKK on TF1 activation :: Activation of TF1
		omegaE2,	//[4]	Coefficient of the active ERK on TF1 activation :: Activation of TF1
		km2,		//[5]	With IKK and ERK control, the half-maximal activation concentration of the inactive TF1 :: Activation of TF1
		k3,		//[6]	Maximal inactivation rate of TF1 :: Activation of TF1
		km3,		//[7]	Half-maximal inactivation concentration of the active TF1 :: Activation of TF1
		k4,		//[8]	Association constant of TF1 binding to DNA :: Sequential model
		h4,		//[9]	Number of the working TF1 on the promoter :: Sequential model
		k5,		//[10]	Dissociation constant of TF1 leaving DNA :: Sequential model
		k6,		//[11]	Natural maximal activation rate of TF2 :: Activation of TF2
		km6,		//[12]	Half-maximal activation concentration of the inactive TF2 :: Activation of TF2
		k7,		//[13]	Coefficient of TF2 relative level to proteins :: Activation of TF2
		omegaI7,	//[14]	Coefficient of the active IKK on TF2 activation :: Activation of TF2
		omegaE7,	//[15]	Coefficient of the active ERK on TF2 activation :: Activation of TF2
		km7,		//[16]	With IKK and ERK control, the half-maximal activation concentration of the inactive TF2 :: Activation of TF2
		k8,		//[17]	Maximal inactivation rate of TF2 :: Activation of TF2
		km8,		//[18]	Half-maximal inactivation concentration of the active TF2 :: Activation of TF2
		k9,		//[19]	Association constant of TF2 binding to DNA :: Sequential model
		h9,		//[20]	Number of the working TF2 on the promoter :: Sequential model
		k10,		//[21]	Dissociation constant of TF2 leaving DNA :: Sequential model
		xi4,		//[22]	Cooperativity of the bound TF2 on TF1 binding to the TF2_promoter complex; xi4 > 1 means that TF1 becomes easier to bind because of the bound TF2, and vice versa :: Sequential model
		xi9,		//[23]	Cooperativity of the bound TF1 on TF2 binding to the TF1_promoter complex; xi9 > 1 means that TF2 becomes easier to bind because of the bound TF1, and vice versa :: Sequential model
		k11,		//[24]	Dissociation constant of TF1 leaving the TF2_promoter complex :: Sequential model
		k12,		//[25]	Dissociation constant of TF2 leaving the TF1_promoter complex :: Sequential model
		k19,		//[26]	Synthesis rate constant of mRNA by the TF1_TF2_promoter complex :: mRNA production
		k20,		//[27]	Degradation rate constant of mRNA :: mRNA production
		K_n };
//dependent variable
  enum {
		aTF1,		//[0]
		aTF2,		//[1]
		TF1_P,		//[2]
		TF2_P,		//[3]
		TF1_TF2_P,	//[4]
		mRNA,		//[5]
		M_n };

#endif

