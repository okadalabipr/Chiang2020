#include "main_TF.h"

void parameterSetting(double *params, double *getData, int *pSign){

	int i;

	for(i=0;i<(M_n+K_n);i++) params[i]=0.0;//Initialization

  //----------------------------------Parameter values
	//If parameter is search parameter, getData[i] is set.
	//If parameter is constant parameter, value is set.
	params[M_n+k1]		=pSign[0] * exp(getData[0]);
        params[M_n+km1]		=pSign[1] * exp(getData[1]);
	params[M_n+k2]		=pSign[2] * exp(getData[2]);
	params[M_n+omegaI2]	=pSign[3] * exp(getData[3]);
	params[M_n+omegaE2]	=pSign[4] * exp(getData[4]);
	params[M_n+k3]		=pSign[5] * exp(getData[5]);
	params[M_n+km3]		=pSign[6] * exp(getData[6]);

}

