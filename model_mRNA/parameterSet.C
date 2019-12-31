#include "main_mRNA.h"

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
	params[M_n+km2]		=pSign[5] * exp(getData[5]);
	params[M_n+k3]		=pSign[6] * exp(getData[6]);
	params[M_n+km3]		=pSign[7] * exp(getData[7]);
	params[M_n+k4]		=pSign[8] * exp(getData[8]);
	params[M_n+h4]		=pSign[9] * exp(getData[9]);
	params[M_n+k5]		=pSign[10] * exp(getData[10]);
	params[M_n+k6]		=pSign[11] * exp(getData[11]);
	params[M_n+km6]		=pSign[12] * exp(getData[12]);
	params[M_n+k7]		=pSign[13] * exp(getData[13]);
	params[M_n+omegaI7]	=pSign[14] * exp(getData[14]);
	params[M_n+omegaE7]	=pSign[15] * exp(getData[15]);
	params[M_n+km7]		=pSign[16] * exp(getData[16]);
	params[M_n+k8]		=pSign[17] * exp(getData[17]);
	params[M_n+km8]		=pSign[18] * exp(getData[18]);
	params[M_n+k9]		=pSign[19] * exp(getData[19]);
	params[M_n+h9]		=pSign[20] * exp(getData[20]);
	params[M_n+k10]		=pSign[21] * exp(getData[21]);
	params[M_n+xi4]		=pSign[22] * exp(getData[22]);
	params[M_n+xi9]		=pSign[23] * exp(getData[23]);
	params[M_n+k11]		=pSign[24] * exp(getData[24]);
	params[M_n+k12]		=pSign[25] * exp(getData[25]);
	params[M_n+k19]		=pSign[26] * exp(getData[26]);
	params[M_n+k20]		=pSign[27] * exp(getData[27]);

}

