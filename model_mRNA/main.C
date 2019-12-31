//**********************************
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "main_mRNA.h"


extern int simulation(int track, double *params, realtype tstep[], int tnum, double **yout);
extern void parameterSetting(double *params, double *getData, int *pSign);

void init_yout(double**);
void deft_yout(double**);

//***************************************************************************
int main(int argc,char *argv[])
{
	int i,j,ret,ind;
	char str1[32],str2[2],str3[2],str4[32],str5[32];
	double getData[K_n];
	int pSign[K_n];
	char **pNota = NULL;
	int nline;
	static double *params = NULL;


	if (!(argc == 3)) {
		fprintf(stderr,
			"2 argument (parameter set# & savefile name) is required!\n");
		goto END;
	}

	pNota = (char**)malloc(sizeof(char*) * (K_n));
	params = (double *)malloc(sizeof(double) * (M_n+K_n));

  //*******************************************
  {
    	char *file;
    	struct stat buf;
    	FILE *fp = NULL;

    	file = argv[1];
    	if (stat(file, &buf) == -1) {
      		fprintf(stderr,"Parameter Setting File(%s) read error!\n",file);
      		goto END;
    	} else if ((fp = fopen(file, "r")) == NULL) {
      		fprintf(stderr,"Parameter Setting File(%s) read error!\n",file);
      		goto END;
    	} else {
      		fprintf(stderr,"Loading a parameter set from %s\n",file);

    		nline=0;
			for (;;) {
    			ret = fscanf(fp,"%s %s %s %s %s",str1,str2,str3,str4,str5);
				nline++;
    			if (ret == EOF) break;
				if (nline<3) continue;
				ind = nline - 3;
				pNota[ind] = (char*)malloc(sizeof(char) * (strlen(str1)+1));
				strcpy(pNota[ind],str1);
				if (strcmp(str3, "N") == 0) {
				    pSign[ind] = -1;
				} else if (strcmp(str3,"0") == 0) {
				    pSign[ind] = 0;
				} else {
				    pSign[ind] = 1;
				}
				getData[ind] = atof(str4);
			}
			parameterSetting(params, getData, pSign);
    	}
    	if (fp != NULL) fclose(fp);
  }

	for(i=0;i<M_n;i++){
		fprintf(stderr,"Y0[%d]=%e\n",i,params[i]);
	}
	for(i=0;i<K_n;i++){
		fprintf(stderr,"%s=%e\n",pNota[i],params[M_n+i]);
	}

  //-------------------------------------------------Main Simulation


  {
	realtype tstep[MEMALLOC];

	int rv_wt,rv_erk,rv_ikk;
	double **yout_wt = NULL;
	double **yout_erk = NULL;
	double **yout_ikk = NULL;

	for(i=0;i<MEMALLOC;i++){
		tstep[i] = (double)i*(double)STEPSIZE;
	}

	if (yout_wt!=NULL) {deft_yout(yout_wt);free(yout_wt);yout_wt = NULL;}
	yout_wt = (double**)malloc(sizeof(double*) * (MEMALLOC));
	init_yout(yout_wt);

	if (yout_erk!=NULL) {deft_yout(yout_erk);free(yout_erk);yout_erk = NULL;}
	yout_erk = (double**)malloc(sizeof(double*) * (MEMALLOC));
	init_yout(yout_erk);

	if (yout_ikk!=NULL) {deft_yout(yout_ikk);free(yout_ikk);yout_ikk = NULL;}
	yout_ikk = (double**)malloc(sizeof(double*) * (MEMALLOC));
	init_yout(yout_ikk);

	rv_wt = simulation(_No_perturbation_, params, tstep, MEMALLOC, yout_wt);
	rv_erk = simulation(_erk_inhibition_, params, tstep, MEMALLOC, yout_erk);
	rv_ikk = simulation(_ikk_inactivation_, params, tstep, MEMALLOC, yout_ikk);

	if((rv_wt==0)||(rv_erk==0)||(rv_ikk==0)){
	    FILE *fo;
	    char fileName[1024];

	    sprintf(fileName,"%s", argv[2]);
	//Write file
	    fo=fopen(fileName,"w");
	    fprintf(fo,"Time\taTF1n\taTF2n\tTF1-DNAn\tTF2-DNAn\tTF1-TF2-DNAn\tmRNAn\taTF1e\taTF2e\tTF1-DNAe\tTF2-DNAe\tTF1-TF2-DNAe\tmRNAe\taTF1i\taTF2i\tTF1-DNAi\tTF2-DNAi\tTF1-TF2-DNAi\tmRNAi\n");
	    for(i=0;i<MEMALLOC;i++){
		fprintf(fo,"%d",(int)tstep[i]);// Output time

		if (rv_wt==0)
		    for(j=0;j<M_n;j++)
			fprintf(fo,"\t%e", yout_wt[i][j]);// Output time courses of active TF, TF-promoter complex, gene expression level
		else
		    for(j=0;j<M_n;j++)
			fprintf(fo,"\tnd");

		if (rv_erk==0)
		    for(j=0;j<M_n;j++)
			fprintf(fo,"\t%e", yout_erk[i][j]);// Output time courses of active TF, TF-promoter complex, gene expression level
		else
		    for(j=0;j<M_n;j++)
			fprintf(fo,"\tnd");

		if (rv_ikk==0)
		    for(j=0;j<M_n;j++)
			fprintf(fo,"\t%e", yout_ikk[i][j]);// Output time courses of active TF, TF-promoter complex, gene expression level
		else
		    for(j=0;j<M_n;j++)
			fprintf(fo,"\tnd");

		fprintf(fo,"\n");
	    }
	    fclose(fo);

	    fprintf(stderr,"***%s\n", argv[2]);
	}else{
	    fprintf(stderr,"Simulation completely error!\n");
	}

	if (yout_wt != NULL) {deft_yout(yout_wt);free(yout_wt);}yout_wt = NULL;
	if (yout_erk != NULL) {deft_yout(yout_erk);free(yout_erk);}yout_erk = NULL;
	if (yout_ikk != NULL) {deft_yout(yout_ikk);free(yout_ikk);}yout_ikk = NULL;
  }

END:
	if(pNota != NULL) {
	    for (i=0;i<K_n;i++)
		if (pNota[i] != NULL) {free(pNota[i]);pNota[i] = NULL;}
	    free(pNota);pNota = NULL;
	}
	if (params != NULL) free(params);params = NULL;
	return 0;
}

void init_yout(double **yout)
{
    int i;

    for(i=0;i<MEMALLOC;i++) {
	yout[i] = NULL;
	yout[i] = (double*)malloc(sizeof(double) * (M_n));
    }
}

void deft_yout(double **yout)
{
    int i;

    for(i=0;i<MEMALLOC;i++)
	if (yout[i] != NULL) {free(yout[i]);yout[i] = NULL;}
}

