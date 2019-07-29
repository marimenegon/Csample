#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#define nMAX 10
#define NfMAX 3*(nMAX+1)*(nMAX+1) 

#define TRUE 1
#define FALSE 0

#define PI 3.141592

#define T 10000

#define bondL 0.5
#define latticeS 3.0
#define Natoms_f 21

#define LLp 0.3
#define bondCoeff 50.0
#define angleInit 180.0

#define Nx 16
#define Ny 18
#define Nz 16

#define tdc 100000
#define tc 10000000
#define ti 0

#define pIniInter 5
#define pFimInter 5

/* variables */

int i, j, k, l, m;

FILE *fileA, *fileB, *fileC, *fileD;
char fileName[100];

int pInter=0;

float Lx, Ly, Lz;

float stdDev=0;

int Natoms=1, Nbonds=1, Nangles=1, Nf=0;

int c1=0, c2=0;

int BOOL, boolErro=TRUE;

char notUseful[100];

struct xyz{
	float x;
	float y;
	float z;
};

typedef struct xyz ponto;


/* Functions */

float randFloat(float a){
	float b;
	b=a*(float)rand()/(float)RAND_MAX; /* Result between 0 and a */
	return b-a/2;      /* Result between -a/2 and a/2 */
}

float checkProj(float x1, float y1, float z1, float x2, float y2, float z2, float x, float y, float z, float *dist){
	ponto vecR, vecS, vecN;
	float mod2R, RdotS, modN;
	ponto P;
	float dist1, dist2;

	vecR.x=x2-x1;
	vecR.y=y2-y1;
	vecR.z=z2-z1;
	mod2R=pow(vecR.x,2)+pow(vecR.y,2)+pow(vecR.z,2);

	vecS.x=x-x1;
	vecS.y=y-y1;
	vecS.z=z-z1;

	RdotS=vecR.x*vecS.x+vecR.y*vecS.y+vecR.z*vecS.z;

	if((RdotS/mod2R>=0) && (RdotS/mod2R<=1)){
		P.x=x1+(RdotS/mod2R)*vecR.x;
		P.y=y1+(RdotS/mod2R)*vecR.y;
		P.z=z1+(RdotS/mod2R)*vecR.z;

		vecN.x=x-P.x;
		vecN.y=y-P.y;
		vecN.z=z-P.z;

		*dist=sqrt(pow(vecN.x,2)+pow(vecN.y,2)+pow(vecN.z,2));

		return 1;	
	}
	else{
		dist1=pow(x-x1,2)+pow(y-y1,2)+pow(z-z1,2);
		dist2=pow(x-x2,2)+pow(y-y2,2)+pow(z-z2,2);

		if(dist1<dist2) *dist=sqrt(dist1);
		else *dist=sqrt(dist2);
		return 0;
	}
}

void msgErro(char *nomeArq){
	printf("Error: file %s inexistent... End.\n", nomeArq);
	exit(0);	
}

void msgErro2Proceed(char *nomeArq){
        printf("Error: file %s inexistent... End.\n", nomeArq);
	boolErro=TRUE;
} 

float Heaviside(float x){
	if(x<0) return 0.0;
	if(x==0) return 0.5;
	if(x>0) return 1.0;
	else return 0.0; 
}

float AvgStdDev(int a, float *vec, float *stdDev){
	int counter;
	float avg=0;
 
	*stdDev=0;
	
	for(counter=0;counter<a;counter++) avg=avg+*(vec+counter);
	avg=avg/(float)a;
	for(counter=0;counter<a;counter++) *stdDev=*stdDev+pow(*(vec+counter)-avg,2);
	*stdDev=sqrt(*stdDev/a);
	return avg;
}
