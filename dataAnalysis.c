/******************************************************************************
 * Data analysis						       	      *
 * Last modification: 29 June 2019					      *
 * Organize the LAMMPS output file: atom.dump                                 *
 * Creates the result files: 01aTime.dat and 01bParticle.dat	              *
 * Creates the initial configuration for other simulations: 01cIC.dat	      *
 * Modified from 01patchy.c: We are not going to consider the guest rods in   *
 * the phase classification						      *
 ******************************************************************************/

#include "mutual.h"

int main(int argc, char *argv[]){
	Natoms = Natoms_f*Nx*Ny*Nz;
	Nbonds=(Natoms_f-1)*Nx*Ny*Nz;
	Nangles=(Natoms_f-2)*Nx*Ny*Nz;

	int Nhost=Nx*Ny*Nz;
	int Nguest=0;
	Nf = Nhost+Nguest; 
	//Nf = (Nx*Ny-1)*Nz+Nz*0.5;	

	//Storing per particle quantities
	ponto *P;
	P=malloc(3*Natoms*sizeof(float));	
	
	int *ID;
	ID=malloc(Natoms*sizeof(int));
	
	int *TYPE;
	TYPE=malloc(Natoms*sizeof(int));

	int *LOCAL;
        LOCAL=malloc(Natoms*sizeof(int));

	int *MOLTYPE;
        MOLTYPE=malloc(Nf*sizeof(int));

	//Storing time dependent quantities
	ponto *boxSize;
	boxSize = malloc(3*(tc/tdc+1)*sizeof(int));

	float *fractionUP;
	fractionUP = malloc((tc/tdc+1)*sizeof(float));

	//Storing particle related quatities
	//Computing each particle direction, length, and position of the center of mass
	ponto *CM;
	CM=malloc(3*Nf*sizeof(float));
	for(i=0;i<Nf;i++)
	{
 		CM[i].x = 0.0;
		CM[i].y = 0.0;
		CM[i].z = 0.0;
	}

	ponto *U;
        U=malloc(3*Nf*sizeof(float));
        for(i=0;i<Nf;i++)
        {
                U[i].x = 0.0;
                U[i].y = 0.0;
                U[i].z = 0.0;
        }

	float *bondSum;
	bondSum = malloc(Nf*sizeof(float));
	for(i=0;i<Nf;i++) bondSum[i] = 0.0;

	int recomecar=TRUE, chance=TRUE;

	char timek[] = "ITEM: NUMBER OF ATOMS";
	char inits[] = "ITEM: TIMESTEP";
	char *word;
	char insINT[100];
	float number1, number2;
	int charCom;

	ponto beadRef, bondSize;
	float bondModulus, direcModulus;

	//Computing the nematic OP
	ponto diagA={0.0,0.0,0.0}, diagB={0.0,0.0,0.0}, restA={0.0,0.0,0.0}, restB={0.0,0.0,0.0};

	ponto *eigenVec;
	eigenVec = malloc(3*(tc/tdc+1)*sizeof(float));
		
	float *S;
	S = malloc((tc/tdc+1)*sizeof(float));
	for(i=0;i<tc/tdc+1;i++) S[i] = 0.0;

	float p, q, detB, phi, eigenVecModulus, lambda;

	//Computing the smectic OP
	ponto deltaCM, direcREF, REF={0.0,1.0,0.0}, projPart;
	
	float direcREFmodulus, direcDotDeltaCM, projCM, projPartModulus;
	float cosTheta, sumNeigCos, sumNeigSen, sumCos, sumSen;

	float partSizeAvg=0.0;
		
	int neigCount=0, hexaCount=0;

	float *PSI;
        PSI = malloc((tc/tdc+1)*sizeof(float));
        for(i=0;i<tc/tdc+1;i++) PSI[i] = 0.0;

	//Computing the ratio of patches up and down
	int directorUP, counterUP;

	float eigenVecDotU;

	//Computing the volume fraction
	int compAtom=Natoms_f;
	float L0, v0, volFrac;

        if(compAtom==13) L0=6.4;
        else L0=10.7;
	printf("ZERO PRESSURE LENGTH %.1f\n", L0);

	v0=(PI/4)*L0+(PI/6);

	//counters for the atoms of a given particle... long guest in the smectic phase project :-)
	int cNatoms_f, maxNatoms_f;
	int cTESTE=0;

	pInter=atoi(argv[1]);

	printf("Analysis from p=%d.%d to %d.%d...\n", pInter, pIniInter, pInter, pFimInter);

	for(m=pIniInter;m<pFimInter+1;m++)
	{
		if(pInter==0 && m==0) m++;

		printf("\nFolder p%d_%d...\n", pInter, m);

		sprintf(fileName,"../p%d_%d/01aTime.dat", pInter, m);
                printf("\nCreating %s file...\n",fileName);
                fileB=fopen(fileName,"w");

		sprintf(fileName,"../p%d_%d/01bParticle.dat", pInter, m);
                printf("\nCreating %s file...\n",fileName);
                fileC=fopen(fileName,"w");

		sprintf(fileName,"../p%d_%d/01cIC.dat", pInter, m);
                printf("\nCreating %s file...\n",fileName);
                fileD=fopen(fileName,"w");

		fprintf(fileD,"\n%d atoms\n%d bonds\n%d angles\n\n",Natoms,Nbonds,Nangles);	
		fprintf(fileD,"2 atom types\n1 bond types\n2 angle types\n\n");

		sprintf(fileName,"../p%d_%d/dump.atom", pInter, m);
		printf("\nReading %s file...\n",fileName);
		fileA=fopen(fileName,"r");
		if(fileA==NULL) msgErro(fileName);

		i=0;
		c1 = 0;	
		while(!feof(fileA))
		{
			if(fscanf(fileA, "%d %d %f %f %f", &ID[i], &TYPE[i], &P[i].x, &P[i].y, &P[i].z)==5)
			{
				LOCAL[ID[i]-1]=i;
				i++;
			}
			else
			{
	 			fgets(notUseful,100,fileA);

				charCom = strncmp(notUseful,timek,21);
				if(charCom == 0)
				{
		 			printf("TIME STEP: %d\n", ID[i]);
					k = ID[i];
				}

				charCom = strncmp(notUseful,inits,14);
				if(i > Nf && (charCom == 0 || feof(fileA))) //No particular reason, just a large number 
				{
					fprintf(fileC, "TIME STEP: %d\n", k);					
					
					cNatoms_f=0;
					c2=1;

					for(j=0;j<Natoms;j++)
					{
						cNatoms_f++;

						//Print the IC file from the last frame
						if(k==20000000) fprintf(fileD,"%d %d %d %f %f %f\n",ID[LOCAL[j]],c2,TYPE[LOCAL[j]],Lx*P[LOCAL[j]].x,Ly*P[LOCAL[j]].y,Lz*P[LOCAL[j]].z);

						if(CM[c2-1].x == 0.0 && CM[c2-1].y == 0.0 && CM[c2-1].z == 0.0)
						{ 
							beadRef.x = Lx*P[LOCAL[j]].x;
							beadRef.y = Ly*P[LOCAL[j]].y;
							beadRef.z = Lz*P[LOCAL[j]].z;

							bondModulus = 0.0;
						}
						else
						{
							bondSize.x = Lx*P[LOCAL[j]].x-beadRef.x;
							bondSize.y = Ly*P[LOCAL[j]].y-beadRef.y;
							bondSize.z = Lz*P[LOCAL[j]].z-beadRef.z;
		
							if(bondSize.x > Lx/2) P[LOCAL[j]].x = P[LOCAL[j]].x-1;
							else if(bondSize.x < -Lx/2) P[LOCAL[j]].x = P[LOCAL[j]].x+1;
							bondSize.x = Lx*P[LOCAL[j]].x-beadRef.x;

							if(bondSize.y > Ly/2) P[LOCAL[j]].y = P[LOCAL[j]].y-1;
							else if(bondSize.y < -Ly/2) P[LOCAL[j]].y = P[LOCAL[j]].y+1;
							bondSize.y = Ly*P[LOCAL[j]].y-beadRef.y;

							if(bondSize.z > Lz/2) P[LOCAL[j]].z = P[LOCAL[j]].z-1;
							else if(bondSize.z < -Lz/2) P[LOCAL[j]].z = P[LOCAL[j]].z+1;
							bondSize.z = Lz*P[LOCAL[j]].z-beadRef.z;

							bondModulus = pow(bondSize.x,2)+pow(bondSize.y,2)+pow(bondSize.z,2);
							bondModulus = sqrt(bondModulus);
							if(bondModulus > bondL+1.0) printf("eita! %.2f \n", bondModulus);
						
							beadRef.x = Lx*P[LOCAL[j]].x;
							beadRef.y = Ly*P[LOCAL[j]].y;
							beadRef.z = Lz*P[LOCAL[j]].z;
						}
	
					 	CM[c2-1].x = CM[c2-1].x + Lx*P[LOCAL[j]].x;
						CM[c2-1].y = CM[c2-1].y + Ly*P[LOCAL[j]].y;
						CM[c2-1].z = CM[c2-1].z + Lz*P[LOCAL[j]].z;

						bondSum[c2-1] = bondSum[c2-1]+bondModulus;
						
						if(bondModulus != 0.0)
						{
							U[c2-1].x = U[c2-1].x+bondSize.x/bondModulus;
							U[c2-1].y = U[c2-1].y+bondSize.y/bondModulus;
							U[c2-1].z = U[c2-1].z+bondSize.z/bondModulus;
						}

                                                if(cNatoms_f == Natoms_f) 
						{
                                                        
							MOLTYPE[c2-1]=TYPE[LOCAL[j]];
                                                        c2++;
                                                        cNatoms_f=0;
                                                }

					}			

					partSizeAvg = 0.0;

					for(j=0;j<Nf;j++) 
					{	
						if(MOLTYPE[c2-1] == 1) cNatoms_f=Natoms_f;		

						//Position of the of CM for each particle
						CM[j].x = CM[j].x/Natoms_f;
						CM[j].y = CM[j].y/Natoms_f;
						CM[j].z = CM[j].z/Natoms_f;

						U[j].x = U[j].x/(Natoms_f-1);
                                                U[j].y = U[j].y/(Natoms_f-1);
                                                U[j].z = U[j].z/(Natoms_f-1);

						direcModulus = sqrt(pow(U[j].x,2)+pow(U[j].y,2)+pow(U[j].z,2));

						U[j].x = U[j].x/direcModulus;
						U[j].y = U[j].y/direcModulus;
						U[j].z = U[j].z/direcModulus;

						//Place CM inside box
						if(CM[j].x > Lx) CM[j].x = CM[j].x-Lx;
						else if(CM[j].x < 0.0) CM[j].x = CM[j].x+Lx;
						if(CM[j].y > Ly) CM[j].y = CM[j].y-Ly;
                                                else if(CM[j].y < 0.0) CM[j].y = CM[j].y+Ly; 
						if(CM[j].z > Lz) CM[j].z = CM[j].z-Lz;
                                                else if(CM[j].z < 0.0) CM[j].z = CM[j].z+Lz;

						//Compute the nematic OP
						diagA.x = diagA.x+1.5*pow(U[j].x,2)-0.5;
						diagA.y = diagA.y+1.5*pow(U[j].y,2)-0.5;
						diagA.z = diagA.z+1.5*pow(U[j].z,2)-0.5;

						restA.x = restA.x+1.5*U[j].y*U[j].z;
						restA.y = restA.y+1.5*U[j].x*U[j].z;
						restA.z = restA.z+1.5*U[j].x*U[j].y;

						partSizeAvg = partSizeAvg+bondSum[j];
					}

					partSizeAvg = partSizeAvg/Nhost;
					printf("PARTICLE SIZE %f\n", partSizeAvg+1.0);

					diagA.x = diagA.x/(Nx*Ny*Nz);
					diagA.y = diagA.y/(Nx*Ny*Nz);
					diagA.z = diagA.z/(Nx*Ny*Nz);

					restA.x = restA.x/(Nx*Ny*Nz);
					restA.y = restA.y/(Nx*Ny*Nz);
					restA.z = restA.z/(Nx*Ny*Nz);

					q = (diagA.x+diagA.y+diagA.z)/3.0;

					p = sqrt((2*pow(restA.x,2)+2*pow(restA.y,2)+2*pow(restA.z,2)+pow(diagA.x-q,2)+pow(diagA.y-q,2)+pow(diagA.z-q,2))/6.0);					

					diagB.x = (diagA.x-q)/p;
					diagB.y = (diagA.y-q)/p;
					diagB.z = (diagA.z-q)/p;

					restB.x = restA.x/p;
					restB.y = restA.y/p;
					restB.z = restA.z/p;

					detB = (diagB.x*diagB.y*diagB.z+2*restB.x*restB.y*restB.z-pow(restB.x,2)*diagB.x-pow(restB.y,2)*diagB.y-pow(restB.z,2)*diagB.z)/2.0;

					if(detB <= -1.0) phi = PI/3.0;
					else if(detB >= 1.0) phi = 0.0;
					else phi = acos(detB)/3.0;

					lambda = q+2*p*cos(phi);

					eigenVec[k/tdc].x = 1.0;
					eigenVec[k/tdc].z = (pow(restA.z,2)-(diagA.x-lambda)*(diagA.y-lambda))/((diagA.y-lambda)*restA.y-restA.z*restA.x);
					eigenVec[k/tdc].y = -1.0*(diagA.x-lambda+restA.y*eigenVec[k/tdc].z)/restA.z;

					eigenVecModulus = sqrt(pow(eigenVec[k/tdc].x,2)+pow(eigenVec[k/tdc].y,2)+pow(eigenVec[k/tdc].z,2));

					eigenVec[k/tdc].x = eigenVec[k/tdc].x/eigenVecModulus;
					eigenVec[k/tdc].y = eigenVec[k/tdc].y/eigenVecModulus;
					eigenVec[k/tdc].z = eigenVec[k/tdc].z/eigenVecModulus;

					REF.x = 0.0;
					REF.y = 1.0;
					REF.z = 0.0;
					
					if(sqrt(eigenVec[k/tdc].y) > 0.9)
					{
						REF.x = 0.0;
						REF.y = 0.0;
						REF.z = 1.0;
					}

					direcREF.x = REF.y*eigenVec[k/tdc].z-REF.z*eigenVec[k/tdc].y;
					direcREF.y = REF.z*eigenVec[k/tdc].x-REF.x*eigenVec[k/tdc].z;
					direcREF.z = REF.x*eigenVec[k/tdc].y-REF.y*eigenVec[k/tdc].x;
					
					direcREFmodulus = sqrt(pow(direcREF.x,2)+pow(direcREF.y,2)+pow(direcREF.z,2));
					direcREF.x = direcREF.x/direcREFmodulus;
					direcREF.y = direcREF.y/direcREFmodulus;
					direcREF.z = direcREF.z/direcREFmodulus;

					hexaCount = 0;
					sumCos = 0.0;
					sumSen = 0.0;

					counterUP = 0;

					for(j=0;j<Nf;j++) //if(MOLTYPE[j] == 1)
					{
						eigenVecDotU = eigenVec[k/tdc].x*U[j].x+eigenVec[k/tdc].y*U[j].y+eigenVec[k/tdc].z*U[j].z;
	
	 					S[k/tdc] = S[k/tdc]+pow(eigenVecDotU,2);

						neigCount = 0;
						sumNeigCos = 0.0;
						sumNeigSen = 0.0;
	
						for(l=0;l<Nf;l++)
						{	
							if(l != j)
							{
								deltaCM.x = CM[l].x-CM[j].x;
								deltaCM.y = CM[l].y-CM[j].y;
								deltaCM.z = CM[l].z-CM[j].z;

								//Due to the periodic boundary conditions the particles are never more distant than half box side
								if(deltaCM.x > boxSize[k/tdc].x/2.0) deltaCM.x = deltaCM.x-boxSize[k/tdc].x;
								else if(deltaCM.x < -boxSize[k/tdc].x/2.0) deltaCM.x = deltaCM.x+boxSize[k/tdc].x;
								if(deltaCM.y > boxSize[k/tdc].y/2.0) deltaCM.y = deltaCM.y-boxSize[k/tdc].y;
								else if(deltaCM.y < -boxSize[k/tdc].y/2.0) deltaCM.y = deltaCM.y+boxSize[k/tdc].y;
								if(deltaCM.z > boxSize[k/tdc].z/2.0) deltaCM.z = deltaCM.z-boxSize[k/tdc].z;
								else if(deltaCM.z < -boxSize[k/tdc].z/2.0) deltaCM.z = deltaCM.z+boxSize[k/tdc].z;

 								direcDotDeltaCM = eigenVec[k/tdc].x*deltaCM.x+eigenVec[k/tdc].y*deltaCM.y+eigenVec[k/tdc].z*deltaCM.z;
					
								projCM = sqrt(pow(deltaCM.x,2)+pow(deltaCM.y,2)+pow(deltaCM.z,2)-pow(direcDotDeltaCM,2));

								if(direcDotDeltaCM <= (partSizeAvg+1.0)*0.5 && direcDotDeltaCM >= -(partSizeAvg+1.0)*0.5 && projCM < 1.7)
								{
									neigCount++;

									projPart.x = deltaCM.x-direcDotDeltaCM*eigenVec[k/tdc].x; 
									projPart.y = deltaCM.y-direcDotDeltaCM*eigenVec[k/tdc].y;
									projPart.z = deltaCM.z-direcDotDeltaCM*eigenVec[k/tdc].z;

									projPartModulus = sqrt(pow(projPart.x,2)+pow(projPart.y,2)+pow(projPart.z,2));

									direcREFmodulus = sqrt(pow(direcREF.x,2)+pow(direcREF.y,2)+pow(direcREF.z,2));

									cosTheta = (projPart.x*direcREF.x+projPart.y*direcREF.y+projPart.z*direcREF.z)/(projPartModulus*direcREFmodulus);

									if(cosTheta > 1.0) cosTheta = 1.0;
									else if(cosTheta < -1.0) cosTheta = -1.0;

									sumNeigCos = sumNeigCos+cos(6*acos(cosTheta));
									sumNeigSen = sumNeigSen+sin(6*acos(cosTheta));			
								}		
							}

						}

						if(neigCount >= 1)
						{
							sumNeigCos = sumNeigCos/(float)neigCount;
							sumNeigSen = sumNeigSen/(float)neigCount;
						
							hexaCount++;

							sumCos = sumCos+sumNeigCos;
							sumSen = sumSen+sumNeigSen;
						}

						fprintf(fileC, "%d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %f %d\n", j+1, U[j].x, U[j].y, U[j].z, bondSum[j], CM[j].x, CM[j].y, CM[j].z, sqrt(pow(sumNeigCos,2)+pow(sumNeigSen,2)), neigCount);

						//Clear variables for next particle
						CM[j].x = 0.0; 
                                                CM[j].y = 0.0;
		                                CM[j].z = 0.0;                                                      
                          
                             			bondSum[j] = 0.0;           
                            
                                                U[j].x = 0.0;
                                                U[j].y = 0.0;        
                                                U[j].z = 0.0;
					}

					S[k/tdc] = S[k/tdc]/Nhost;
					S[k/tdc] = 1.5*S[k/tdc]-0.5;	
	
					printf("NEMATIC OP: %f\n", S[k/tdc]);

					sumCos = sumCos/(float)hexaCount;
                                        sumSen = sumSen/(float)hexaCount;

                                        PSI[k/tdc] = sqrt(pow(sumCos,2)+pow(sumSen,2));

					printf("HEXAGONAL OP: %f\n", PSI[k/tdc]);

					i = 0;
				}
				sprintf(insINT, "%d", ID[i]);
				
				strcat(insINT,notUseful);
				
				number1 = strtof(insINT, &word);	
				number2 = strtof(word, NULL);

				if(number1 != 0.0 && number2 != 0.0)
				{
					if(c1 == 0)
		 			{
						Lx = number2-number1;
						if(k==20000000) fprintf(fileD,"%.1f %.1f xlo xhi\n",number1,number2);
					}
					if(c1 == 1)
					{
						Ly = number2-number1;
						if(k==20000000) fprintf(fileD,"%.1f %.1f ylo yhi\n",number1,number2);
					}
					if(c1 == 2)
					{
 						Lz = number2-number1;
						if(k==20000000)
						{ 
							fprintf(fileD,"%.1f %.1f zlo zhi\n\n",number1,number2);
							fprintf(fileD,"Masses\n\n1 1\n2 1\n\n");
							fprintf(fileD,"Bond Coeffs\n\n1 50.0 0.5\n\n");
							fprintf(fileD,"Angle Coeffs\n\n1 100.0 180.0\n2 45.0 180.0\n\n");
							fprintf(fileD,"Atoms\n\n");
						}
					}

					c1++;
					if(c1 > 2)
					{ 
						c1 = 0;

						printf("BOX SIZE: %f %f %f\n", Lx, Ly, Lz);

						boxSize[k/tdc].x = Lx;
                                        	boxSize[k/tdc].y = Ly;
                                        	boxSize[k/tdc].z = Lz;
					} 
				}

				
			} 
	        }

		for(j=0;j<=tc/tdc;j++)
		{
			volFrac = v0*Nx*Ny*Nz/(boxSize[j].x*boxSize[j].y*boxSize[j].z);

			if(j==0) fprintf(fileB, "%d %f 1.0 0.0 0.0 1.0 1.0 %f %f %f 0.5\n", j*tdc, volFrac, boxSize[j].x, boxSize[j].y, boxSize[j].z); 
			else fprintf(fileB, "%d %f %f %f %f %f %f %f %f %f %f\n", j*tdc, volFrac, S[j], eigenVec[j].x, eigenVec[j].y, eigenVec[j].z, PSI[j], boxSize[j].x, boxSize[j].y, boxSize[j].z, fractionUP[j]);
			S[j] = 0.0;
		}
		
		fclose(fileA);
		fclose(fileB);
		fclose(fileC);

		fclose(fileD);

	}

	free(P);
	P = NULL;
 
	free(ID);
	ID = NULL;

	free(TYPE);
	TYPE = NULL;

	free(MOLTYPE);
        MOLTYPE = NULL;

	free(LOCAL);
	LOCAL = NULL;

	free(CM);
	CM = NULL;

	free(U);
        U = NULL;

	free(bondSum);
        bondSum = NULL;

	free(boxSize);
	boxSize = NULL;

	free(fractionUP);
	fractionUP = NULL;

	free(eigenVec);
	eigenVec = NULL;

	free(S);
	S = NULL;

	free(PSI);
        PSI = NULL;

	printf("\nEnd.\n\n");
	return 0;
}
