/*
###############################################################
#
#
#
#
#
###############################################################
*/
// Complete code to execute all the functions for the 
// metropolis hastings MCMC for dynoNEM
// the wrapper function for r is provided separately as wrapper.c file


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netlearn.h>
#include "functionsMCMC.h"
#define INF 999999
#define CONVERGENCE 0.05
#define ACCEPT 0.25
#define SAMPLE 10000
#define BURNIN 10000



// copy the initial network 
// modify the copied network
// check update factor if acceptable 
// calculate the MI with the original network 
// count the number of accepted nework 
// if the number in greater than BURNIN and Mutual Information is less than acceptable
// start modifying the last network
// calculate edge probability for all the networks in the sample



//########## Function to calculate Hastings ratio /// checked and changes suggested 


double updateFactor (double **netOld, double **netNew, int nsgenes, int negenes, int T, double ***D, double **Egene_prior, int type, int nrep, double alpha, double beta, double hRatio){
	
	double neighOld = counterNeighbor(netOld, nsgenes, T);  //Neighbors of the old network
	double neighNew = counterNeighbor(netNew, nsgenes, T);  //Neighbors of the new network
	double likLogOld = network_likelihood(netOld, nsgenes, nsgenes, T, D, Egene_prior, type, nrep, alpha, beta);     //Likelihood for the old network
	double likLogNew = network_likelihood(netNew, nsgenes, nsgenes, T, D, Egene_prior, type, nrep, alpha, beta);     //Likelihood for the new network
	double ratio = (likLogNew - likLogOld + neighNew - neighOld);   // hastings ratio
	
	if (ratio < 1)
		hRatio = ratio;
	else 
		hRatio = 1;
	return(hRatio);  // minimum of 1 and hastings ratio
	
}

// ######### Function to find the number of neighbors of a network // checked and changes suggested

double counterNeighbor(double **adjMat, int nsgenes, int T)// Gives reversible and decreased neighbors
{
	int i, j, dec_count, inc_count, rev_count; 
	dec_count = 0;
	inc_count = 0;
	rev_count = 0;
	
	for (i=0; i<=nsgenes; i++)
    {
		for (j = 0; j <=nsgenes; j++) 
		{
			if(adjMat[i][j] > 0)
				dec_count++;
			if((adjMat[i][j] > 0) && (adjMat[i][j] < T))
				inc_count++;
			if((adjMat[i][j] != 0) && (adjMat[j][i] == 0))
				rev_count++;
			
		}
	}
	return(log(dec_count+inc_count+rev_count));
}

//####################### log likelihood calculation (Plugged in from dynoNEM) #################### //Function reviewed

double** getPerturbProb(double** Psi, int T, int nsgenes, int k){
	double** perturb_prob = (double**) calloc(nsgenes, sizeof(double*));
	int parent_perturb_prob;
	int s, t, p;
	for(s = 0; s < nsgenes; s++){
		perturb_prob[s] = (double*) calloc(T, sizeof(double));
	}
	for(t = 0; t < T; t++){
		for(s = 0; s < nsgenes; s++){
			perturb_prob[s][t] = 0; // default: s is unperturbed
			perturb_prob[k][t] = 1; // the perturbed gene is always inactive
			if(s != k){
				for(p = 0; p < nsgenes; p++){
					if(Psi[p][s] != 0 && abs(Psi[p][s]) <= t){ // p is a parent
						if(t > 0)
							parent_perturb_prob = perturb_prob[p][t-1];
						else
							parent_perturb_prob = (p == k);
						if(parent_perturb_prob){
							perturb_prob[s][t] = 1;
							break;
						}
					}
				}
			}
		}
	}
	
	return(perturb_prob);
}

double network_likelihood(double** Psi, int nsgenes, int negenes, int T, double*** D, double** egene_prior, int type, int nrep, double alpha, double beta){
	double*** perturb_prob = (double***) calloc(nsgenes, sizeof(double**));
	int s, k, t, i;
	for(k = 0; k < nsgenes; k++)
		perturb_prob[k] = (double**) getPerturbProb(Psi, T, nsgenes, k);
	
	double tmp;
	double loglik0;
	double loglik=0;
	double loglik_tmp;
	for (i=0; i<negenes; i++) {
		loglik_tmp=0;
		for (s=0; s<nsgenes; s++) {
			tmp=0;
			for (k=0; k<nsgenes; k++) {
				for (t=0; t<T; t++) {
					if(egene_prior[i][s] > 0)
						if(type == PVAL_DENS) // for p-value densities:
							tmp += log((D[t][i][k]*perturb_prob[k][s][t] + (1-perturb_prob[k][s][t])*1) * egene_prior[i][s]);
						else if(type == EFFECT_PROB) // for effect probabilities
							tmp += log((D[t][i][k]*perturb_prob[k][s][t] + (1-perturb_prob[k][s][t])*(1 - D[t][i][k])) * egene_prior[i][s]);
						else // for count data:
							tmp += log((pow(1-beta, D[t][i][k]*perturb_prob[k][s][t])*pow(beta, (nrep-D[t][i][k])*perturb_prob[k][s][t]) +
										pow(alpha, D[t][i][k]*(1-perturb_prob[k][s][t]))*pow(1-alpha, (nrep-D[t][i][k])*(1-perturb_prob[k][s][t]))) * egene_prior[i][s]);
				}
			}
			if (s==0) {
				loglik0 = tmp;
			}
			else {
				loglik_tmp += exp(tmp - loglik0);
			}
			
		}
		loglik += log(1 + loglik_tmp) + loglik0;
	}
	for(k = 0; k < nsgenes; k++){
		for(s = 0; s < nsgenes; s++)
			free(perturb_prob[k][s]);
		free(perturb_prob[k]);
	}
	free(perturb_prob);
	return(loglik);
}

//############ Function to change the network // ***** To be reviewed

double alterNet(double** net, int nsgenes, int negenes, int T, int number, double hRatio, double** newNet){
	//int step = 0;
	double* minfo = (double*) calloc(1, sizeof(double));
	double** temp1 = (double**) calloc(nsgenes, sizeof(double*));
	double** temp2 = (double**) calloc(nsgenes, sizeof(double*));
	int i, j, k;
	for (i = 0; i < nsgenes; i++) {
		temp1[i]= (double*) calloc(nsgenes, sizeof(double));
		temp2[i]= (double*) calloc(nsgenes, sizeof(double));
	}
	copyNet(nsgenes, net, temp1);
	
	int iter = 0;
	//double hRatio =0; // ***** Initial declaration of hRatio 
	double hratioInc, hratioDecr, hratioDel, hratioRev, type, nrep, alpha, beta;
	double*** D;
	double** Egene_Prior;
	while ((*minfo > CONVERGENCE) && (iter > number)) {
		iter++;
		for (i = 0; i < nsgenes; i++) {
			for (j = 0; j < nsgenes; j++) {
				if (i != j) {
					
					if ( temp1[i][j] < T) { //  to increment the edge weight
						copyNet(nsgenes, temp1, temp2);
						temp2[i][j] += 1;
						hratioInc = updateFactor(net, temp2, nsgenes, negenes, T, D, Egene_Prior, type, nrep, alpha, beta, hRatio);
						if (hratioInc > hRatio) {
							hRatio = hratioInc;
							
							 mutInfo(net, temp2, nsgenes, minfo);
							//addToList(temp2, iter, nsgenes, list);// not needed here, only for the second call
						}
					}
					
					if (temp1[i][j] > 0) {  //  to decrement the edge weight 
						copyNet(nsgenes, temp1, temp2);
						temp2[i][j] -= 1;
						hratioDecr = updateFactor(net, temp2, nsgenes, negenes, T, D, Egene_Prior, type, nrep, alpha, beta, hRatio);
						if (hratioDecr > hRatio) {
							hRatio = hratioDecr;
							mutInfo(net, temp2, nsgenes, minfo);
							//addToList(temp2, iter, nsgenes, list); // not needed here, only for the second call
						}
						
					}
					if ((temp1[i][j] != 0) && temp1[j][i] == 0) {//  to reverse the edge
						copyNet(nsgenes, temp1, temp2);
						temp2[j][i] = temp1[i][j];
						hratioRev = updateFactor(net, temp2, nsgenes, negenes, T, D, Egene_Prior, type, nrep, alpha, beta, hRatio);
						if (hratioRev > hRatio) {
							hRatio = hratioRev;
							mutInfo(net, temp2, nsgenes, minfo);
							//addToList(temp2, iter, nsgenes, list);// not needed here only for the second call
						}
					}
					if (temp1[i][j] != 0) { //  to delete an edge
						copyNet(nsgenes, temp1, temp2);
						temp2[i][j] = 0;
						hratioDel = updateFactor(net, temp2, nsgenes, negenes, T, D, Egene_Prior, type, nrep, alpha, beta, hRatio);
						if (hratioDel > hRatio) {
							hRatio = hratioDel;
							mutInfo(net, temp2, nsgenes, minfo);
							//addToList(temp2, iter, nsgenes, list);// not needed only for the second call
						}
					}
					
				}
			}
		}
	}
	return (temp2, hRatio);
	for(i = 0; i < nsgenes; i++){
		free(temp1[i]);
		free(temp2[i]);
	}
	free(temp1);
	free(temp2);
	
}


// ##################### Function for Mutual Information of vectors  //  Function reviewed last time // 

double mutInfo(double *v1, double *v2, int *ns, double *MI)  //function declaration to calculate mutual information
{  
	double *Pa, *Pb, *Pab, n; // pointer declaration to probability of 'a', 'b' and their joint probability repectively
	// Memory allocation to the pointers
	Pa = (double*) calloc(3, sizeof(int));
	Pb = (double*) calloc(3, sizeof(int));
	Pab = (double*) calloc(9, sizeof(int));
	//Assigning size of the vector
	int i,j,k;
	
	for (k=0; k<9; k++) 
		Pab[k] = 0;
	for (k = 0; k < 3; k++) {
		Pb[k] = 0;
		for (j = 0; j < *ns; j++)  
			if (v2[j] == k+1) 
				Pb[k]++;
	}
	for (k = 0; k < 3; k++) {
		Pa[k] = 0;
		for (j = 0; j < *ns; j++) {
			if (v1[j] == k+1) {
				Pa[k]++;
				for (i = 0; i < 3; i++) 
					if (v2[j] == i+1) 
						Pab[i*2 +k]++;
			}
		}
	}
	
	*MI = n * log(*ns);
	for (k = 0; k < 9; k++) 
		if (Pab[k]>0) 
			*MI = *MI + (Pab[k] * log(Pab[k])) ;
	
	for (k = 0; k < 3; k++) {
		if (Pb[k]>0) 
			*MI = *MI - (Pb[k] *log(Pb[k])) ;
		if (Pa[k]>0) 
			*MI = *MI - (Pa[k] *log(Pa[k])) ;
	}
	
	*MI = *MI / (double) *ns;
	*MI = *MI/ log(2);
	printf ("%f\n\n",*MI);
	free (Pa);
	free (Pb);
	free (Pab);
	
}

// ################# Converting 2D arrays to 1D vectors for MI // Function reviewed

void matConversion(double** A, int nsgenes, double* B){
int i, j, k;
B = (double*) calloc((nsgenes*nsgenes), sizeof(double));
k = 0;
for(i = 0;i < nsgenes; i++)
{
	for(j = 0;j<nsgenes;j++)
	{
		B[k++] = A[i][j];
	}
	
}
	return(B);
}

//############//Function to copy network  //  Function checked last time

void copyNet(int nsgenes, double** net, double** netCopy)
{
	int i, j;
	for (i = 0; i < nsgenes; i++){
		for (j = 0; j < nsgenes; j++){
			netCopy[i][j] = net[i][j];
		}
	}
}
// ############## the fuction to run the entire set of subroutines // ***** to be reviewed

double mcmcMH ( double** net, int nsgenes, int negenes, double* number, int T, double*** D, double **Egene_prior, double* type, int nrep, double* alpha, double* beta )  { 
  copyNet(nsGenes, net, netCopy)
  int counter = 0;
  int counter2 = 0;
  updateFactor(net, netCopy, nsgenes, negenes, T, D, Egene_prior, type, nrep, alpha, beta, hRatio);

  while ((counter < BURNIN) & (mutInfo(net, netCopy) > CONVERGENCE)){ // ***** Should the 2D net matrix be converted to 1D matrix for MI or the data is a 1D matrix already

	alterNet(netCopy, nsgenes, negenes, T, number, hRatio, newNet);
	updateFactor(net, newNet, nsgenes, negenes, T, D, Egene_prior, type, nrep, alpha, beta, hRatio_new);
		if(hRatio_new > hRatio){
		  hRatio = hRatio_new;
		  netCopy = newNet;
		  counter++;
      }
  
  }

  while(counter2 < SAMPLE){

	alterNet(netCopy, nsgenes, negenes, T, number, hRatio, newNet);
	updateFactor(net, newNet, nsgenes, negenes, T, D, Egene_prior, type, nrep, alpha, beta, hRatio_new);
		if(hRatio_new > hRatio){
		hRatio = hRatio_new;
		netCopy = newNet;
		addTolist(netCopy, list);// ***** doubts: How to add the networks to create the required sample just by appending to an empty matrix or some other ways
		counter2++;
    
	}
  return list;// if required return the list of the sampled network
  }

}

// calculating frequency in sampled networks // ***** To be reviewed
//############################ ###################
//
double edgeprob( double** net, int nsgenes, int negenes, double* number, int T, double*** D, double **Egene_prior, double* type, int nrep, double* alpha, double* beta){
double** matrix = (double**) calloc(nsgenes, sizeof(double*))
for(int i = 0; i<SAMPLE; i++){
  matrix[i] = (double*) calloc (nsgenes, sizeof(double))
}
int counter2 = 0;
while(counter2 < SAMPLE){
  
  alterNet(netCopy, nsgenes, negenes, T, number, hRatio, newNet);
	updateFactor(net, newNet, nsgenes, negenes, T, D, Egene_prior, type, nrep, alpha, beta, hRatio_new);
		if(hRatio_new > hRatio){
		hRatio = hRatio_new;
		netCopy = newNet;
		for(i= 0; i < nsgenes; i++){
			for(j = 0; j < nsgenes; j++){
			if(netCopy[i ][j ] > 0){
				matrix[i][j]++;
				}
			}
		 }
	}	
	counter2++;
 }
 
 for(i = 0; i < SAMPLE; i++){
	for(j = 0; j < SAMPLE; j++){
		prob_matrix[i][j] = matrix[i][j]/SAMPLE;
	}
 }
   return prob_matrix;
   for(i = 0; i < SAMPLE; i++){
		free(matrix[i]);
	}
	free(matrix);
}
    

