#include "netlearn.h"

double** getPerturbProb(double** Psi, int T, int nsgenes, int k){
        double** perturb_prob = (double**) Calloc(nsgenes, double*);
        int parent_perturb_prob;
        int s, t, p;
        for(s = 0; s < nsgenes; s++){
                perturb_prob[s] = (double*) Calloc(T, double);
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
        double*** perturb_prob = (double***) Calloc(nsgenes, double**);
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
                        Free(perturb_prob[k][s]);
                Free(perturb_prob[k]);
        }
        Free(perturb_prob);
        return(loglik);
}


void array_copy(int nsgenes, double** net_matrix, double** net_copy){
	int s, k;
	for (s=0; s<nsgenes; s++) {
		for (k=0; k<nsgenes; k++) {
			net_copy[s][k]=net_matrix[s][k];
		}
	}
}

double logPrior(int nsgenes, double** net, double** prior, double nu){
        int s;
        int k;

        if (prior==NULL) {
                return 0;
        }
        else {
                double p=0;
                for (s=0; s<nsgenes; s++) {
                        for (k=0; k<nsgenes; k++) {
                                p -= abs(net[s][k]-prior[s][k])/nu;
                        }
                }
                p -= nsgenes*nsgenes*log(2*nu);
                return p;
        }
        return 0;
}

void print_network(double** net, int nsgenes){
    int s, k;
    for(s = 0; s < nsgenes; s++){
        for(k = 0; k < nsgenes; k++)
            Rprintf("%g ", net[s][k]);
        Rprintf("\n");
    }
    Rprintf("============================\n");
}


double learn_network(int T, int nsgenes,int negenes, double*** D, double** initial, double** network_prior, double** Egene_prior, double prior_scale, double **net, int type, int nrep, double alpha, double beta){
        double** tmp = (double**) Calloc(nsgenes, double*);
        double** tmp2 = (double**) Calloc(nsgenes, double*);
        int i, s, k;
	for(i = 0; i < nsgenes; i++){
		tmp[i] = (double*) Calloc(nsgenes, double);
		tmp2[i] = (double*) Calloc(nsgenes, double);
	}
	double lik_incr;
	double lik_decr;
	double lik_rev;
	double lik_switch;
	int converged=0;
	if (Egene_prior==NULL) {		
		for (i=0; negenes; i++) {
			for (s=0; s< nsgenes; s++) {
				Egene_prior[i][s]=1/(nsgenes+1);
			}
		}
	}
	double loglik = network_likelihood(initial, nsgenes, negenes, T, D, Egene_prior, type, nrep, alpha, beta) + logPrior(nsgenes, initial, network_prior, prior_scale);
	array_copy(nsgenes, initial, net);
	print_network(net, nsgenes);
	//hill climbing algorithm implemented
	Rprintf("initial log-likelihood = %g\n", loglik);
	int iterator = 0;
	int improved = 0;
	while (converged == 0) {
	        improved = 0;
		iterator++;
		for (s=0; s<nsgenes; s++) {
			for (k=0; k<nsgenes; k++) {
			        if(s != k){
                                        //Check for increase in edge weight
                                        if (net[s][k]<T) {
                                                array_copy(nsgenes, net, tmp);
                                                tmp[s][k] += 1.0;
                                                lik_incr = network_likelihood(tmp, nsgenes, negenes, T, D, Egene_prior, type, nrep, alpha, beta) + logPrior(nsgenes, tmp, network_prior, prior_scale);
                                                if (lik_incr > loglik) {
                                                        loglik = lik_incr;
                                                        array_copy(nsgenes, tmp, tmp2);
                                                        improved = 1;
                                                }
                                        }

                                        //check for decrement in edge weight
                                        if (net[s][k]>0){
                                                array_copy(nsgenes, net, tmp);
                                                tmp[s][k] -= 1.0;
                                                lik_decr = network_likelihood(tmp, nsgenes, negenes, T, D, Egene_prior, type, nrep, alpha, beta) + logPrior(nsgenes, tmp, network_prior, prior_scale);
                                                if (lik_decr > loglik) { // loglik is now the likelihood after adding an edge!
                                                        loglik = lik_decr;
                                                        array_copy(nsgenes, tmp, tmp2);
                                                        improved = 1;
                                                }
                                        }

                                        //Check for edge reversal
                                        if ((net[s][k] != 0) && (net[k][s] == 0)) {
                                                array_copy(nsgenes, net, tmp);
                                                tmp[k][s] = tmp[s][k];
                                                tmp[s][k] = 0.0;
                                                lik_rev = network_likelihood(tmp, nsgenes, negenes, T, D, Egene_prior,  type, nrep, alpha, beta) + logPrior(nsgenes, tmp, network_prior, prior_scale);
                                                if (lik_rev > loglik) { // loglik is now max{likelihood after adding an edge, likelihood after removing an edge}!
                                                        loglik = lik_rev;
                                                        array_copy(nsgenes, tmp, tmp2);
                                                        improved = 1;
                                                }
                                        }
			        }
			}
		}
		if (!improved)
			break;
		array_copy(nsgenes, tmp2, net);
		print_network(net, nsgenes);
		Rprintf("\n\niteration %d: log-likelihood = %g\n", iterator, loglik);
	}

	for(i = 0; i < nsgenes; i++){
		Free(tmp[i]);
		Free(tmp2[i]);
	}
	Free(tmp);
	Free(tmp2);
	loglik = network_likelihood(net, nsgenes, negenes, T, D, Egene_prior,  type, nrep, alpha, beta);
	return(loglik);
}

