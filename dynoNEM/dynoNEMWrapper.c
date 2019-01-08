#include "netlearn.h"


void dynoNEMWrapper(int* T, int* nsgenes, int* negenes, double* D_R, double* initial_R, double* network_prior_R, double* Egene_prior_R, double* prior_scale, int* type, int* nrep, double* alpha, double* beta, double *res_network, double* res_likelihood);
void dynoNEM_getPerturbProb(double* Psi_R, int* T, int* nsgenes, int* k, double* res);

void dynoNEMWrapper(int* T, int* nsgenes, int* negenes, double* D_R, double* initial_R, double* network_prior_R, double* Egene_prior_R, double* prior_scale, int* type, int* nrep, double* alpha, double* beta, double *res_network, double* res_likelihood){
    double*** D = (double***) R_alloc(*T, sizeof(double**));
    double** initial = (double**) R_alloc(*nsgenes, sizeof(double*));
    double** network_prior = (double**) R_alloc(*nsgenes, sizeof(double*));
    double** network = (double**) R_alloc(*nsgenes, sizeof(double*));
    double** Egene_prior = (double**) R_alloc(*negenes, sizeof(double*));
    int t, i, s, k;
    double likelihood;
    for(t = 0; t < *T; t++){
        D[t] = (double**) R_alloc(*negenes, sizeof(double*));
        for(i = 0; i < *negenes; i++){
            D[t][i] = (double*) R_alloc(*nsgenes, sizeof(double));
            Egene_prior[i] = (double*) R_alloc(*nsgenes, sizeof(double));
            for(s = 0; s < *nsgenes; s++){
                D[t][i][s] = D_R[s * *T * *negenes + i * *T + t];
                Egene_prior[i][s] = Egene_prior_R[s * *nsgenes + i];
            }
        }
    }
    for(s = 0; s < *nsgenes; s++){
        initial[s] = (double*) R_alloc(*nsgenes, sizeof(double));
        network_prior[s] = (double*) R_alloc(*nsgenes, sizeof(double));
        network[s] = (double*) R_alloc(*nsgenes, sizeof(double));
        for(k = 0; k < *nsgenes; k++){
            initial[s][k] = initial_R[k * *nsgenes + s];
            network_prior[s][k] = network_prior_R[k * *nsgenes + s];
            network[s][k] = 0.0;
        }
    }
    likelihood = learn_network(*T, *nsgenes, *negenes, D, initial, network_prior, Egene_prior, *prior_scale, network, *type, *nrep, *alpha, *beta);
    Rprintf("final likelihood = %g\n", likelihood);
    res_likelihood[0] = likelihood;
    for(s = 0; s < *nsgenes; s++){
        for(k = 0; k < *nsgenes; k++){
            res_network[k * *nsgenes + s] = network[s][k];
        }
    }
}

void dynoNEM_getPerturbProb(double* Psi_R, int* T, int* nsgenes, int* k, double* res){
    double** Psi = (double**) R_alloc(*nsgenes, sizeof(double*));
    int s, q, t;
    for(s = 0; s < *nsgenes; s++){
        Psi[s] = (double*) R_alloc(*nsgenes, sizeof(double));
        for(q = 0; q < *nsgenes; q++)
            Psi[s][q] = Psi_R[q* *nsgenes + s];
    }
    double** perturb_prob = getPerturbProb(Psi, *T, *nsgenes, *k);
    for(s = 0; s < *nsgenes; s++){
        for(t = 0; t < *T; t++){
            res[t * *nsgenes + s] = perturb_prob[s][t];
        }
    }
}
