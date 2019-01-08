/*
 * netlearn.h
 *
 *  Created on: Aug 2, 2010
 *      Author: frohlich
 */

#ifndef NETLEARN_H_
#define NETLEARN_H_

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <float.h>

#define PVAL_DENS 0
#define EFFECT_PROB 1
#define DISCRETE 2

double** getPerturbProb(double** Psi, int T, int nsgenes, int k);
double learn_network(int T, int nsgenes,int negenes, double*** D, double** initial, double** network_prior, double** Egene_prior, double prior_scale, double **net, int type, int nrep, double alpha, double beta);
double network_likelihood(double** Psi, int nsgenes, int negenes, int T, double*** D, double** egene_prior, int type, int nrep, double alpha, double beta);
double logPrior(int nsgenes, double** net, double** prior, double nu);
void array_copy(int nsgenes,double** net_matrix, double** net_copy);

#endif /* NETLEARN_H_ */
