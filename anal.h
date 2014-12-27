#pragma once

#include <vector>
#include "readplink.h"
void wrap(const plink *p,const std::vector<double> &phe,const std::vector<double> &ad, const std::vector<double> &sex,const std::vector<double> &freq);

/*
  Below is the pars.
  This will contain all data needed for analyzing one site.
  The idea is that we allocate pars[nthreads].
  The internal arrays are nsites long. And will therefore not require reallocating.
  
  We assume no missing data (samples are filteret in the set_pars)

 */
typedef struct{
  int len; //length
  double **covs;//covs is a design matrix.
  int ncov; //number of covariates <ncov_max
  int ncov_max;//this is the maximum nr of covariates. Idea is to avoid realloc.

  char *gs;//genotypes gs \in {0,1,2,3},counts of allele2, 3=NA/missing
  double *ys;//pheno
  double *qvec;//admix prop
  double *mafs;//freqs
  int model;
  int quant;
  double *start;
  int maxIter;
  double *sex;
  double tol;
  int doOptim;
  int retPrior;
  int doLean;
  double *res;//should contain the persite results. Takes from Line this is a double[47]
}pars;
pars *init_pars(int l,int ncov);
void set_pars(pars*p,char *g,const std::vector<double> &phe,const std::vector<double> &ad, const std::vector<double> &sex,const std::vector<double> &freq);
void kill_pars(pars *p,int l);
