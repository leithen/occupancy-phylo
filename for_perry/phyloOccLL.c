#include <Rmath.h>

inline double expit(double x) {
  return(1./(1.+exp(-x)));
}

void sample_phylo_occupancy_ll(double *expit_p_0_samples, 
			       double *psi_beta_samples, 
			       double *psi_site_samples, 
			       double *psi_year_samples, 
			       double *env, 
			       double *psi_0_vec, 
			       int *boolXzero,
			       int *boolZzero,
			       int *pnumSamples, 
			       int *pnumSites, 
			       int *pnumYears, 
			       int *pnumReps, 
			       int *pnumSpecies, 
			       double *LL ) {
  int iSample, iSite, iYear, iRep, iSpecies, iX;
  double sumQQ, psi_sp_int, psi_sp_slope, psi, p;
  int numSamples = *pnumSamples;
  int numSites = *pnumSites;
  int numYears = *pnumYears;
  int numReps = *pnumReps;
  int numSpecies = *pnumSpecies;
  double C1, C2;
  
  for(iSample = 0; iSample < numSamples; ++iSample) {
    sumQQ = 0;
    iX = 0;
    for(iSpecies = 0; iSpecies < numSpecies; ++iSpecies) {
      p = expit_p_0_samples[iSample + iSpecies * numSamples];
      C2 = psi_beta_samples[iSample + iSpecies * numSamples];
      for(iRep = 0; iRep < numReps; ++iRep) {
	for(iYear = 0; iYear < numYears; ++iYear) {
	  C1 = psi_year_samples[iSample + iYear * numSamples] + psi_0_vec[iSpecies];
	  for(iSite = 0; iSite < numSites; ++iSite) {
	    psi_sp_int = psi_site_samples[iSample + iSite * numSamples] + C1;
	    psi_sp_slope = env[iSite] * C2;
	    psi = expit(psi_sp_int + psi_sp_slope);

	    /* p.mat[boolXzero] <- 1-p.mat[boolXzero] */
	    /* QQ <- log(psi.mat*p.mat + (1-psi.mat)*boolZzero) */

	    if(boolXzero[iX]==0) {
	      sumQQ += log(psi * p);
	    }
	    if(boolXzero[iX]==1 && boolZzero[iX]==0) {
	      sumQQ += log(psi * (1.-p));
	    }
	    if(boolXzero[iX]==1 && boolZzero[iX]==1) {
	      sumQQ += log(psi * (1.-p) + (1.-psi));
	    }
	    
	    /* if(boolXzero[iX]) { */
	    /*   sumQQ += log(psi * (1.-p) + (1.-psi)); */
	    /* } else { */
	    /*   sumQQ += log(psi * p); */
	    /* } */

	    ++iX;
	  }
	}
      }
    }
    LL[iSample] = sumQQ;
  }
};

