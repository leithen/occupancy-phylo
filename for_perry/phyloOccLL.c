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
  int iSample, iSite, iYear, iRep, iSpecies, iX, iZ;
  double prodPP, sumQQ, psi_sp_int, psi_sp_slope, psi, p;
  int numSamples = *pnumSamples;
  int numSites = *pnumSites;
  int numYears = *pnumYears;
  int numReps = *pnumReps;
  int numSpecies = *pnumSpecies;
  double C1, C2;
  
  for(iSample = 0; iSample < numSamples; iSample++) {
    sumQQ = 0;
    iZ = 0;

    for(iSpecies = 0; iSpecies < numSpecies; iSpecies++) {
      C2 = psi_beta_samples[iSample + iSpecies * numSamples];
      for(iYear = 0; iYear < numYears; iYear++) {
	C1 = psi_0_vec[iSpecies] + 
	  psi_year_samples[iSample + iYear * numSamples];	
	for(iSite = 0; iSite < numSites; iSite++) {

	  /* detectability */
	  prodPP = 1;
	  p = expit_p_0_samples[iSpecies * numSamples + iSample];
	  for(iRep = 0; iRep < numReps; iRep++) {
	    if(boolXzero[iSite + 
			 numSites*iYear + 
			 numSites*numYears*iRep +
			 numSites*numYears*numReps*iSpecies]==0) { prodPP *= p; }
	    if(boolXzero[iSite + 
			 numSites*iYear + 
			 numSites*numYears*iRep +
			 numSites*numYears*numReps*iSpecies]==1) { prodPP *= (1.-p); }
	  }

	  /* occupancy */
	  psi_sp_int = psi_site_samples[iSample + iSite * numSamples] + C1;
	  psi_sp_slope = env[iSite] * C2;
	  psi = expit(psi_sp_int + psi_sp_slope);

	  sumQQ += log(psi * prodPP + (1.-psi)*boolZzero[iZ++]);

	}
      }
    }
    LL[iSample] = sumQQ;
  }
};
