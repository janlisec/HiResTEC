#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "R.h"
#include "Rdefines.h"
#include <R_ext/Rdynload.h>

#undef TRUE
#define TRUE    1
#undef FALSE
#define FALSE    0

//#define DEBUG

struct scanStruct
{
   double  mz;
   double  intensity;
};

int lowerBound(double val,double *mzval,int first, int length){
  int half,mid;
  while (length > 0) {
    half = length >> 1;
    mid = first;
    mid += half;
    if ( mzval[mid] < val){
      first = mid;
      first ++;
      length = length - half -1;
    }
    else length = half;
  }
  return(first);
}

int upperBound(double val,double *mzval,int first, int length){
  int half,mid;
  while (length > 0) {
    half = length >> 1;
    mid = first;
    mid += half;
    if (val < mzval[mid]){
      length = half;
    }
    else {
      first = mid;
      first ++;
      length = length - half -1;
    }
  }
  return(first);
}

struct scanStruct getScanBPC(int scan,
			     double from,
			     double to,
			     double *pmz,
			     double *pintensity,
			     int *pscanindex,
			     int nmz,
			     int lastScan,
			     int returneic) {
  int idx,idx1,idx2;
  double max=0.0;
  double mzatmax=0.0;
  struct scanStruct res;

  idx1 =  pscanindex[scan -1] +1;
  if (scan == lastScan)  idx2 =  nmz-1;
  else idx2 =  pscanindex[scan];

  int idx1b = lowerBound(from, pmz, idx1-1, idx2-idx1);
  int idx2b = upperBound(to, pmz, idx1b, idx2-idx1b);

  for (idx=idx1b; idx <= idx2b; idx++)
    {
      double mzval = pmz[idx-1];
    if(!returneic) {
      if ((mzval <= to) && (mzval >= from) && (pintensity[idx-1] > max)) {
        max = pintensity[idx-1];
        mzatmax = mzval;
      }
    } else {
      if ((mzval <= to) && (mzval >= from)) {
        max += pintensity[idx-1];
        mzatmax = mzval;
      }
    }
    }
  res.mz = mzatmax;
  res.intensity = max;
  return(res);
}

int CMA(double **values, double *averages, int size, int periods)
  // simple moving average, borrowed from https://stackoverflow.com/questions/12884600/how-to-calculate-simple-moving-average-faster-in-c
  // modified to include centering (lag correction)
{
  int i;
  double sum = 0;
  int lag = (int)ceil((periods - 1) / 2);
  for (i = 0; i < size; i++) {
    if (i < periods) {
      sum += *values[i];
      if (i - lag >= 0) {
        averages[i-lag] = (i == periods - 1) ? sum / (double)periods : 0;
      }
    } else {
      sum = sum - *values[i - periods] + *values[i];
      averages[i-lag] = sum / (double)periods;
    }
  }
  for (i = size - 1 - lag; i < size; i++) {
    averages[i] = 0; // fill the rest with 0's
  }
  return (size - periods + 1 > 0) ? size - periods + 1 : 0; // number of calculated averages
}


SEXP getMultipleBPC_C_(SEXP mz,
		       SEXP intensity,
		       SEXP scanindex,
		       SEXP tgtmz,
		       SEXP mzdev_lower,
		       SEXP mzdev_upper,
		       SEXP scanrange,
		       SEXP lastscan,
		       SEXP smooth,
		       SEXP returnEIC) {
  double *pmz, *pintensity, *p_vint, *ptgtmz, *pmzdev_lower, *pmzdev_upper, *p_vmz;
  int i,j,k,l, *pscanindex, *p_scan, scanrangeFrom, scanrangeTo, ilastScan, nmz, ctScan, buflength, returneic;
  int ntgtmz, nscans;
  int *psmooth;
  SEXP list_names, reslist, vscan, vint, vmz;
  pmz = REAL(mz);
  nmz = GET_LENGTH(mz);
  pintensity = REAL(intensity);
  pscanindex = INTEGER(scanindex);
  psmooth = INTEGER(smooth);
  int firstScan = 1;   // is always 1
  ilastScan = INTEGER(lastscan)[0];
  ptgtmz = REAL(tgtmz);
  ntgtmz = GET_LENGTH(tgtmz);
  pmzdev_lower = REAL(mzdev_lower);
  pmzdev_upper = REAL(mzdev_upper);
  double mzrangeFrom[ntgtmz];
  double mzrangeTo[ntgtmz];
  for(i=0; i<ntgtmz; i++) {
    mzrangeFrom[i] = ptgtmz[i] - pmzdev_lower[i];
    mzrangeTo[i] = ptgtmz[i] + pmzdev_upper[i];
  }
  scanrangeFrom = INTEGER(scanrange)[0];
  scanrangeTo = INTEGER(scanrange)[1];
  nscans = scanrangeTo - scanrangeFrom + 1;
  if ((scanrangeFrom <  firstScan) || (scanrangeFrom > ilastScan) || (scanrangeTo < firstScan) || (scanrangeTo > ilastScan))
    error("Error in scanrange \n");
  char *names[3] = {"scan", "intensity", "mz"};
  PROTECT(list_names = allocVector(STRSXP, 3));
  for(i = 0; i < 3; i++)
    SET_STRING_ELT(list_names, i,  mkChar(names[i]));
  returneic = INTEGER(returnEIC)[0];

  /* buflength = scanrangeTo - scanrangeFrom +1; */
  buflength = nscans * ntgtmz;
  PROTECT(reslist = allocVector(VECSXP, 3));
  PROTECT(vscan = NEW_INTEGER(buflength));
  p_scan = INTEGER_POINTER(vscan);
  PROTECT(vint = NEW_NUMERIC(buflength));
  p_vint = NUMERIC_POINTER(vint);
  PROTECT(vmz = NEW_NUMERIC(buflength));
  p_vmz = NUMERIC_POINTER(vmz);

  // get data for requested scans and m/z ranges and save into one long vector
  i=0;
  for (ctScan=scanrangeFrom;ctScan<=scanrangeTo;ctScan++) {
    for (j=0; j<ntgtmz; j++) {
      p_scan[i] = ctScan;
      struct scanStruct tmp = getScanBPC(ctScan,mzrangeFrom[j],mzrangeTo[j],pmz,pintensity,pscanindex,nmz,ilastScan,returneic);
      p_vint[i] = tmp.intensity;
      p_vmz[i] = tmp.mz;
      i++;
    }
  }

  // smooth if requested
  if (*psmooth) {
    double *pp_vint[nscans]; // array of pointers to p_vint (effectively pointer to pointer)
    double averages[nscans], *paverages; // buffer for averages
    paverages = averages;
    for (i=0; i<ntgtmz; i++) { // for each target mass ...
      k=0;
      for (j=i; j<buflength; j=j+ntgtmz) { // .. pass index of addresses within linear intensity vector to smoother ...
        pp_vint[k++] = &p_vint[j];
      }
      CMA(pp_vint, paverages, nscans, *psmooth);
      for(l=0; l<nscans; l++) {
        *pp_vint[l] = paverages[l]; // ... and finally replace original intensities by averaged ones.
      }
    }
  }

  SET_VECTOR_ELT(reslist, 0, vscan); // attaching integer vector scan to list
  SET_VECTOR_ELT(reslist, 1, vint); // attaching double vector intentsity to list
  SET_VECTOR_ELT(reslist, 2, vmz); // attaching double vector m/z to list
  setAttrib(reslist, R_NamesSymbol, list_names); //and attaching the vector names

  UNPROTECT(5);
  return(reslist);
}

// register native routines (i.e. .C, .Call, etc.)
static const R_CallMethodDef callMethods[]  = {
  {"getMultipleBPC_C", (DL_FUNC) &getMultipleBPC_C_, 10},
  {NULL, NULL, 0}
};

void R_init_HiResTEC(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
