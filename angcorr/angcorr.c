#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_spline.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.1415926535897
#define ch 2998.0 //speed of light/(100 km/s/Mpc); or c*h Mpc

/*Get the window function squared.

  Args:
      l: angular wavenumber of the power spectrum
      k: wavenumber of the power spectrum; h/Mpc

  Returns:
      W^2(l/k): window function squared
 */
double get_W2(double l, double k){
  double x = l/k/ch; //unitless
  double sinx = sin(x);
  double cosx = cos(x);
  return (sinx-x*cosx)*(sinx-x*cosx)/x/x/x;
}

/*Integrand of the kernel.

  Note: the arguments have to be passed through the F_params struct.

  Args:
      l: angular wavenumber of the power spectrum
      theta: angle on the sky
      k: wavenumber of the power spectrum; h/Mpc
 */
typedef struct F_params{
  double theta;
  double k;
}F_params;
double F_integrand(double l, void*params){
  F_params pars = *(F_params*)params;
  double theta = pars.theta;
  double k = pars.k;
  return gsl_sf_bessel_J0(theta*l)*get_W2(k,l);
}

/*Angular correlation kernel, F(k,theta).

  Args:
      theta: angle on the sky
      k: wavenumber of the power spectrum; h/Mpc

  Returns:
      F(theta,k): angular correlation function kernel
 */
double F(double theta, double k){
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc(1000);

  double result,error;

  F_params * params = malloc(sizeof(F_params));
  params->theta = theta;
  params->k = k;
  
  gsl_function F;
  F.function = &F_integrand;
  F.params = &params;
  
  gsl_integration_qagiu(&F, 0, 1e-7, 1e-7, 1000, w, &result, &error);

  return result/k;
}

/*Evaluate the matter power spectrum at k.

  Args:
      ki: wavenumber to calculate P at; h/Mpc
      k: array of wavenumbers; h/Mpc
      P: matter power spectrum; (h/Mpc)^3
      Nk: length of k and P arrays
      alpha: P power law indices
      A: P power law amplitudes
      Pspl: spline for P
      acc: accelerator for the spline

  Returns:
     P(ki): matter power spectrum at ki
 */
double get_P(double ki,double*k,double*P,int Nk,
	     double*alpha,double*A,
	     gsl_spline*Pspl,gsl_interp_accel*acc){
  double kmin = k[0];
  double kmax = k[Nk-1];
  if (ki < kmin){
    return A[0]*pow(ki,alpha[0]);
  }else if (ki > kmax){
    return A[1]*pow(ki,alpha[1]);
  }// Assume power laws at ends
  return gsl_spline_eval(Pspl,ki,acc);
}

/*
  Compute the power law behavior of the power spectrum
  at its ends.

  Args:
      k: wavenumber - units in either h/Mpc
      P: power spectrum - units in either (h/Mpc)^3
      Nk: number of k and P points
      alpha: power law indices - first element is low-k, second element is high-k (output)
      A: power law amplitude - first element is low-k, second element is high-k (output)
 */
int calc_power_law(double*k,double*P,int Nk,double*alpha,double*A){
  alpha[0] = log(P[1]/P[0])/log(k[1]/k[0]);
  A[0] = P[0]/pow(k[0],alpha[0]);
  alpha[1] = log(P[Nk-1]/P[Nk-2])/log(k[Nk-1]/k[Nk-2]);
  A[1] = P[Nk-1]/pow(k[Nk-1],alpha[1]);
  return 0;
}

/*Integrand of w(theta).

  Note: the arguments have to be passed through the w_params struct.

  Args:
      ki: wavenumber of the power spectrum; h/Mpc
      theta: angle on the sky
      k: array of wavenumbers; h/Mpc
      P: matter power spectrum; (h/Mpc)^3
      Nk: length of k and P arrays
      alpha: P power law indices
      A: P power law amplitudes
      Pspl: spline for P
      acc: accelerator for the spline
 */
typedef struct w_params{
  double theta;
  double*k;
  double*P;
  int Nk;
  double*alpha;
  double*A;
  gsl_spline*Pspl;
  gsl_interp_accel*acc;
}w_params;
double w_theta_integrand(double ki, void*params){
  w_params pars = *(w_params*)params;
  double theta = pars.theta;
  double*k = pars.k;
  double*P = pars.P;
  int Nk = pars.Nk;
  double*alpha = pars.alpha;
  double*A = pars.A;
  gsl_spline*Pspl = pars.Pspl;
  gsl_interp_accel*acc = pars.acc;
  return ki*get_P(ki, k, P, Nk, alpha, A, Pspl, acc)*F(theta,ki);
}

/*Angular correlation function, w(theta) at all angles.

  Args:
      k: array of wavenumbers; h/Mpc
      P: matter power spectrum; (h/Mpc)^3
      Nk: length of k and P arrays
      theta: angles on the sky; degrees
      w; angular correlation function
      Nt; length of theta and w arrays
 */
int ang_corr(double*k, double*P, int Nk, double*theta, double*w, int Nt){
  int i, status=0;
  gsl_integration_workspace * ws
    = gsl_integration_workspace_alloc(1000);

  gsl_spline*Pspl = gsl_spline_alloc(gsl_interp_cspline,Nk);
  gsl_spline_init(Pspl,k,P,Nk);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  double*alpha=(double*)malloc(2*sizeof(double));
  double*A=(double*)malloc(2*sizeof(double));
  status |= calc_power_law(k,P,Nk,alpha,A);
  
  w_params * params = malloc(sizeof(w_params));
  params->k = k;
  params->P = P;
  params->Nk = Nk;
  params->alpha = alpha;
  params->A = A;
  params->Pspl = Pspl;
  params->acc = acc;

  gsl_function F;
  F.function = &w_theta_integrand;
  F.params = &params;
  double result,error;

  for(i=0; i<Nt; i++){
    result = 0;
    error = 0;
    params->theta = theta[i];
    gsl_integration_qagui(&F, 0, 1e-7, 1e-7, 1000, ws, &result, &error);
    w[i] = result*9./(2*PI);
  }

  free(alpha),free(A);
}
