#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#ifdef _MSC_VER
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// Interface to a general H(B) relation with a variable number of material
//  model parameters and a variable number of states.
EXPORT int eval(double *oldB,         // Magnetic field, previous converged step, 3-vector, input
                double *B,            // Magnetic field, 3-vector, input
                double *H,            // Magnetic flux density, 3-vector, state
				double *Jac,          // Jacobian of B with respect to H, 3-by-3 matrix in row-major order, output
		        int *nPar,            // Number of material model parameters, scalar, input
                double *par,          // Material model parameters, nPar-vector, input
		        int *nStates,         // Number of states, scalar, input
		        double *states) {     // States, nStates-vector
  const double eps = 3e-16;
  const double u0 = 3.141692*4e-7;    // Permeability of vacuum
  double normH2, normH, normEpsH2, normEpsH, normB, dNormB, s;
  double Br0, Hc, jHc, ur, tol;                           // Input parameters
  double oldBy, oldHy, oldBr;                             // State variables
  double Br, Bx, By, Bz, xi, Bknee, Hknee, k1, k2, Br2;      // Calculation variables
  double Hy;
  
  // Check inputs
  if (nPar[0]!=5)
    return 1;	              // error code 1 = "Wrong number of parameters"
  if (nStates[0]!=3)
	return 2;                 // error code 2 = "Wrong number of states"

  // Read input parameters from parameters vector, call convention:
  // { Br, mur, jHc,  }
  Br0 = par[0];
  Hc = par[1];
  jHc = par[2];
  ur = par[3];
  tol = par[4];

  // Read state variables  
  // {oldBr, oldBy, oldHy}
  if ( states[0]==0 && states[1]==0 && states[2]==0 ) oldBr=Br0; 
  else oldBr = states[0];
  // By
  oldBy = states[1];
  // Hy
  oldHy = states[2];

  // Read input
  Bx = B[0];
  By = B[1];
  Bz = B[2];

  // Calculation of the BH-curve
  k1 = u0*ur;
  k2 = (0-(k1*jHc))/(Hc-jHc);
  Br2 = -k2 * Hc;
  Hknee = (Br2 - oldBr) / (k1-k2);
  Bknee = k1*Hknee + oldBr;

  // H(B), depending on if By is above or below the knee point. 
  if (By>Bknee) Hy = (By-oldBr)/k1;
  else Hy = (By-Br2)/k2;

  // New Br?
  Br = MIN(oldBr, By-(u0*ur*Hy));
  // Br = oldBr;
  
  // Hy = (B[1]/u0 - oldMy)/ur;

  // // Calculate dB and M
  // Ms = Br0/u0;
  // xi = ur - 1;
  // dH = Hy - oldHy;

  // if ( Hy<oldHknee ) {
  //   // Lowering the magnetization
  //   // M = oldMy + xi * dH; // Old version
  //   M = oldMy + 0.75*(Hy-oldHknee); // New version - steeper demagnetization

  //   Hknee = Hy;
  // } else {
  //   M = oldMy;
  //   Hknee = oldHknee;
  // }
    if ( abs(Br)<Br && abs(oldBr)<(Br-tol) ) {
    // B[1] = u0 * (H[1] + M);
    // H[1] = B[1]/u0 - M;
    Jac[0] = 1/(u0*ur);
    Jac[1] = 0;
    Jac[2] = 0;
    Jac[3] = 0;
    Jac[4] = 1/(u0*ur);
    Jac[5] = 0;
    Jac[6] = 0;
    Jac[7] = 0;
    Jac[8] = 0;
  } else {
    if ( abs(Br)>0 ) states[0] = Br / abs(Br) * Br;
    else states[0] = 0.0;
    Jac[0] = 1/(u0*ur) + B[0]*B[0]/(Br*Br);
    Jac[1] = + B[0]*B[1]/(Br*Br);
    Jac[2] = + B[0]*B[2]/(Br*Br);
    Jac[3] = + B[1]*B[0]/(Br*Br);
    Jac[4] = 1/(u0*ur) + B[1]*B[1]/(Br*Br);
    Jac[5] = + B[1]*B[2]/(Br*Br);
    Jac[6] = + B[2]*B[0]/(Br*Br);
    Jac[7] = + B[2]*B[1]/(Br*Br);
    Jac[8] = 1/(u0*ur) + B[2]*B[2]/(Br*Br);
  }

  H[0] = B[0]/(u0 * ur);
  H[1] = (B[1]-Br)/(u0*ur);
  H[2] = B[2]/(u0 * ur);

  // Write new states variables
  states[0] = Br;
  states[1] = B[1]; 
  states[2] = H[1];

  return 0;  // Return value 0 is success, any other value triggers an exception
}