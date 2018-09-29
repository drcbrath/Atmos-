// COESA (aka ISA) atmosphere model; standard and otherwise

#include <limits>
#include <math.h>
#include <vector>
#include "Atmos.h"

// Atm functions
// For simple use of the atmosphere model
// These functions use the standard lapse rates and optionally
// a constant temperature deviation, dT, applied to the profile.

// Atmos class
// This class provides more consistency with the math model; i.e. less able to 
// arbtrarily set temperatures and violate the physics and math underlying the model
// This class constructs the full temperature profile at construction/instantiation 
// of an Atmos object variable which then stays with that object. This costs some 
// extra computation at construction time, but less at use time; soa simulation
// that might make thousands of calls to evaluate the atmosphere properties 
// should see an overall reduction of computation time.

// helpers

// layer local pressure ratio (from bottom of layer at k up to altitude h)
inline double pr(double h, double Hk, double T, double Tk, double Tgradk, double GMR)
{
	if (Tgradk != 0)
		return(pow((T / Tk), (-GMR / Tgradk)));
	else
		return(exp(-GMR*(h-Hk)/Tk));
}

// layer dh from pressure ratio
inline double dh_pr(double PR, double PRk, double T, double Tk, double Tgradk, double GMR)
{
   double dh;
   // local, layer pressure ratio from pressure ratio w.r.t. datum at h & bottom of layer k
   double pr = PR / PRk;

   if (Tgradk != 0)
      dh = (pow(pr, (-Tgradk / GMR))-1.0)*Tk/Tgradk;
   else
      dh = -(Tk / GMR)*log(pr);

   return(dh);
}

// layer dh from density ratio
inline double dh_sgm(double sgm, double sgmk, double T, double Tk, double Tgradk, double GMR)
{
   double dh;

   // local, layer density ratio from density ratio w.r.t. datum at h & bottom of layer k
   double sr = sgm / sgmk;

   if (Tgradk != 0)
      dh = (Tk/Tgradk)*( pow(sr,(1/(1+GMR/Tgradk)))-1.0);
   else
      dh = -(Tk / GMR)*log(sr);

   return(dh);
}

//------- AtmProfile class, to hold temperature vs altitude profile -------

AtmProfile::AtmProfile() {};

AtmProfile::AtmProfile(std::vector<double> H, std::vector<double> T)
{
	nLayers = H.size() - 1;
	Hk = H;
	Tk = T;
   PRk.resize(H.size());
   DRk.resize(H.size());
   Tgradk.resize(nLayers);

	// should test gradient, if below a threshold, then truncate to zero to guard 
	// against round-off error inappropriately generating linear thermal layers

	PRk[0] = 1.0;
	DRk[0] = 1.0;
	for (int k = 0; k < nLayers; k++)
	{
      Tgradk[k] = (Tk[k + 1] - Tk[k]) / (Hk[k + 1] - Hk[k]);
      PRk[k + 1] = PRk[k] * pr(Hk[k + 1], Hk[k], Tk[k + 1], Tk[k], Tgradk[k], GMRearth);
		DRk[k + 1] = PRk[k + 1] * T0/Tk[k + 1];
	}

};

AtmProfile::AtmProfile(std::vector<double> H, std::vector<double> Tgrad, double T_1st)
{
	// assemble profile using gradient, starting from first point in H,Tgrad profile
   nLayers = H.size() - 1;
   Hk = H;
	Tgradk = Tgrad;
   Tk.resize(H.size());
   PRk.resize(H.size());
   DRk.resize(H.size());
 
	Tk[0] = T_1st;
   PRk[0] = 1.0;
   DRk[0] = 1.0;
   for (int k = 0; k < nLayers; k++)
   {
      Tk[k + 1] = Tk[k] + Tgradk[k] * (Hk[k + 1] - Hk[k]);
      PRk[k + 1] = PRk[k] * pr(Hk[k + 1], Hk[k], Tk[k + 1], Tk[k], Tgradk[k], GMRearth);
      DRk[k + 1] = PRk[k + 1] * T0 / Tk[k + 1];
   }

};

AtmProfile::AtmProfile(std::vector<double> H, std::vector<double> Tgrad, double Hic, double Tic)
{
	// assemble profile using gradient
	// find Hic in H, then build up & down from there
};

//------- Atm functions -------

int AtmSIRatios(double Hgp, double &theta, double &sigma, double &delta, double &kappa)
{
   // find layer n
   int n; for (n = 1; n <= StdDayProfile.nLayers && StdDayProfile.Hk[n] < Hgp; n++); n--;
   if (n >= StdDayProfile.nLayers) n = StdDayProfile.nLayers - 1;

   // compute properties from bottom of layer n up to given geopotential altitude hgp
   double T = StdDayProfile.Tk[n] + StdDayProfile.Tgradk[n] * (Hgp - StdDayProfile.Hk[n]);
   theta = T / T0;
   delta = StdDayProfile.PRk[n] * pr(Hgp, StdDayProfile.Hk[n], T, StdDayProfile.Tk[n], StdDayProfile.Tgradk[n], GMRearth);
   sigma = delta / theta;
   kappa = sqrt(theta);

   return(0);
};

int AtmSIRatiosdT(double Hgp, double &theta, double &sigma, double &delta, double &kappa, double dT)
{
   // find layer n
   int n; for (n = 1; n <= StdDayProfile.nLayers && StdDayProfile.Hk[n] < Hgp; n++); n--;
   if (n >= StdDayProfile.nLayers) n = StdDayProfile.nLayers - 1;

   // compute profile properies at breakpoints up to n, including dT effect
   std::vector<double> PRk(n + 1);  // allocate vector for pressure ratio breakpoints in dT modified profile to be made below
   PRk[0] = 1.0;
   for (int k = 0; k < n; k++)
      PRk[k + 1] = PRk[k] * pr(StdDayProfile.Hk[k + 1], StdDayProfile.Hk[k], StdDayProfile.Tk[k + 1] + dT, StdDayProfile.Tk[k] + dT, StdDayProfile.Tgradk[k], GMRearth);

   // compute properties from bottom of layer n up to given geopotential altitude hgp
   double T = StdDayProfile.Tk[n] + StdDayProfile.Tgradk[n] * (Hgp - StdDayProfile.Hk[n]) + dT;
   theta = T / T0;
   delta = PRk[n] * pr(Hgp, StdDayProfile.Hk[n], T, StdDayProfile.Tk[n] + dT, StdDayProfile.Tgradk[n], GMRearth);
   sigma = delta / theta;
   kappa = sqrt(theta);

   return(0);
};

int AtmSIRatiosT(double Hgp, double &theta, double &sigma, double &delta, double &kappa, double T)
{
   // compute ambient T from std day
   double theta_std;
   AtmSIRatios(Hgp, theta_std, sigma, delta, kappa);

   // compute dT
   double dT = T - theta_std*T0;

   // now compute desired results
   AtmSIRatiosdT(Hgp, theta, sigma, delta, kappa, dT);

   return(0);
};

int AtmSIRatiosCP(double Hgp, double &theta, double &sigma, double &delta, double &kappa, AtmProfile atm)
{
   // find layer n
   int n; for (n = 1; n <= atm.nLayers && atm.Hk[n] < Hgp; n++); n--;
   if (n >= atm.nLayers) n = atm.nLayers - 1;

   // compute properties from bottom of layer n up to given geopotential altitude hgp
   double T = atm.Tk[n] + atm.Tgradk[n] * (Hgp - atm.Hk[n]);
   theta = T / T0;
   delta = atm.PRk[n] * pr(Hgp, atm.Hk[n], T, atm.Tk[n], atm.Tgradk[n], GMRearth);
   sigma = delta / theta;
   kappa = sqrt(theta);

   return(0);
};


int AtmSI(double Hgp, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa)
{
   AtmSIRatios(Hgp, theta, sigma, delta, kappa);
   T = theta*T0;
   rho = sigma*rho0;
   P = delta*P0;
   a = kappa*a0;

   return(0);
};

int AtmSIdT(double Hgp, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT)
{
   AtmSIRatiosdT(Hgp, theta, sigma, delta, kappa, dT);
   T = theta*T0;
   rho = sigma*rho0;
   P = delta*P0;
   a = kappa*a0;

   return(0);
};

int AtmSIT(double Hgp, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT)
{
   // compute ambient T from std day
   double T_std;
   AtmSI(Hgp, T_std, rho, P, a, theta, sigma, delta, kappa);

   // compute dT
   double dT = T - T_std;

   // now compute desired results
   AtmSIdT(Hgp, T, rho, P, a, theta, sigma, delta, kappa,dT);

   return(0);
};

int AtmSICP(double Hgp, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, AtmProfile atm)
{
   AtmSIRatiosCP(Hgp, theta, sigma, delta, kappa, atm);
   T = theta*T0;
   rho = sigma*rho0;
   P = delta*P0;
   a = kappa*a0;

   return(0);
};


inline int AtmSIHgpdT(double &Hgm, double Hgp, double &Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT)
{
	// compute properties & ratios at (Hgp,dT)
		AtmSI(Hgp, dT, T, rho, P, a, theta, sigma, delta, kappa);

	// Hgm, geometric altitude
		Hgm = Hgp*Rearth / (Rearth - Hgp);

	// Hpa, pressure altitude
		// find layer n for delta in std day profile
      int n; for (n = 1; n <= StdDayProfile.nLayers && delta < StdDayProfile.PRk[n]; n++);
      n--;   // decrement so n points to bottom of layer
      if (n >= StdDayProfile.nLayers) n = StdDayProfile.nLayers - 1;   // max n is nLayers, if all loop tests fail, point is above table so must limit n here

		// from delta, compute dh in layer std day profile, and Hpa
		Hpa = StdDayProfile.Hk[n] + dh_pr(delta, StdDayProfile.PRk[n], T, StdDayProfile.Tk[n], StdDayProfile.Tgradk[n], GMRearth);

	// Hda, density altitude
		// find layer n for sigma in std day profile
      for (n = 1; n < StdDayProfile.DRk.size() && sigma < StdDayProfile.DRk[n]; n++);
      n--;   // decrement so n points to bottom of layer
      if (n >= StdDayProfile.nLayers) n = StdDayProfile.nLayers - 1;

		// from sigma, compute dh in layer std day profile, and Hda
		Hda = StdDayProfile.Hk[n] + dh_sgm(sigma, StdDayProfile.DRk[n], T, StdDayProfile.Tk[n], StdDayProfile.Tgradk[n], GMRearth);
	
	return(0);
};

inline int AtmSIHgpT(double &Hgm, double Hgp, double &Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT)
{
	double Tstd;

	// compute ambient T from std day
	AtmSI(Hgp, 0.0, Tstd, rho, P, a, theta, sigma, delta, kappa);

	// compute dT
	dT = T - Tstd;

	// now compute desired results
	AtmSIHgpdT(Hgm, Hgp, Hpa, Hda, dT, T, rho, P, a, theta, sigma, delta, kappa);

	return(0);
};

inline int AtmSIHgmdT(double Hgm, double &Hgp, double &Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT)
{
   // Hgp, geopotential altitude
   Hgp = Hgm*Rearth / (Rearth + Hgm);

   // compute properties & ratios at (Hgp,dT), Hpa, Hda
   AtmSIHgpdT(Hgm, Hgp, Hpa, Hda, dT, T, rho, P, a, theta, sigma, delta, kappa);

   return(0);
};

inline int AtmSIHgmT(double Hgm, double &Hgp, double &Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT)
{
   double Tstd;

   Hgp = Hgm*Rearth / (Rearth + Hgm);

   // compute ambient T for std day
   AtmSI(Hgp, 0.0, Tstd, rho, P, a, theta, sigma, delta, kappa);

   // compute dT
   dT = T - Tstd;

   // now compute desired results
   AtmSIHgpdT(Hgm, Hgp, Hpa, Hda, dT, T, rho, P, a, theta, sigma, delta, kappa);

   return(0);
};

inline int AtmSIHpadT(double &Hgm, double &Hgp, double Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT)
{
   // default, use StdDayProfile as base for dT mods

   // given pressure altitude, Hpa, & temperature deviation, dT
   // find geopotential, geometric, density altitudes, and atmosphere properties
   // compute properties & ratios at (Hgp,dT)

   // find pressure ratio, delta, at Hpa
   AtmSIRatios(Hpa, 0.0, theta, sigma, delta, kappa);   // dT=0.0 for std day per pressure altitude definiton

   // find layer, n, in dT modified profile
   std::vector<double> Tk(StdDayProfile.nLayers + 1);   // allocate vector for temperature breakpoints in dT modified profile to be made below
   std::vector<double> PRk(StdDayProfile.nLayers + 1);  // allocate vector for pressure ratio breakpoints in dT modified profile to be made below
   int n;
   Tk[0] = T0 + dT;
   PRk[0] = 1.0;
   for (n = 0; n < StdDayProfile.nLayers && delta < PRk[n]; n++)
   {
      Tk[n + 1] = Tk[n] + StdDayProfile.Tgradk[n] * (StdDayProfile.Hk[n + 1] - StdDayProfile.Hk[n]);
      PRk[n + 1] = PRk[n] * pr(StdDayProfile.Hk[n + 1], StdDayProfile.Hk[n], Tk[n + 1], Tk[n], StdDayProfile.Tgradk[n], GMRearth);
   }
   if(n != 0) n--;   // decrement so n points to bottom of layer
   if (n >= StdDayProfile.nLayers) n = StdDayProfile.nLayers - 1;   // max n is nLayers, if all loop tests fail, point is above table so must limit n here

   // from delta find dh to Hgp in layer n of dT modified profile, thence Hgp
   Hgp = StdDayProfile.Hk[n] + dh_pr(delta, PRk[n], T, Tk[n], StdDayProfile.Tgradk[n], GMRearth);

   // now with (Hgp,dT), compute other altitudes and properties

   // properties
   AtmSI(Hgp, dT, T, rho, P, a, theta, sigma, delta, kappa);

   // Hgm, geometric altitude
   Hgm = Hgp*Rearth / (Rearth - Hgp);

   // Hda, density altitude
   // find layer n for sigma in std day profile
   for (n = 1; n < StdDayProfile.DRk.size() && sigma < StdDayProfile.DRk[n]; n++);
   n--;   // decrement so n points to bottom of layer
   if (n >= StdDayProfile.nLayers) n = StdDayProfile.nLayers - 1;   // max n is nLayers, if all loop tests fail, point is above table so must limit n here

   // from sigma, compute dh in layer std day profile, and Hda
   Hda = StdDayProfile.Hk[n] + dh_sgm(sigma, StdDayProfile.DRk[n], T, StdDayProfile.Tk[n], StdDayProfile.Tgradk[n], GMRearth);
   
   return(0);
};

// !#!#!#!#!
inline int AtmSIHpaT(double &Hgm, double &Hgp, double Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT)
{
   double Tstd;

   // input T is ambient temperature at current Hgp (i.e. Hgp corresponding to input Hpa and T)
   // Tstd is T at Hgp in std day atm
   // dT is T - Tstd
   
   // consider scenario: flight test measured pressure altitude (Hpa), hence measured pressure P, and measured OAT (T)
   // this occurs at some pareticular, but as yet unknown, Hgp and temperature deviation from std day atm
   // now, the question is: what (Hgp,dT) generates matching pressure and temperature (P,T) to those measured?
   // once answered, this (Hgp,dT) can be used to find the desired results

   // temperature cannot uniquely determine an altitude, not even temperature altitude, because of isothermal layers where T = constant over a significant altitude range
   // yet dT must be uniquely determined from inputs---it must fall to Hpa to assist in finding Tstd hence dT

   // in actually, temperature affects pressure altitude, so Hpa = Hpa(T) --> input is AtmSIHpaT( ..., Hpa(T), T, ... ), hence it is no correct to determine Tstd from Hpa alone
   // dT must be determined using both Hpa and T, simultaneously
   // i.e. must solve simultaneously, (Hpa,T) --> (Hgp,dT)

   // given: Hpa,T, std profile Hk, Tk, Tgradk, T0, GMR
   // Hpa --> delta
   // known: delta,T
   // unknown: Hgp,dT
   // known or unknown? n
   // equations:
   // T = Tk[n] + Tgradk[n]*(Hgp-Hk[n]+dT, (eqn 1) temperature-altitude relation
   // delta = (T/T0)^(-GMR/Tgradk[n]),     (eqn 2) pressure-temperature relation

   // estimate n
   // compute Hgp,dT
   // compute n at Hgp,dT, check vs estimate
   // revise n if needed (should be relatively rarely, will only differ if dT sufficient move Hpa across Hk[n] breakpoint)


	return(0);
};
// !#!#!#!#!

inline int AtmSIHdadT(double &Hgm, double &Hgp, double &Hpa, double Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT)
{
   //// default, use StdDayProfile as base for dT mods
   //const AtmProfile &Profile = StdDayProfile;
   //AtmProfile prof = StdDayProfile;

   // given pressure altitude, Hpa, & temperature deviation, dT
   // find geopotential, geometric, density altitudes, and atmosphere properties
   // compute properties & ratios at (Hgp,dT)

   // find density ratio, sigma, at Hda
   AtmSIRatios(Hda, 0.0, theta, sigma, delta, kappa);   // dT=0.0 for std day per density altitude definiton

   // find layer, n, in dT modified profile
   std::vector<double> Tk(StdDayProfile.nLayers + 1);   // allocate vector for temperature breakpoints in dT modified profile to be made below
   std::vector<double> PRk(StdDayProfile.nLayers + 1);  // allocate vector for pressure ratio breakpoints in dT modified profile to be made below
   std::vector<double> DRk(StdDayProfile.nLayers + 1);  // allocate vector for density ratio breakpoints in dT modified profile to be made below
   int n;
   Tk[0] = T0 + dT;
   PRk[0] = 1.0;
   DRk[0] = 1.0;
   for (n = 0; n < StdDayProfile.nLayers && sigma < DRk[n]; n++)
   {
      Tk[n + 1] = Tk[n] + StdDayProfile.Tgradk[n] * (StdDayProfile.Hk[n + 1] - StdDayProfile.Hk[n]);
      PRk[n + 1] = PRk[n] * pr(StdDayProfile.Hk[n + 1], StdDayProfile.Hk[n], Tk[n + 1], Tk[n], StdDayProfile.Tgradk[n], GMRearth);
      DRk[n + 1] = PRk[n + 1] / (Tk[n + 1] / T0);
   }
   if (n != 0) n--;   // decrement so n points to bottom of layer
   if (n >= StdDayProfile.nLayers) n = StdDayProfile.nLayers - 1;   // max n is nLayers, if all loop tests fail, point is above table so must limit n here

   //double Tkn, Tknm1 = T0, PRkn = 1.0, DRkn = 1.0;
   //for (n = 1; n < StdDayProfile.nLayers; n++)
   //{
   //   Tkn = Tknm1 + StdDayProfile.Tgradk[n-1] * (StdDayProfile.Hk[n + 1] - StdDayProfile.Hk[n-1]);
   //   PRkn = PRkn * pr(StdDayProfile.Hk[n], StdDayProfile.Hk[n-1], Tkn, Tknm1, StdDayProfile.Tgradk[n-1], GMRearth);
   //   DRkn = PRkn / (Tkn / Tknm1);
   //   if ( sigma < DRkn) break;
   //}
   //n--;   // decrement so n points to bottom of layer
   //if (n >= StdDayProfile.nLayers) n = StdDayProfile.nLayers - 1;


   // from delta find dh to Hgp in layer n of dT modified profile, thence Hgp
   Hgp = StdDayProfile.Hk[n] + dh_sgm(sigma, DRk[n], T, Tk[n], StdDayProfile.Tgradk[n], GMRearth);

   // now with (Hgp,dT), compute other altitudes and properties

   // properties
   AtmSI(Hgp, dT, T, rho, P, a, theta, sigma, delta, kappa);

   // Hgm, geometric altitude
   Hgm = Hgp*Rearth / (Rearth - Hgp);

   // Hpa, pressure altitude
   // find layer n for delta in std day profile
   for (n = 1; n < StdDayProfile.PRk.size() && sigma < StdDayProfile.PRk[n]; n++);
   n--;   // decrement so n points to bottom of layer
   if (n >= StdDayProfile.nLayers) n = StdDayProfile.nLayers - 1;   // max n is nLayers, if all loop tests fail, point is above table so must limit n here

   // from delta, compute dh in layer std day profile, and Hda
   Hpa = StdDayProfile.Hk[n] + dh_pr(delta, StdDayProfile.PRk[n], T, StdDayProfile.Tk[n], StdDayProfile.Tgradk[n], GMRearth);

	return(0);
};

// !#!#!#!#!
inline int AtmSIHdaT(double &Hgm, double &Hgp, double &Hpa, double Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT)
{
   double Tstd;
   return(0);
};
// !#!#!#!#!

int AtmSI_HgmdT(double Hgm, double &Hgp, double &Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT)
{
   return(AtmSIHgmdT(Hgm, Hgp, Hpa, Hda, dT, T, rho, P, a, theta, sigma, delta, kappa));
};

int AtmSI_HgmT(double Hgm, double &Hgp, double &Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT)
{
   return(AtmSIHgmT(Hgm, Hgp, Hpa, Hda, dT, T, rho, P, a, theta, sigma, delta, kappa));
};

int AtmSI_HgpdT(double &Hgm, double Hgp, double &Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT)
{
	return(AtmSIHgpdT(Hgm, Hgp, Hpa, Hda, dT, T, rho, P, a, theta, sigma, delta, kappa));
};

int AtmSI_HgpT( double &Hgm, double Hgp, double &Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT)
{
	return(AtmSIHgpT(Hgm, Hgp, Hpa, Hda, dT, T, rho, P, a, theta, sigma, delta, kappa));
};

int AtmSI_HpadT(double &Hgm, double &Hgp, double Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT)
{
	return(AtmSIHpadT(Hgm, Hgp, Hpa, Hda, dT, T, rho, P, a, theta, sigma, delta, kappa));
};

int AtmSI_HpaT( double &Hgm, double &Hgp, double Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT)
{
	return(AtmSIHpaT(Hgm, Hgp, Hpa, Hda, dT, T, rho, P, a, theta, sigma, delta, kappa));
};

int AtmSI_HdadT(double &Hgm, double &Hgp, double &Hpa, double Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT)
{
   return(AtmSIHdadT(Hgm, Hgp, Hpa, Hda, dT, T, rho, P, a, theta, sigma, delta, kappa));
};

int AtmSI_HdaT( double &Hgm, double &Hgp, double &Hpa, double Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT)
{
	return(AtmSIHdaT(Hgm, Hgp, Hpa, Hda, dT, T, rho, P, a, theta, sigma, delta, kappa));
};


Atm AtmSI_Hgp(double Hgp)
{
   Atm atm;
   atm.Hgp = Hgp;
   atm.dT = 0.0;
   AtmSIHgpdT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);
   return(atm);
};
Atm AtmSI_HgpdT(double Hgp, double dT)
{
	Atm atm;
	atm.Hgp = Hgp;
	atm.dT = dT;
	AtmSIHgpdT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);
	return(atm);
};
Atm AtmSI_HgpT(double Hgp, double T)
{
	Atm atm;
	atm.Hgp = Hgp;
   atm.T = T;
   AtmSIHgpT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);

   return(atm);
};

Atm AtmSI_Hgm(double Hgm)
{
   Atm atm;
   atm.Hgm = Hgm;
   atm.dT = 0.0;
   AtmSIHgmdT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);
   return(atm);
};
Atm AtmSI_HgmdT(double Hgm, double dT)
{
   Atm atm;
   atm.Hgm = Hgm;
   atm.dT = dT;
   AtmSIHgmdT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);
   return(atm);
};
Atm AtmSI_HgmT(double Hgm, double T)
{
	Atm atm;
	atm.Hgm = Hgm;
   atm.T = T;
   AtmSIHgmT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);

   return(atm);
};

Atm AtmSI_Hpa(double Hpa)
{
   Atm atm;
   atm.Hpa = Hpa;
   atm.dT = 0.0;
   AtmSIHpadT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);
   return(atm);
};
Atm AtmSI_HpadT(double Hpa, double dT)
{
   Atm atm;
   atm.Hpa = Hpa;
   atm.dT = dT;
   AtmSIHpadT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);
   return(atm);
};
Atm AtmSI_HpaT(double Hpa, double T)
{
	Atm atm;
	atm.Hpa = Hpa;
   atm.T = T;
   AtmSIHpaT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);

   return(atm);
};

Atm AtmSI_Hda(double Hda)
{
   Atm atm;
   atm.Hda = Hda;
   atm.dT = 0.0;
   AtmSIHdadT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);
   return(atm);
};
Atm AtmSI_HdadT(double Hda, double dT)
{
   Atm atm;
   atm.Hda = Hda;
   atm.dT = dT;
   AtmSIHdadT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);
   return(atm);
};
Atm AtmSI_HdaT(double Hda, double T)
{
	Atm atm;
	atm.Hda = Hda;
   atm.T = T;
   AtmSIHdaT(atm.Hgm, atm.Hgp, atm.Hpa, atm.Hda, atm.dT, atm.T, atm.rho, atm.P, atm.a, atm.theta, atm.sigma, atm.delta, atm.kappa);

   return(atm);
};


//------- Atmos class implementation -------

AtmosSI::AtmosSI()   // construct standard atmosphere
{
};

AtmosSI::AtmosSI(double dT)   // construct standard + dT deviation added to temperature profile vs altitude
{
};

AtmosSI::AtmosSI(double Hic, double Tic)   // construct standard lapse rate profile used to construct temperature profile through (Hic,Tic)
{
};

AtmosSI::AtmosSI(std::vector<double> H, std::vector<double> T)   // construct custom temperature profile vs altitude
{
	Hk = H;
	Tk = T;
	nLayers = Hk.size() - 1;
};

AtmosSI(std::vector<double> Hprofile, std::vector<double> Tprofile)   // custom temperature profile vs altitude
{

};

AtmosSI(std::vector<double> Hprofile, std::vector<double> Tgradprofile, double T_1st)   // profile constructed from T gradient & T at 1st point
{

};

AtmosSI(std::vector<double> Hprofile, std::vector<double> Tgradprofile, double Hic, double Tic)   // profile constructed from T gradient &  specified initial condition (Hic,Tic)
{

};

AtmosSI::~AtmosSI(){};

// evaluate properties at gepotential altitude
// but not at other altitude definitions;  in fact, invalidate others to avoid 
// performance penalty and avoid incorrect/out-of-sync values
int AtmosSI::at(double hgp)
{
	return(0);
};

int AtmosSI::operator()(double hgp)
{
	return(at(hgp));
};

// evaluate atmosphere properties, including alternative altitude definitions, at 
int AtmosSI::atHgp(double hgp)   // given geoptential altitude
{
	return(0);
};
int AtmosSI::atHgm(double hgm)   // given geometric altitude
{
	return(0);
};
int AtmosSI::atHpa(double Hpa)   // given pressure altitude
{
	return(0);
};
int AtmosSI::atHda(double hda)   // given density altitude
{
	return(0);
};

// return altitude, evaluate if necessary first to provide valid result in case
// hgm, hpa, or hda does not already correspond to current hgp

double AtmosSI::Hgp(void)   // geoptential altitude
{
	return(hgp);
};

double AtmosSI::Hgm(void)	// geometric altitude
{
	if (isnan(hgm))   // then hpa is not valid, does not correspond to current hgp
	{
		atHgp(hgp);   // update all altitudes by evaluating properties including all altitudes
	}
	return(hgm);
};

double AtmosSI::Hpa(void)	// pressure altitude
{
	if (isnan(hpa))   // then hpa is not valid, does not correspond to current hgp
	{
		atHgp(hgp);   // update all altitudes by evaluating properties including all altitudes
	}
	return(hpa);
};

double AtmosSI::Hda(void)	// density altitude
{
	if (isnan(hda))   // then hpa is not valid, does not correspond to current hgp
	{
		atHgp(hgp);   // update all altitudes by evaluating properties including all altitudes
	}
	return(hda);
};
