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

// find n, layer number which contains hgp_in
inline int findLayer(int nLayers, double hgp_in, std::vector<double> Hk)
{
   int n; for (n = 1; n <= nLayers && Hk[n] < hgp_in; n++); n--;
   if (n >= nLayers) n = nLayers - 1;

   return(n);
}

// find n, layer number which contains pressure or density ratio
inline int findLayer_PD(int nLayers, double pd, std::vector<double> PDk)
{
   // find layer n for delta or sigma (pd) in profile PRk or DRk (PDk)
   int n;
   for (n = 1; n <= nLayers && pd < PDk[n]; n++);
   n--;   // decrement so n points to bottom of layer
   if (n > nLayers) n = nLayers;   // max n is nLayers, if all loop tests fail, point is above table so must limit n here

   return(n);
}
// layer local pressure ratio (from bottom of layer at k up to altitude h)
inline double pr(double h, double Hk, double T, double Tk, double Tgradk, double GMR)
{
   if (Tgradk != 0)
      return(pow((T / Tk), (-GMR / Tgradk)));
   else
      return(exp(-GMR*(h - Hk) / Tk));
}

// layer dh from pressure ratio
inline double dh_pr(double PR, double PRk, double T, double Tk, double Tgradk, double GMR)
{
   double dh;
   // local, layer pressure ratio from pressure ratio w.r.t. datum at h & bottom of layer k
   double pr = PR / PRk;

   if (Tgradk != 0)
      dh = (pow(pr, (-Tgradk / GMR)) - 1.0)*Tk / Tgradk;
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
      dh = (Tk / Tgradk)*(pow(sr, (1 / (1 + GMR / Tgradk))) - 1.0);
   else
      dh = -(Tk / GMR)*log(sr);

   return(dh);
}

// find pressure altitude from pressure ratio
inline double fHpa(double delta, std::vector<double> StdDayHk, std::vector<double> StdDayTk, std::vector<double> StdDayTgradk, std::vector<double> StdDayPRk, double T, double GMR)
{
   double hpa_out;

   int nLayers = StdDayHk.size()-1;

   int n = findLayer_PD(nLayers, delta, StdDayPRk);

   // from delta, compute dh in layer of std day profile, and Hpa
   hpa_out = StdDayHk[n] + dh_pr(delta, StdDayPRk[n], T, StdDayTk[n], StdDayTgradk[n], GMR);

   return(hpa_out);
}

// find density altitude from density ratio
inline double fHda(double sigma, std::vector<double> StdDayHk, std::vector<double> StdDayTk, std::vector<double> StdDayTgradk, std::vector<double> StdDayDRk, double T, double GMR)
{
   double hda_out;
   int n;
   int nLayers = StdDayHk.size();

   // find layer n containing sigma in std day profile
   int n = findLayer_PD(nLayers, sigma, StdDayDRk);

   // from sigma, compute dh in layer std day profile, and Hda
   hda_out = StdDayHk[n] + dh_sgm(sigma, StdDayDRk[n], T, StdDayTk[n], StdDayTgradk[n], GMR);

   return(hda_out);
}

inline void initializeProfile(double T0, double GMR, std::vector<double> Hk, std::vector<double> Tk, int &nLayers, std::vector<double> &Tgradk, std::vector<double> &PRk, std::vector<double> &DRk)
{
   nLayers = Hk.size() - 1;

   Tgradk.resize(nLayers);
   PRk.resize(nLayers + 1);
   DRk.resize(nLayers + 1);

   // should test gradient, if below a threshold, then truncate to zero to guard 
   // against round-off error inappropriately generating linear thermal layers

   PRk[0] = 1.0;
   DRk[0] = 1.0;
   for (int k = 0; k < nLayers; k++)
   {
      Tgradk[k] = (Tk[k + 1] - Tk[k]) / (Hk[k + 1] - Hk[k]);
      PRk[k + 1] = PRk[k] * pr(Hk[k + 1], Hk[k], Tk[k + 1], Tk[k], Tgradk[k], GMR);
      DRk[k + 1] = PRk[k + 1] * T0 / Tk[k + 1];
   }

}

inline void initializePRkDRk(double T0, double GMR, int n, double Hic, double Tic, double Pic, std::vector<double> Hk, std::vector<double> Tk, std::vector<double> Tgradk,
   std::vector<double> PRk, std::vector<double> DRk)
{
   int nLayers = Hk.size() - 1;

   double PicbyPRk = pr(Hic, Hk[n], Tic, Tk[n], Tgradk[n], GMR);
   PRk[n] = Pic / PicbyPRk;   // P at bottom of layer holding initial condition point

   int k;
   // for each 0 <= k < n, compute Tk & Pk below initial condition point
   for (k = n - 1; k >= 0; k--)
   {
      PRk[k] = PRk[k + 1] / pr(Hk[k + 1], Hk[n], Tk[k + 1], Tk[n], Tgradk[n], GMR);
   }

   // for each n < k < nLayers+1, compute Tk & Pk above initial condition point
   for (k = n; n <= nLayers; k++)
   {
      PRk[k + 1] = PRk[k] * pr(Hk[k + 1], Hk[n], Tk[k + 1], Tk[n], Tgradk[n], GMR);
   }

   // construct density ratio profile
   DRk[0] = 1.0;
   for (int k = 0; k <= nLayers; k++)
   {
      DRk[k + 1] = PRk[k + 1] * T0 / Tk[k + 1];
   }

}

//------- Atm functions -------

int AtmRatios(double hgp,  double T0, double GMR, std::vector<double> Hk, std::vector<double> Tk, std::vector<double> Tgradk, std::vector<double> PRk, std::vector<double> DRk,
   double &theta, double &sigma, double &delta, double &kappa)
{
   int nLayers = Hk.size() - 1;

   // find layer n
   int n = findLayer(nLayers, hgp, Hk);

   // compute properties from bottom of layer n up to given geopotential altitude hgp
   double T = Tk[n] + Tgradk[n] * (hgp - Hk[n]);
   theta = T / T0;
   delta = PRk[n] * pr(hgp, Hk[n], T, Tk[n], Tgradk[n], GMR);
   sigma = delta / theta;
   kappa = sqrt(theta);

   return(0);
};

//------- Atmos class implementation -------

// construct standard atmosphere, SI units
Atmos::Atmos(AtmosParameters AtmPrms = AtmosParameters_si)
{
   Re = AtmPrms.Re;               // planet radius
   GMR = AtmPrms.GMR;             // combined gravity and gas constant of planet atmosphere
   T0 = AtmPrms.T0;               // temperature
   rho0 = AtmPrms.rho0;           // density
   P0 = AtmPrms.P0;               // pressure
   a0 = AtmPrms.a0;               // sonic speed
   StdDayHk = AtmPrms.StdDayHk;
   StdDayTk = AtmPrms.StdDayTk;
   nStdLayers = StdDayHk.size() - 1;
   initializeProfile(T0, GMR, StdDayHk, StdDayTk, nStdLayers, StdDayTgradk, StdDayPRk, StdDayDRk);

   Hk = StdDayHk_si;
   Tk = StdDayTk_si;
   nLayers = Hk.size() - 1;

   Tgradk.resize(nLayers);
   PRk.resize(nLayers + 1);
   DRk.resize(nLayers + 1);

   initializeProfile(T0, GMR, Hk, Tk, nLayers, Tgradk, PRk, DRk);

   this->at(0.0);   // initialize object to datum (i.e. sea level)

};

// construct standard + dT deviation added to temperature profile vs altitude, SI units
Atmos::Atmos(double dT, AtmosParameters AtmPrms = AtmosParameters_si)
{
   Re = AtmPrms.Re;               // planet radius
   GMR = AtmPrms.GMR;             // combined gravity and gas constant of planet atmosphere
   T0 = AtmPrms.T0;               // temperature
   rho0 = AtmPrms.rho0;           // density
   P0 = AtmPrms.P0;               // pressure
   a0 = AtmPrms.a0;               // sonic speed
   StdDayHk = AtmPrms.StdDayHk;
   StdDayTk = AtmPrms.StdDayTk;
   nStdLayers = StdDayHk.size() - 1;
   initializeProfile(T0, GMR, StdDayHk, StdDayTk, nStdLayers, StdDayTgradk, StdDayPRk, StdDayDRk);

   Hk = StdDayHk_si;
   Tk = StdDayTk_si;
   nLayers = Hk.size() - 1;

   for (int m = 0; m < Hk.size; m++)
      Tk[m] += dT;

   initializeProfile(T0, GMR, Hk, Tk, nLayers, Tgradk, PRk, DRk);

   this->at(0.0);   // initialize object to datum (i.e. sea level)

};

// construct using standard lapse rate profile used to derive temperature profile through (Hic,Tic,Pic), SI units
Atmos::Atmos(double Hic, double Tic, double Pic, AtmosParameters AtmPrms = AtmosParameters_si)
{
   Re = AtmPrms.Re;               // planet radius
   GMR = AtmPrms.GMR;             // combined gravity and gas constant of planet atmosphere
   T0 = AtmPrms.T0;               // temperature
   rho0 = AtmPrms.rho0;           // density
   P0 = AtmPrms.P0;               // pressure
   a0 = AtmPrms.a0;               // sonic speed
   StdDayHk = AtmPrms.StdDayHk;
   StdDayTk = AtmPrms.StdDayTk;
   nStdLayers = StdDayHk.size() - 1;
   initializeProfile(T0, GMR, StdDayHk, StdDayTk, nStdLayers, StdDayTgradk, StdDayPRk, StdDayDRk);

   Hk = StdDayHk_si;
   Tgradk = StdDayTgradk_si;

   nLayers = Hk.size() - 1;
   Tk.resize(nLayers + 1);
   PRk.resize(nLayers + 1);
   DRk.resize(nLayers + 1);

   // construct temperature profile
   // find n such that Hk[n] < Hic < Hk[n + 1];
   int n;
   for (n = 1; n <= nLayers && StdDayHk_si[n] < Hic; n++); n--;
   if (n >= nLayers) n = nLayers - 1;   // ???

   Tk[n] = Tic - Tgradk[n] * (Hic - Hk[n]);   // T at bottom of layer holding initial condition point

   int k;
   // for each 0 <= k < n, compute Tk & Pk below initial condition point
   for (k = n - 1; k >= 0; k--)
   {
      Tk[k] = Tk[k + 1] - Tgradk[k] * (Hk[k + 1] - Hk[k]);
   }

   // for each n < k < nLayers+1, compute Tk & Pk above initial condition point
   for (k = n; n <= nLayers; k++)
   {
      Tk[k + 1] = Tk[k] + Tgradk[k] * (Hk[k + 1] - Hk[k]);
   }

   initializePRkDRk(T0, GMR, n, Hic, Tic, Pic, Hk, Tk, Tgradk, PRk, DRk);

   this->at(0.0);   // initialize object to datum (i.e. sea level)

};

// custom from initial conditions Hic, Pic, breakpoints, and temperature profile, units determined by input
Atmos::Atmos(double Hic, double Pic, std::vector<double> Hj, std::vector<double> Tj, AtmosParameters AtmPrms = AtmosParameters_si)
{
   Re = AtmPrms.Re;               // planet radius
   GMR = AtmPrms.GMR;             // combined gravity and gas constant of planet atmosphere
   T0 = AtmPrms.T0;               // temperature
   rho0 = AtmPrms.rho0;           // density
   P0 = AtmPrms.P0;               // pressure
   a0 = AtmPrms.a0;               // sonic speed
   StdDayHk = AtmPrms.StdDayHk;
   StdDayTk = AtmPrms.StdDayTk;
   nStdLayers = StdDayHk.size() - 1;
   initializeProfile(T0, GMR, StdDayHk, StdDayTk, nStdLayers, StdDayTgradk, StdDayPRk, StdDayDRk);

   Hk = Hj;
   Tk = Tj;

   nLayers = Hk.size() - 1;
   Tk.resize(nLayers + 1);
   PRk.resize(nLayers + 1);
   DRk.resize(nLayers + 1);
   
// construct temperature gradient profile
   for (int k = 0; k < nLayers; k++)
   {
      Tgradk[k] = (Tk[k + 1] - Tk[k]) / (Hk[k + 1] - Hk[k]);
   }

// construct pressure ratio profile
   // find n such that Hk[n] < Hic < Hk[n + 1];
   int n;
   for (n = 1; n <= nLayers && StdDayHk_si[n] < Hic; n++); n--;
   if (n >= nLayers) n = nLayers - 1;   // ???

   double Tic = Tk[n] + Tgradk[n] * (Hic - Hk[n]);   // temperature at initial condition point from given temperature profile

   initializePRkDRk(T0, GMR, n, Hic, Tic, Pic, Hk, Tk, Tgradk, PRk, DRk);

   this->at(0.0);   // initialize object to datum (i.e. sea level)

};

// custom from initial conditions (Hic, Tic, Pic), breakpoints, lapse rates, units determined by input
Atmos::Atmos(double Hic, double Tic, double Pic, std::vector<double> Hj, std::vector<double> Tgradj, AtmosParameters AtmPrms = AtmosParameters_si)
{
   Re = AtmPrms.Re;               // planet radius
   GMR = AtmPrms.GMR;             // combined gravity and gas constant of planet atmosphere
   T0 = AtmPrms.T0;               // temperature
   rho0 = AtmPrms.rho0;           // density
   P0 = AtmPrms.P0;               // pressure
   a0 = AtmPrms.a0;               // sonic speed
   StdDayHk = AtmPrms.StdDayHk;
   StdDayTk = AtmPrms.StdDayTk;
   nStdLayers = StdDayHk.size() - 1;
   initializeProfile(T0, GMR, StdDayHk, StdDayTk, nStdLayers, StdDayTgradk, StdDayPRk, StdDayDRk);

   Hk = Hj;
   Tgradk = Tgradj;

   nLayers = Hk.size() - 1;
   Tk.resize(nLayers + 1);
   PRk.resize(nLayers + 1);
   DRk.resize(nLayers + 1);

// construct temperature profile
  // find n such that Hk[n] < Hic < Hk[n + 1];
   int n;
   for (n = 1; n <= nLayers && StdDayHk_si[n] < Hic; n++); n--;
   if (n >= nLayers) n = nLayers - 1;   // ???

   Tk[n] = Tic - Tgradk[n] * (Hic - Hk[n]);   // T at bottom of layer holding initial condition point

   int k;
   // for each 0 <= k < n, compute Tk & Pk below initial condition point
   for (k = n - 1; k >= 0; k--)
   {
      Tk[k] = Tk[k + 1] - Tgradk[k] * (Hk[k + 1] - Hk[k]);
   }

   // for each n < k < nLayers+1, compute Tk & Pk above initial condition point
   for (k = n; n <= nLayers; k++)
   {
      Tk[k + 1] = Tk[k] + Tgradk[k] * (Hk[k + 1] - Hk[k]);
   }

   initializePRkDRk(T0, GMR, n, Hic, Tic, Pic, Hk, Tk, Tgradk, PRk, DRk);

   this->at(0.0);   // initialize object to datum (i.e. sea level)

};

// destructor
Atmos::~Atmos() {};

// evaluate properties at geopotential altitude
// but not at other altitude definitions; in fact, invalidate others to avoid 
// performance penalty and avoid incorrect/out-of-sync values

int Atmos::at(double hgp_in)
{

   hgp = hgp_in;

   int n = findLayer(nLayers, hgp, Hk);

   // compute properties from bottom of layer n up to given geopotential altitude hgp
   TT = Tk[n] + Tgradk[n] * (hgp_in - Hk[n]);
   Theta = TT / T0;
   Delta = PRk[n] * pr(hgp_in, Hk[n], TT, Tk[n], Tgradk[n], GMR);
   Sigma = Delta / Theta;
   Kappa = sqrt(Theta);

   hgm = hpa = hda = NAN;   // invalidate other altitudes

   return(0);
};

int Atmos::operator()(double hgp_in)
{
   return(at(hgp_in));
};

// evaluate atmosphere properties, including alternative altitude definitions
int Atmos::atHgp(double hgp_in)   // given geoptential altitude
{

   this->at(hgp_in);

   hgm = hgp_in * Re / (Re - hgp_in);                                         // geometric altitude
   hpa = fHpa(Delta, StdDayHk, StdDayTk, StdDayTgradk, StdDayPRk, TT, GMR);   // pressure altitude
   hda = fHda(Sigma, StdDayHk, StdDayTk, StdDayTgradk, StdDayPRk, TT, GMR);   // density altitude

   return(0);
};

int Atmos::atHgm(double hgm_in)   // given geometric altitude
{
   hgm = hgm_in;
   hgp = hgm*Re / (Re + hgm);                                                 // geopotential altitude

   this->at(hgp);

   hgm = hgm_in;                                                              // at() invalidates hgm, so it must be reset
   hpa = fHpa(Delta, StdDayHk, StdDayTk, StdDayTgradk, StdDayPRk, TT, GMR);   // pressure altitude
   hda = fHda(Sigma, StdDayHk, StdDayTk, StdDayTgradk, StdDayPRk, TT, GMR);   // density altitude

   return(0);
};

int Atmos::atHpa(double hpa_in)   // given pressure altitude
{
   //  find pressure ratio, delta, at Hpa
   double thta, sgma, dlta, kppa;

   AtmRatios(hpa_in, T0, GMR, StdDayHk, StdDayTk, StdDayTgradk, StdDayPRk, StdDayDRk, thta, sgma, dlta, kppa);

   // find layer, n, in this.profile containing pressure ratio, dlta
   int n = findLayer_PD(nLayers, dlta, PRk);

   // from delta, find dh to Hgp in layer n of this.profile, thence hgp
   hgp = StdDayHk[n] + dh_pr(dlta, PRk[n], thta*T0, Tk[n], StdDayTgradk[n], GMR);

   // now with hgp, compute other altitudes and properties

   this->atHgp(hgp);

   return(0);
};

int Atmos::atHda(double hda_in)   // given density altitude
{
   //  find pressure ratio, sigma, at Hda
   double thta, sgma, dlta, kppa;

   AtmRatios(hda_in, T0, GMR, StdDayHk, StdDayTk, StdDayTgradk, StdDayPRk, StdDayDRk, thta, sgma, dlta, kppa);

   // find layer, n, in this.profile containing density ratio, sgma
   int n = findLayer_PD(nLayers, sgma, DRk);

   // from sigma, find dh to Hgp in layer n of this.profile, thence hgp
   hgp = StdDayHk[n] + dh_sgm(sgma, DRk[n], thta*T0, Tk[n], StdDayTgradk[n], GMR);

   // now with hgp, compute other altitudes and properties

   this->atHgp(hgp);

   return(0);
};

// return altitude, evaluate if necessary first to provide valid result in case
// hgm, hpa, or hda does not already correspond to current hgp

double Atmos::Hgp(void)   // geoptential altitude
{
   return(hgp);
};

double Atmos::Hgm(void)	// geometric altitude
{
   if (isnan(hgm))   // then hpa is not valid, does not correspond to current hgp
   {
      atHgp(hgp);   // update all altitudes by evaluating properties including all altitudes
   }
   return(hgm);
};

double Atmos::Hpa(void)	// pressure altitude
{
   if (isnan(hpa))   // then hpa is not valid, does not correspond to current hgp
   {
      atHgp(hgp);   // update all altitudes by evaluating properties including all altitudes
   }
   return(hpa);
};

double Atmos::Hda(void)	// density altitude
{
   if (isnan(hda))   // then hpa is not valid, does not correspond to current hgp
   {
      atHgp(hgp);   // update all altitudes by evaluating properties including all altitudes
   }
   return(hda);
};

// return properties at current altitude

double Atmos::T(void)
{
   return(TT);
}

double Atmos::rho(void)
{
   return(Rho);
}
double Atmos::P(void)
{
   return(p);
}

double Atmos::a(void)
{
   return(Sonic);
}

double Atmos::visc(void)
{
   return(Visc);
}

double Atmos::theta(void)
{
   return(Theta);
}

double Atmos::sigma(void)
{
   return(Sigma);
}

double Atmos::delta(void)
{
   return(Delta);
}

double Atmos::kappa(void)
{
   return(Kappa);
}