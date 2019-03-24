// Atmosphere properties model; standard and otherwise

#include <limits>
#include <math.h>
#include <vector>
#include "Atmos.h"

// Atmos class
// This class constructs the full temperature profile at construction/instantiation 
// of an Atmos object variable which then stays with that object. This costs some 
// extra computation at construction time, but less at use time; so a simulation
// that might make thousands of calls to evaluate the atmosphere properties 
// should see an overall reduction of computation time.

//------- constants available externally -------

// to derive from fundamentals, use these formulas instead of constants below
//const double g0 = 9.90665;               // (m/s^2) standard referencence gravity
//const double Ru = ?;                     // () universal gas constant
//const double M = ?;                      // () dry air molal mass
//const double gma = 1.4;                  // dry air ratio of specific heats
//const double a0 = sqrt(gma*Ru / M*T0);   // (m/s), speed of sound at SL std temperature
//const double GMR = g0*M / Rearth;        // (degK/m) combined gravity and gas constant of dry air on earth

// SI units (default), use with atmos constructor to construct SI units based atmos object

const double Re_si = 6369000;                    // (m) radius of the earth
const double GMR_si = 0.034163195;               // (degK/m) combined gravity and gas constant of dry air on earth
const double H0_si = 0.0;                        // (m) datum, sea level
const double T0_si = 288.15;                     // (K) SL std temp
const double rho0_si = 1.225;                    // (kg/m^3) SL std density
const double P0_si = 101325;                     // (N/m^2) SL std pressure
const double a0_si = 340.2686;                   // (m/s) speed of sound at SL std temperature
const double visc0_si = 1.789380278077583e-05;   // (N/m^2) air dynamic viscosity at SL std
const double S_si = 110.4;                       // (K) air Sutherland temperature for viscosity computation

// defined atmosphere profiles <need to revise alternate day profiles! base on ref Mil 3013? references !>
const std::vector<double> StdDayHk_si({ 0.0,  11000.0,  20000.0,  32000.0,  47000.0,  51000.0,  71000.0,  84852.0 });
const std::vector<double> StdDayTk_si({ 288.15,   216.65,   216.65,   228.65,   270.65,   270.65,   214.65,   186.95 });
const std::vector<double> StdDayTgradk_si({ -0.0065, 0.0000, 0.0010, 0.0028, 0.0000, -0.0028, -0.0019997, 0.0 });

const std::vector<double> HotDayHk_si({ 0.0, 11000., 20000. });
const std::vector<double> HotDayTk_si({ 308.15,269.65,237.65 });

const std::vector<double> ColdDayHk_si({ 0.0, 11000., 20000. });
const std::vector<double> ColdDayTk_si({ 308.15,269.65,237.65 });

const std::vector<double> TropicalDayHk_si({ 0.0, 11000., 20000. });
const std::vector<double> TropicalDayTk_si({ 308.15,269.65,237.65 });

const std::vector<double> PolarDayHk_si({ 0.0, 11000., 20000. });
const std::vector<double> PolarDayTk_si({ 308.15,269.65,237.65 });

const AtmosParameters AtmosParameters_si = { Re_si, GMR_si, H0_si, T0_si, rho0_si, P0_si, a0_si, visc0_si, S_si, StdDayHk_si, StdDayTk_si };

// US units, use with atmos constructor to construct US units based atmos object
const double Re_us = Re_si / 0.3048;                                   // (m) radius of the earth
const double GMR_us = GMR_si * (1.8*0.3048);                           // (degR/ft) combined gravity and gas constant of dry air on earth
const double H0_us = 0.0;                                              // (ft) datum, sea level
const double T0_us = T0_si * 1.8;                                      // (R), SL std temp
const double rho0_us = rho0_si * (0.068521766*0.3048*0.3048*0.3048);   // (sl/ft^3), SL std density
const double P0_us = P0_si * (0.3048*0.3048/4.4482216152605);          // (lbf/ft^2), SL std pressure
const double a0_us = a0_si / 0.3048;                                   // (ft/s), speed of sound at SL std temperature
const double visc0_us = visc0_si * (0.3048*0.3048/4.4482216152605);    // (N/m^2) air dynamic viscosity at SL std
const double S_us = S_si*1.8;                                          // (R) air Sutherland temperature for viscosity computation

   // defined atmosphere profiles <need to revise alternate day profiles! base on ref Mil 3013? references !>
const std::vector<double> StdDayHk_us({ 0.000, 3352.800, 6096.000, 9753.600, 14325.600, 15544.800, 21640.800, 25862.890 });
const std::vector<double> StdDayTk_us({ 518.670, 389.970, 389.970, 411.570, 487.170, 487.170, 386.370, 336.510 });
const std::vector<double> StdDayTgradk_us({ 0 });   // definitely wrong

const std::vector<double> HotDayHk_us({ 0.000, 3352.800, 6096.000 });
const std::vector<double> HotDayTk_us({ 554.670, 485.370, 427.770 });

const std::vector<double> ColdDayHk_us({ 0.000, 3352.800, 6096.000 });
const std::vector<double> ColdDayTk_us({ 554.670, 485.370, 427.770 });

const std::vector<double> TropicalDayHk_us({ 0.000, 3352.800, 6096.000 });
const std::vector<double> TropicalDayTk_us({ 554.670, 485.370, 427.770 });

const std::vector<double> PolarDayHk_us({ 0.000, 3352.800, 6096.000 });
const std::vector<double> PolarDayTk_us({ 554.670, 485.370, 427.770 });

const AtmosParameters AtmosParameters_us = { Re_us, GMR_us, H0_us, T0_us, rho0_us, P0_us, a0_us, visc0_us, S_us, StdDayHk_us, StdDayTk_us };

// defined atmosphere models
Atmos StdDay_si(AtmosParameters_si);
Atmos HotDay_si(H0_si, P0_si, HotDayHk_si, HotDayTk_si, AtmosParameters_si);
Atmos ColdDay_si(H0_si, P0_si, ColdDayHk_si, ColdDayTk_si, AtmosParameters_si);
Atmos TropicalDay_si(H0_si, P0_si, TropicalDayHk_si, TropicalDayTk_si, AtmosParameters_si);
Atmos PolarDay_si(H0_si, P0_si, PolarDayHk_si, PolarDayTk_si, AtmosParameters_si);

Atmos StdDay_us(AtmosParameters_us);
Atmos HotDay_us(H0_us, P0_us, HotDayHk_us, HotDayTk_us, AtmosParameters_us);
Atmos ColdDay_us(H0_us, P0_us, ColdDayHk_us, ColdDayTk_us, AtmosParameters_us);
Atmos TropicalDay_us(H0_us, P0_us, TropicalDayHk_us, TropicalDayTk_us, AtmosParameters_us);
Atmos PolarDay_us(H0_us, P0_us, PolarDayHk_us, PolarDayTk_us, AtmosParameters_us);

//------- helper functions -------

// find n, layer number which contains hgp_in
inline int findLayer(int nLayers, double hgp_in, std::vector<double> Hk)
{
   int n; for (n = 1; n <= nLayers && Hk[n] < hgp_in; n++); n--;
   //if (n >= nLayers) n = nLayers - 1;

   return(n);
}

// find n, layer number which contains pressure or density ratio
inline int findLayer_PD(int nLayers, double pd, std::vector<double> PDk)
{
   // find layer n for delta or sigma (pd) in profile PRk or DRk (PDk)
   int n;
   for (n = 1; n <= nLayers && pd < PDk[n]; n++);
   n--;   // decrement so n points to bottom of layer
   if (n >= nLayers) n = nLayers-1;   // max n (index) is nLayers-1, if all loop tests fail, point is above table so must limit n here

   return(n);
}
// layer local pressure ratio (from bottom of layer at k up to altitude h)
inline double pr(double h, double Hk, double T, double Tk, double Tgradk, double GMR)
{
   if (Tgradk != 0)                             // linear thermal layer
      return(pow((T / Tk), (-GMR / Tgradk)));
   else                                         // isothermal layer
      return(exp(-GMR*(h - Hk) / Tk));
}

// layer dh from pressure ratio
inline double dh_pr(double PR, double PRk, double Tk, double Tgradk, double GMR)
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
inline double dh_sgm(double sgm, double sgmk, double Tk, double Tgradk, double GMR)
{
   double dh;

   if (Tgradk != 0)
      dh = Tk * (pow(sgmk/sgm, 1/(1+GMR/Tgradk)) - 1) / Tgradk;
   else
      dh = -(Tk / GMR)*log(sgm/sgmk);

   return(dh);
}

// find pressure altitude from pressure ratio
inline double fHpa(double delta, std::vector<double> StdDayHk, std::vector<double> StdDayTk, std::vector<double> StdDayTgradk, std::vector<double> StdDayPRk, double GMR)
{
   double hpa_out;

   int nLayers = StdDayHk.size()-1;

   int n = findLayer_PD(nLayers, delta, StdDayPRk);

   // from delta, compute dh in layer of std day profile, and Hpa
   hpa_out = StdDayHk[n] + dh_pr(delta, StdDayPRk[n], StdDayTk[n], StdDayTgradk[n], GMR);

   return(hpa_out);
}

// find density altitude from density ratio
inline double fHda(double sigma, std::vector<double> StdDayHk, std::vector<double> StdDayTk, std::vector<double> StdDayTgradk, std::vector<double> StdDayDRk, double GMR)
{
   double hda_out;
   
   int nLayers = StdDayHk.size()-1;

   // find layer n containing sigma in std day profile
   int n = findLayer_PD(nLayers, sigma, StdDayDRk);

   // from sigma, compute dh in layer std day profile, and Hda
   hda_out = StdDayHk[n] + dh_sgm(sigma, StdDayDRk[n], StdDayTk[n], StdDayTgradk[n], GMR);

   return(hda_out);
}

inline void initializeProfile(double T0, double GMR, std::vector<double> Hk, std::vector<double> Tk, int &nLayers, std::vector<double> &Tgradk, std::vector<double> &PRk, std::vector<double> &DRk)
{
   nLayers = Hk.size() - 1;

   Tgradk.resize(nLayers + 1);
   PRk.resize(nLayers + 1);
   DRk.resize(nLayers + 1);

   // for gradient, if below a threshold, then truncate to zero to guard 
   // against round-off error inappropriately generating linear thermal layers
   double eps = std::numeric_limits<double>::epsilon();

   for (int k = 0; k < nLayers; k++)
   {
      Tgradk[k] = (Tk[k + 1] - Tk[k]) / (Hk[k + 1] - Hk[k]);
      if( Tgradk[k] < 10 * eps)
         Tgradk[k] = 0.0;
   }
   Tgradk[nLayers] = 0.0;

   PRk[0] = 1.0;
   DRk[0] = PRk[0]*T0/Tk[0];
   for (int k = 0; k < nLayers; k++)
   {
      PRk[k + 1] = PRk[k] * pr(Hk[k + 1], Hk[k], Tk[k + 1], Tk[k], Tgradk[k], GMR);
      DRk[k + 1] = PRk[k + 1] * T0 / Tk[k + 1];
   }

}

inline void initializePRkDRk(double T0, double P0, double GMR, int n, double Hic, double Tic, double Pic, std::vector<double> Hk, std::vector<double> Tk, std::vector<double> Tgradk,
   std::vector<double> &PRk, std::vector<double> &DRk)
{
   int nLayers = Hk.size() - 1;

   double PicbyPn = pr(Hic, Hk[n], Tic, Tk[n], Tgradk[n], GMR);
   PRk[n] = (Pic / PicbyPn) / P0;   // P at bottom of layer holding initial condition point

   int k;
   // for each 0 <= k < n, compute Tk & Pk below initial condition point
   for (k = n - 1; k >= 0; k--)
   {
      PRk[k] = PRk[k + 1] / pr(Hk[k + 1], Hk[k], Tk[k + 1], Tk[k], Tgradk[k], GMR);
   }

   // for each n < k < nLayers+1, compute Tk & Pk above initial condition point
   for (k = n; k < nLayers; k++)
   {
      PRk[k + 1] = PRk[k] * pr(Hk[k + 1], Hk[k], Tk[k + 1], Tk[k], Tgradk[k], GMR);
   }

   // construct density ratio profile
   for (int k = 0; k < nLayers + 1; k++)
   {
      DRk[k] = PRk[k] * T0 / Tk[k];
   }

}


//inline void initializeTgPDrk(double T0, double P0, double GMR, double Hic, double Tic, double Pic, std::vector<double> Hk, std::vector<double> &Tk, std::vector<double> Tgradk,
//   std::vector<double> &PRk, std::vector<double> &DRk)
//{
//   int nLayers = Hk.size() - 1;
//
//   // construct temperature profile from gradient profile, (Hk,Tgradk)
// // find n such that Hk[n] < Hic < Hk[n + 1];
//   int n;
//   for (n = 1; n <= nLayers && Hk[n] < Hic; n++); n--;
//
//   Tk[n] = Tic - Tgradk[n] * (Hic - Hk[n]);   // T at bottom of layer holding initial condition point
//
//   // for each 0 <= k < n, compute Tk & Pk below initial condition point
//   for (int k = n - 1; k >= 0; k--)
//      Tk[k] = Tk[k + 1] - Tgradk[k] * (Hk[k + 1] - Hk[k]);
//
//   // for each n < k < nLayers+1, compute Tk & Pk above initial condition point
//   for (int k = n; k < nLayers; k++)
//      Tk[k + 1] = Tk[k] + Tgradk[k] * (Hk[k + 1] - Hk[k]);
//
//
//   double PicbyPn = pr(Hic, Hk[n], Tic, Tk[n], Tgradk[n], GMR);
//   PRk[n] = (Pic / PicbyPn) / P0;   // P at bottom of layer holding initial condition point
//
//   int k;
//   // for each 0 <= k < n, compute Tk & Pk below initial condition point
//   for (k = n - 1; k >= 0; k--)
//   {
//      PRk[k] = PRk[k + 1] / pr(Hk[k + 1], Hk[k], Tk[k + 1], Tk[k], Tgradk[k], GMR);
//   }
//
//   // for each n < k < nLayers+1, compute Tk & Pk above initial condition point
//   for (k = n; k < nLayers; k++)
//   {
//      PRk[k + 1] = PRk[k] * pr(Hk[k + 1], Hk[k], Tk[k + 1], Tk[k], Tgradk[k], GMR);
//   }
//
//   // construct density ratio profile
//   for (int k = 0; k < nLayers + 1; k++)
//   {
//      DRk[k] = PRk[k] * T0 / Tk[k];
//   }
//
//}

inline int AtmRatios(double hgp,  double T0, double GMR, std::vector<double> Hk, std::vector<double> Tk, std::vector<double> Tgradk, std::vector<double> PRk, std::vector<double> DRk,
   double &theta, double &sigma, double &delta, double &kappa)
{
   int nLayers = Hk.size() - 1;

   // find layer index, n
   int n = findLayer(nLayers, hgp, Hk);

   // compute properties from bottom of layer n up to given geopotential altitude hgp
   double T = Tk[n] + Tgradk[n] * (hgp - Hk[n]);
   theta = T / T0;
   delta = PRk[n] * pr(hgp, Hk[n], T, Tk[n], Tgradk[n], GMR);
   sigma = delta / theta;
   kappa = sqrt(theta);

   return(n);
};

inline double AirViscRatio(double theta, double Sb0)
{
   // Dyanmic(absolute) viscosity of air as a function of temperature for moderate pressures,
   // like atmospheric pressure, from Sutherland's law. This nondimensional form has been derived
   // to ease use of alternatives units. Form and constants derived from 1976 US Standard Atmosphere,
   // page 19 (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539.pdf), which is 
   // a widely used, primary atmosphere property model.
   //
   // Input
   //    theta === temperature ratio to sea level standard conditions, T/T0
   //    Sb0   === Sutherland temperature ratio by sea level std, S/T0
   //
   // Output
   //    vr    === absolute viscosity ratio to sea level standard, visc0

   double vr = pow(theta, 1.5) *(1 + Sb0) / (theta + Sb0);

   return(vr);
}


//------- Atmos class implementation -------

//Atmos::Atmos(AtmosParameters AtmPrms)
Atmos::Atmos(AtmosParameters AtmPrms)
{
   Re = AtmPrms.Re;               // planet radius
   GMR = AtmPrms.GMR;             // combined gravity and gas constant of planet atmosphere
   T0 = AtmPrms.T0;               // temperature at sea level std
   rho0 = AtmPrms.rho0;           // density at sea level std
   P0 = AtmPrms.P0;               // pressure at sea level std
   a0 = AtmPrms.a0;               // sonic speed at sea level std
   visc0 = AtmPrms.visc0;         // air viscosity at sea level std
   Sb0 = AtmPrms.S / T0;          // Sutherland temperature ratio at sea level std

   StdDayHk = AtmPrms.StdDayHk;
   StdDayTk = AtmPrms.StdDayTk;
   nStdLayers = StdDayHk.size() - 1;
   initializeProfile(T0, GMR, StdDayHk, StdDayTk, nStdLayers, StdDayTgradk, StdDayPRk, StdDayDRk);

   Hk = StdDayHk;
   Tk = StdDayTk;
   nLayers = Hk.size() - 1;
   initializeProfile(T0, GMR, Hk, Tk, nLayers, Tgradk, PRk, DRk);

   this->at(0.0);   // initialize object to datum (i.e. sea level)

};

// construct standard + dT deviation added to temperature profile vs altitude, SI units
Atmos::Atmos(double dT, AtmosParameters AtmPrms)
{
   Re = AtmPrms.Re;               // planet radius
   GMR = AtmPrms.GMR;             // combined gravity and gas constant of planet atmosphere
   T0 = AtmPrms.T0;               // temperature at sea level std
   rho0 = AtmPrms.rho0;           // density at sea level std
   P0 = AtmPrms.P0;               // pressure at sea level std
   a0 = AtmPrms.a0;               // sonic speed at sea level std
   visc0 = AtmPrms.visc0;         // air viscosity at sea level std
   Sb0 = AtmPrms.S / T0;          // Sutherland temperature ratio at sea level std

   StdDayHk = AtmPrms.StdDayHk;
   StdDayTk = AtmPrms.StdDayTk;
   nStdLayers = StdDayHk.size() - 1;
   initializeProfile(T0, GMR, StdDayHk, StdDayTk, nStdLayers, StdDayTgradk, StdDayPRk, StdDayDRk);

   Hk = StdDayHk;
   Tk = StdDayTk;
   nLayers = Hk.size() - 1;

   for (int m = 0; m < Hk.size(); m++)
      Tk[m] += dT;

   initializeProfile(T0, GMR, Hk, Tk, nLayers, Tgradk, PRk, DRk);

   this->at(0.0);   // initialize object to datum (i.e. sea level)

};

// construct using standard lapse rate profile used to derive temperature profile through (Hic,Tic,Pic), SI units

Atmos::Atmos(double Hic, double Tic, double Pic, AtmosParameters AtmPrms)
{
   Re = AtmPrms.Re;               // planet radius
   GMR = AtmPrms.GMR;             // combined gravity and gas constant of planet atmosphere
   T0 = AtmPrms.T0;               // temperature at sea level std
   rho0 = AtmPrms.rho0;           // density at sea level std
   P0 = AtmPrms.P0;               // pressure at sea level std
   a0 = AtmPrms.a0;               // sonic speed at sea level std
   visc0 = AtmPrms.visc0;         // air viscosity at sea level std
   Sb0 = AtmPrms.S / T0;          // Sutherland temperature ratio at sea level std

   StdDayHk = AtmPrms.StdDayHk;
   StdDayTk = AtmPrms.StdDayTk;
   nStdLayers = StdDayHk.size() - 1;
   initializeProfile(T0, GMR, StdDayHk, StdDayTk, nStdLayers, StdDayTgradk, StdDayPRk, StdDayDRk);

   Hk = StdDayHk;
   Tgradk = StdDayTgradk;

   nLayers = Hk.size() - 1;
   Tk.resize(nLayers + 1);
   PRk.resize(nLayers + 1);
   DRk.resize(nLayers + 1);

   // construct temperature profile from gradient profile, (Hk,Tgradk)
   // find n such that Hk[n] < Hic < Hk[n + 1];
   int n;
   for (n = 1; n <= nLayers && Hk[n] < Hic; n++); n--;

   Tk[n] = Tic - Tgradk[n] * (Hic - Hk[n]);   // T at bottom of layer holding initial condition point

   // for each 0 <= k < n, compute Tk & Pk below initial condition point
   for (int k = n - 1; k >= 0; k--)
      Tk[k] = Tk[k + 1] - Tgradk[k] * (Hk[k + 1] - Hk[k]);

   // for each n < k < nLayers+1, compute Tk & Pk above initial condition point
   for (int k = n; k < nLayers; k++)
      Tk[k + 1] = Tk[k] + Tgradk[k] * (Hk[k + 1] - Hk[k]);

   initializePRkDRk(T0, P0, GMR, n, Hic, Tic, Pic, Hk, Tk, Tgradk, PRk, DRk);

   this->at(0.0);   // initialize object to datum (i.e. sea level)

};

// custom from initial conditions Hic, Pic, breakpoints, and temperature profile, units determined by input
Atmos::Atmos(double Hic, double Pic, std::vector<double> Hj, std::vector<double> Tj, AtmosParameters AtmPrms)
{
   Re = AtmPrms.Re;               // planet radius
   GMR = AtmPrms.GMR;             // combined gravity and gas constant of planet atmosphere
   T0 = AtmPrms.T0;               // temperature at sea level std
   rho0 = AtmPrms.rho0;           // density at sea level std
   P0 = AtmPrms.P0;               // pressure at sea level std
   a0 = AtmPrms.a0;               // sonic speed at sea level std
   visc0 = AtmPrms.visc0;         // air viscosity at sea level std
   Sb0 = AtmPrms.S / T0;          // Sutherland temperature ratio at sea level std

   StdDayHk = AtmPrms.StdDayHk;
   StdDayTk = AtmPrms.StdDayTk;
   nStdLayers = StdDayHk.size() - 1;
   initializeProfile(T0, GMR, StdDayHk, StdDayTk, nStdLayers, StdDayTgradk, StdDayPRk, StdDayDRk);

   Hk = Hj;
   Tk = Tj;

   nLayers = Hk.size() - 1;
   Tk.resize(nLayers + 1);
   Tgradk.resize(nLayers + 1);
   PRk.resize(nLayers + 1);
   DRk.resize(nLayers + 1);

   // construct temperature gradient profile
   // for gradient, if below a threshold, then truncate to zero to guard 
   // against round-off error inappropriately generating linear thermal layers
   double eps = std::numeric_limits<double>::epsilon();

   for (int k = 0; k < nLayers; k++)
   {
      Tgradk[k] = (Tk[k + 1] - Tk[k]) / (Hk[k + 1] - Hk[k]);
      if (Tgradk[k] < 10 * eps)
         Tgradk[k] = 0.0;
   }
   Tgradk[nLayers] = 0.0;

// construct pressure ratio profile
   // find n such that Hk[n] < Hic < Hk[n + 1];
   int n;
   for (n = 1; n <= nLayers && StdDayHk_si[n] < Hic; n++); n--;

   double Tic = Tk[n] + Tgradk[n] * (Hic - Hk[n]);   // temperature at initial condition point from given temperature profile

   initializePRkDRk(T0, P0, GMR, n, Hic, Tic, Pic, Hk, Tk, Tgradk, PRk, DRk);


   this->at(0.0);   // initialize object to datum (i.e. sea level)

};

// custom from initial conditions (Hic, Tic, Pic), breakpoints, lapse rates, units determined by input
Atmos::Atmos(double Hic, double Tic, double Pic, std::vector<double> Hj, std::vector<double> Tgradj, AtmosParameters AtmPrms)
{
   Re = AtmPrms.Re;               // planet radius
   GMR = AtmPrms.GMR;             // combined gravity and gas constant of planet atmosphere
   T0 = AtmPrms.T0;               // temperature at sea level std
   rho0 = AtmPrms.rho0;           // density at sea level std
   P0 = AtmPrms.P0;               // pressure at sea level std
   a0 = AtmPrms.a0;               // sonic speed at sea level std
   visc0 = AtmPrms.visc0;         // air viscosity at sea level std
   Sb0 = AtmPrms.S / T0;          // Sutherland temperature ratio at sea level std

   StdDayHk = AtmPrms.StdDayHk;
   StdDayTk = AtmPrms.StdDayTk;
   nStdLayers = StdDayHk.size() - 1;
   initializeProfile(T0, GMR, StdDayHk, StdDayTk, nStdLayers, StdDayTgradk, StdDayPRk, StdDayDRk);

   // should test Hj & Tgradj have same number of elements
   // should impose Tgradk[last] = 0.0 to safe extrapolation above last layer
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

   Tk[n] = Tic - Tgradk[n] * (Hic - Hk[n]);   // T at bottom of layer holding initial condition point

   // for each 0 <= k < n, compute Tk & Pk below initial condition point
   for (int k = n - 1; k >= 0; k--)
      Tk[k] = Tk[k + 1] - Tgradk[k] * (Hk[k + 1] - Hk[k]);

   // for each n < k < nLayers+1, compute Tk & Pk above initial condition point
   for (int k = n; k < nLayers; k++)
      Tk[k + 1] = Tk[k] + Tgradk[k] * (Hk[k + 1] - Hk[k]);

   initializePRkDRk(T0, P0, GMR, n, Hic, Tic, Pic, Hk, Tk, Tgradk, PRk, DRk);

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

   AtmRatios(hgp_in, T0, GMR, Hk, Tk, Tgradk, PRk, DRk, Theta, Sigma, Delta, Kappa);

   TT = T0 * Theta;
   Rho = rho0 * Sigma;
   p = P0 * Delta;
   Sonic = a0 * Kappa;

   Visc = visc0*AirViscRatio(Theta, Sb0);

   hgm = hpa = hda = NAN;   // invalidate other altitudes

   return(0);
};

int Atmos::operator()(double hgp_in)
{
   hgp = hgp_in;

   AtmRatios(hgp_in, T0, GMR, Hk, Tk, Tgradk, PRk, DRk, Theta, Sigma, Delta, Kappa);

   TT = T0 * Theta;
   Rho = rho0 * Sigma;
   p = P0 * Delta;
   Sonic = a0 * Kappa;

   Visc = visc0 * AirViscRatio(Theta, Sb0);

   hgm = hpa = hda = NAN;   // invalidate other altitudes

   return(0);
};

// evaluate atmosphere properties, including alternative altitude definitions
int Atmos::atHgp(double hgp_in)   // given geoptential altitude
{
   this->at(hgp_in);

   hgm = hgp_in * Re / (Re - hgp_in);                                         // geometric altitude
   hpa = fHpa(Delta, StdDayHk, StdDayTk, StdDayTgradk, StdDayPRk, GMR);   // pressure altitude
   hda = fHda(Sigma, StdDayHk, StdDayTk, StdDayTgradk, StdDayDRk, GMR);   // density altitude

   return(0);
};

int Atmos::atHgm(double hgm_in)   // given geometric altitude
{
   hgp = hgm_in*Re / (Re + hgm_in);                                           // geopotential altitude

   this->at(hgp);

   hgm = hgm_in;                                                              // at() invalidates hgm, so it must be reset
   hpa = fHpa(Delta, StdDayHk, StdDayTk, StdDayTgradk, StdDayPRk, GMR);   // pressure altitude
   hda = fHda(Sigma, StdDayHk, StdDayTk, StdDayTgradk, StdDayPRk, GMR);   // density altitude

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
   hgp = Hk[n] + dh_pr(dlta, PRk[n], Tk[n], Tgradk[n], GMR);

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
   hgp = Hk[n] + dh_sgm(sgma, DRk[n], Tk[n], Tgradk[n], GMR);

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
   if (isnan(hgm))   // then hpm is not valid, does not correspond to current hgp
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