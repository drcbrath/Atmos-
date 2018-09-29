#pragma once

// US Standard Atmosphere 1976 (aka ISA) atmosphere model
// standard, non-standard, custom

// The atmosphere model is based on the hydrostatic equation and ideal gas law.
// When combined with empirical values for temperature profile versus altitude 
// and a few physical constants, a consistent model is obtained of atmospheric 
// properties from below sea level to the top of the given temperature profile.
// Physical constants and standard temperature profiles are supplied from the US
// Standard Atmosphere of 1976. Additionally, the temperature profile can be 
// modified by addition of a constant delta deviation, dT, added to the entire
// profile; this gives a simple, effective, and commonly used means of 
// representing hot or cold days or approximating specific flight conditions.
// Furthermore, a completely custom temperature profile may be used.

// Atmos class
// This class provides consistency with the math model; i.e. less able to 
// arbirarily set temperatures and violate the physics and math underlying the model.
// This class constructs the full temperature profile at construction/instantiation 
// of an Atmos object variable which then stays with that object. 

#include <vector>

// to derive from fundamentals, use these formulas instead of constants below
//const double g0 = 9.90665;               // (m/s^2) standard referencence gravity
//const double Ru = ?;                     // () universal gas constant
//const double M = ?;                      // () dry air molal mass
//const double gma = 1.4;                  // dry air ratio of specific heats
//const double a0 = sqrt(gma*Ru / M*T0);   // (m/s), speed of sound at SL std temperature
//const double GMR = g0*M / Rearth;        // (degK/m) combined gravity and gas constant of dry air on earth

// to define and hold constants necessary for the construction and evaluation of
// an atmos object, using a typedef & structure allows this to be passed to a
// constructor to override default values, hence allowing the Atmos class itself
// to be neutral regarding units, it can default to SI, but over-ride to any
// other unit system according to the units of the constants given in this structure

typedef struct AtmParms {
   const double Re;     // radius of the earth
   const double GMR;    // combined gravity and gas constant of dry air on earth
   const double H0;     // datum, sea level
   const double T0;     // SL std temp
   const double rho0;   // SL std density
   const double P0;     // SL std pressure
   const double a0;     // speed of sound at SL std temperature
   // std day profile (for use computing Hpa & Hda, without units specific code in the class nor its helper functions)
   std::vector<double> StdDayHk;
   std::vector<double> StdDayTk;

} AtmosParameters;

//------- constants -------

// SI units (default), use with atmos constructor to construct SI units based atmos object

const double Re_si = 6369000;          // (m) radius of the earth
const double GMR_si = 0.034163195;     // (degK/m) combined gravity and gas constant of dry air on earth
const double H0_si = 0.0;              // (m) datum, sea level
const double T0_si = 288.15;           // (K), SL std temp
const double rho0_si = 1.225;          // (kg/m^3), SL std density
const double P0_si = 101325;           // (N/m^2), SL std pressure
const double a0_si = 340.2686;         // (m/s), speed of sound at SL std temperature

   // defined atmosphere profiles <need to revise alternate day profiles! base on ref Mil 3013? references !>
const std::vector<double> StdDayHk_si({ 0.0,  11000.0,  20000.0,  32000.0,  47000.0,  51000.0,  71000.0,  84852.0 });
const std::vector<double> StdDayTk_si({ 288.15,   216.65,   216.65,   228.65,   270.65,   270.65,   214.65,   186.95 });
const std::vector<double> StdDayTgradk_si({ 288.15,   216.65,   216.65,   228.65,   270.65,   270.65,   214.65,   186.95 });

const std::vector<double> HotDayHk_si({ 0.0, 11000., 20000. });
const std::vector<double> HotDayTk_si({ 308.15,269.65,237.65 });

const std::vector<double> ColdDayHk_si({ 0.0, 11000., 20000. });
const std::vector<double> ColdDayTk_si({ 308.15,269.65,237.65 });

const std::vector<double> TropicalDayHk_si({ 0.0, 11000., 20000. });
const std::vector<double> TropicalDayTk_si({ 308.15,269.65,237.65 });

const std::vector<double> PolarDayHk_si({ 0.0, 11000., 20000. });
const std::vector<double> PolarDayTk_si({ 308.15,269.65,237.65 });

const AtmosParameters AtmosParameters_si = { Re_si, H0_si, T0_si, rho0_si, P0_si, a0_si, GMR_si, StdDayHk_si, StdDayTk_si };

// US units, use with atmos constructor to construct US units based atmos object
const double Re_us = 20895669.;        // 6369000 / 0.3048;                             // (m) radius of the earth
const double GMR_us = 0.01874329530;   // 0.034163195*(1.8*0.3048);                     // (degR/ft) combined gravity and gas constant of dry air on earth
const double H0_us = 0.0;              //                                               // (ft) datum, sea level
const double T0_us = 518.67;           // 288.15*1.8;                                   // (R), SL std temp
const double rho0_us = 0.00237689;     // 1.225*(0.068521766*0.3048*0.3048*0.3048);     // (sl/ft^3), SL std density
const double P0_us = 2116.2166;        // 101325 * (0.3048*0.3048 / 4.4482216152605);   // (lbf/ft^2), SL std pressure
const double a0_us = 1116.36680;       // 340.2686 / 0.3048;                            // (ft/s), speed of sound at SL std temperature

   // defined atmosphere profiles <need to revise alternate day profiles! base on ref Mil 3013? references !>
const std::vector<double> StdDayHk_us({ 0.000, 3352.800, 6096.000, 9753.600, 14325.600, 15544.800, 21640.800, 25862.890 });
const std::vector<double> StdDayTk_us({ 518.670, 389.970, 389.970, 411.570, 487.170, 487.170, 386.370, 336.510 });
const std::vector<double> StdDayTgradk_us({ 1701.673, 1279.429, 1279.429, 1350.295, 1598.327, 1598.327, 1267.618, 1104.035 });

const std::vector<double> HotDayHk_us({ 0.000, 3352.800, 6096.000 });
const std::vector<double> HotDayTk_us({ 554.670, 485.370, 427.770 });

const std::vector<double> ColdDayHk_us({ 0.000, 3352.800, 6096.000 });
const std::vector<double> ColdDayTk_us({ 554.670, 485.370, 427.770 });

const std::vector<double> TropicalDayHk_us({ 0.000, 3352.800, 6096.000 });
const std::vector<double> TropicalDayTk_us({ 554.670, 485.370, 427.770 });

const std::vector<double> PolarDayHk_us({ 0.000, 3352.800, 6096.000 });
const std::vector<double> PolarDayTk_us({ 554.670, 485.370, 427.770 });

const AtmosParameters AtmosParameters_us = { Re_us, H0_us, T0_us, rho0_us, P0_us, a0_us, GMR_us, StdDayHk_us, StdDayTk_us };

//------- Atmos class definitions -------

// profile defined from:
// (1) given Hic, Tic, Pic, Hk, Tgradk; derived Tk and Pk
// (2) given Hic, Pic, Hk, Tk; derived TgradK, Pk, and Tic=T(Hic), i.e. Tic evaluated at Hic based on (Hk,Tk)
// variations from (1) & (2) are actually specific selections opr means of specifying initial conditions or profiles

class Atmos
{
public:
   Atmos(AtmosParameters AtmPrms);                             // standard atmosphere
   Atmos(double dT, AtmosParameters AtmPrms);                            // standard + dT deviation added to temperature profile vs altitude
   Atmos(double Hic, double Tic, double Pic, AtmosParameters AtmPrms);   // standard lapse rate profile used to construct temperature profile and associated pressure breakpoints through initial condition point (Hic,Tic,Pic)
   Atmos(double Hic, double Pic, std::vector<double> Hj, std::vector<double> Tj, AtmosParameters AtmPrms);   // (2) custom temperature profile vs altitude used to construct profile with given initial condition Pic at Hic
   Atmos(double Hic, double Tic, double Pic, std::vector<double> Hj, std::vector<double> Tgradj, AtmosParameters AtmPrms);   // (1) custom from breakpoints and lapse rates, and initial condition (Hic, Tic, Pic)

	~Atmos();

   // evaluate properties at geopotential altitude, but not other altitude definitions to avoid performance penalty
	// in fact, invalidate others to avoid possibility of incorrect/out-of-sync values
	int at(double hgp);
	int operator()(double hgp);

	// evaluate atmosphere properties, including alternative altitude definitions, at given altitude
	int atHgp(double hgp);   // given geopotential altitude
	int atHgm(double hgm);   // given geometric altitude
	int atHpa(double Hpa);   // given pressure altitude
	int atHda(double hda);   // given density altitude

	// return altitude, evaluate if necessary first to provide valid result in case hgm, hpa, or hda does not already correspond to current hgp
	double Hgp(void);   // geoptential altitude
	double Hgm(void);   // geometric altitude
	double Hpa(void);   // pressure altitude
	double Hda(void);   // density altitude

   // return individual properties for current altitude
   double T(void);
   double rho(void);
   double P(void);
   double a(void);
   double visc(void);
   double theta(void);
   double sigma(void);
   double delta(void);
   double kappa(void);

private:
	// temperature vs geopotential altitude profile table, (Hk,Tk)
	std::vector<double> Hk;       // altitude breakpoint vector
   std::vector<double> Tk;       // temperature breakpoint vector
   std::vector<double> Tgradk;   // temperature gradient breakpoint vectors, dT/dH in each layer
   std::vector<double> PRk;      // pressure ratio breakpoint vector, P/P0 in each layer
   std::vector<double> DRk;      // density ratio breakpoint vector, Rho/Rho0 in each layer
   int nLayers;                  // number of layers in atmosphere

   // std day profile (for use computing Hpa & Hda, without units specific code in the class itself nor in its helper functions)
   std::vector<double> StdDayHk;
   std::vector<double> StdDayTk;
   std::vector<double> StdDayTgradk;
   std::vector<double> StdDayPRk;
   std::vector<double> StdDayDRk;
   int nStdLayers;

   // atmosphere model parameters, values at datum, Hgp=0
	double Re;     // planet radius
	double GMR;    // combined gravity and gas constant of planet atmosphere
   double T0;     // temperature
   double rho0;   // density
   double P0;     // pressure
   double a0;     // sonic speed

	// altitudes
	double hgp;   // geopotential altitude
	double hgm;   // geometric altitude
	double hpa;	  // pressure altitude
	double hda;	  // density altitude

   // properties at the most recently set hgp (either set directly as input or as result of evaluation using any other altitude definition)
   double TT;      // temperature
   double Rho;     // density
   double p;       // pressure
   double Sonic;   // sonic speed
   double Visc;    // viscosity
   double Theta;   // temperature ratio, T/T0
   double Sigma;   // density ratio, rho/rho0
   double Delta;   // pressure ratio, P/P0
   double Kappa;   // sonic speed ratio, a/a0

};

// defined atmosphere models
Atmos StdDay_si;
Atmos HotDay_si(H0_si, P0_si, HotDayHk_si, HotDayTk_si);
Atmos ColdDay_si(H0_si, P0_si, ColdDayHk_si, ColdDayTk_si);
Atmos TropicalDay_si(H0_si, P0_si, TropicalDayHk_si, TropicalDayTk_si);
Atmos PolarDay_si(H0_si, P0_si, PolarDayHk_si, PolarDayTk_si);

Atmos StdDay_us;
Atmos HotDay_us(H0_us, P0_us, HotDayHk_us, HotDayTk_us, AtmosParameters_us);
Atmos ColdDay_us(H0_us, P0_us, ColdDayHk_us, ColdDayTk_us, AtmosParameters_us);
Atmos TropicalDay_us(H0_us, P0_us, TropicalDayHk_us, TropicalDayTk_us, AtmosParameters_us);
Atmos PolarDay_us(H0_us, P0_us, PolarDayHk_us, PolarDayTk_us, AtmosParameters_us);
