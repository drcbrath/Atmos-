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
// This class constructs the full temperature profile at construction/instantiation 
// of an Atmos object variable which then stays with that object. 

#include <vector>

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
   const double visc0;  // air dynamic (absolute) viscosity at SL std
   const double S;      // air Sutherland temperature, for computation of viscosity
   // std day profile (for computing Hpa & Hda, without units specific code in class nor helper functions)
   const std::vector<double> StdDayHk;
   const std::vector<double> StdDayTk;

} AtmosParameters;

//------- constants available externally -------

// SI units (default), use with atmos constructor to construct SI units based atmos object

extern const double Re_si;      // (m) radius of the earth
extern const double GMR_si;     // (degK/m) combined gravity and gas constant of dry air on earth
extern const double H0_si;      // (m) datum, sea level
extern const double T0_si;      // (K) SL std temp
extern const double rho0_si;    // (kg/m^3) SL std density
extern const double P0_si;      // (N/m^2) SL std pressure
extern const double a0_si;      // (m/s) speed of sound at SL std temperature
extern const double visc0_si;   // (N/m^2) air dynamic viscosity at SL std
extern const double S_si;       // (K) air Sutherland temperature for viscosity computation

// defined atmosphere profiles <need to revise alternate day profiles! base on ref Mil 3013? references !>
extern const std::vector<double> StdDayHk_si;
extern const std::vector<double> StdDayTk_si;
extern const std::vector<double> StdDayTgradk_si;

extern const std::vector<double> HotDayHk_si;
extern const std::vector<double> HotDayTk_si;

extern const std::vector<double> ColdDayHk_si;
extern const std::vector<double> ColdDayTk_si;

extern const std::vector<double> TropicalDayHk_si;
extern const std::vector<double> TropicalDayTk_si;

extern const std::vector<double> PolarDayHk_si;
extern const std::vector<double> PolarDayTk_si;

extern const AtmosParameters AtmosParameters_si;

// US units, use with atmos constructor to construct US units based atmos object
extern const double Re_us;      // (m) radius of the earth
extern const double GMR_us;     // (degK/m) combined gravity and gas constant of dry air on earth
extern const double H0_us;      // (m) datum, sea level
extern const double T0_us;      // (K) SL std temp
extern const double rho0_us;    // (kg/m^3) SL std density
extern const double P0_us;      // (N/m^2) SL std pressure
extern const double a0_us;      // (m/s) speed of sound at SL std temperature
extern const double visc0_us;   // (N/m^2) air dynamic viscosity at SL std
extern const double S_us;       // (K) air Sutherland temperature for viscosity computation

//------- Atmos class definitions -------

// profile defined from:
// (1) given Hic, Tic, Pic, Hk, Tgradk; derived Tk and Pk
// (2) given Hic, Pic, Hk, Tk; derived TgradK, Pk, and Tic=T(Hic), i.e. Tic evaluated at Hic based on (Hk,Tk)
// variations from (1) & (2) are actually specific selections or means of specifying initial conditions or profiles

class Atmos
{
public:
   // note: for all constructors, last parameter, AtmPrms, is optional;
   // if omitted, SI units AtmosParameters_si will be used. Alternatively,
   // one may use AtmosParameters_us for US customary units, or assemble a
   // custom AtmosParameter to suit.

   Atmos(AtmosParameters AtmPrms = AtmosParameters_si);                             // standard atmosphere
   Atmos(double dT, AtmosParameters AtmPrms = AtmosParameters_si);                            // standard + dT deviation added to temperature profile vs altitude
   Atmos(double Hic, double Tic, double Pic, AtmosParameters AtmPrms = AtmosParameters_si);   // standard lapse rate profile used to construct temperature profile and associated pressure breakpoints through initial condition point (Hic,Tic,Pic)
   Atmos(double Hic, double Pic, std::vector<double> Hj, std::vector<double> Tj, AtmosParameters AtmPrms = AtmosParameters_si);   // (2) custom temperature profile vs altitude used to construct profile with given initial condition Pic at Hic
   Atmos(double Hic, double Tic, double Pic, std::vector<double> Hj, std::vector<double> Tgradj, AtmosParameters AtmPrms = AtmosParameters_si);   // (1) custom temperature gradient profile (i.e. lapse rates) and initial condition (Hic, Tic, Pic)

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

   // atmosphere model parameters, values at datum on standard day (typically sea level std)
	double Re;     // planet radius
	double GMR;    // combined gravity and gas constant of planet atmosphere
   double T0;     // temperature
   double rho0;   // density
   double P0;     // pressure
   double a0;     // sonic speed
   double visc0;  // air viscosity
   double Sb0;    // Sutherland temperature ratio by T0

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


extern const AtmosParameters AtmosParameters_us;

// defined atmosphere models
extern Atmos StdDay_si;
extern Atmos HotDay_si;
extern Atmos ColdDay_si;
extern Atmos TropicalDay_si;
extern Atmos PolarDay_si;

extern Atmos StdDay_us;
extern Atmos HotDay_us;
extern Atmos ColdDay_us;
extern Atmos TropicalDay_us;
extern Atmos PolarDay_us;
