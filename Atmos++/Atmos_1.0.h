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

// Atm functions
// For simple use of the atmosphere model
// These functions use the standard lapse rates and optionally
// a constant temperature deviation, dT, applied to the profile.

// Atmos class
// This class provides more consistency with the math model; i.e. less able to 
// arbirarily set temperatures and violate the physics and math underlying the model.
// This class constructs the full temperature profile at construction/instantiation 
// of an Atmos object variable which then stays with that object. This costs some 
// extra computation at construction time, but less at use time; so a simulation
// that might make many calls to evaluate the atmosphere properties should see
// an overall reduction in computation time.

#include <vector>

//------- constants -------
//const double g0 = 9.90665;                          // (m/s^2) standard referencence gravity
//const double Ru = ?;                             // () universal gas constant
//const double M = ?;                              // () dry air molal mass
//const double gma = 1.4;                             // dry air ratio of specific heats

//------- Atmos parameters -------
const double Rearth = 6369000;         // (m) radius of the earth
const double T0 = 288.15;              // (K), SL std temp
const double rho0 = 1.225;             // (kg/m^3), SL std density
const double P0 = 101325;              // (N/m^2), SL std pressure
const double a0 = 340.2686;            // (m/s), speed of sound at SL std temperature
const double GMRearth = 0.034163195;   // (degK/m) combined gravity and gas constant of dry air on earth

//const double a0 = sqrt(gma*Ru / M*T0);   // (m/s), speed of sound at SL std temperature
//const double GMRearth = g0*M / Rearth;   // (degK/m) combined gravity and gas constant of dry air on earth

//------- AtmProfile class, to hold temperature vs altitude profile -------
   // should add parameters & reference values to AtmProfile properties since these are necessary to define the profile
class AtmProfile {
public:
	AtmProfile();
	AtmProfile(std::vector<double> Hk, std::vector<double> Tk);
	AtmProfile(std::vector<double> Hk, std::vector<double> Tgradk, double T_1st);
	AtmProfile(std::vector<double> Hk, std::vector<double> Tgradk, double Hic, double Tic);

	std::vector<double> Hk;
	std::vector<double> Tk;
	std::vector<double> PRk;
	std::vector<double> DRk;
	std::vector<double> Tgradk;
	size_t nLayers;
};

// defined atmosphere profiles
const AtmProfile StdDayProfile(     {    0.0,  11000.0,  20000.0,  32000.0,  47000.0,  51000.0,  71000.0,  84852.0},
                                    { 288.15,   216.65,   216.65,   228.65,   270.65,   270.65,   214.65,   186.95} );
const AtmProfile HotDayProfile(     { 0.0, 11000., 20000. },
                                    { 308.15,269.65,237.65 });
const AtmProfile ColdDayProfile(    { 0.0, 11000., 20000. },
                                    { 308.15,269.65,237.65 });
const AtmProfile TropicalDayProfile({ 0.0, 11000., 20000. },
                                    { 308.15,269.65,237.65 });
const AtmProfile PolarDayProfile(   { 0.0, 11000., 20000. },
                                    { 308.15,269.65,237.65 });


//------- Atm functions -------

int AtmSIRatios(double Hgp, double &theta, double &sigma, double &delta, double &kappa);
int AtmSIRatiosdT(double Hgp, double &theta, double &sigma, double &delta, double &kappa, double dT);
int AtmSIRatiosT(double Hgp, double &theta, double &sigma, double &delta, double &kappa, double T);
int AtmSIRatiosCP(double Hgp, double &theta, double &sigma, double &delta, double &kappa, AtmProfile atm);


int AtmSI_Hgp(double Hgp, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa);
int AtmSI_HgpdT(double Hgp, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT);
int AtmSI_HgpT(double Hgp, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT);
int AtmSI_HgpCP(double Hgp, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, AtmProfile atm);

int AtmSI_Hgm(double Hgm, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa);
int AtmSI_HgmdT(double Hgm, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT);
int AtmSI_HgmT(double Hgm, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT);
int AtmSI_HgmCP(double Hgm, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, AtmProfile atm);

int AtmSI_Hpa(double Hpa, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa);
int AtmSI_HpadT(double Hpa, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT);
int AtmSI_HpaT(double Hpa, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT);
int AtmSI_HpaCP(double Hpa, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, AtmProfile atm);

int AtmSI_Hda(double Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa);
int AtmSI_HdadT(double Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT);
int AtmSI_HdaT(double Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT);
int AtmSI_HdaCP(double Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, AtmProfile atm);


int AtmSIall_Hgm(double Hgm, double &Hgp, double &Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa);
int AtmSIall_HgmdT(double Hgm, double &Hgp, double &Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT);
int AtmSIall_HgmT(double Hgm, double &Hgp, double &Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT);
int AtmSIall_HgmCP(double Hgm, double &Hgp, double &Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, AtmProfile atm);

int AtmSIall_Hgp(double &Hgm, double Hgp, double &Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa);
int AtmSIall_HgpdT(double &Hgm, double Hgp, double &Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT);
int AtmSIall_HgpT(double &Hgm, double Hgp, double &Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT);
int AtmSIall_HgpCP(double &Hgm, double Hgp, double &Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, AtmProfile atm);

int AtmSIall_Hpa(double &Hgm, double &Hgp, double Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa);
int AtmSIall_HpadT(double &Hgm, double &Hgp, double Hpa, double &Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT);
int AtmSIall_HpaT(double &Hgm, double &Hgp, double Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT);
int AtmSIall_HpaCP(double &Hgm, double &Hgp, double Hpa, double &Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, AtmProfile atm);

int AtmSIall_Hda(double &Hgm, double &Hgp, double &Hpa, double Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa);
int AtmSIall_HdadT(double &Hgm, double &Hgp, double &Hpa, double Hda, double &T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double dT);
int AtmSIall_HdaT(double &Hgm, double &Hgp, double &Hpa, double Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, double &dT);
int AtmSIall_HdaCP(double &Hgm, double &Hgp, double &Hpa, double Hda, double T, double &rho, double &P, double &a, double &theta, double &sigma, double &delta, double &kappa, AtmProfile atm);


typedef struct Atm { double Hgp; double Hgm; double Hpa; double Hda; double T; double rho; double P; double a; double theta; double sigma; double delta; double kappa; double dT; } Atm;

Atm AtmSI_Hgp(double Hgp);
Atm AtmSI_HgpdT(double Hgp, double dT);
Atm AtmSI_HgpT(double Hgp, double T);
Atm AtmSI_HgpCP(double Hgp, AtmProfile atm);

Atm AtmSI_Hgm(double Hgm);
Atm AtmSI_HgmdT(double Hgm, double dT);
Atm AtmSI_HgmT(double Hgm, double T);
Atm AtmSI_HgmCP(double Hgm, AtmProfile atm);

Atm AtmSI_Hpa(double Hpa);
Atm AtmSI_HpadT(double Hpa, double dT);
Atm AtmSI_HpaT(double Hpa, double T);
Atm AtmSI_HpaCP(double Hpa, AtmProfile atm);

Atm AtmSI_Hda(double Hda);
Atm AtmSI_HdadT(double Hda, double dT);
Atm AtmSI_HdaT(double Hda, double T);
Atm AtmSI_HdaCP(double Hda, AtmProfile atm);


Atm AtmSIall_Hgp(double Hgp);
Atm AtmSIall_HgpdT(double Hgp, double dT);
Atm AtmSIall_HgpT(double Hgp, double T);
Atm AtmSIall_HgpCP(double Hgp, AtmProfile atm);

Atm AtmSIall_Hgm(double Hgm);
Atm AtmSIall_HgmdT(double Hgm, double dT);
Atm AtmSIall_HgmT(double Hgm, double T);
Atm AtmSIall_HgmCP(double Hgm, AtmProfile atm);

Atm AtmSIall_Hpa(double Hpa);
Atm AtmSIall_HpadT(double Hpa, double dT);
Atm AtmSIall_HpaT(double Hpa, double T);
Atm AtmSIall_HpaCP(double Hpa, AtmProfile atm);

Atm AtmSIall_Hda(double Hda);
Atm AtmSIall_HdadT(double Hda, double dT);
Atm AtmSIall_HdaT(double Hda, double T);
Atm AtmSIall_HdaCP(double Hda, AtmProfile atm);


//------- Atmos class definitions -------

class AtmosSI
{
public:
	AtmosSI();                                // standard atmosphere
   AtmosSI(double dT);                       // standard + dT deviation added to temperature profile vs altitude
   AtmosSI(double Hic, double Tic);          // standard lapse rate profile used to construct temperature profile through (Hic,Tic)
   AtmosSI(std::vector<double> Hk, std::vector<double> Tk);   // custom temperature profile vs altitude
   AtmosSI(std::vector<double> Hk, std::vector<double> Tk, double dT);   // custom temperature profile vs altitude, modified by additional dT
   AtmosSI(std::vector<double> Hk, std::vector<double> Tgradk, double Hic, double Tic);   // custom temperature lapse rate vs altitude used to construct temperature profile through (Hic,Tic)

	~AtmosSI();

   // evaluate properties at geopotential altitude, but not other altitude definitions
	// in fact, invalidate others to avoid performance penalty and avoid incorrect/out-of-sync values
	int at(double hgp);
	int operator()(double hgp);

	// evaluate atmosphere properties, including alternative altitude definitions, at given altitude
	int atHgp(double hgp);   // given geopotential altitude
	int atHgm(double hgm);   // given geometric altitude
	int atHpa(double Hpa);   // given pressure altitude
	int atHda(double hda);   // given density altitude

	// return altitude, evaluate if necessary first to provide valid result in case hgm, hpa, or hda does not already correspond to current hgp
	double Hgp(void); // geoptential altitude
	double Hgm(void);	// geometric altitude
	double Hpa(void);	// pressure altitude
	double Hda(void);	// density altitude

   // return properties for current altitude
   double T(void);
   double rho(void);
   double P(void);
   double a(void);
   double theta(void);
   double sigma(void);
   double delta(void);
   double kappa(void);

private:
	// temperature vs geopotential altitude profile table, (Hk,Tk)
	std::vector<double> Hk;     // altitude breakpoint vector
   std::vector<double> Tk;     // temperature breakpoint vector
   std::vector<double> Tgradk;   // temperature gradient breakpoint vectors, dT/dH in each layer
   std::vector<double> PRk;     // pressure ratio breakpoint vector, P/P0 in each layer
   std::vector<double> DRk;     // density ratio breakpoint vector, Rho/Rho0 in each layer
   size_t nLayers;             // number of layers in atmosphere
	double Re  = Rearth;        // (m) planet radius; default is earth
	double GMR = GMRearth;      // (degK/m) combined gravity and gas constant of planet atmosphere; default is dry air on earth


	// altitudes
	double hgp;   // geopotential altitude
	double hgm;   // geometric altitude
	double hpa;	  // pressure altitude
	double hda;	  // density altitude

   // properties at the most recently set Hgp (either set directly as input or as result of evaluation using any other altitude definition)
   double T;
   double rho;
   double P;
   double a;
   double theta;
   double sigma;
   double delta;
   double kappa;

};

// defined atmosphere models
// these do not need to be const; since the profile member data is private, it cannot be changed later by accident
AtmosSI StdDay(StdDayProfile);
AtmosSI HotDay(HotDayProfile);
AtmosSI ColdDay(ColdDayProfile);
AtmosSI TropicalDay(TropicalDayProfile);
AtmosSI PolarDay(PolarDayProfile);
