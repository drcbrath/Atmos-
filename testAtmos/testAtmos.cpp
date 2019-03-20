// testAtmos.cpp : Test Atmos++ library, Atmos class
//
// (1) test constructors and atmosphere property computation over full domain of altitudes
// (2) then test alternative altitude computations
// (3) test method operator(), should be the same as method at(Hgp), and at(Hgp) method is used by the other at methods: atHgp, atHgm, ... so at(HgP) will       have been tested sufficiently by now

#include "stdafx.h"
#include <vector>
#include "Atmos.h"

int writeTestLine(double Hgm, double Hgp, double Hpa, double Hda, double T, double rho, double P, double sonic, double theta, double sigma, double delta, double sonicratio, double visc)
{
   std::cout << std::fixed << std::setw(7) << std::setprecision(1) << Hgm
      << " " << std::fixed << std::setw(7) << std::setprecision(1) << Hgp
      << " " << std::fixed << std::setw(7) << std::setprecision(1) << Hpa
      << " " << std::fixed << std::setw(7) << std::setprecision(1) << Hda
      << " " << std::fixed << std::setw(7) << std::setprecision(3) << T
      << " " << std::scientific << std::setprecision(6) << std::setw(13) << rho
      << " " << std::scientific << std::setprecision(6) << std::setw(13) << P
      << " " << std::fixed << std::setw(7) << std::setprecision(2) << sonic
      << " " << std::scientific << std::setw(13) << std::setprecision(6) << theta
      << " " << std::scientific << std::setw(13) << std::setprecision(6) << sigma
      << " " << std::scientific << std::setw(13) << std::setprecision(6) << delta
      << " " << std::scientific << std::setw(13) << std::setprecision(6) << sonicratio
      << " " << std::scientific << std::setw(13) << std::setprecision(6) << visc
      << std::endl;

   return(0);
}

int main()
{
   // test altitudes
      // test at: altitudes below table, at all breakpoints, between all breakpoints, above table
   std::vector<double> Hgp_in = { -100, 0.0, 5500, 11000.0, 15500, 20000.0, 26000, 32000.0, 39500, 47000.0, 49000, 51000.0, 61000, 71000.0, 80000, 84852.0, 90000 };
   std::vector<double> Hgm_in(Hgp_in.size());
   std::vector<double> Hpa_in(Hgp_in.size());
   std::vector<double> Hda_in(Hgp_in.size());
   std::vector<double> T_in(Hgp_in.size());

   // test initial conditions
   double Hic, Tic, Pic;
   // test dT offset
   double dT = 30.0;

   // test profiles
   std::vector<double> Hj;
   std::vector<double> Tj;
   std::vector<double> Tgradj;

   // test results
   std::vector<double> Hgp_out(Hgp_in.size());
   std::vector<double> Hgm_out(Hgp_in.size());
   std::vector<double> Hpa_out(Hgp_in.size());
   std::vector<double> Hda_out(Hgp_in.size());
   std::vector<double> T_out(Hgp_in.size());

   
   // SI
   // (1) standard atmosphere, test omitting AtmPrms, should generate std day object and results
   Atmos Atmos1(AtmosParameters_si);   // standard atmosphere, SI units

   // point tests
   // troposphere
   Atmos1(5500.0); std::cout << "Atmos1() --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.at(5500.0); std::cout << "Atmos1.at --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHgp(5500.0); std::cout << "Atmos1.atHgp --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHgm(Atmos1.Hgm()); std::cout << "Atmos1.atHgm --- Hgm, P, rho:" << Atmos1.Hgm() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHpa(Atmos1.Hpa()); std::cout << "Atmos1.atHpa --- Hpa, P, rho:" << Atmos1.Hpa() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHda(Atmos1.Hda()); std::cout << "Atmos1.atHda --- Hda, P, rho:" << Atmos1.Hda() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   // stratosphere
   Atmos1(15500.0); std::cout << "Atmos1() --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.at(15500.0); std::cout << "Atmos1.at --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHgp(15500.0); std::cout << "Atmos1.atHgp --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHgm(Atmos1.Hgm()); std::cout << "Atmos1.atHgm --- Hgm, P, rho:" << Atmos1.Hgm() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHpa(Atmos1.Hpa()); std::cout << "Atmos1.atHpa --- Hpa, P, rho:" << Atmos1.Hpa() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHda(Atmos1.Hda()); std::cout << "Atmos1.atHda --- Hda, P, rho:" << Atmos1.Hda() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;

   // full range tests

   std::cout << "Atmosphere Properties Model Test Data-- - Validated against \"US Standard Atmosphere, 1976\"" << std::endl << std::endl << std::endl;
   std::cout << "(1) Comparison of : std day, evaluated at test altitudes" << std::endl;
   std::cout << "    Hgm     Hgp     Hpa     Hda       T           rho             P       a         theta         sigma         delta         kappa          visc" << std::endl;
   std::cout << "    (m)     (m)     (m)     (m)     (K)      (kg/m^3)       (N/m^2)   (m/s)           (1)           (1)           (1)           (1)     (N*s/m^2)" << std::endl;

   for (size_t i = 0; i < Hgp_in.size(); i++)
   {
      Atmos1.atHgp(Hgp_in[i]);   // compute all properties, including alternate altitude definitions, at test point
      writeTestLine(Atmos1.Hgm(), Atmos1.Hgp(), Atmos1.Hpa(), Atmos1.Hda(), Atmos1.T(), Atmos1.rho(), Atmos1.P(), Atmos1.a(), Atmos1.theta(), Atmos1.sigma(), Atmos1.delta(), Atmos1.kappa(), Atmos1.visc());
   }

   // (2) dT offset, use dT from case (3) of test data, should match results of test data
   Atmos Atmos2(dT, AtmosParameters_si);   // standard + dT deviation added to temperature profile vs altitude
   
   // point tests
   // troposphere
   Atmos2.atHgp(5500.0); std::cout << "Atmos2.atHgp --- Hgp, P, rho:" << Atmos2.Hgp() << " " << Atmos2.P() << " " << Atmos2.rho() << std::endl;
   Atmos2.atHgm(Atmos2.Hgm()); std::cout << "Atmos2.atHgm --- Hgm, P, rho:" << Atmos2.Hgm() << " " << Atmos2.P() << " " << Atmos2.rho() << std::endl;
   Atmos2.atHpa(Atmos2.Hpa()); std::cout << "Atmos2.atHpa --- Hpa, P, rho:" << Atmos2.Hpa() << " " << Atmos2.P() << " " << Atmos2.rho() << std::endl;
   Atmos2.atHda(Atmos2.Hda()); std::cout << "Atmos2.atHda --- Hda, P, rho:" << Atmos2.Hda() << " " << Atmos2.P() << " " << Atmos2.rho() << std::endl;
   // stratosphere
   Atmos2.atHgp(15500.0); std::cout << "Atmos2.atHgp --- Hgp, P, rho:" << Atmos2.Hgp() << " " << Atmos2.P() << " " << Atmos2.rho() << std::endl;
   Atmos2.atHgm(Atmos2.Hgm()); std::cout << "Atmos2.atHgm --- Hgm, P, rho:" << Atmos2.Hgm() << " " << Atmos2.P() << " " << Atmos2.rho() << std::endl;
   Atmos2.atHpa(Atmos2.Hpa()); std::cout << "Atmos2.atHpa --- Hpa, P, rho:" << Atmos2.Hpa() << " " << Atmos2.P() << " " << Atmos2.rho() << std::endl;
   Atmos2.atHda(Atmos2.Hda()); std::cout << "Atmos2.atHda --- Hda, P, rho:" << Atmos2.Hda() << " " << Atmos2.P() << " " << Atmos2.rho() << std::endl;
   // full range tests

   std::cout << "Atmosphere Properties Model Test Data-- - Validated against \"US Standard Atmosphere, 1976\"" << std::endl << std::endl << std::endl;
   std::cout << "(2) dT offset, evaluated at test altitudes, should match results of test data case (3) results, " << std::endl;
   std::cout << "    Hgm     Hgp     Hpa     Hda       T           rho             P       a         theta         sigma         delta         kappa          visc" << std::endl;
   std::cout << "    (m)     (m)     (m)     (m)     (K)      (kg/m^3)       (N/m^2)   (m/s)           (1)           (1)           (1)           (1)     (N*s/m^2)" << std::endl;

   for (size_t i = 0; i < Hgp_in.size(); i++)
   {
      Atmos2.atHgp(Hgp_in[i]);   // compute all properties, including alternate altitude definitions, at test point
      writeTestLine(Atmos2.Hgm(), Atmos2.Hgp(), Atmos2.Hpa(), Atmos2.Hda(), Atmos2.T(), Atmos2.rho(), Atmos2.P(), Atmos2.a(), Atmos2.theta(), Atmos2.sigma(), Atmos2.delta(), Atmos2.kappa(), Atmos2.visc());
   }

   //// (3) use standard lapse rates and SI units from AtmosParameters_si with an initial condition selected to match an output from (2), though not at an existing breakpoint (for generality). Since the intial condition does not change the temperature profile, the results should match those obtained from test (2) above.

//   Hic = 0.0;
//   Tic = 288.15;
//   Pic = 101325.0; 
//   Atmos(Hic, Tic, Pic, AtmosParameters_si);   // standard lapse rates used to construct profiles through initial condition point
//
//   // (4) custom initial conditions for H & P and custom temperature profile (T initial condition is then interpolated from custom profile)
//   // construct a custom profile that is equivalent in layer extents and lapse rates to one from above, though not at the same breakpoints, then the results must match.
//
//   Hic = 0.0;
//   Pic = 101325.0;
//   Atmos(Hic, Pic, Hj, Tj, AtmosParameters_si);
//
//   // (5) fully custom, initial conditions and temperature gradient (lapse rates). How to test? (1) construct conditions and profile equivalent to case (4); then results must match. Plus, to test that the custom input really do alter the results and the match in (1) is not some kind of bug, create another and different custom profile and conditions to show that the results do change with input.
//
//
//   Hic = 0.0;        // set same as (4) above
//   Tic = 288.15;     // set from output of (4) at Hic
//   Pic = 101325.0;   // set same as (4) above
//   // Tgradj = asdfasdf;  // compute using (Hj,Tj) from (4)
//   Atmos(Hic, Tic, Pic, Hj, Tgradj, AtmosParameters_si);   // (1) custom temperature gradient profile (i.e. lapse rates) and initial condition (Hic, Tic, Pic)
//
//   // repeat with alternative inputs to get different object
//
////------- tests -------

}
