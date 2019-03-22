// testAtmos.cpp : Test Atmos++ library, Atmos class
//
// (1) test constructors and atmosphere property computation over full domain of altitudes
// (2) then test alternative altitude computations
// (3) test method operator(), should be the same as method at(Hgp), and at(Hgp) method is used by the other at methods: atHgp, atHgm, ... so at(HgP) will have been tested sufficiently by now

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

void writeTestTable(Atmos atmos, std::vector<double> Hgp_in)
{
   std::cout << "    Hgm     Hgp     Hpa     Hda       T           rho             P       a         theta         sigma         delta         kappa          visc" << std::endl;
   std::cout << "    (m)     (m)     (m)     (m)     (K)      (kg/m^3)       (N/m^2)   (m/s)           (1)           (1)           (1)           (1)     (N*s/m^2)" << std::endl;

   for (size_t i = 0; i < Hgp_in.size(); i++)
   {
      atmos.atHgp(Hgp_in[i]);   // compute all properties, including alternate altitude definitions, at test point
      writeTestLine(atmos.Hgm(), atmos.Hgp(), atmos.Hpa(), atmos.Hda(), atmos.T(), atmos.rho(), atmos.P(), atmos.a(), atmos.theta(), atmos.sigma(), atmos.delta(), atmos.kappa(), atmos.visc());
   }
}

void test_atH(Atmos atmos, double hgp_in)
{
   atmos.atHgp(hgp_in);  writeTestLine(atmos.Hgm(), atmos.Hgp(), atmos.Hpa(), atmos.Hda(), atmos.T(), atmos.rho(), atmos.P(), atmos.a(), atmos.theta(), atmos.sigma(), atmos.delta(), atmos.kappa(), atmos.visc());
   atmos.atHgm(atmos.Hgm());  writeTestLine(atmos.Hgm(), atmos.Hgp(), atmos.Hpa(), atmos.Hda(), atmos.T(), atmos.rho(), atmos.P(), atmos.a(), atmos.theta(), atmos.sigma(), atmos.delta(), atmos.kappa(), atmos.visc());
   atmos.atHpa(atmos.Hpa());  writeTestLine(atmos.Hgm(), atmos.Hgp(), atmos.Hpa(), atmos.Hda(), atmos.T(), atmos.rho(), atmos.P(), atmos.a(), atmos.theta(), atmos.sigma(), atmos.delta(), atmos.kappa(), atmos.visc());
   atmos.atHda(atmos.Hda());  writeTestLine(atmos.Hgm(), atmos.Hgp(), atmos.Hpa(), atmos.Hda(), atmos.T(), atmos.rho(), atmos.P(), atmos.a(), atmos.theta(), atmos.sigma(), atmos.delta(), atmos.kappa(), atmos.visc());
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
   std::vector<double> Hj = { 0.0,  8250.0, 11000.0, 17750.0, 20000.0, 29000.0, 32000.0, 43250.0, 47000.0, 50000.0, 51000.0, 66000.0, 71000.0, 81389.0, 84852.0 };
   std::vector<double> Tj = { 288.15,234.53,216.65,216.65,216.65,225.65,228.65,260.15,270.65,270.65,270.65,228.65,214.65,193.88,186.95 };
   std::vector<double> T30j = { 318.15,234.53,216.65,216.65,216.65,225.65,228.65,260.15,270.65,270.65,270.65,228.65,214.65,193.88,186.95 };
   std::vector<double> Tgradj = { -0.0065,-0.0065, 0.0000, 0.0000, 0.0010, 0.0010, 0.0028, 0.0028, 0.0000, 0.0000,-0.0028,-0.0028,-0.0020,-0.0020, 0.0 };

   // test results
   std::vector<double> Hgp_out(Hgp_in.size());
   std::vector<double> Hgm_out(Hgp_in.size());
   std::vector<double> Hpa_out(Hgp_in.size());
   std::vector<double> Hda_out(Hgp_in.size());
   std::vector<double> T_out(Hgp_in.size());

   std::cout << std::endl << "Atmosphere Properties Model Test" << std::endl << std::endl;

    // SI
   // (1) standard atmosphere, test omitting AtmPrms, should generate std day object and results
   std::cout << std::endl << "(1) standard day" << std::endl;

   Atmos Atmos1(AtmosParameters_si);   // standard atmosphere, SI units

   // test operator(), at(), atHgp() --- should give same results
   Atmos1(5500.0); std::cout << "Atmos1() --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.at(5500.0); std::cout << "Atmos1.at --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHgp(5500.0); std::cout << "Atmos1.atHgp --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;

   //point tests
   // troposphere
   test_atH(Atmos1, 5500.0);
   // stratosphere
   test_atH(Atmos1, 15500.0);

   // full range test
   writeTestTable(Atmos1, Hgp_in);

   // (2) dT offset, use dT from case (3) of test data, should match results of test data
   std::cout << std::endl << "(2) dT offset" << std::endl;

   Atmos Atmos2(dT, AtmosParameters_si);   // standard + dT deviation added to temperature profile vs altitude
   
   //point tests
   // troposphere
   test_atH(Atmos2, 5500.0);
   // stratosphere
   test_atH(Atmos2, 15500.0);

   // full range test
   writeTestTable(Atmos2, Hgp_in);

   // (3) use standard lapse rates and SI units from AtmosParameters_si with an initial condition selected to match an output from (2), though not at an existing breakpoint (for generality). Since the intial condition does not change the temperature profile, the results should match those obtained from test (2) above.
   std::cout << std::endl << "(3) custom intial conditions Hic, Tic, Pic; test with dT offset and standard lapse rates, should match results of test (2) above" << std::endl;

   // get intial condition point, say at 500 m
   Atmos2.atHgp(500);

   Hic = Atmos2.Hgp();
   Tic = Atmos2.T();
   Pic = Atmos2.P();

   Atmos Atmos3(Hic, Tic, Pic, AtmosParameters_si);   // standard lapse rates used to construct profiles through initial condition point

   //point tests
   // troposphere
   test_atH(Atmos3, 5500.0);
   // stratosphere
   test_atH(Atmos3, 15500.0);

   // full range test
   writeTestTable(Atmos3, Hgp_in);

//   // (4) custom initial conditions for H & P and custom temperature profile (T initial condition is then interpolated from custom profile)
//   // construct a custom profile that is equivalent in layer extents and lapse rates to one from above, though not at the same breakpoints, then the results must match.
   std::cout << std::endl << "(4) custom intial conditions, altitude and temperature profile; input chosen such that results should match standard day test (1) above" << std::endl;

   Atmos1.atHgp(500);

   Hic = Atmos1.Hgp();
   Pic = Atmos1.P();

   Atmos Atmos4(Hic, Pic, Hj, Tj, AtmosParameters_si);

   //point tests
   // troposphere
   test_atH(Atmos4, 5500.0);
   // stratosphere
   test_atH(Atmos4, 15500.0);

   // full range test
   writeTestTable(Atmos4, Hgp_in);

   // (5) fully custom, initial conditions and temperature gradient (lapse rates). How to test? (1) construct conditions and profile equivalent to case (4); then results must match. Plus, to test that the custom input really do alter the results and the match in (1) is not some kind of bug, create another and different custom profile and conditions to show that the results do change with input.
   std::cout << std::endl << "(5) fully custom, intial conditions, altitude breakpoints, and lapse rates; input chosen such that results should match dT offset results of test (2) above" << std::endl;

   // get intial condition point, say at 500 m
   Atmos2.atHgp(500);

   Hic = Atmos2.Hgp();
   Tic = Atmos2.T();
   Pic = Atmos2.P();

   Atmos Atmos5(Hic, Tic, Pic, Hj, Tgradj, AtmosParameters_si);

   //point tests
   // troposphere
   test_atH(Atmos5,  5500.0);
   // stratosphere
   test_atH(Atmos5, 15500.0);

   // full range test
   writeTestTable(Atmos5, Hgp_in);

}
