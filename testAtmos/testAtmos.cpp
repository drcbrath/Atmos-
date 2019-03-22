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
   std::vector<double> Tgradj = { -0.0065,-0.0065, 0.0000, 0.0000, 0.0010, 0.0010, 0.0028, 0.0028, 0.0000, 0.0000,-0.0028,-0.0028,-0.0020,-0.0020,-0.0020 };

   // test results
   std::vector<double> Hgp_out(Hgp_in.size());
   std::vector<double> Hgm_out(Hgp_in.size());
   std::vector<double> Hpa_out(Hgp_in.size());
   std::vector<double> Hda_out(Hgp_in.size());
   std::vector<double> T_out(Hgp_in.size());

    // SI
   // (1) standard atmosphere, test omitting AtmPrms, should generate std day object and results
   std::cout << std::endl << "Atmosphere Properties Model Test Data" << std::endl << std::endl;
   std::cout << "(1) standard day" << std::endl;

   Atmos Atmos1(AtmosParameters_si);   // standard atmosphere, SI units

   // test operator(), at(), atHgp() --- should give same results
   Atmos1(5500.0); std::cout << "Atmos1() --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.at(5500.0); std::cout << "Atmos1.at --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHgp(5500.0); std::cout << "Atmos1.atHgp --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;

   // point tests, call with all altitude types, walk through types as input such that all results should be the same
   // troposphere
   Atmos1.atHgp(5500.0); std::cout << "Atmos1.atHgp --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHgm(Atmos1.Hgm()); std::cout << "Atmos1.atHgm --- Hgm, P, rho:" << Atmos1.Hgm() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHpa(Atmos1.Hpa()); std::cout << "Atmos1.atHpa --- Hpa, P, rho:" << Atmos1.Hpa() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHda(Atmos1.Hda()); std::cout << "Atmos1.atHda --- Hda, P, rho:" << Atmos1.Hda() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   // stratosphere
   Atmos1.atHgp(15500.0); std::cout << "Atmos1.atHgp --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHgm(Atmos1.Hgm()); std::cout << "Atmos1.atHgm --- Hgm, P, rho:" << Atmos1.Hgm() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHpa(Atmos1.Hpa()); std::cout << "Atmos1.atHpa --- Hpa, P, rho:" << Atmos1.Hpa() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHda(Atmos1.Hda()); std::cout << "Atmos1.atHda --- Hda, P, rho:" << Atmos1.Hda() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;

   // full range tests

   std::cout << "    Hgm     Hgp     Hpa     Hda       T           rho             P       a         theta         sigma         delta         kappa          visc" << std::endl;
   std::cout << "    (m)     (m)     (m)     (m)     (K)      (kg/m^3)       (N/m^2)   (m/s)           (1)           (1)           (1)           (1)     (N*s/m^2)" << std::endl;

   for (size_t i = 0; i < Hgp_in.size(); i++)
   {
      Atmos1.atHgp(Hgp_in[i]);   // compute all properties, including alternate altitude definitions, at test point
      writeTestLine(Atmos1.Hgm(), Atmos1.Hgp(), Atmos1.Hpa(), Atmos1.Hda(), Atmos1.T(), Atmos1.rho(), Atmos1.P(), Atmos1.a(), Atmos1.theta(), Atmos1.sigma(), Atmos1.delta(), Atmos1.kappa(), Atmos1.visc());
   }

   // (2) dT offset, use dT from case (3) of test data, should match results of test data
   std::cout << std::endl << "Atmosphere Properties Model Test Data" << std::endl << std::endl;
   std::cout << "(2) dT offset" << std::endl;

   Atmos Atmos2(dT, AtmosParameters_si);   // standard + dT deviation added to temperature profile vs altitude
   
   // point tests, call with all altitude types, walk through types as input such that all results should be the same
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

   std::cout << "    Hgm     Hgp     Hpa     Hda       T           rho             P       a         theta         sigma         delta         kappa          visc" << std::endl;
   std::cout << "    (m)     (m)     (m)     (m)     (K)      (kg/m^3)       (N/m^2)   (m/s)           (1)           (1)           (1)           (1)     (N*s/m^2)" << std::endl;

   for (size_t i = 0; i < Hgp_in.size(); i++)
   {
      Atmos2.atHgp(Hgp_in[i]);   // compute all properties, including alternate altitude definitions, at test point
      writeTestLine(Atmos2.Hgm(), Atmos2.Hgp(), Atmos2.Hpa(), Atmos2.Hda(), Atmos2.T(), Atmos2.rho(), Atmos2.P(), Atmos2.a(), Atmos2.theta(), Atmos2.sigma(), Atmos2.delta(), Atmos2.kappa(), Atmos2.visc());
   }

   // (3) use standard lapse rates and SI units from AtmosParameters_si with an initial condition selected to match an output from (2), though not at an existing breakpoint (for generality). Since the intial condition does not change the temperature profile, the results should match those obtained from test (2) above.
   std::cout << std::endl << "Atmosphere Properties Model Test Data" << std::endl << std::endl;
   std::cout << "(3) custom intial conditions Hic, Tic, Pic; test with dT offset and standard lapse rates, should match results of test (2) above" << std::endl;

   // get intial condition point, say at 500 m
   Atmos2.atHgp(500);

   Hic = Atmos2.Hgp();
   Tic = Atmos2.T();
   Pic = Atmos2.P();

   Atmos Atmos3(Hic, Tic, Pic, AtmosParameters_si);   // standard lapse rates used to construct profiles through initial condition point

      // point tests, call with all altitude types, walk through types as input such that all results should be the same
   // troposphere
   Atmos3.atHgp(5500.0); std::cout << "Atmos3.atHgp --- Hgp, P, rho:" << Atmos3.Hgp() << " " << Atmos3.P() << " " << Atmos3.rho() << std::endl;
   Atmos3.atHgm(Atmos3.Hgm()); std::cout << "Atmos3.atHgm --- Hgm, P, rho:" << Atmos3.Hgm() << " " << Atmos3.P() << " " << Atmos3.rho() << std::endl;
   Atmos3.atHpa(Atmos3.Hpa()); std::cout << "Atmos3.atHpa --- Hpa, P, rho:" << Atmos3.Hpa() << " " << Atmos3.P() << " " << Atmos3.rho() << std::endl;
   Atmos3.atHda(Atmos3.Hda()); std::cout << "Atmos3.atHda --- Hda, P, rho:" << Atmos3.Hda() << " " << Atmos3.P() << " " << Atmos3.rho() << std::endl;
   // stratosphere
   Atmos3.atHgp(15500.0); std::cout << "Atmos3.atHgp --- Hgp, P, rho:" << Atmos3.Hgp() << " " << Atmos3.P() << " " << Atmos3.rho() << std::endl;
   Atmos3.atHgm(Atmos3.Hgm()); std::cout << "Atmos3.atHgm --- Hgm, P, rho:" << Atmos3.Hgm() << " " << Atmos3.P() << " " << Atmos3.rho() << std::endl;
   Atmos3.atHpa(Atmos3.Hpa()); std::cout << "Atmos3.atHpa --- Hpa, P, rho:" << Atmos3.Hpa() << " " << Atmos3.P() << " " << Atmos3.rho() << std::endl;
   Atmos3.atHda(Atmos3.Hda()); std::cout << "Atmos3.atHda --- Hda, P, rho:" << Atmos3.Hda() << " " << Atmos3.P() << " " << Atmos3.rho() << std::endl;
   // full range tests

   std::cout << "    Hgm     Hgp     Hpa     Hda       T           rho             P       a         theta         sigma         delta         kappa          visc" << std::endl;
   std::cout << "    (m)     (m)     (m)     (m)     (K)      (kg/m^3)       (N/m^2)   (m/s)           (1)           (1)           (1)           (1)     (N*s/m^2)" << std::endl;

   for (size_t i = 0; i < Hgp_in.size(); i++)
   {
      Atmos3.atHgp(Hgp_in[i]);   // compute all properties, including alternate altitude definitions, at test point
      writeTestLine(Atmos3.Hgm(), Atmos3.Hgp(), Atmos3.Hpa(), Atmos3.Hda(), Atmos3.T(), Atmos3.rho(), Atmos3.P(), Atmos3.a(), Atmos3.theta(), Atmos3.sigma(), Atmos3.delta(), Atmos3.kappa(), Atmos3.visc());
   }

//   // (4) custom initial conditions for H & P and custom temperature profile (T initial condition is then interpolated from custom profile)
//   // construct a custom profile that is equivalent in layer extents and lapse rates to one from above, though not at the same breakpoints, then the results must match.
   std::cout << std::endl << "Atmosphere Properties Model Test Data" << std::endl << std::endl;
   std::cout << "(4) custom intial conditions, altitude and temperature profile; input chosen such that results should match standard day test (1) above" << std::endl;

   std::vector<double> H_custom = { 0,5500,11000,15500,20000,26000,32000,39500,47000,49000,51000,61000,71000,77926,84852};
   std::vector<double> T_custom = {288.15,252.40,216.65,216.65,216.65,222.65,228.65,249.65,270.65,270.65,270.65,242.65,214.65,200.80,186.95};

   Atmos1.atHgp(500);

   Hic = Atmos1.Hgp();
   Pic = Atmos1.P();

   Atmos Atoms4(Hic, Pic, H_custom, T_custom, AtmosParameters_si);

      // point tests, call with all altitude types, walk through types as input such that all results should be the same
   // troposphere
   Atoms4.atHgp(5500.0); std::cout << "Atoms4.atHgp --- Hgp, P, rho:" << Atoms4.Hgp() << " " << Atoms4.P() << " " << Atoms4.rho() << std::endl;
   Atoms4.atHgm(Atoms4.Hgm()); std::cout << "Atoms4.atHgm --- Hgm, P, rho:" << Atoms4.Hgm() << " " << Atoms4.P() << " " << Atoms4.rho() << std::endl;
   Atoms4.atHpa(Atoms4.Hpa()); std::cout << "Atoms4.atHpa --- Hpa, P, rho:" << Atoms4.Hpa() << " " << Atoms4.P() << " " << Atoms4.rho() << std::endl;
   Atoms4.atHda(Atoms4.Hda()); std::cout << "Atoms4.atHda --- Hda, P, rho:" << Atoms4.Hda() << " " << Atoms4.P() << " " << Atoms4.rho() << std::endl;
   // stratosphere
   Atoms4.atHgp(15500.0); std::cout << "Atoms4.atHgp --- Hgp, P, rho:" << Atoms4.Hgp() << " " << Atoms4.P() << " " << Atoms4.rho() << std::endl;
   Atoms4.atHgm(Atoms4.Hgm()); std::cout << "Atoms4.atHgm --- Hgm, P, rho:" << Atoms4.Hgm() << " " << Atoms4.P() << " " << Atoms4.rho() << std::endl;
   Atoms4.atHpa(Atoms4.Hpa()); std::cout << "Atoms4.atHpa --- Hpa, P, rho:" << Atoms4.Hpa() << " " << Atoms4.P() << " " << Atoms4.rho() << std::endl;
   Atoms4.atHda(Atoms4.Hda()); std::cout << "Atoms4.atHda --- Hda, P, rho:" << Atoms4.Hda() << " " << Atoms4.P() << " " << Atoms4.rho() << std::endl;
   // full range tests

   std::cout << "    Hgm     Hgp     Hpa     Hda       T           rho             P       a         theta         sigma         delta         kappa          visc" << std::endl;
   std::cout << "    (m)     (m)     (m)     (m)     (K)      (kg/m^3)       (N/m^2)   (m/s)           (1)           (1)           (1)           (1)     (N*s/m^2)" << std::endl;

   for (size_t i = 0; i < Hgp_in.size(); i++)
   {
      Atoms4.atHgp(Hgp_in[i]);   // compute all properties, including alternate altitude definitions, at test point
      writeTestLine(Atoms4.Hgm(), Atoms4.Hgp(), Atoms4.Hpa(), Atoms4.Hda(), Atoms4.T(), Atoms4.rho(), Atoms4.P(), Atoms4.a(), Atoms4.theta(), Atoms4.sigma(), Atoms4.delta(), Atoms4.kappa(), Atoms4.visc());
   }

//   // (5) fully custom, initial conditions and temperature gradient (lapse rates). How to test? (1) construct conditions and profile equivalent to case (4); then results must match. Plus, to test that the custom input really do alter the results and the match in (1) is not some kind of bug, create another and different custom profile and conditions to show that the results do change with input.
   std::cout << std::endl << "Atmosphere Properties Model Test Data" << std::endl << std::endl;
   std::cout << "(5) fully custom, intial conditions, altitude breakpoints, and lapse rates; input chosen such that results should match dT offset results of test (2) above" << std::endl;

   // get intial condition point, say at 500 m
   Atmos2.atHgp(500);

   Hic = Atmos2.Hgp();
   Tic = Atmos2.T();
   Pic = Atmos2.P();

   Atmos Atmos5(Hic, Tic, Pic, Hj, Tgradj, AtmosParameters_si);
//
   // troposphere
   test_atH(Atmos5,  5500.0);
   //stratosphere
   test_atH(Atmos5, 15500.0);

   // full range test
   writeTestTable(Atmos5, Hgp_in);

}
