// testAtmos.cpp : Test Atmos++ library, Atmos class
//
// (1) test constructors and atmosphere property computation over full domain of altitudes
// (2) then test alternative altitude computations
// (3) test method operator(), should be the same as method at(Hgp), and at(Hgp) method is used by the other at methods: atHgp, atHgm, ... so at(HgP) will have been tested sufficiently by now

#include "stdafx.h"
#include <string>
#include <vector>
#include "Atmos.h"

int writeTestLine(double Hgm, double Hgp, double Hpa, double Hda, double T, double rho, double P, double sonic, double visc, double theta, double sigma, double delta, double sonicratio)
{
   std::cout << std::fixed << std::setw(8) << std::setprecision(1) << Hgm
      << " " << std::fixed << std::setw(8) << std::setprecision(1) << Hgp
      << " " << std::fixed << std::setw(8) << std::setprecision(1) << Hpa
      << " " << std::fixed << std::setw(8) << std::setprecision(1) << Hda
      << " " << std::scientific << std::setw(10) << std::setprecision(4) << T
      << " " << std::scientific << std::setw(10) << std::setprecision(4) << rho
      << " " << std::scientific << std::setw(10) << std::setprecision(4) << P
      << " " << std::scientific << std::setw(10) << std::setprecision(4) << sonic
      << " " << std::scientific << std::setw(10) << std::setprecision(4) << visc
      << " " << std::scientific << std::setw(10) << std::setprecision(4) << theta
      << " " << std::scientific << std::setw(10) << std::setprecision(4) << sigma
      << " " << std::scientific << std::setw(10) << std::setprecision(4) << delta
      << " " << std::scientific << std::setw(10) << std::setprecision(4) << sonicratio
      << std::endl;

   return(0);
}

void writeTestTable(std::string HdrUnits,Atmos atmos, std::vector<double> Hgp_in)
{
   std::cout << "     Hgm      Hgp      Hpa      Hda          T        rho          P          a       visc      theta      sigma      delta        kappa" << std::endl;
   std::cout << HdrUnits << std::endl;

   for (size_t i = 0; i < Hgp_in.size(); i++)
   {
      atmos.atHgp(Hgp_in[i]);   // compute all properties, including alternate altitude definitions, at test point
      writeTestLine(atmos.Hgm(), atmos.Hgp(), atmos.Hpa(), atmos.Hda(), atmos.T(), atmos.rho(), atmos.P(), atmos.a(), atmos.visc(), atmos.theta(), atmos.sigma(), atmos.delta(), atmos.kappa());
   }
}

void test_atH(std::string HdrUnits, Atmos atmos, double hgp_in)
{
   std::cout << "            Hgm      Hgp      Hpa      Hda          T        rho          P          a       visc      theta      sigma      delta      kappa" << std::endl;
   std::cout << "       " << HdrUnits << std::endl;
   atmos.atHgp(hgp_in); std::cout << "atHgp: "; writeTestLine(atmos.Hgm(), atmos.Hgp(), atmos.Hpa(), atmos.Hda(), atmos.T(), atmos.rho(), atmos.P(), atmos.a(), atmos.visc(), atmos.theta(), atmos.sigma(), atmos.delta(), atmos.kappa());
   atmos.atHgm(atmos.Hgm()); std::cout << "atHgm: ";  writeTestLine(atmos.Hgm(), atmos.Hgp(), atmos.Hpa(), atmos.Hda(), atmos.T(), atmos.rho(), atmos.P(), atmos.a(), atmos.visc(), atmos.theta(), atmos.sigma(), atmos.delta(), atmos.kappa());
   atmos.atHpa(atmos.Hpa()); std::cout << "atHpa: ";  writeTestLine(atmos.Hgm(), atmos.Hgp(), atmos.Hpa(), atmos.Hda(), atmos.T(), atmos.rho(), atmos.P(), atmos.a(), atmos.visc(), atmos.theta(), atmos.sigma(), atmos.delta(), atmos.kappa());
   atmos.atHda(atmos.Hda()); std::cout << "atHda: ";  writeTestLine(atmos.Hgm(), atmos.Hgp(), atmos.Hpa(), atmos.Hda(), atmos.T(), atmos.rho(), atmos.P(), atmos.a(), atmos.visc(), atmos.theta(), atmos.sigma(), atmos.delta(), atmos.kappa());
}


void testAtmos_si()
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
   std::vector<double> Hj = { 0.0,  2750.0, 11000.0, 13250.0, 20000.0, 23000.0, 32000.0, 35750.0, 47000.0, 48000.0, 51000.0, 56000.0, 71000.0, 74463.0, 84852.0 };
   std::vector<double> Tj = { 288.15,270.27,216.65,216.65,216.65,219.65,228.65,239.15,270.65,270.65,270.65,256.65,214.65,207.72,186.95 };
   std::vector<double> T30j = { 318.15,300.27,246.65,246.65,246.65,249.65,258.65,269.15,300.65,300.65,300.65,286.65,244.65,237.72,216.95 };
   std::vector<double> Tgradj = { -0.0065,-0.0065, 0.0000, 0.0000, 0.0010, 0.0010, 0.0028, 0.0028, 0.0000, 0.0000,-0.0028,-0.0028,-0.0020,-0.0020, 0.0 };

   // test results
   std::vector<double> Hgp_out(Hgp_in.size());
   std::vector<double> Hgm_out(Hgp_in.size());
   std::vector<double> Hpa_out(Hgp_in.size());
   std::vector<double> Hda_out(Hgp_in.size());
   std::vector<double> T_out(Hgp_in.size());

   std::string HdrUnits = "     (m)      (m)      (m)      (m)        (K)   (kg/m^3)    (N/m^2)      (m/s)  (N*s/m^2)        (1)        (1)        (1)          (1)";

   //---------------------------------------------------------------------------------------------------------------------------------------------------------
   std::cout << "Atmos++ Test Data - SI Units --- Atmosphere Model Properties validated against \"US Standard Atmosphere, 1976\"" << std::endl;

   // (B.1) Atmos(), standard atmosphere, test omitting AtmPrms, should generate std day object and results

   Atmos Atmos1(AtmosParameters_si);   // standard atmosphere, SI units

   std::cout << std::endl << "(B.1) Atmos()--- standard atmosphere; compare to 1976 US Standard Atmosphere, authoritative reference work" << std::endl;

   writeTestTable(HdrUnits,Atmos1, Hgp_in);

   // (B.2) Atmos(dT), dT offset

   Atmos Atmos2(dT, AtmosParameters_si);   // standard + dT deviation added to temperature profile vs altitude

   // for comparison to hand calc in linear thermal layer
   double H_B2a_L = 11000.0;
   Atmos2.at(H_B2a_L);
   double T_B2a_L = Atmos2.T();
   double P_B2a_L = Atmos2.P();

   // for comparison to had clac in isothermal layer
   double H_B2a_I = 15500.0;
   Atmos2.at(H_B2a_I);
   double T_B2a_I = Atmos2.T();
   double P_B2a_I = Atmos2.P();

   std::cout << std::endl;
   std::cout << "(B.2) Atmos(dT) --- standard + dT deviation added to temperature profile vs altitude" << std::endl;
   std::cout << std::endl;
   std::cout << "(B.2a) Comparison to hand calculations, one in linear layer, one in isothermal layer" << std::endl;
   std::cout << "Hand calculations followed equations and constants from \"US Standard Atmosphere, 1976\"" << std::endl;
   std::cout << "deviating only by application of temperature offset dT to temperature profile" << std::endl;
   std::cout << "" << std::endl;
   std::cout << "H, T, P only, since other results follow from these, as is verified in (B.1) above" << std::endl;
   std::cout << "       H          T          P" << std::endl;
   std::cout << "     (m)        (K)    (N/m^2)" << std::endl;
   std::cout << H_B2a_L << " " << T_B2a_L << " " << P_B2a_L << " --- linear layer with dT = " << dT << std::endl;
   std::cout << H_B2a_I << " " << T_B2a_I << " " << P_B2a_I << " --- isothermal layer with dT = " << dT << std::endl;
   std::cout << std::endl;

   std::cout << "(B.2b) Atmos(dT), dT offset, std+(" << dT << ") day, evaluated at test altitudes" << std::endl;
   std::cout << "Compare at H=11000 & 15500 m, T & P should be as above in (B.2a) \"hand calcs\"." << std::endl;
   writeTestTable(HdrUnits,Atmos2, Hgp_in);

   // (B.3) use standard lapse rates and SI units with an initial condition selected to match an output from (B.2b), though not at an existing breakpoint (for generality). Since the intial condition does not change the temperature profile, the results should match those obtained from test (B.2b) above.
   std::cout << std::endl;
   std::cout << "(B.3) Atmos(Hic, Tic, Pic), standard lapse rates used to construct profile through given initial conditions." << std::endl;
   std::cout << "      Chosen to create profile equivalent to dT offset case in (B.2b) above for comparison." << std::endl;

   // get intial condition point, say at 500 m
   Atmos2.atHgp(1000);

   Hic = Atmos2.Hgp();
   Tic = Atmos2.T();
   Pic = Atmos2.P();

   Atmos Atmos3(Hic, Tic, Pic);   // standard lapse rates used to construct profiles through initial condition point

   writeTestTable(HdrUnits,Atmos3, Hgp_in);

   // (B.4) custom initial conditions for H & P and custom temperature profile (T initial condition is then interpolated from custom profile)
   // construct a custom profile that is equivalent in layer extents and lapse rates to one from above, though not at the same breakpoints, then the results must match.
   std::cout << std::endl;
   std::cout << "(B.4) Atmos(Hic,Pic, Hj,Tj)" << std::endl;
   std::cout << "      Custom temperature profile constructed with given initial condition Hic, Pic with Tic interpolated from input (Hj,Tj)." << std::endl;
   std::cout << "      Initial conditions and temperature profile chosen to produce a profile equivalent to (B.2b) for comparison." << std::endl;

   Atmos2.atHgp(1000);

   Hic = Atmos2.Hgp();
   Pic = Atmos2.P();

   Atmos Atmos4(Hic, Pic, Hj, T30j, AtmosParameters_si);

   writeTestTable(HdrUnits,Atmos4, Hgp_in);

   // (B.5) fully custom, initial conditions and temperature gradient (lapse rates). How to test? (1) construct conditions and profile equivalent to case (4); then results must match. Plus, to test that the custom input really do alter the results and the match in (1) is not some kind of bug, create another and different custom profile and conditions to show that the results do change with input.
   std::cout << std::endl;
   std::cout << "(B.5) Atmos(Hic,Tic,Pic, Hj, Tgradj)" << std::endl;
   std::cout << "      custom temperature gradient profile (i.e. lapse rates) and initial conditions" << std::endl;
   std::cout << "      setup to produce profile equivalent to that of (B.2b) for comparison" << std::endl;

   // get intial condition point, say at 500 m
   Atmos2.atHgp(1000);

   Hic = Atmos2.Hgp();
   Tic = Atmos2.T();
   Pic = Atmos2.P();

   Atmos Atmos5(Hic, Tic, Pic, Hj, Tgradj, AtmosParameters_si);

   writeTestTable(HdrUnits,Atmos5, Hgp_in);

   // (C) ------------ Method Function Tests ---------------
   // check computation and access to results through all the method functions

   // use standard atmosphere for comparison to authoritative reference

   // since values over full domain have been checked above in (B.1), it suffices here to check at only one altitude
   double Hgp = 5500;

   // (C.1), (C.2), (C.3), evaluate atmosphere properties at given Hgp, these should produce the same results
      // test operator(), at(), atHgp() --- should give same results
   std::cout << std::endl << "Test Method Functions" << std::endl;
   Atmos1(Hgp); std::cout << "(C.1) Atmos1() --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.at(Hgp); std::cout << "(C.2) Atmos1.at --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHgp(Hgp); std::cout << "(C.3) Atmos1.atHgp --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;

   std::cout << "(C.4) Atmos1.Hgm(): " << Atmos1.Hgm() << std::endl;
   std::cout << "(C.5) Atmos1.Hpa(): " << Atmos1.Hpa() << std::endl;
   std::cout << "(C.6) Atmos1.Hda(): " << Atmos1.Hda() << std::endl;
   std::cout << "(C.7) Atmos1.T(): " << Atmos1.T() << std::endl;
   std::cout << "(C.8) Atmos1.rho(): " << Atmos1.rho() << std::endl;
   std::cout << "(C.9) Atmos1.P(): " << Atmos1.P() << std::endl;
   std::cout << "(C.10) Atmos1.a(): " << Atmos1.a() << std::endl;
   std::cout << "(C.11) Atmos1.visc(): " << Atmos1.visc() << std::endl;
   std::cout << "(C.12) Atmos1.theta(): " << Atmos1.theta() << std::endl;
   std::cout << "(C.13) Atmos1.sigma(): " << Atmos1.sigma() << std::endl;
   std::cout << "(C.14) Atmos1.delta(): " << Atmos1.delta() << std::endl;
   std::cout << "(C.15) Atmos1.kappa(): " << Atmos1.kappa() << std::endl;

   // test computation of altitude types and method functions atHgm, atHpa, atHda
   // start with given Hgp, the object should compute the other types, using these successively should generate then same results
   // use object with dT offset so that Hpa and Hda will differ from Hgp (they are the same in std atmosphere due to their definition)
   std::cout << std::endl << "(C.16), (C.17), (C.18), (C.19) evaluate at altitude of each type, using atmos with dT = " << dT << std::endl;
   std::cout << "with non-zero dT, altitudes Hpa and Hda differ from standard atmosphere Hgp." << std::endl;
   std::wcout << "types computed and evaluated successively (Hgp, Hgm, Hpa, Hda) such that four lines of results should be the same" << std::endl;
   test_atH(HdrUnits,Atmos2, Hgp);

}

void testAtmos_us()
{
   // test altitudes
   // test at: altitudes below table, at all breakpoints, between all breakpoints, above table
   std::vector<double> Hgp_in = { -328.1, 0.0, 18044.6, 36089.2, 50853.0,  65616.8,  85301.8, 104986.9, 129593.2, 154199.5, 160761.2, 167322.8, 200131.2, 232939.6, 262467.2, 278385.8, 295275.6 };
   std::vector<double> Hgm_in(Hgp_in.size());
   std::vector<double> Hpa_in(Hgp_in.size());
   std::vector<double> Hda_in(Hgp_in.size());
   std::vector<double> T_in(Hgp_in.size());

   // test initial conditions
   double Hic, Tic, Pic;
   // test dT offset
   double dT = 30.0*1.8;

   // test profiles
   std::vector<double> Hj = { 0.0,  9022.3, 36089.2, 43471.1, 65616.8, 75459.3,104986.9,117290.0,154199.5,157480.3,167322.8,183727.0,232939.6,244301.2,278385.8 };
   std::vector<double> Tj = { 518.67,486.50,389.97,389.97,389.97,395.37,411.57,430.47,487.17,487.17,487.17,461.97,386.37,373.90,336.51 };
   std::vector<double> T30j = { 572.67,540.50,443.97,443.97,443.97,449.37,465.57,484.47,541.17,541.17,541.17,515.97,440.37,427.90,390.51 };
   std::vector<double> Tgradj = { -0.00356616, -0.00356616,  0.00000000,  0.00000000,  0.00054864,  0.00054864,  0.00153619,  0.00153619,  0.00000000,  0.00000000, -0.00153619, -0.00153619, -0.00109712, -0.00109712,  0.00000000 };

   // test results
   std::vector<double> Hgp_out(Hgp_in.size());
   std::vector<double> Hgm_out(Hgp_in.size());
   std::vector<double> Hpa_out(Hgp_in.size());
   std::vector<double> Hda_out(Hgp_in.size());
   std::vector<double> T_out(Hgp_in.size());

   std::string HdrUnits = "    (ft)     (ft)     (ft)     (ft)        (R)  (sl/ft^3) (lbf/ft^2)     (ft/s) (lbf*s/ft^2)      (1)        (1)        (1)          (1)";

   //---------------------------------------------------------------------------------------------------------------------------------------------------------
   std::cout << "Atmos++ Test Data - US Customary Units --- Atmosphere Model Properties validated against \"US Standard Atmosphere, 1976\"" << std::endl;

   // (B.1) Atmos(), standard atmosphere, test omitting AtmPrms, should generate std day object and results

   Atmos Atmos1(AtmosParameters_us);   // standard atmosphere, SI units

   std::cout << std::endl << "(B.1) Atmos()--- standard atmosphere; compare to 1976 US Standard Atmosphere, authoritative reference work" << std::endl;

   writeTestTable(HdrUnits,Atmos1, Hgp_in);

   // (B.2) Atmos(dT), dT offset

   Atmos Atmos2(dT, AtmosParameters_us);   // standard + dT deviation added to temperature profile vs altitude

   // for comparison to hand calc in linear thermal layer
   double H_B2a_L = 36089.2;
   Atmos2.at(H_B2a_L);
   double T_B2a_L = Atmos2.T();
   double P_B2a_L = Atmos2.P();

   // for comparison to had clac in isothermal layer
   double H_B2a_I = 50853.0;
   Atmos2.at(H_B2a_I);
   double T_B2a_I = Atmos2.T();
   double P_B2a_I = Atmos2.P();

   std::cout << std::endl;
   std::cout << "(B.2) Atmos(dT) --- standard + dT deviation added to temperature profile vs altitude" << std::endl;
   std::cout << std::endl;
   std::cout << "(B.2a) Comparison to hand calculations, one in linear layer, one in isothermal layer" << std::endl;
   std::cout << "Hand calculations followed equations and constants from \"US Standard Atmosphere, 1976\"" << std::endl;
   std::cout << "deviating only by application of temperature offset dT to temperature profile" << std::endl;
   std::cout << "" << std::endl;
   std::cout << "H, T, P only, since other results follow from these, as is verified in (B.1) above" << std::endl;
   std::cout << "       H          T          P" << std::endl;
   std::cout << "     (m)        (K)    (N/m^2)" << std::endl;
   std::cout << H_B2a_L << " " << T_B2a_L << " " << P_B2a_L << " --- linear layer with dT = " << dT << std::endl;
   std::cout << H_B2a_I << " " << T_B2a_I << " " << P_B2a_I << " --- isothermal layer with dT = " << dT << std::endl;
   std::cout << std::endl;

   std::cout << "(B.2b) Atmos(dT), dT offset, std+(" << dT << ") day, evaluated at test altitudes" << std::endl;
   std::cout << "Compare at H=" << H_B2a_L << " & " << H_B2a_I << " ft, T & P should be as above in (B.2a) \"hand calcs\"." << std::endl;
   writeTestTable(HdrUnits,Atmos2, Hgp_in);

   // (B.3) use standard lapse rates and SI units with an initial condition selected to match an output from (B.2b), though not at an existing breakpoint (for generality). Since the intial condition does not change the temperature profile, the results should match those obtained from test (B.2b) above.
   std::cout << std::endl;
   std::cout << "(B.3) Atmos(Hic, Tic, Pic), standard lapse rates used to construct profile through given initial conditions." << std::endl;
   std::cout << "      Chosen to create profile equivalent to dT offset case in (B.2b) above for comparison." << std::endl;

   // get intial condition point, say at 1000 m
   Atmos2.atHgp(1000/0.3048);

   Hic = Atmos2.Hgp();
   Tic = Atmos2.T();
   Pic = Atmos2.P();

   Atmos Atmos3(Hic, Tic, Pic, AtmosParameters_us);   // standard lapse rates used to construct profiles through initial condition point

   writeTestTable(HdrUnits,Atmos3, Hgp_in);

   // (B.4) custom initial conditions for H & P and custom temperature profile (T initial condition is then interpolated from custom profile)
   // construct a custom profile that is equivalent in layer extents and lapse rates to one from above, though not at the same breakpoints, then the results must match.
   std::cout << std::endl;
   std::cout << "(B.4) Atmos(Hic,Pic, Hj,Tj)" << std::endl;
   std::cout << "      Custom temperature profile constructed with given initial condition Hic, Pic with Tic interpolated from input (Hj,Tj)." << std::endl;
   std::cout << "      Initial conditions and temperature profile chosen to produce a profile equivalent to (B.2b) for comparison." << std::endl;

   Atmos2.atHgp(1000 / 0.3048);

   Hic = Atmos2.Hgp();
   Pic = Atmos2.P();

   Atmos Atmos4(Hic, Pic, Hj, T30j, AtmosParameters_us);

   writeTestTable(HdrUnits,Atmos4, Hgp_in);

   // (B.5) fully custom, initial conditions and temperature gradient (lapse rates). How to test? (1) construct conditions and profile equivalent to case (4); then results must match. Plus, to test that the custom input really do alter the results and the match in (1) is not some kind of bug, create another and different custom profile and conditions to show that the results do change with input.
   std::cout << std::endl;
   std::cout << "(B.5) Atmos(Hic,Tic,Pic, Hj, Tgradj)" << std::endl;
   std::cout << "      custom temperature gradient profile (i.e. lapse rates) and initial conditions" << std::endl;
   std::cout << "      setup to produce profile equivalent to that of (B.2b) for comparison" << std::endl;

   // get intial condition point, say at 500 m
   Atmos2.atHgp(1000 / 0.3048);

   Hic = Atmos2.Hgp();
   Tic = Atmos2.T();
   Pic = Atmos2.P();

   Atmos Atmos5(Hic, Tic, Pic, Hj, Tgradj, AtmosParameters_us);

   writeTestTable(HdrUnits,Atmos5, Hgp_in);

   // (C) ------------ Method Function Tests ---------------
   // check computation and access to results through all the method functions

   // use standard atmosphere for comparison to authoritative reference

   // since values over full domain have been checked above in (B.1), it suffices here to check at only one altitude
   double Hgp = 5500 / 0.3048;

   // (C.1), (C.2), (C.3), evaluate atmosphere properties at given Hgp, these should produce the same results
      // test operator(), at(), atHgp() --- should give same results
   std::cout << std::endl << "Test Method Functions" << std::endl;
   Atmos1(Hgp); std::cout << "(C.1) Atmos1() --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.at(Hgp); std::cout << "(C.2) Atmos1.at --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;
   Atmos1.atHgp(Hgp); std::cout << "(C.3) Atmos1.atHgp --- Hgp, P, rho:" << Atmos1.Hgp() << " " << Atmos1.P() << " " << Atmos1.rho() << std::endl;

   std::cout << "(C.4) Atmos1.Hgm(): " << Atmos1.Hgm() << std::endl;
   std::cout << "(C.5) Atmos1.Hpa(): " << Atmos1.Hpa() << std::endl;
   std::cout << "(C.6) Atmos1.Hda(): " << Atmos1.Hda() << std::endl;
   std::cout << "(C.7) Atmos1.T(): " << Atmos1.T() << std::endl;
   std::cout << "(C.8) Atmos1.rho(): " << Atmos1.rho() << std::endl;
   std::cout << "(C.9) Atmos1.P(): " << Atmos1.P() << std::endl;
   std::cout << "(C.10) Atmos1.a(): " << Atmos1.a() << std::endl;
   std::cout << "(C.11) Atmos1.visc(): " << Atmos1.visc() << std::endl;
   std::cout << "(C.12) Atmos1.theta(): " << Atmos1.theta() << std::endl;
   std::cout << "(C.13) Atmos1.sigma(): " << Atmos1.sigma() << std::endl;
   std::cout << "(C.14) Atmos1.delta(): " << Atmos1.delta() << std::endl;
   std::cout << "(C.15) Atmos1.kappa(): " << Atmos1.kappa() << std::endl;

   // test computation of altitude types and method functions atHgm, atHpa, atHda
   // start with given Hgp, the object should compute the other types, using these successively should generate then same results
   // use object with dT offset so that Hpa and Hda will differ from Hgp (they are the same in std atmosphere due to their definition)
   std::cout << std::endl << "(C.16), (C.17), (C.18), (C.19) evaluate at altitude of each type, using atmos with dT = " << dT << std::endl;
   std::cout << "with non-zero dT, altitudes Hpa and Hda differ from standard atmosphere Hgp." << std::endl;
   std::wcout << "types computed and evaluated successively (Hgp, Hgm, Hpa, Hda) such that four lines of results should be the same" << std::endl;
   test_atH(HdrUnits,Atmos2, Hgp);

}

int main()
{

   testAtmos_si();

   testAtmos_us();

}
