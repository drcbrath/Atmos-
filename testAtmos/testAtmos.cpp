// testAtmos.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include "Atmos.h"

int writeTestLine(double Hgm, double Hgp, double Hpa, double Hda, double dT, double T, double rho, double P, double sonic, double theta, double sigma, double delta, double sonicratio)
{
   std::cout << std::scientific << std::setprecision(6) << std::setw(13) << Hgm
      << std::setprecision(6) << std::setw(14) << Hgp
      << std::setprecision(6) << std::setw(14) << Hpa
      << std::setprecision(6) << std::setw(14) << Hda
      << std::setprecision(6) << std::setw(14) << dT
      << std::setprecision(6) << std::setw(14) << T
      << std::setprecision(6) << std::setw(14) << rho
      << std::setprecision(6) << std::setw(14) << P
      << std::setprecision(6) << std::setw(14) << sonic
      << std::setprecision(6) << std::setw(14) << theta
      << std::setprecision(6) << std::setw(14) << sigma
      << std::setprecision(6) << std::setw(14) << delta
      << std::setprecision(6) << std::setw(14) << sonicratio
      << std::endl;

   return(0);
}

int main()
{
   // 1. Validate basis of other functions, AtmSIRatios, on std day by comparison to independent authoritative reference, US Std Atm
   // 2. Verify & validate dT modified off std day by comparison to hand calculations in the two types of layers, linear thermal and isothermal. All other logic is unchanged, so dT=/=0 does not affect other aspects; hence demonstration at just two points suffices to validate dT=/=0 computation
   // 3. Validate dimensional function AtmSI results for T, rho, P, a, by comparison to US Std Atm std day. Compare off standard day to AtmSIRatios results multiplied by reference values. IN fact, this multiplication is all that AtmSI adds to AtmSIRatios---so it is somewhat trivial validation.
   // 4. The following eight functions use the above tested AtmSI for the property results, the added code to test is that which computes the various altitude types (geometric, pressure, density altitudes; Hgm, Hpa, Hda) corresponding to the geopotential (Hgp) and temperature. Using same Hgp input, the above results serve as the "should be" comparison data for the property (T,rho,P, ...) computation tests below. All eight functions do the same computations but with different altitude type and temperature type (dT constant or T at alt) as inputs, so given equivalent input they should give equivalent results. Geopotential altitude is the basic altitude type for the atmosphere models, so:
      // 4.1 Test AtmSI_HgpdT. Compare altitudes Hgm, Hpa, Hda to hand calculations (or other authoritative or carefully verified source). Once proven these outputs become the "should be" results for the remaining seven function tests.
      // Now, the altitude results Hgm, Hpa, Hda output from AtmSI_HgpdT are saved for use as input to the remaining seven functions, and the temperature is also saved for use as input to those functions with T input rather that dT (i.e. AtmSI_H??T())---these constitute the equivalent input, since they correspond to the Hgp and T used in the first tests.
      // 4.2 Test AtmSI_HgpT using T output above in (4.1), compare to (4.1)
      // 4.3 Test AtmSI_HgmdT using Hgm output above in (4.1), compare to (4.1)
      // 4.4 Test AtmSI_HgmT using Hgm and T output above in (4.1), compare to (4.1)
      // 4.5 Test AtmSI_HpadT using Hpa output above in (4.1), compare to (4.1)
      // 4.6 Test AtmSI_HpaT using Hpa and T output above in (4.1), compare to (4.1)
      // 4.7 Test AtmSI_HdadT using Hda output above in (4.1), compare to (4.1)
      // 4.8 Test AtmSI_HdaT using Hda and T output above in (4.1), compare to (4.1)
   
   // test point inputs
   double dTcold = -20;
   double dTstd = 0;
   double dThot = 20;

   // std, hot, cold
   std::vector<double> dTs = { 0.0, 20.0, -20.0};

   // test altitudes
      // test at: altitudes below table, at all breakpoints, between all breakpoints, above table
   std::vector<double> Hgp_in = { -305.0, 0.0,  5500.0, 11000.0, 15500.0,  20000.0, 26000.0, 32000.0, 39500.0, 47000.0, 49000.0, 51000.0, 61000.0, 71000.0, 77927.0, 84852.0, 90000.0 };
   std::vector<double> Hgm_in(Hgp_in.size());
   std::vector<double> Hpa_in(Hgp_in.size());
   std::vector<double> Hda_in(Hgp_in.size());
   std::vector<double> T_in(Hgp_in.size());

   // test results
   double Hpa, Hda, dT_out, T, rho, P, a, theta, sigma, delta, kappa;
   std::vector<double> Hgp_out(Hgp_in.size());
   std::vector<double> Hgm_out(Hgp_in.size());
   std::vector<double> Hpa_out(Hgp_in.size());
   std::vector<double> Hda_out(Hgp_in.size());
   std::vector<double> T_out(Hgp_in.size());

//------- tests -------

// (1) & (2) validate AtmSIRatios

   for (int m = 0; m < dTs.size(); m++)
   {
      std::cout << "AtmSIRatios(Hgp,dT=" << dTs[m] << ")" << std::endl;
      std::cout << "        Hgp(m)          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
      for (int j = 0; j <= Hgp_in.size(); j++)
      {
         AtmSIRatios(Hgp_in[m], dTstd, theta, sigma, delta, kappa);
         for (int m = 0; m < Hgp_in.size(); m++)
         {
            AtmSIRatios(Hgp_in[m], dTs[m], theta, sigma, delta, kappa);
            std::cout << std::scientific << std::setprecision(6) << std::setw(14) << Hgp_in[m]
               << std::setprecision(6) << std::setw(14) << theta
               << std::setprecision(6) << std::setw(14) << sigma
               << std::setprecision(6) << std::setw(14) << delta
               << std::setprecision(6) << std::setw(14) << kappa
               << std::endl;
         }
         std::cout << std::endl << std::endl;
      }
   }

// (3) validate AtmSI
// calls AtmSIRatio (tested above), results differ only by multiplication by reference values, so it would  suffice to validate only at std day conditions, but other dTs are also computed to test coverage of inputs

   for (int m = 0; m < dTs.size(); m++)
   {
      std::cout << "AtmSI(Hgp,dT=" << dTs[m] << ")" << std::endl;
      std::cout << "        Hgp(m)             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
      for (int j = 0; j < Hgp_in.size(); j++)
      {
         AtmSI(Hgp_in[m], dTs[m], T, rho, P, a, theta, sigma, delta, kappa);
         std::cout << std::scientific << std::setprecision(6) << std::setw(14) << Hgp_in[m]
            << std::setprecision(6) << std::setw(14) << T
            << std::setprecision(6) << std::setw(14) << rho
            << std::setprecision(6) << std::setw(14) << P
            << std::setprecision(6) << std::setw(14) << a
            << std::setprecision(6) << std::setw(14) << theta
            << std::setprecision(6) << std::setw(14) << sigma
            << std::setprecision(6) << std::setw(14) << delta
            << std::setprecision(6) << std::setw(14) << kappa
            << std::endl;
      }
   }

   //------- tests -------

	// AtmSIRatio
	// dT for std, hot, cold
   std::cout << "AtmSIRatio(Hgp,dT=" << dTstd << "), std" << std::endl;
   std::cout << "        Hgp(m)          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hgp_in.size(); m++)
   {
      AtmSIRatios(Hgp_in[m], dTstd, theta, sigma, delta, kappa);
      std::cout << std::scientific << std::setprecision(6) << std::setw(14) << Hgp_in[m]
         << std::setprecision(6) << std::setw(14) << theta
         << std::setprecision(6) << std::setw(14) << sigma
         << std::setprecision(6) << std::setw(14) << delta
         << std::setprecision(6) << std::setw(14) << kappa
         << std::endl;
   }
   std::cout << std::endl << std::endl;

   std::cout << "AtmSIRatio(Hgp,dT=" << dThot << "), hot" << std::endl;
   std::cout << "        Hgp(m)          T/T0      rho/rho0         P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hgp_in.size(); m++)
   {
      AtmSIRatios(Hgp_in[m], dThot, theta, sigma, delta, kappa);
      std::cout << std::scientific << std::setprecision(6) << std::setw(14) << Hgp_in[m]
         << std::setprecision(6) << std::setw(14) << theta
         << std::setprecision(6) << std::setw(14) << sigma
         << std::setprecision(6) << std::setw(14) << delta
         << std::setprecision(6) << std::setw(14) << kappa
         << std::endl;
   }
   std::cout << std::endl << std::endl;

   std::cout << "AtmSIRatio(Hgp,dT=" << dTcold << "), cold" << std::endl;
   std::cout << "        Hgp(m)          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hgp_in.size(); m++)
   {
      AtmSIRatios(Hgp_in[m], dTcold, theta, sigma, delta, kappa);
      std::cout << std::scientific << std::setprecision(6) << std::setw(14) << Hgp_in[m]
         << std::setprecision(6) << std::setw(14) << theta
         << std::setprecision(6) << std::setw(14) << sigma
         << std::setprecision(6) << std::setw(14) << delta
         << std::setprecision(6) << std::setw(14) << kappa
         << std::endl;
   }
   std::cout << std::endl;

	// AtmSI
	// calls AtmSIRatio (tested above), results differ only by multiplication by reference values, so it suffices to test only at std day conditions
   std::cout << "AtmSI(Hgp,dT=" << dTstd << "), std" << std::endl;
   std::cout << "        Hgp(m)             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hgp_in.size(); m++)
   {
      AtmSI(Hgp_in[m], dTstd, T, rho, P, a, theta, sigma, delta, kappa);
      std::cout << std::scientific << std::setprecision(6) << std::setw(14) << Hgp_in[m]
         << std::setprecision(6) << std::setw(14) << T
         << std::setprecision(6) << std::setw(14) << rho
         << std::setprecision(6) << std::setw(14) << P
         << std::setprecision(6) << std::setw(14) << a
         << std::setprecision(6) << std::setw(14) << theta
         << std::setprecision(6) << std::setw(14) << sigma
         << std::setprecision(6) << std::setw(14) << delta
         << std::setprecision(6) << std::setw(14) << kappa
         << std::endl;
   }
   std::cout << std::endl << std::endl;

// (4) validate the eight AtmSI_H... functions
  // The following eight functions use the above tested AtmSI for the property results, the added code to test is that which computes the various altitude types (geometric, pressure, density altitudes; Hgm, Hpa, Hda) corresponding to the geopotential (Hgp) and temperature. Using same Hgp input, the above results can serve as the "should be" comparison data for the property (T,rho,P, ...) computation tests below.

   // All eight functions do the same computations but with different altitude type and temperature type (dT constant or T at alt) as inputs, so given equivalent input they should give equivalent results. Geopotential altitude is the basic altitude type for the atmosphere models, so test AtmSI_HgpdT & AtmSI_HgpT first. Compare altitudes Hgm, Hpa, Hda to hand calculations (or other authoritative or carefully verified source). Once proven these outputs become the "should be" results for the remaining six function tests.
   // So, the altitude results Hgm, Hpa, Hda output from AtmSI_HgpdT are saved for use as input to the remaining six functions---these constitute the equivalent input, since they correspond to the Hgp used in the first two tests.

// loop through dTs, testing the eight functions in sequence at each---all eight should generate the same results for the same dT (i.e. all std should match, then all hot, then all cold)

   for (int n = 0; n < dTs.size(); n++)
   {
      double dT = dTs[n];
      // 4.1 Test AtmSI_HgpdT. Compare altitudes Hgm, Hpa, Hda to hand calculations (or other authoritative or carefully verified source). Once proven these outputs become the "should be" results for the remaining seven function tests.
         // AtmSI_HgpdT -- Hgp,dT input
         // for geometric altitude, Hgm=Hgp=0; but above alt 0, then Hgm > Hgp
         // for std day, Hgp=Hpa=Hda
         // for hot day Hda > Hgp, Hpa < Hgp
         // for cold day Hda < Hgp, Hpa > Hgp
      std::cout << "AtmSI_HgpdT(Hgp,dT=" << dT << ")" << std::endl;
      std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
      for (int m = 0; m < Hgp_in.size(); m++)
      {
         AtmSI_HgpdT(Hgm_out[m], Hgp_in[m], Hpa_out[m], Hda_out[m], dT, T_out[m], rho, P, a, theta, sigma, delta, kappa);
         writeTestLine(Hgm_out[m], Hgp_in[m], Hpa_out[m], Hda_out[m], dT, T_out[m], rho, P, a, theta, sigma, delta, kappa);
      }
      std::cout << std::endl << std::endl;

      // save output for use as inputs to the tests which follow, as described above
      Hgm_in = Hgm_out;
      Hpa_in = Hpa_out;
      Hda_in = Hda_out;
      T_in = T_out;

      // 4.2 Test AtmSI_HgpT using T output above in (4.1), compare to (4.1)
         // AtmSI_HgpT  -- Hgp,T input
         // differs from AtmSI_HgpdT only by input of T instead of dT; dT is computed from subtraction std day T (which is computed by AtmSI already tested). The test then uses the T output from AtmSI_HgpdT tested above, so the inputs are equivalent and results should be the same.
      std::cout << "AtmSI_HgpT(Hgp,T)" << std::endl;
      std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
      for (int m = 0; m < Hgp_in.size(); m++)
      {
         AtmSI_HgpT(Hgm_out[m], Hgp_in[m], Hpa_out[m], Hda_out[m], dT_out, T_in[m], rho, P, a, theta, sigma, delta, kappa);
         writeTestLine(Hgm_out[m], Hgp_in[m], Hpa_out[m], Hda_out[m], dT, T_in[m], rho, P, a, theta, sigma, delta, kappa);
      }
      std::cout << std::endl << std::endl;

      // 4.3 Test AtmSI_HgmdT using Hgm output above in (4.1), compare to (4.1)
         // AtmSI_HgmdT -- Hgm,dT input
         // very, very similar to AtmSI_HgpdT above, differs by input of Hgm instead of Hgp, order of computation changed from AtmSI_HgpdT, but otherwise is the same
         // using the geometric altitude output from test of AtmSI_HgpdT as input here should then return all the same results, including Hgp_out == Hgp_in
      std::cout << "AtmSI_HgmdT(Hgm,dT=" << dT << ")" << std::endl;
      std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
      for (int m = 0; m < Hgm_in.size(); m++)
      {
         AtmSI_HgmdT(Hgm_in[m], Hgp_out[m], Hpa_out[m], Hda_out[m], dT, T_out[m], rho, P, a, theta, sigma, delta, kappa);
         writeTestLine(Hgm_in[m], Hgp_out[m], Hpa_out[m], Hda_out[m], dT, T_out[m], rho, P, a, theta, sigma, delta, kappa);
      }
      std::cout << std::endl << std::endl;

      // 4.4 Test AtmSI_HgmT using Hgm and T output above in (4.1), compare to (4.1)
         // AtmSI_HgmT  -- Hgm,dT input
         // differs from AtmSI_HgmdT only by input of T instead of d; dT is computed from subtraction std day T (which is computed by AtmSI already tested). The test then uses the T output from AtmSI_HgmdT tested above, so the inputs are equivalent and results should be the same.
      std::cout << "AtmSI_HgmT(Hgm,T)" << std::endl;
      std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
      for (int m = 0; m < Hgm_in.size(); m++)
      {
         AtmSI_HgmT(Hgm_in[m], Hgp_out[m], Hpa_out[m], Hda_out[m], dT_out, T_in[m], rho, P, a, theta, sigma, delta, kappa);
         writeTestLine(Hgm_in[m], Hgp_out[m], Hpa_out[m], Hda_out[m], dT_out, T_in[m], rho, P, a, theta, sigma, delta, kappa);
      }
      std::cout << std::endl << std::endl;

      // 4.5 Test AtmSI_HpadT using Hpa output above in (4.1), compare to (4.1)
         // AtmSI_HpadT -- Hpa,dT input
         // adds computation of altitudes, all other results from AtmSI, need only add test altitude computations
      std::cout << "AtmSI_HpadT(Hpa,dT=" << dT << ")" << std::endl;
      std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
      for (int m = 0; m < Hpa_in.size(); m++)
      {
         AtmSI_HpadT(Hgm_out[m], Hgp_out[m], Hpa_in[m], Hda_out[m], dT, T_out[m], rho, P, a, theta, sigma, delta, kappa);
         writeTestLine(Hgm_out[m], Hgp_out[m], Hpa_in[m], Hda_out[m], dT, T_out[m], rho, P, a, theta, sigma, delta, kappa);
      }
      std::cout << std::endl << std::endl;

   // 4.6 Test AtmSI_HpaT using Hpa and T output above in (4.1), compare to (4.1)
      // AtmSI_HpaT  -- Hpa,T input
      // differs from AtmSI_HpaT only by input of T instead of d; dT is computed from subtraction std day T (which is computed by AtmSI already tested). The test then uses the T output from AtmSI_HpadT tested above, so the inputs are equivalent and results should be the same.
      std::cout << "AtmSI_HpaT(Hgm,T)" << std::endl;
      std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
      for (int m = 0; m < Hpa_in.size(); m++)
      {
         AtmSI_HpaT(Hgm_out[m], Hgp_out[m], Hpa_in[m], Hda_out[m], dT_out, T_in[m], rho, P, a, theta, sigma, delta, kappa);
         writeTestLine(Hgm_out[m], Hgp_out[m], Hpa_in[m], Hda_out[m], dT_out, T_in[m], rho, P, a, theta, sigma, delta, kappa);
      }
      std::cout << std::endl << std::endl;

      // 4.7 Test AtmSI_HdadT using Hda output above in (4.1), compare to (4.1)
         // AtmSI_HdadT -- Hds,dT input
         // adds computation of altitudes, all other results from AtmSI, need only add test altitude computations
      std::cout << "AtmSI_HdadT(Hpa,dT=" << dT << ")" << std::endl;
      std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
      for (int m = 0; m < Hda_in.size(); m++)
      {
         AtmSI_HdadT(Hgm_out[m], Hgp_out[m], Hpa_out[m], Hda_in[m], dT, T_out[m], rho, P, a, theta, sigma, delta, kappa);
         writeTestLine(Hgm_out[m], Hgp_out[m], Hpa_out[m], Hda_in[m], dT, T_out[m], rho, P, a, theta, sigma, delta, kappa);
      }
      std::cout << std::endl << std::endl;

      // 4.8 Test AtmSI_HdaT using Hda and T output above in (4.1), compare to (4.1)
         // AtmSI_HdaT  -- Hda,dT input
         // differs from AtmSI_HdadT only by input of T instead of d; dT is computed from subtraction std day T (which is computed by AtmSI already tested). The test then uses the T output from AtmSI_HdadT tested above, so the inputs are equivalent and results should be the same.
      std::cout << "AtmSI_HpaT(Hgm,T)" << std::endl;
      std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
      for (int m = 0; m < Hda_in.size(); m++)
      {
         AtmSI_HpaT(Hgm_out[m], Hgp_out[m], Hpa_out[m], Hda_in[m], dT_out, T_in[m], rho, P, a, theta, sigma, delta, kappa);
         writeTestLine(Hgm_out[m], Hgp_out[m], Hpa_out[m], Hda_in[m], dT_out, T_in[m], rho, P, a, theta, sigma, delta, kappa);
      }
      std::cout << std::endl << std::endl;

   } // end of loops through dTs comparing results of tests (4.1) through (4.8)

//---------------------------------------------------------------------------------

   // These functions are overloaded versions of the above. The return results are supplied in an Atm type variable, so the argument list is shorter. Otherwise the functions, their operation, and results are the same as those above. In fact, internally each of these functions works by calling the above functions and placing the results into the return Atm variable. So once the above functions are proven, these functions results' can be directly compared with the above of the same name for validation.

   // Atm AtmSI_HgpdT(double Hgp, double dT = 0.0)
   std::cout << "AtmSI_HgpdT(Hgp,dT=" << dTstd << "), std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hgp_in.size(); m++)
   {
      Atm atmsiHgpdT = AtmSI_HgpdT(Hgp_in[m], dTstd);
      writeTestLine(atmsiHgpdT.Hgm, atmsiHgpdT.Hgp, atmsiHgpdT.Hpa, atmsiHgpdT.Hda, atmsiHgpdT.dT, atmsiHgpdT.T,
         atmsiHgpdT.rho, atmsiHgpdT.P, atmsiHgpdT.a, atmsiHgpdT.theta, atmsiHgpdT.sigma, atmsiHgpdT.delta, atmsiHgpdT.kappa);
   }
   std::cout << std::endl << std::endl;
   // without dT
   std::cout << "AtmSI_HgpdT(Hgp), default to std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hgp_in.size(); m++)
   {
      Atm atmsiHgp = AtmSI_Hgp(Hgp_in[m]);
      writeTestLine(atmsiHgp.Hgm, atmsiHgp.Hgp, atmsiHgp.Hpa, atmsiHgp.Hda, atmsiHgp.dT, atmsiHgp.T,
         atmsiHgp.rho, atmsiHgp.P, atmsiHgp.a, atmsiHgp.theta, atmsiHgp.sigma, atmsiHgp.delta, atmsiHgp.kappa);
   }
   std::cout << std::endl << std::endl;

   // Atm AtmSI_HgpT(double Hgp, double T = 288.15);
   // with T
   std::cout << "AtmSI_HgpT(Hgp,T=" << T0 << "), std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hgp_in.size(); m++)
   {
      Atm atmsiHgpT = AtmSI_HgpT(Hgp_in[m], T_in[m]);
      writeTestLine(atmsiHgpT.Hgm, atmsiHgpT.Hgp, atmsiHgpT.Hpa, atmsiHgpT.Hda, atmsiHgpT.dT, atmsiHgpT.T,
         atmsiHgpT.rho, atmsiHgpT.P, atmsiHgpT.a, atmsiHgpT.theta, atmsiHgpT.sigma, atmsiHgpT.delta, atmsiHgpT.kappa);
   }
   std::cout << std::endl << std::endl;

   // Atm AtmSI_HgmdT(double Hgm, double dT = 0.0);
   // with dT
   std::cout << "AtmSI_HgmdT(Hgp,dT=" << dTstd << "), std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hgm_in.size(); m++)
   {
      Atm atmsiHgmdT = AtmSI_HgmdT(Hgm_in[m], dTstd);
      writeTestLine(atmsiHgmdT.Hgm, atmsiHgmdT.Hgp, atmsiHgmdT.Hpa, atmsiHgmdT.Hda, atmsiHgmdT.dT, atmsiHgmdT.T,
         atmsiHgmdT.rho, atmsiHgmdT.P, atmsiHgmdT.a, atmsiHgmdT.theta, atmsiHgmdT.sigma, atmsiHgmdT.delta, atmsiHgmdT.kappa);
   }
   std::cout << std::endl << std::endl;
   // without dT
   std::cout << "AtmSI_HgmdT(Hgp), default to std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hgm_in.size(); m++)
   {
      Atm atmsiHgmdT = AtmSI_HgmdT(Hgp_in[m]);
      writeTestLine(atmsiHgmdT.Hgm, atmsiHgmdT.Hgp, atmsiHgmdT.Hpa, atmsiHgmdT.Hda, atmsiHgmdT.dT, atmsiHgmdT.T,
         atmsiHgmdT.rho, atmsiHgmdT.P, atmsiHgmdT.a, atmsiHgmdT.theta, atmsiHgmdT.sigma, atmsiHgmdT.delta, atmsiHgmdT.kappa);
   }
   std::cout << std::endl << std::endl;

   // Atm AtmSI_HgmT(double Hgm, double T = 288.15);
   // with T
   std::cout << "AtmSI_HgmT(Hgm,T=" << T0 << ")" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hgm_in.size(); m++)
   {
      Atm atmsiHgmdT = AtmSI_HgmT(Hgp_in[m], T_in[m]);
      writeTestLine(atmsiHgmdT.Hgm, atmsiHgmdT.Hgp, atmsiHgmdT.Hpa, atmsiHgmdT.Hda, atmsiHgmdT.dT, atmsiHgmdT.T,
         atmsiHgmdT.rho, atmsiHgmdT.P, atmsiHgmdT.a, atmsiHgmdT.theta, atmsiHgmdT.sigma, atmsiHgmdT.delta, atmsiHgmdT.kappa);
   }
   std::cout << std::endl << std::endl;
   // without T
   std::cout << "AtmSI_HgmT(Hgm), default to std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hgm_in.size(); m++)
   {
      Atm atmsiHgmdT = AtmSI_HgmT(Hgp_in[m]);
      writeTestLine(atmsiHgmdT.Hgm, atmsiHgmdT.Hgp, atmsiHgmdT.Hpa, atmsiHgmdT.Hda, atmsiHgmdT.dT, atmsiHgmdT.T,
         atmsiHgmdT.rho, atmsiHgmdT.P, atmsiHgmdT.a, atmsiHgmdT.theta, atmsiHgmdT.sigma, atmsiHgmdT.delta, atmsiHgmdT.kappa);
   }
   std::cout << std::endl << std::endl;

   // Atm AtmSI_HpadT(double Hpa, double dT = 0.0);
   // with dT
   std::cout << "AtmSI_HpadT(Hpa,dT=" << dTstd << "), std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hpa_in.size(); m++)
   {
      Atm atmsiHpadT = AtmSI_HpadT(Hpa_in[m], dTstd);
      writeTestLine(atmsiHpadT.Hgm, atmsiHpadT.Hgp, atmsiHpadT.Hpa, atmsiHpadT.Hda, atmsiHpadT.dT, atmsiHpadT.T,
         atmsiHpadT.rho, atmsiHpadT.P, atmsiHpadT.a, atmsiHpadT.theta, atmsiHpadT.sigma, atmsiHpadT.delta, atmsiHpadT.kappa);
   }
   std::cout << std::endl << std::endl;
   // without dT
   std::cout << "AtmSI_HpadT(Hpa), default to std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hpa_in.size(); m++)
   {
      Atm atmsiHpadT = AtmSI_HpadT(Hpa_in[m]);
      writeTestLine(atmsiHpadT.Hgm, atmsiHpadT.Hgp, atmsiHpadT.Hpa, atmsiHpadT.Hda, atmsiHpadT.dT, atmsiHpadT.T,
         atmsiHpadT.rho, atmsiHpadT.P, atmsiHpadT.a, atmsiHpadT.theta, atmsiHpadT.sigma, atmsiHpadT.delta, atmsiHpadT.kappa);
   }
   std::cout << std::endl << std::endl;

   // Atm AtmSI_HpaT(double Hpa, double T = 288.15);
   // with T
   std::cout << "AtmSI_HpaT(Hgp,T=" << T0 << "), std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hpa_in.size(); m++)
   {
      Atm atmsiHpaT = AtmSI_HpaT(Hgp_in[m], T_in[m]);
      writeTestLine(atmsiHpaT.Hgm, atmsiHpaT.Hgp, atmsiHpaT.Hpa, atmsiHpaT.Hda, atmsiHpaT.dT, atmsiHpaT.T,
         atmsiHpaT.rho, atmsiHpaT.P, atmsiHpaT.a, atmsiHpaT.theta, atmsiHpaT.sigma, atmsiHpaT.delta, atmsiHpaT.kappa);
   }
   std::cout << std::endl << std::endl;
   // without T
   std::cout << "AtmSI_HpaT(Hgp), default to std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hpa_in.size(); m++)
   {
      Atm atmsiHpaT = AtmSI_HpaT(Hgp_in[m]);
      writeTestLine(atmsiHpaT.Hgm, atmsiHpaT.Hgp, atmsiHpaT.Hpa, atmsiHpaT.Hda, atmsiHpaT.dT, atmsiHpaT.T,
         atmsiHpaT.rho, atmsiHpaT.P, atmsiHpaT.a, atmsiHpaT.theta, atmsiHpaT.sigma, atmsiHpaT.delta, atmsiHpaT.kappa);
   }
   std::cout << std::endl << std::endl;

   // Atm AtmSI_HdadT(double Hda, double dT = 0.0);
   // with dT
   std::cout << "AtmSI_HdadT(Hgp,dT=" << dTstd << "), std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hda_in.size(); m++)
   {
      Atm atmsiHdadT = AtmSI_HdadT(Hgp_in[m], dTstd);
      writeTestLine(atmsiHdadT.Hgm, atmsiHdadT.Hgp, atmsiHdadT.Hpa, atmsiHdadT.Hda, atmsiHdadT.dT, atmsiHdadT.T,
         atmsiHdadT.rho, atmsiHdadT.P, atmsiHdadT.a, atmsiHdadT.theta, atmsiHdadT.sigma, atmsiHdadT.delta, atmsiHdadT.kappa);
   }
   std::cout << std::endl << std::endl;
   // without dT
   std::cout << "AtmSI_HdadT(Hgp), default to std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hda_in.size(); m++)
   {
      Atm atmsiHdadT = AtmSI_HdadT(Hgp_in[m]);
      writeTestLine(atmsiHdadT.Hgm, atmsiHdadT.Hgp, atmsiHdadT.Hpa, atmsiHdadT.Hda, atmsiHdadT.dT, atmsiHdadT.T,
         atmsiHdadT.rho, atmsiHdadT.P, atmsiHdadT.a, atmsiHdadT.theta, atmsiHdadT.sigma, atmsiHdadT.delta, atmsiHdadT.kappa);
   }
   std::cout << std::endl << std::endl; 

   // Atm AtmSI_HdaT(double Hda, double T = 288.15);
   // with T
   std::cout << "AtmSI_HdaT(Hgp,T=" << T0 << "), std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hda_in.size(); m++)
   {
      Atm atmsiHdaT = AtmSI_HdaT(Hgp_in[m], T_in[m]);
      writeTestLine(atmsiHdaT.Hgm, atmsiHdaT.Hgp, atmsiHdaT.Hpa, atmsiHdaT.Hda, atmsiHdaT.dT, atmsiHdaT.T,
         atmsiHdaT.rho, atmsiHdaT.P, atmsiHdaT.a, atmsiHdaT.theta, atmsiHdaT.sigma, atmsiHdaT.delta, atmsiHdaT.kappa);
   }
   std::cout << std::endl << std::endl;
   // without T
   std::cout << "AtmSI_HdaT(Hgp), default to std day" << std::endl;
   std::cout << "       Hgm(m)        Hgp(m)        Hpa(m)        Hda(m)            dT             T           rho             P             a          T/T0      rho/rho0          P/P0          a/a0" << std::endl;
   for (int m = 0; m < Hda_in.size(); m++)
   {
      Atm atmsiHdaT = AtmSI_HdaT(Hgp_in[m]);
      writeTestLine(atmsiHdaT.Hgm, atmsiHdaT.Hgp, atmsiHdaT.Hpa, atmsiHdaT.Hda, atmsiHdaT.dT, atmsiHdaT.T,
         atmsiHdaT.rho, atmsiHdaT.P, atmsiHdaT.a, atmsiHdaT.theta, atmsiHdaT.sigma, atmsiHdaT.delta, atmsiHdaT.kappa);
   }
   std::cout << std::endl << std::endl;


	return(0);
}
