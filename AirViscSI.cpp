inline double AirViscSI(double TdegK)
{
   // dyanmic(absolute) viscosity of air at T(deg K)
   // for moderate pressures, like atmospheric
   // from Sutherland's law as found in 1976 US Standard Atmosphere, page 19
   // https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539.pdf
   //
   // Note: this is the two coefficient form of Sutherland's law, its results will 
   // differ slightly from the three coefficient form found elsewhere. Its use is chosen
   // here for consistency with the US Standard Atmosphere, which is a widely used,
   // primary atmosphere property model.
   //
   // Input
   //    TdegK == = (K) air temperature
   //
   // Output
   //    visc == = (N*s / m ^ 2) absolute viscosity

   double bta = 1.458e-6;   // (kg / (s*m*K^0.5))
   double S = 110.4;        // (deg K) Sutherland temperature

   double visc = bta * pow(TdegK, 1.5) / (TdegK + S);

   return(visc);
}
